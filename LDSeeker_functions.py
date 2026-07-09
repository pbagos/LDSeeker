#!/usr/bin/env python
import os
import re
import gc
import math
import shutil
import polars as pl
import numpy as np
from numba import njit
from statistics import NormalDist

pl.Config.set_engine_affinity("streaming")


# -------------------------------------------------------------------
# Helper: Unify target SNPs entirely natively (No Python sets)
# -------------------------------------------------------------------
def _get_snp_series(rs_series: pl.Series, imp_snp_list: list | None) -> pl.Series:
    """Combines GWAS SNPs and Imputed SNPs natively in Rust/Polars."""
    if imp_snp_list:
        imp_series = pl.Series(imp_snp_list, dtype=rs_series.dtype).drop_nulls()
        return pl.concat([rs_series, imp_series]).unique(maintain_order=False)
    return rs_series.unique(maintain_order=False)


# -------------------------------------------------------------------
# Significance filtering of the INPUT variants (Wald z = BETA/SE,
# two-tailed p from z). Lets the user keep either the significant or
# the non-significant variants of the study file before LD retrieval.
# -------------------------------------------------------------------
_STD_NORM = NormalDist(0.0, 1.0)


def two_tailed_p_from_z(z: float) -> float:
    """
    Two-tailed Wald p-value from a z-score, numerically stable in the tail.

        p = 2 * (1 - Phi(|z|)) = erfc(|z| / sqrt(2))

    Using erfc avoids the catastrophic underflow of ``1 - Phi(|z|)`` (which
    rounds to 0 for |z| beyond ~8 in float64); erfc stays accurate down to
    ~1e-300. Scalar utility for annotating variants, not used in the hot loop.
    """
    if z is None or (isinstance(z, float) and math.isnan(z)):
        return float("nan")
    return math.erfc(abs(z) / math.sqrt(2.0))


def z_critical_from_p(p_threshold: float) -> float:
    """
    Two-tailed critical |z| for a p-value threshold:  z_c = Phi^{-1}(1 - p/2).

    Filtering ``two-tailed p < alpha`` is *exactly* equivalent to ``|z| > z_c``.
    We convert the threshold once (scalar, via the stdlib inverse-normal) and
    then filter |z| natively in Polars, which is both faster and free of any
    tail underflow that per-row p-value computation would suffer.
    """
    if not (0.0 < p_threshold < 1.0):
        raise ValueError("p_threshold must be in the open interval (0, 1).")
    return _STD_NORM.inv_cdf(1.0 - p_threshold / 2.0)


def _z_expr(cols, z_col: str = "Z", beta_col: str = "BETA", se_col: str = "SE"):
    """
    Polars expression for the Wald z-score, or None if the columns are absent.
    Priority: an explicit Z column, otherwise BETA / SE. SE <= 0 (or non-finite)
    yields a null z so the variant is treated as undecidable rather than dividing
    by zero.
    """
    if z_col and z_col in cols:
        return pl.col(z_col).cast(pl.Float64)
    if beta_col in cols and se_col in cols:
        se = pl.col(se_col).cast(pl.Float64)
        return (
            pl.when(se > 0.0)
            .then(pl.col(beta_col).cast(pl.Float64) / se)
            .otherwise(None)
        )
    return None


def _significance_filter_expr(cols, keep: str, metric: str = "P", threshold: float = 5e-8,
                              z_col: str = "Z", beta_col: str = "BETA", se_col: str = "SE",
                              p_col: str = "P"):
    """
    Build a boolean Polars expression selecting the rows to KEEP.

        keep      : "significant"  -> keep associated variants
                    "nonsignificant" -> keep the complement
        metric    : "P" -> two-tailed p derived from z (z = Z col, else BETA/SE)
                    "Z" -> filter directly on |z| = |BETA/SE|
        threshold : p-value cutoff (metric "P") or |z| cutoff (metric "Z")

    Returns (expr, description). expr is None (with a reason in description) when
    the required columns are missing, so the caller can skip filtering gracefully.
    Variants with a null/undecidable z (or null P in the fallback) are excluded
    from the kept set under BOTH keep-directions.
    """
    keep = keep.lower()
    if keep not in ("significant", "nonsignificant"):
        raise ValueError("keep must be 'significant' or 'nonsignificant'.")
    metric = metric.upper()

    z = _z_expr(cols, z_col, beta_col, se_col)

    if metric == "Z":
        if z is None:
            return None, "Z filtering needs a 'Z' column or both 'BETA' and 'SE'."
        zc = float(threshold)
        sig = z.abs() > pl.lit(zc, dtype=pl.Float64)
        desc = f"|z| > {zc:g}"
    elif metric == "P":
        if z is not None:
            zc = z_critical_from_p(float(threshold))
            sig = z.abs() > pl.lit(zc, dtype=pl.Float64)   # p < alpha  <=>  |z| > z_c
            desc = f"two-tailed p < {threshold:g}  (|z| > {zc:.4f})"
        elif p_col in cols:
            sig = pl.col(p_col).cast(pl.Float64) < pl.lit(float(threshold), dtype=pl.Float64)
            desc = f"{p_col} < {threshold:g}"
        else:
            return None, "P filtering needs 'Z' or 'BETA'/'SE' (to derive p), or a 'P' column."
    else:
        raise ValueError("metric must be 'P' or 'Z'.")

    keep_expr = sig if keep == "significant" else ~sig
    # Null predicate -> row dropped by .filter(); make the intent explicit so an
    # undecidable variant is never silently retained as "non-significant".
    keep_expr = keep_expr & sig.is_not_null()
    return keep_expr, desc


def _semi_join_filter(lazy_df: pl.LazyFrame, col_name: str, series_to_match: pl.Series) -> pl.LazyFrame:
    """
    Performance Trick: Replaces massive `.is_in()` calls.
    Building an AST with millions of elements in .is_in() slows the Polars
    optimizer to a crawl. Doing a semi-join allows Rust to use high-speed hash tables.
    """
    match_df = series_to_match.to_frame("MATCH_COL").lazy()
    return lazy_df.join(match_df, left_on=col_name, right_on="MATCH_COL", how="semi")


def _merge_tmp_files(output_file: str, tmp_files: list[str]):
    """
    Merges temporary chromosome files sequentially using a 16MB block buffer.
    Ensures O(1) memory usage regardless of dataset size.
    """
    if not tmp_files:
        return False

    with open(output_file, 'wb') as outfile:
        for i, f in enumerate(tmp_files):
            with open(f, 'rb') as infile:
                if i != 0:
                    # Skip header line for all but the first file
                    infile.readline()
                # 16MB chunked copy to strictly bound memory
                shutil.copyfileobj(infile, outfile, length=1024 * 1024 * 16)
            os.remove(f)  # Clean up the temp file immediately
    return True


# -------------------------------------------------------------------
# SNP -> Gene annotation (--snp_to_gene YES)
# -------------------------------------------------------------------
# The gene reference lives in ``snps_genes_ref/`` in the working directory,
# split per chromosome as ``snp_gene_map_chr{N}.parquet`` and holding the
# columns ``Chromosome``, ``SNP`` and ``Gene``. Files are read through the
# Arrow library (pyarrow), projecting only ``SNP``/``Gene`` so the memory
# footprint stays small, and collapsed to a single Gene value per SNP
# (multiple overlapping genes are joined with ";"). The resulting map is
# left-joined onto every rsID/SNP identifier column of the final results
# just before it is streamed to disk, so no LD row is ever duplicated.

# Any of: rsid, snp, rsID1/rsid1/rsID_1, rsID2, rsID_A/rsID_B, SNP1/SNP2,
# SNP_A/SNP_B ... (case-insensitive). Position/allele columns such as POS1,
# REF1, MAF1, Uniq_ID_1 do NOT start with rsid/snp and are therefore ignored.
_ID_COL_RE = re.compile(r"^(rsid|snp)(?:[_\-]?([0-9]+|[abAB]))?$", re.IGNORECASE)


def _gene_col_for(id_col: str) -> str:
    """Name of the Gene column paired with a given identifier column.

    rsID1 / rsid1 / rsID_1 -> Gene1 ; rsID2 -> Gene2 ;
    rsID_A / SNP_A -> Gene_A ; bare SNP / rsid -> Gene.
    """
    m = _ID_COL_RE.match(id_col)
    suffix = m.group(2) if (m and m.group(2)) else ""
    if suffix == "":
        return "Gene"
    if suffix.isdigit():
        return f"Gene{suffix}"
    return f"Gene_{suffix.upper()}"


def build_snp_gene_map(chroms, snps_genes_dir: str = "snps_genes_ref"):
    """Build a unique SNP -> Gene lookup for the requested chromosomes.

    Reads ``snps_genes_dir/snp_gene_map_chr{N}.parquet`` with pyarrow (only the
    ``SNP`` and ``Gene`` columns), concatenates the needed chromosomes and
    collapses to one row per SNP, with overlapping genes joined by ";".
    Returns a small collected ``pl.DataFrame`` (columns ``SNP``, ``Gene``) that
    is reused for every output batch, or ``None`` if annotation cannot proceed
    (missing folder / files / pyarrow). Returning ``None`` makes annotation a
    graceful no-op rather than an error.
    """
    import glob

    if not os.path.isdir(snps_genes_dir):
        print(f"[snp_to_gene] Reference folder '{snps_genes_dir}' not found - skipping gene annotation.")
        return None

    # Prefer only the chromosome files actually needed by this run.
    files = []
    if chroms:
        for c in chroms:
            f = os.path.join(snps_genes_dir, f"snp_gene_map_chr{c}.parquet")
            if os.path.exists(f):
                files.append(f)
    if not files:  # fall back to every available map file
        files = sorted(glob.glob(os.path.join(snps_genes_dir, "snp_gene_map_chr*.parquet")))
    if not files:
        print(f"[snp_to_gene] No 'snp_gene_map_chr*.parquet' files in '{snps_genes_dir}' - skipping.")
        return None

    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError:
        print("[snp_to_gene] pyarrow (arrow library) not available - skipping gene annotation.")
        return None

    tables = []
    for f in files:
        try:
            # Column projection: read only what we need for a light footprint.
            tables.append(pq.read_table(f, columns=["SNP", "Gene"]))
        except Exception as e:
            print(f"[snp_to_gene] Could not read '{f}': {e}")
    if not tables:
        print("[snp_to_gene] No readable gene-map tables - skipping gene annotation.")
        return None

    combined = pa.concat_tables(tables)

    # Hand the Arrow table to Polars, then reduce to a unique SNP key so the
    # downstream left-join is strictly 1:1 and can never inflate LD row counts.
    gene_map = (
        pl.from_arrow(combined)
        .lazy()
        .select([pl.col("SNP").cast(pl.Utf8), pl.col("Gene").cast(pl.Utf8)])
        .drop_nulls(subset=["SNP"])
        .group_by("SNP")
        .agg(pl.col("Gene").drop_nulls().unique(maintain_order=True).alias("Gene"))
        .with_columns(pl.col("Gene").list.join(";"))
        .collect()
    )
    print(f"[snp_to_gene] Loaded {gene_map.height} unique SNP->Gene entries from {len(files)} file(s).")
    return gene_map


def annotate_genes_lazy(lazy_frame: pl.LazyFrame, gene_map) -> pl.LazyFrame:
    """Left-join a ``Gene`` column beside every rsID/SNP identifier column.

    Detects identifier columns via ``_ID_COL_RE`` on the (lazy) schema and, for
    each, joins the pre-built ``gene_map`` on that column, placing the resulting
    gene column immediately to the right of its source. Unknown SNPs get a null
    gene. The whole operation stays lazy so it fuses into the existing
    ``sink_csv`` streaming plan. If ``gene_map`` is ``None`` or no identifier
    column is present, the frame is returned unchanged.
    """
    if gene_map is None:
        return lazy_frame

    schema_names = lazy_frame.collect_schema().names()
    id_cols = [c for c in schema_names if _ID_COL_RE.match(c)]
    if not id_cols:
        return lazy_frame

    gene_map_lazy = gene_map.lazy() if isinstance(gene_map, pl.DataFrame) else gene_map

    used = set(schema_names)
    final_order = list(schema_names)
    lf = lazy_frame
    for id_col in id_cols:
        gcol = _gene_col_for(id_col)
        if gcol in used:                      # avoid clobbering an existing / duplicate name
            gcol = f"Gene_{id_col}"
        used.add(gcol)
        gm = gene_map_lazy.rename({"SNP": id_col, "Gene": gcol})
        # Cast the key to Utf8 on both sides so the join never fails on dtype.
        lf = lf.with_columns(pl.col(id_col).cast(pl.Utf8)).join(gm, on=id_col, how="left")
        idx = final_order.index(id_col)
        final_order.insert(idx + 1, gcol)     # gene column sits right after its SNP
    return lf.select(final_order)


# -------------------------------------------------------------------
# High-Performance Numba Graph Pruning (Zero-Allocation Streaming)
# -------------------------------------------------------------------

@njit(fastmath=True)
def _numba_greedy_prune(n_snps: int, edges_u: np.ndarray, edges_v: np.ndarray) -> np.ndarray:
    """
    O(E + V) JIT-compiled Greedy Graph Pruning.
    Zero-Allocation strategy: Assumes edges are directed (u < v) and sorted by u ascending.
    Streams through the edges directly without building an adjacency matrix.
    Returns a boolean mask of SNPs to keep.
    """
    keep = np.zeros(n_snps, dtype=np.bool_)
    excluded = np.zeros(n_snps, dtype=np.bool_)

    n_edges = len(edges_u)
    edge_idx = 0

    for i in range(n_snps):
        is_kept = not excluded[i]
        keep[i] = is_kept

        # Fast-forward edge_idx if it somehow fell behind
        while edge_idx < n_edges and edges_u[edge_idx] < i:
            edge_idx += 1

        # Process all outgoing edges from node i
        while edge_idx < n_edges and edges_u[edge_idx] == i:
            if is_kept:
                v = edges_v[edge_idx]
                excluded[v] = True
            edge_idx += 1

    return keep


def ld_prune_pairs(
        study_lazy: pl.LazyFrame,
        ld_pairs_lazy: pl.LazyFrame,
        snp_col: str = "SNP",
        p_col: str | None = "P",
        p_threshold: float | None = None,
        p_threshold_mode: str = "below",
        z_col: str | None = "Z",
        z_threshold: float | None = 2.0,
        prune_metric: str = "P",  # Can be "P" or "Z" to dictate which SNP to prioritize
        snp1_col: str = "SNP1",
        snp2_col: str = "SNP2"
) -> pl.LazyFrame:
    """
    Massively optimized Polars + Numba LD pruning.
    Supports filtering and sorting by both P-values and Z-values.
    """
    active_lazy = study_lazy
    study_cols = active_lazy.collect_schema().names()

    # 1) Pre-filter based on P-value threshold if provided
    if (p_threshold is not None) and (p_col is not None) and (p_col in study_cols):
        p_expr = pl.col(p_col).cast(pl.Float32)
        if p_threshold_mode == "above":
            active_lazy = active_lazy.filter(p_expr > pl.lit(p_threshold, dtype=pl.Float32))
        else:
            active_lazy = active_lazy.filter(p_expr < pl.lit(p_threshold, dtype=pl.Float32))

    # 2) Pre-filter based on Z-value threshold if provided
    if (z_threshold is not None) and (z_col is not None) and (z_col in study_cols):
        # We take the absolute value of Z, expecting it to be > threshold (e.g. > 2.0)
        z_expr = pl.col(z_col).cast(pl.Float32).abs()
        active_lazy = active_lazy.filter(z_expr > pl.lit(z_threshold, dtype=pl.Float32))

    # 3) Setup sorting expression so the "most significant" SNPs are processed first by Numba
    sort_expr = None
    sort_desc = False

    if prune_metric.upper() == "P" and (p_col is not None) and (p_col in study_cols):
        sort_expr = pl.col(p_col).cast(pl.Float32)
        sort_desc = (p_threshold_mode == "above")
    elif prune_metric.upper() == "Z" and (z_col is not None) and (z_col in study_cols):
        sort_expr = pl.col(z_col).cast(pl.Float32).abs()
        sort_desc = True  # We want largest Z-values first

    # 4) Extract ordered SNPs natively
    select_cols = [snp_col]
    if (p_col is not None) and (p_col in study_cols): select_cols.append(p_col)
    if (z_col is not None) and (z_col in study_cols): select_cols.append(z_col)

    ordered_snps_lazy = (
        active_lazy.select(select_cols)
        .drop_nulls(subset=[snp_col])
        .unique(subset=[snp_col], maintain_order=False)
    )

    # Use streaming=True to ensure out-of-core chunk processing for massive sorting
    if sort_expr is not None:
        ordered_snps_df = ordered_snps_lazy.sort(sort_expr, descending=sort_desc).collect(streaming=True)
    else:
        ordered_snps_df = ordered_snps_lazy.collect(streaming=True)

    # Add 0 to N-1 integer IDs mapping to feed Numba natively
    n_snps = len(ordered_snps_df)
    ordered_snps_df = ordered_snps_df.with_columns(
        pl.int_range(0, n_snps, dtype=pl.Int32).alias("_numba_id")
    )

    # 5) Memory-Safe Graph Construction:
    # Drops the massive String columns and R2 values at the Rust layer.
    ld_mapped = (
        ld_pairs_lazy
        .select([snp1_col, snp2_col])
        .join(ordered_snps_df.lazy(), left_on=snp1_col, right_on=snp_col, how="inner")
        .rename({"_numba_id": "u"})
        .join(ordered_snps_df.lazy(), left_on=snp2_col, right_on=snp_col, how="inner")
        .rename({"_numba_id": "v"})
        .with_columns([
            pl.min_horizontal("u", "v").alias("u_directed"),
            pl.max_horizontal("u", "v").alias("v_directed")
        ])
        .filter(pl.col("u_directed") != pl.col("v_directed"))  # Remove self-loops
        .select([
            pl.col("u_directed").alias("u"),
            pl.col("v_directed").alias("v")
        ])
        .unique()  # Removes duplicated undirected edges (A->B and B->A become just A->B)
        .sort("u")  # Numba requires edges to be sorted ascending by origin node
        .collect()  # Safe to collect here: integer arrays are >90% smaller than string pairs
    )

    if ld_mapped.is_empty():
        return active_lazy

    # 6) Run extreme speed graph algorithm in C via Numba
    edges_u = ld_mapped["u"].to_numpy()
    edges_v = ld_mapped["v"].to_numpy()

    keep_mask = _numba_greedy_prune(n_snps, edges_u, edges_v)

    # Aggressively return internal array memory to OS before continuing
    del ld_mapped, edges_u, edges_v

    # 7) Filter by the boolean array and return Semi-Join AST node
    kept_snps = ordered_snps_df.filter(pl.Series(keep_mask))[snp_col]
    return _semi_join_filter(active_lazy, snp_col, kept_snps)


# -------------------------------------------------------------------
# TOP-LD (pairwise)
# -------------------------------------------------------------------

def TOP_LD_info_pairwise(rs_series, chrom, population, maf_threshold, R2_threshold, imp_snp_list=None, ld_prune=False,
                         pairwise: bool = True):
    print(f"Building query plan: TOP-LD files for {population} chr{chrom}...")

    snp_series = _get_snp_series(rs_series, imp_snp_list)

    maf_path = f"D:/ref/TOP_LD/{population}/SNV/{population}_chr{chrom}_no_filter_0.2_1000000_info_annotation.parquet"
    ld_path = f"D:/ref/TOP_LD/{population}/SNV/{population}_chr{chrom}_no_filter_0.2_1000000_LD.parquet"

    if ld_prune:
        maf_lazy = pl.scan_parquet(maf_path).select(["Uniq_ID", "rsID"])
        maf_lazy1 = _semi_join_filter(maf_lazy, "rsID", rs_series).rename({"Uniq_ID": "Uniq_ID_1", "rsID": "SNP1"})
        maf_lazy2 = _semi_join_filter(maf_lazy, "rsID", snp_series).rename({"Uniq_ID": "Uniq_ID_2", "rsID": "SNP2"})

        ld_lazy = pl.scan_parquet(ld_path).select(["Uniq_ID_1", "Uniq_ID_2", "R2"])
        if R2_threshold is not None:
            ld_lazy = ld_lazy.filter(pl.col("R2").cast(pl.Float32) >= pl.lit(R2_threshold, dtype=pl.Float32))

        return ld_lazy.join(maf_lazy1, on="Uniq_ID_1", how="inner").join(maf_lazy2, on="Uniq_ID_2", how="inner").select(
            ["SNP1", "SNP2", "R2"])

    maf_lazy = (
        pl.scan_parquet(maf_path)
        .select(["Uniq_ID", "rsID", "MAF", "REF", "ALT"])
        .filter(pl.col("MAF").cast(pl.Float32) >= pl.lit(maf_threshold, dtype=pl.Float32))
    )

    # Anchor (SNP1) is always the inserted rsIDs. The target (SNP2) is restricted
    # to the inserted set only in pairwise mode; in global mode it spans the whole
    # (MAF-filtered) panel so every LD partner of each inserted SNP is returned.
    if pairwise:
        maf_target = _semi_join_filter(maf_lazy, "rsID", snp_series)
    else:
        maf_target = maf_lazy
    valid_ids1 = _semi_join_filter(maf_lazy, "rsID", rs_series).select("Uniq_ID")
    valid_ids2 = maf_target.select("Uniq_ID")

    ld_lazy = pl.scan_parquet(ld_path)
    if R2_threshold is not None:
        ld_lazy = ld_lazy.filter(pl.col("R2").cast(pl.Float32) >= pl.lit(R2_threshold, dtype=pl.Float32))

    ld_lazy = (
        ld_lazy
        .join(valid_ids1, left_on="Uniq_ID_1", right_on="Uniq_ID", how="semi")
        .join(valid_ids2, left_on="Uniq_ID_2", right_on="Uniq_ID", how="semi")
    )

    maf1 = maf_lazy.rename({"Uniq_ID": "Uniq_ID_1", "rsID": "rsID1", "MAF": "MAF1", "REF": "REF1", "ALT": "ALT1"})
    maf2 = maf_lazy.rename({"Uniq_ID": "Uniq_ID_2", "rsID": "rsID2", "MAF": "MAF2", "REF": "REF2", "ALT": "ALT2"})

    return ld_lazy.join(maf1, on="Uniq_ID_1", how="inner").join(maf2, on="Uniq_ID_2", how="inner").select([
        pl.col("Uniq_ID_1").alias("pos1"), pl.col("Uniq_ID_2").alias("pos2"),
        "R2", "+/-corr", "Dprime", "rsID1", "rsID2", "MAF1", "MAF2", "REF1", "ALT1", "REF2", "ALT2"
    ])


def TOP_LD_process_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list, ld_prune=False):
    return TOP_LD_info_pairwise(rs_series, chromosome, population, maf_input, r2threshold, imp_snp_list, ld_prune)


# -------------------------------------------------------------------
# HapMap (pairwise)
# -------------------------------------------------------------------

def Hap_Map_LD_info_dask_pairwise(rs_series, chrom, population, maf_threshold, R2_threshold, imp_snp_list,
                                  ld_prune=False, pairwise: bool = True):
    print(f"Building query plan: Hap Map files ({population}) chr{chrom}...")
    snp_series = _get_snp_series(rs_series, imp_snp_list)

    maf_file = f'D:/ref/Hap_Map/allele_freqs_chr{chrom}_{population}_phase3.2_nr.b36_fwd.parquet'
    ld_file = f'D:/ref/Hap_Map/ld_chr{chrom}_{population}.parquet'

    ld_raw = pl.scan_parquet(ld_file)
    cols = ld_raw.collect_schema().names()

    if ld_prune:
        ld_lazy = ld_raw.select([cols[3], cols[4], cols[6]])
        ld_lazy = _semi_join_filter(ld_lazy, cols[3], rs_series)
        ld_lazy = _semi_join_filter(ld_lazy, cols[4], snp_series)

        ld_lazy = ld_lazy.rename({cols[3]: "SNP1", cols[4]: "SNP2", cols[6]: "R2"})
        if R2_threshold is not None:
            ld_lazy = ld_lazy.filter(pl.col("R2").cast(pl.Float32) >= pl.lit(R2_threshold, dtype=pl.Float32))
        return ld_lazy

    maf_df = (
        pl.scan_parquet(maf_file)
        .select(["rs#", "refallele", "otherallele", "otherallele_freq"])
        .rename({"rs#": "rsID", "refallele": "REF", "otherallele": "ALT", "otherallele_freq": "MAF"})
        .with_columns(pl.col("MAF").cast(pl.Float32))
        .filter(pl.col("MAF") >= pl.lit(maf_threshold, dtype=pl.Float32))
    )
    # In pairwise mode restrict the MAF table to the inserted set; in global mode
    # keep the whole (MAF-filtered) panel so any SNP2 partner can be annotated.
    if pairwise:
        maf_df = _semi_join_filter(maf_df, "rsID", snp_series)

    ld_df = ld_raw.filter(pl.col(cols[3]) != pl.col(cols[4]))
    ld_df = _semi_join_filter(ld_df, cols[3], rs_series)          # SNP1 anchor: inserted rsIDs
    if pairwise:
        ld_df = _semi_join_filter(ld_df, cols[4], snp_series)     # SNP2: inserted set only (pairwise)

    if R2_threshold is not None:
        ld_df = ld_df.filter(pl.col(cols[6]).cast(pl.Float32) >= pl.lit(R2_threshold, dtype=pl.Float32))

    ld_df = ld_df.rename({
        cols[0]: "pos1", cols[1]: "pos2", cols[2]: "pop",
        cols[3]: "rsID1", cols[4]: "rsID2", cols[5]: "Dprime", cols[6]: "R2"
    }).select(["pos1", "pos2", "pop", "rsID1", "rsID2", "Dprime", "R2"])

    return (
        ld_df.join(maf_df, left_on="rsID1", right_on="rsID", how="inner")
        .join(maf_df, left_on="rsID2", right_on="rsID", how="inner", suffix="_2")
        .select([
            pl.col("pos1"), pl.col("pos2"), pl.col("rsID1"), pl.col("rsID2"),
            pl.col("MAF").alias("MAF1"), pl.col("MAF_2").alias("MAF2"),
            pl.col("REF").alias("REF1"), pl.col("REF_2").alias("REF2"),
            pl.col("ALT").alias("ALT1"), pl.col("ALT_2").alias("ALT2"),
            pl.col("R2"), pl.col("Dprime")
        ])
    )


def Hap_Map_process_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list, ld_prune=False):
    return Hap_Map_LD_info_dask_pairwise(rs_series, chromosome, population, maf_input, r2threshold, imp_snp_list,
                                         ld_prune)


# -------------------------------------------------------------------
# PhenoScanner (pairwise)
# -------------------------------------------------------------------

def pheno_Scanner_LD_info_dask_pairwise(rs_series, chrom, population, maf_threshold, R2_threshold, imp_snp_list,
                                        ld_prune=False, pairwise: bool = True):
    if R2_threshold is not None and R2_threshold < 0.8:
        print("Pheno Scanner includes data with R2 >= 0.8. Setting R2_threshold = 0.8")
        R2_threshold = 0.8

    print(f"Building query plan: Pheno Scanner files chr{chrom}...")
    snp_series = _get_snp_series(rs_series, imp_snp_list)

    maf_file = "D:/ref/Pheno_Scanner/1000G.parquet"
    ld_file = f"D:/ref/Pheno_Scanner/1000G_{population}/1000G_{population}_chr{chrom}.parquet"
    pop_map = {"EUR": "eur", "EAS": "eas", "AFR": "afr", "AMR": "amr", "SAS": "sas"}
    maf_pop = pop_map.get(population)

    if ld_prune:
        ld_lazy = pl.scan_parquet(ld_file).select(["ref_rsid", "rsid", "r2"])
        ld_lazy = _semi_join_filter(ld_lazy, "ref_rsid", rs_series)
        ld_lazy = _semi_join_filter(ld_lazy, "rsid", snp_series)

        ld_lazy = ld_lazy.rename({"ref_rsid": "SNP1", "rsid": "SNP2", "r2": "R2"})
        if R2_threshold is not None:
            ld_lazy = ld_lazy.filter(pl.col("R2").cast(pl.Float32) >= pl.lit(R2_threshold, dtype=pl.Float32))
        return ld_lazy

    maf1_lazy = (
        pl.scan_parquet(maf_file)
        .select(["hg19_coordinates", "chr", "rsid", maf_pop, "a1", "a2"])
        .filter((pl.col("chr") == chrom) & (pl.col(maf_pop) != "-"))
        .with_columns(pl.col(maf_pop).cast(pl.Float32))
        .filter(pl.col(maf_pop) >= pl.lit(maf_threshold, dtype=pl.Float32))
    )

    # SNP2 annotation: restrict to inserted set in pairwise mode, whole panel in global.
    maf2_base = _semi_join_filter(maf1_lazy, "rsid", snp_series) if pairwise else maf1_lazy
    maf2_lazy = maf2_base.rename({
        "hg19_coordinates": "pos2(hg19)", maf_pop: "MAF2", "a1": "ALT2", "a2": "REF2"
    }).select(["pos2(hg19)", "rsid", "MAF2", "ALT2", "REF2"])

    maf1_lazy = _semi_join_filter(maf1_lazy, "rsid", rs_series).rename({
        "hg19_coordinates": "ref_hg19_coordinates", "rsid": "ref_rsid", maf_pop: "MAF1", "a1": "ALT1", "a2": "REF1"
    }).select(["ref_hg19_coordinates", "ref_rsid", "MAF1", "ALT1", "REF1"])

    ld_lazy = pl.scan_parquet(ld_file).select(["ref_hg19_coordinates", "ref_rsid", "rsid", "r2", "r", "dprime"]).filter(
        pl.col("ref_rsid") != pl.col("rsid"))
    ld_lazy = _semi_join_filter(ld_lazy, "ref_rsid", rs_series)   # SNP1 anchor
    if pairwise:
        ld_lazy = _semi_join_filter(ld_lazy, "rsid", snp_series)  # SNP2: pairwise only

    if R2_threshold is not None:
        ld_lazy = ld_lazy.filter(pl.col("r2").cast(pl.Float32) >= pl.lit(R2_threshold, dtype=pl.Float32))

    return (
        ld_lazy.join(maf1_lazy, on="ref_rsid", how="inner").join(maf2_lazy, on="rsid", how="inner")
        .rename({"ref_rsid": "rsID1", "rsid": "rsID2", "ref_hg19_coordinates": "pos1(hg19)", "r2": "R2"})
        .with_columns([
            pl.col("pos1(hg19)").str.split(":").list.last().cast(pl.Int32).alias("pos1(hg19)"),
            pl.col("pos2(hg19)").str.split(":").list.last().cast(pl.Int32).alias("pos2(hg19)"),
            pl.lit(chrom).cast(pl.Int16).alias("CHR")
        ])
        .select(
            ["rsID1", "pos1(hg19)", "rsID2", "dprime", "pos2(hg19)", "R2", "r", "MAF1", "MAF2", "ALT1", "REF1", "ALT2",
             "REF2", "CHR"])
    )


def pheno_Scanner_process_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list,
                                   ld_prune=False):
    return pheno_Scanner_LD_info_dask_pairwise(rs_series, chromosome, population, maf_input, r2threshold, imp_snp_list,
                                               ld_prune)


# -------------------------------------------------------------------
# 1000G hg38 (pairwise)
# -------------------------------------------------------------------

def hg38_1kg_LD_info_pairwise(rs_series, chrom, population, maf_threshold, R2_threshold, imp_snp_list, ld_prune=False,
                              pairwise: bool = True):
    print(f"Building query plan: 1000 Genomes Project (hg38) files ({population}) chr{chrom}...")
    snp_series = _get_snp_series(rs_series, imp_snp_list)

    maf_file = f'D:/ref/1000G_hg38/1000G_{population}_0_01.parquet'
    ld_file = f'D:/ref/1000G_hg38/{population}/chr{chrom}_merged.parquet'

    ld_lazy = pl.scan_parquet(ld_file)

    if ld_prune:
        ld_lazy = ld_lazy.select(["SNP_A", "SNP_B", "R"])
        ld_lazy = _semi_join_filter(ld_lazy, "SNP_A", rs_series)
        ld_lazy = _semi_join_filter(ld_lazy, "SNP_B", snp_series)

        ld_lazy = ld_lazy.with_columns(pl.col("R").cast(pl.Float32).pow(2).alias("R2")).rename(
            {"SNP_A": "SNP1", "SNP_B": "SNP2"}).drop("R")
        if R2_threshold is not None:
            ld_lazy = ld_lazy.filter(pl.col("R2") >= pl.lit(R2_threshold, dtype=pl.Float32))
        return ld_lazy

    tmp = pl.scan_parquet(maf_file)
    cols = tmp.collect_schema().names()
    maf_lazy = tmp.rename({old: new for old, new in zip(cols[:6], ['CHR', 'SNP', 'MAF', 'POS', 'REF', 'ALT'])})
    maf_lazy = maf_lazy.with_columns(pl.col("MAF").cast(pl.Float32))
    # Restrict the MAF table to the inserted set only in pairwise mode; global mode
    # keeps the whole (MAF-filtered) panel so any SNP_B partner can be annotated.
    if pairwise:
        maf_lazy = _semi_join_filter(maf_lazy, "SNP", snp_series)

    if maf_threshold is not None:
        maf_lazy = maf_lazy.filter(pl.col("MAF") >= pl.lit(maf_threshold, dtype=pl.Float32))

    ld_lazy = ld_lazy.select(["CHR_A", "SNP_A", "CHR_B", "SNP_B", "R"])
    ld_lazy = _semi_join_filter(ld_lazy, "SNP_A", rs_series)          # SNP_A anchor
    if pairwise:
        ld_lazy = _semi_join_filter(ld_lazy, "SNP_B", snp_series)     # SNP_B: pairwise only
    ld_lazy = ld_lazy.with_columns(pl.col("R").cast(pl.Float32).pow(2).alias("R2"))

    if R2_threshold is not None:
        ld_lazy = ld_lazy.filter(pl.col("R2") >= pl.lit(R2_threshold, dtype=pl.Float32))

    return (
        ld_lazy.join(maf_lazy, left_on=["CHR_A", "SNP_A"], right_on=["CHR", "SNP"], how="inner")
        .rename({"CHR_A": "CHR", "MAF": "MAF_A", "POS": "POS_A", "REF": "REF_A", "ALT": "ALT_A"})
        .join(maf_lazy, left_on=["CHR_B", "SNP_B"], right_on=["CHR", "SNP"], how="inner")
        .rename({"MAF": "MAF_B", "POS": "POS_B", "REF": "REF_B", "ALT": "ALT_B"})
        .drop(["CHR_B"])
    )


def hg38_1kg_process_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list, ld_prune=False):
    return hg38_1kg_LD_info_pairwise(rs_series, chromosome, population, maf_input, r2threshold, imp_snp_list, ld_prune)


def _generic_vcor_pairwise(rs_series, R2_threshold, maf_threshold, ld_file, imp_snp_list, ld_prune, is_phased=False,
                           pairwise: bool = True):
    """Shared fast logic for HGDP and 1KGP high coverage using new fast semi joins."""
    snp_series = _get_snp_series(rs_series, imp_snp_list)
    ld_lazy = pl.scan_parquet(ld_file)
    r_col = "PHASED_R" if is_phased else "UNPHASED_R"

    if ld_prune:
        ld_lazy = ld_lazy.select(["ID_A", "ID_B", r_col])
        ld_lazy = _semi_join_filter(ld_lazy, "ID_A", rs_series)
        ld_lazy = _semi_join_filter(ld_lazy, "ID_B", snp_series)

        ld_lazy = ld_lazy.with_columns(pl.col(r_col).cast(pl.Float32).pow(2).alias("R2")).rename(
            {"ID_A": "SNP1", "ID_B": "SNP2"}).drop(r_col)
        if R2_threshold is not None:
            ld_lazy = ld_lazy.filter(pl.col("R2") >= pl.lit(R2_threshold, dtype=pl.Float32))
        return ld_lazy

    cols = ld_lazy.collect_schema().names()
    ld_lazy = _semi_join_filter(ld_lazy, "ID_A", rs_series)          # ID_A anchor
    if pairwise:
        ld_lazy = _semi_join_filter(ld_lazy, "ID_B", snp_series)     # ID_B: pairwise only

    maf1_col = "MAF_1" if "MAF_1" in cols else ("MAF1" if "MAF1" in cols else "NONMAJ_FREQ_A")
    ld_lazy = ld_lazy.filter(pl.col(maf1_col).cast(pl.Float32) >= pl.lit(maf_threshold, dtype=pl.Float32))

    rename_map = {"#CHROM_A": "CHR", "CHROM_B": "CHR2", 'POS_A': 'POS1', 'POS_B': 'POS2', 'MAJ_A': 'MAJ1',
                  'MAJ_B': 'MAJ2', 'NONMAJ_A': 'NONMAJ1', 'NONMAJ_B': 'NONMAJ2', "ID_A": "rsID1", "ID_B": "rsID2",
                  "NONMAJ_FREQ_A": "MAF1", "NONMAJ_FREQ_B": "MAF2", "DPRIME": "Dprime", r_col: "R"}
    rename_map = {k: v for k, v in rename_map.items() if k in cols}

    ld_lazy = ld_lazy.rename(rename_map).with_columns(pl.col("R").cast(pl.Float32).pow(2).alias("R2"))

    if R2_threshold is not None:
        ld_lazy = ld_lazy.filter(pl.col("R2") >= pl.lit(R2_threshold, dtype=pl.Float32))

    if "CHR2" in ld_lazy.collect_schema().names():
        ld_lazy = ld_lazy.drop("CHR2")

    return ld_lazy


def hg38_1kg_LD_info_high_coverage_pairwise(rs_series, R2_threshold, population, maf_threshold, chrom, imp_snp_list,
                                            ld_prune=False, pairwise: bool = True):
    print(f"Building query plan: 1KGP high-cov files ({population}) chr{chrom}...")
    ld_file = f'D:/ref/1KGP_high_coverage/LD_{population}_r_unphased/{population}_chr{chrom}_r_unphased.vcor.parquet'
    return _generic_vcor_pairwise(rs_series, R2_threshold, maf_threshold, ld_file, imp_snp_list, ld_prune,
                                  is_phased=False, pairwise=pairwise)


def hg38_1kg_process_high_coverage_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list,
                                            ld_prune=False):
    return hg38_1kg_LD_info_high_coverage_pairwise(rs_series, r2threshold, population, maf_input, chromosome,
                                                   imp_snp_list, ld_prune)


def HGDP_LD_info_pairwise(rs_series, R2_threshold, population, maf_threshold, chrom, imp_snp_list, ld_prune=False,
                          pairwise: bool = True):
    pop_map = {
        "EUR": "EUROPE",
        "EAS": "EAST_ASIA",
        "AFR": "AFRICA",
        "AMR": "AMERICA",
        "SAS": "CENTRAL_SOUTH_ASIA",
        "MID": "MIDDLE_EAST",
        "OCE": "OCEANIA"
    }

    # Safely convert to uppercase if input is abbreviation, otherwise fall back to string as-is
    hgdp_pop = pop_map.get(str(population).upper(), population)

    print(f"Building query plan: HGDP files ({hgdp_pop}) chr{chrom}...")
    ld_file = f'D:/ref/HGDP/LD_{hgdp_pop}_r_phased/{hgdp_pop}_chr{chrom}_r_phased.vcor.parquet'
    return _generic_vcor_pairwise(rs_series, R2_threshold, maf_threshold, ld_file, imp_snp_list, ld_prune,
                                  is_phased=True, pairwise=pairwise)


def HGDP_process_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list, ld_prune=False):
    return HGDP_LD_info_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list, ld_prune)


def UKBB_LD_info_pairwise(rs_series, R2_threshold, population, maf_threshold, chrom, imp_snp_list, ld_prune=False,
                          pairwise: bool = True):
    print(f"Building query plan: UK Biobank ({population}) chr{chrom}...")
    snp_series = _get_snp_series(rs_series, imp_snp_list)
    ld_file = f'D:/ref/UKBB/{population}/chr_{chrom}_ld.parquet'
    ld_lazy = pl.scan_parquet(ld_file)

    if ld_prune:
        ld_lazy = ld_lazy.select(["snp1", "snp2", "r"])
        ld_lazy = _semi_join_filter(ld_lazy, "snp1", rs_series)
        ld_lazy = _semi_join_filter(ld_lazy, "snp2", snp_series)

        ld_lazy = ld_lazy.with_columns(pl.col("r").cast(pl.Float32).pow(2).alias("R2")).rename(
            {"snp1": "SNP1", "snp2": "SNP2"}).drop("r")
        if R2_threshold is not None:
            ld_lazy = ld_lazy.filter(pl.col("R2") >= pl.lit(R2_threshold, dtype=pl.Float32))
        return ld_lazy

    cols = ld_lazy.collect_schema().names()
    ld_lazy = _semi_join_filter(ld_lazy, "snp1", rs_series)          # snp1 anchor
    if pairwise:
        ld_lazy = _semi_join_filter(ld_lazy, "snp2", snp_series)     # snp2: pairwise only

    maf1_col = "maf1" if "maf1" in cols else ("MAF1" if "MAF1" in cols else "maf1")
    ld_lazy = ld_lazy.filter(pl.col(maf1_col).cast(pl.Float32) >= pl.lit(maf_threshold, dtype=pl.Float32))

    rename_map = {"chr": "CHR", 'bp1': 'POS1', 'bp2': 'POS2', 'a1_1': 'A1_1', 'a1_2': 'A1_2', 'a2_1': 'A2_1',
                  'a2_2': 'A2_2', "snp1": "rsID1", "snp2": "rsID2", "maf1": "MAF1", "maf2": "MAF2", "r2": "R2",
                  "r": "R"}
    rename_map = {k: v for k, v in rename_map.items() if k in cols}

    ld_lazy = ld_lazy.rename(rename_map).with_columns(pl.col("R").cast(pl.Float32).pow(2).alias("R2"))

    if R2_threshold is not None:
        ld_lazy = ld_lazy.filter(pl.col("R2") >= pl.lit(R2_threshold, dtype=pl.Float32))

    return ld_lazy


def UKBB_process_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list, ld_prune=False):
    return UKBB_LD_info_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list, ld_prune)

# -------------------------------------------------------------------
# LASI-DAD (pairwise)
# -------------------------------------------------------------------

def LASI_DAD_LD_info_pairwise(rs_series, R2_threshold, population, maf_threshold, chrom, imp_snp_list, ld_prune=False,
                              pairwise: bool = True):
    print(f"Building query plan: LASI-DAD files chr{chrom}...")
    snp_series = _get_snp_series(rs_series, imp_snp_list)

    # LASI-DAD is a single-population (South Asian) panel, so `population` is
    # accepted for interface compatibility but not used to build the path.
    # >>> Adjust this path/pattern to match your actual layout. <
    ld_file = f'D:/ref/LASI_DAD/chr{chrom}_ld_with_rsids.parquet'

    ld_lazy = pl.scan_parquet(ld_file)

    # ---- LD-pruning mode: only rsIDs + R2 are needed ----
    if ld_prune:
        ld_lazy = ld_lazy.select(["rsid1", "rsid2", "R2"])
        ld_lazy = _semi_join_filter(ld_lazy, "rsid1", rs_series)
        ld_lazy = _semi_join_filter(ld_lazy, "rsid2", snp_series)

        ld_lazy = ld_lazy.rename({"rsid1": "SNP1", "rsid2": "SNP2"})
        # R2 is stored directly in the panel (no need to square r).
        if R2_threshold is not None:
            ld_lazy = ld_lazy.filter(pl.col("R2").cast(pl.Float32) >= pl.lit(R2_threshold, dtype=pl.Float32))
        return ld_lazy

    # ---- Full retrieval mode ----
    ld_lazy = _semi_join_filter(ld_lazy, "rsid1", rs_series)          # rsid1 anchor
    if pairwise:
        ld_lazy = _semi_join_filter(ld_lazy, "rsid2", snp_series)     # rsid2: pairwise only

    # AF1/AF2 are ALT-allele frequencies (can be > 0.5), so derive the true
    # minor-allele frequency as min(AF, 1 - AF) before applying the MAF filter.
    ld_lazy = ld_lazy.with_columns([
        pl.min_horizontal(
            pl.col("AF1").cast(pl.Float32), 1.0 - pl.col("AF1").cast(pl.Float32)
        ).alias("MAF1"),
        pl.min_horizontal(
            pl.col("AF2").cast(pl.Float32), 1.0 - pl.col("AF2").cast(pl.Float32)
        ).alias("MAF2"),
    ])

    if maf_threshold is not None:
        ld_lazy = ld_lazy.filter(pl.col("MAF1") >= pl.lit(maf_threshold, dtype=pl.Float32))

    if R2_threshold is not None:
        ld_lazy = ld_lazy.filter(pl.col("R2").cast(pl.Float32) >= pl.lit(R2_threshold, dtype=pl.Float32))

    # Reconstruct signed r from the stored Sign column: r = sign * sqrt(R2)
    ld_lazy = ld_lazy.with_columns(
        (pl.col("Sign").cast(pl.Float32) * pl.col("R2").cast(pl.Float32).sqrt()).alias("r")
    )

    return (
        ld_lazy
        .rename({
            "CHROM": "CHR",
            "POS1": "pos1", "POS2": "pos2",
            "rsid1": "rsID1", "rsid2": "rsID2",
        })
        .select([
            "CHR", "pos1", "pos2", "rsID1", "rsID2",
            "MAF1", "MAF2", "REF1", "REF2", "ALT1", "ALT2",
            "R2", "r", "Dprime",
        ])
    )


def LASI_DAD_process_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list, ld_prune=False):
    return LASI_DAD_LD_info_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list, ld_prune)
# -------------------------------------------------------------------
# Non-Pairwise Functions (Wrappers for backward compatibility)
# -------------------------------------------------------------------
def hg38_1kg_LD_info_high_coverage(rs_series, chrom, pop, maf, R2, imp): return hg38_1kg_LD_info_high_coverage_pairwise(
    rs_series, R2, pop, maf, chrom, imp, ld_prune=False, pairwise=False)


def hg38_1kg_process_high_coverage(rs_series, r2, pop, maf, chrom, imp): return hg38_1kg_LD_info_high_coverage(
    rs_series, chrom, pop, maf, r2, imp)


#LASI-DAD
def LASI_DAD_LD_info(rs_series, chrom, pop, maf, R2, imp): return LASI_DAD_LD_info_pairwise(rs_series, R2, pop, maf, chrom, imp, ld_prune=False, pairwise=False)

def LASI_DAD_process(rs_series, r2, pop, maf, chrom, imp): return LASI_DAD_LD_info(rs_series, chrom, pop, maf, r2, imp)



def HGDP_LD_info(rs_series, chrom, pop, maf, R2, imp): return HGDP_LD_info_pairwise(rs_series, R2, pop, maf, chrom, imp, ld_prune=False, pairwise=False)


def HGDP_process(rs_series, r2, pop, maf, chrom, imp): return HGDP_LD_info(rs_series, chrom, pop, maf, r2, imp)


def UKBB_LD_info(rs_series, chrom, pop, maf, R2, imp): return UKBB_LD_info_pairwise(rs_series, R2, pop, maf, chrom, imp, ld_prune=False, pairwise=False)


def UKBB_process(rs_series, r2, pop, maf, chrom, imp): return UKBB_LD_info(rs_series, chrom, pop, maf, r2, imp)


def hg38_1kg_LD_info(rs_series, chrom, pop, maf, R2, imp): return hg38_1kg_LD_info_pairwise(rs_series, chrom, pop, maf,
                                                                                            R2, imp, ld_prune=False,
                                                                                            pairwise=False)


def hg38_1kg_process(rs_series, r2, pop, maf, chrom, imp): return hg38_1kg_LD_info(rs_series, chrom, pop, maf, r2, imp)


def Hap_Map_LD_info_dask(rs_series, chrom, pop, maf, R2, imp): return Hap_Map_LD_info_dask_pairwise(rs_series, chrom,
                                                                                                    pop, maf, R2, imp,
                                                                                                    ld_prune=False,
                                                                                                    pairwise=False)


def Hap_Map_process(rs_series, r2, pop, maf, chrom, imp): return Hap_Map_LD_info_dask(rs_series, chrom, pop, maf, r2,
                                                                                      imp)


def pheno_Scanner_LD_info_dask(rs_series, chrom, pop, maf, R2, imp): return pheno_Scanner_LD_info_dask_pairwise(
    rs_series, chrom, pop, maf, R2, imp, ld_prune=False, pairwise=False)


def pheno_Scanner_process(rs_series, r2, pop, maf, chrom, imp): return pheno_Scanner_LD_info_dask(rs_series, chrom, pop,
                                                                                                  maf, r2, imp)


def TOP_LD_info(rs_series, chrom, pop, maf, R2, imp=None): return TOP_LD_info_pairwise(rs_series, chrom, pop, maf, R2,
                                                                                       imp, ld_prune=False,
                                                                                       pairwise=False)


def TOP_LD_process(rs_series, r2, pop, maf, chrom, imp): return TOP_LD_info(rs_series, chrom, pop, maf, r2, imp)


# -------------------------------------------------------------------
# Data Pipeline execution (Massively Parallel Streaming)
# -------------------------------------------------------------------

def _normalize_chr_values(values):
    out = set()
    for v in values:
        if v is None: continue
        m = re.search(r"(\d+)", str(v).strip())
        if m and 1 <= int(m.group(1)) <= 22:
            out.add(int(m.group(1)))
    return sorted(out)


def _normalize_autosomal_chroms(values) -> list[int]:
    out = set()
    for v in values:
        if v is None: continue
        s = re.sub(r"^(chr|CHR|Chr)", "", str(v).strip()).split(".")[0]
        if s.isdigit() and 1 <= int(s) <= 22:
            out.add(int(s))
    return sorted(out)


def _get_study_info_chunked(file_path: str, is_pairwise: bool = False,
                            sig_keep: str | None = None, sig_metric: str = "P",
                            sig_threshold: float = 5e-8,
                            z_col: str = "Z", beta_col: str = "BETA", se_col: str = "SE",
                            p_col: str = "P"):
    """Loads the GWAS input file in chunks to extract SNPs and Chromosomes without memory spikes.

    When ``sig_keep`` is "significant" or "nonsignificant", each chunk is filtered
    (Wald z = BETA/SE, two-tailed p from z) so only the requested variants survive
    to become LD anchors. Filtering here propagates to every panel and to both the
    global and pairwise retrieval modes, and stays fully chunked / out-of-core.
    """
    print("Scanning input file in chunks to preserve memory...")

    # Using read_csv_batched to iterate over massive files natively chunk-by-chunk
    reader = pl.read_csv_batched(file_path, separator="\t")

    chroms_set = set()
    snps_chunks = []
    study_cols = None
    filt_expr = None  # significance keep-expression, resolved once we see the header

    while True:
        batches = reader.next_batches(5)  # Process 5 chunks into memory at a time
        if not batches:
            break
        for chunk in batches:
            if study_cols is None:
                study_cols = chunk.columns
                if sig_keep is not None:
                    filt_expr, filt_desc = _significance_filter_expr(
                        study_cols, sig_keep, sig_metric, sig_threshold,
                        z_col, beta_col, se_col, p_col)
                    if filt_expr is None:
                        print(f"[significance filter] SKIPPED — {filt_desc}")
                    else:
                        print(f"[significance filter] keeping {sig_keep.upper()} variants: {filt_desc}")

            work = chunk.filter(filt_expr) if filt_expr is not None else chunk

            if "CHR" in study_cols:
                chroms_set.update(work.get_column("CHR").drop_nulls().unique().to_list())
            if "SNP" in study_cols:
                snps_chunks.append(work.get_column("SNP").drop_nulls().unique())

    if snps_chunks:
        rs_series = pl.concat(snps_chunks).unique(maintain_order=False)
    else:
        rs_series = pl.Series("SNP", [], dtype=pl.Utf8)

    if study_cols and "CHR" in study_cols:
        if is_pairwise:
            chroms = _normalize_autosomal_chroms(chroms_set)
        else:
            chroms = _normalize_chr_values(chroms_set)
    else:
        chroms = list(range(1, 23))

    if not chroms:
        chroms = list(range(1, 23))

    return chroms, rs_series, study_cols


def process_data(file_path, r2threshold, population, maf_input, ref_file, imp_snp_list,
                 sig_keep: str | None = None, sig_metric: str = "P", sig_threshold: float = 5e-8,
                 z_col: str = "Z", beta_col: str = "BETA", se_col: str = "SE", p_col: str = "P",
                 snp_to_gene: bool = False, snps_genes_dir: str = "snps_genes_ref"):
    # Determine execution strategy based on file size threshold (500MB)
    file_size_mb = os.path.getsize(file_path) / (1024 * 1024)
    is_large_file = file_size_mb > 500

    with pl.StringCache():
        chroms, rs_series, study_cols = _get_study_info_chunked(
            file_path, is_pairwise=False,
            sig_keep=sig_keep, sig_metric=sig_metric, sig_threshold=sig_threshold,
            z_col=z_col, beta_col=beta_col, se_col=se_col, p_col=p_col)

        # Build the SNP->Gene lookup once (reused for every output batch).
        gene_map = build_snp_gene_map(chroms, snps_genes_dir) if snp_to_gene else None

        processors = {
            "1000G_hg38": hg38_1kg_process, "1000G_hg38_high_cov": hg38_1kg_process_high_coverage,
            "UKBB": UKBB_process, "HGDP": HGDP_process, "TOP_LD": TOP_LD_process,
            "Pheno_Scanner": pheno_Scanner_process, "Hap_Map": Hap_Map_process,
            "LASI_DAD": LASI_DAD_process,
        }

        if ref_file not in processors:
            raise ValueError(f"Unsupported ref_panel: {ref_file}")

        batch_size = 4 if is_large_file else max(1, len(chroms))
        tmp_files = []

        for i in range(0, len(chroms), batch_size):
            batch_chroms = chroms[i:i + batch_size]
            lazy_results_list = []

            for chrom in batch_chroms:
                lazy_data = processors[ref_file](rs_series, r2threshold, population, maf_input, chrom, imp_snp_list)
                if lazy_data is not None:
                    lazy_results_list.append(lazy_data)
                    print(f"Evaluation graph for chr{chrom} constructed.")

            if lazy_results_list:
                batch_lazy = pl.concat(lazy_results_list, how="vertical_relaxed")
                # Attach Gene column(s) beside each rsID/SNP column before writing.
                out_lazy = annotate_genes_lazy(batch_lazy, gene_map)
                if is_large_file and batch_size < len(chroms):
                    tmp_file = f"LD_info_tmp_batch_{i}.txt"
                    print(f"Streaming batch {i // batch_size + 1} (chroms: {batch_chroms}) to disk...")
                    out_lazy.sink_csv(tmp_file, separator="\t")
                    tmp_files.append(tmp_file)
                else:
                    print("Streaming highly parallelized query across all chromosomes...")
                    out_lazy.sink_csv("LD_info_chr_all.txt", separator="\t")
                    print("Check 'LD_info_chr_all.txt' for results")
                    return

            # Aggressively collect garbage to free unreferenced internal memory arrays
            gc.collect()

        if is_large_file and tmp_files:
            print("Merging temporary batch files into final output...")
            _merge_tmp_files("LD_info_chr_all.txt", tmp_files)
            print("Check 'LD_info_chr_all.txt' for results")
        elif not tmp_files and not lazy_results_list:
            print("No matching SNPs were found across any chromosome.")


def process_data_pairwise(
        file_path, r2threshold, population, maf_input, ref_file, imp_snp_list,
        ld_prune: bool = False,
        ld_prune_p_col: str | None = "P",
        ld_prune_p_threshold: float | None = None,
        ld_prune_p_threshold_mode: str = "below",
        ld_prune_z_col: str | None = "Z",
        ld_prune_z_threshold: float | None = 2.0,
        ld_prune_metric: str = "P",
        ld_prune_out_prefix: str = "LD_pruned",
        sig_keep: str | None = None, sig_metric: str = "P", sig_threshold: float = 5e-8,
        z_col: str = "Z", beta_col: str = "BETA", se_col: str = "SE", p_col: str = "P",
        snp_to_gene: bool = False, snps_genes_dir: str = "snps_genes_ref"
):
    # Determine execution strategy based on file size threshold (500MB)
    file_size_mb = os.path.getsize(file_path) / (1024 * 1024)
    is_large_file = file_size_mb > 500

    with pl.StringCache():
        chroms, rs_series, study_cols = _get_study_info_chunked(
            file_path, is_pairwise=True,
            sig_keep=sig_keep, sig_metric=sig_metric, sig_threshold=sig_threshold,
            z_col=z_col, beta_col=beta_col, se_col=se_col, p_col=p_col)

        # Build the SNP->Gene lookup once (reused for every output batch).
        gene_map = build_snp_gene_map(chroms, snps_genes_dir) if snp_to_gene else None

        # We still define the lazy frame to execute sinks/filters downstream in an out-of-core manner.
        study_lazy = pl.scan_csv(file_path, separator="\t")

        # Apply the same significance filter to the lazy scan so the Numba LD-pruner
        # (ld_prune branch below) operates on exactly the kept variant set.
        if sig_keep is not None and study_cols is not None:
            _sig_expr, _ = _significance_filter_expr(
                study_cols, sig_keep, sig_metric, sig_threshold,
                z_col, beta_col, se_col, p_col)
            if _sig_expr is not None:
                study_lazy = study_lazy.filter(_sig_expr)

        panel_cfg = {
            "1000G_hg38": {"fn": hg38_1kg_process_pairwise, "snp1_col": "SNP_A", "snp2_col": "SNP_B"},
            "1000G_hg38_high_cov": {"fn": hg38_1kg_process_high_coverage_pairwise, "snp1_col": "rsID1",
                                    "snp2_col": "rsID2"},
            "HGDP": {"fn": HGDP_process_pairwise, "snp1_col": "rsID1", "snp2_col": "rsID2"},
            "UKBB": {"fn": UKBB_process_pairwise, "snp1_col": "rsID1", "snp2_col": "rsID2"},
            "TOP_LD": {"fn": TOP_LD_process_pairwise, "snp1_col": "rsID1", "snp2_col": "rsID2"},
            "Hap_Map": {"fn": Hap_Map_process_pairwise, "snp1_col": "rsID1", "snp2_col": "rsID2"},
            "Pheno_Scanner": {"fn": pheno_Scanner_process_pairwise, "snp1_col": "rsID1", "snp2_col": "rsID2"},
            "LASI_DAD": {"fn": LASI_DAD_process_pairwise, "snp1_col": "rsID1", "snp2_col": "rsID2"},
        }

        if ref_file not in panel_cfg: raise ValueError(f"Unsupported ref_panel: {ref_file}")

        proc_fn = panel_cfg[ref_file]["fn"]
        snp1_col = "SNP1" if ld_prune else panel_cfg[ref_file]["snp1_col"]
        snp2_col = "SNP2" if ld_prune else panel_cfg[ref_file]["snp2_col"]

        batch_size = 4 if is_large_file else max(1, len(chroms))

        if not ld_prune:
            tmp_files = []

            for i in range(0, len(chroms), batch_size):
                batch_chroms = chroms[i:i + batch_size]
                lazy_results_list = []

                for chrom in batch_chroms:
                    lazy_data = proc_fn(rs_series, r2threshold, population, maf_input, chrom, imp_snp_list,
                                        ld_prune=False)
                    if lazy_data is not None:
                        lazy_results_list.append(lazy_data)
                        print(f"Evaluation graph for chr{chrom} constructed.")

                if lazy_results_list:
                    batch_lazy = pl.concat(lazy_results_list, how="vertical_relaxed")
                    # Attach Gene column(s) beside each rsID/SNP column before writing.
                    out_lazy = annotate_genes_lazy(batch_lazy, gene_map)
                    if is_large_file and batch_size < len(chroms):
                        tmp_file = f"LD_info_pairwise_tmp_batch_{i}.txt"
                        print(f"Streaming batch {i // batch_size + 1} (chroms: {batch_chroms}) to disk...")
                        out_lazy.sink_csv(tmp_file, separator="\t")
                        tmp_files.append(tmp_file)
                    else:
                        print("Streaming highly parallelized pairwise query across all chromosomes...")
                        out_lazy.sink_csv("LD_info_chr_all_pairwise.txt", separator="\t")
                        print("Check 'LD_info_chr_all_pairwise.txt' for results")
                        return pl.DataFrame(), None

                gc.collect()

            if is_large_file and tmp_files:
                print("Merging temporary batch files into final output...")
                _merge_tmp_files("LD_info_chr_all_pairwise.txt", tmp_files)
                print("Check 'LD_info_chr_all_pairwise.txt' for results")
                return pl.DataFrame(), None
            else:
                print("No matching SNPs were found across any chromosome.")
                return pl.DataFrame(), None

        else:  # LD Pruning Mode
            out_name = f"{ld_prune_out_prefix}_kept.txt"
            tmp_files = []

            for i in range(0, len(chroms), batch_size):
                batch_chroms = chroms[i:i + batch_size]
                kept_lazy_list = []

                for chrom in batch_chroms:
                    lazy_data = proc_fn(rs_series, r2threshold, population, maf_input, chrom, imp_snp_list,
                                        ld_prune=True)
                    if lazy_data is None:
                        continue

                    print(f"Evaluating LD graph for Numba pruning engine on chr{chrom}...")

                    # Subset GWAS for the current chromosome to prevent pruning unrelated SNPs
                    if "CHR" in study_cols:
                        chr_study_lazy = study_lazy.filter(
                            pl.col("CHR").cast(pl.Utf8).str.replace(r"^(chr|CHR|Chr)", "").str.split(
                                ".").list.first() == str(chrom)
                        )
                    else:
                        if i == 0:
                            print(
                                f"WARNING: 'CHR' column missing in GWAS data. Pruning against chr{chrom} LD will output ALL chromosomes in this file!")
                        chr_study_lazy = study_lazy

                    kept_gwas_lazy = ld_prune_pairs(
                        study_lazy=chr_study_lazy,
                        ld_pairs_lazy=lazy_data,
                        snp_col="SNP",
                        p_col=ld_prune_p_col,
                        p_threshold=ld_prune_p_threshold,
                        p_threshold_mode=ld_prune_p_threshold_mode,
                        z_col=ld_prune_z_col,
                        z_threshold=ld_prune_z_threshold,
                        prune_metric=ld_prune_metric,
                        snp1_col=snp1_col, snp2_col=snp2_col
                    )
                    kept_lazy_list.append(kept_gwas_lazy)

                if kept_lazy_list:
                    batch_lazy = pl.concat(kept_lazy_list, how="vertical_relaxed")
                    # Attach Gene column(s) beside the kept-GWAS SNP column before writing.
                    out_lazy = annotate_genes_lazy(batch_lazy, gene_map)
                    if is_large_file and batch_size < len(chroms):
                        tmp_file = f"{ld_prune_out_prefix}_tmp_batch_{i}.txt"
                        print(f"Streaming pruning batch {i // batch_size + 1} (chroms: {batch_chroms}) to disk...")
                        out_lazy.sink_csv(tmp_file, separator="\t")
                        tmp_files.append(tmp_file)
                    else:
                        print("Streaming highly parallelized LD pruning across all chromosomes...")
                        out_lazy.sink_csv(out_name, separator="\t")
                        print(f"LD-pruned GWAS successfully saved to '{out_name}'.")
                        return pl.DataFrame(), pl.scan_csv(out_name, separator="\t")

                gc.collect()

            if is_large_file and tmp_files:
                print("Merging pruned results from all batches...")
                _merge_tmp_files(out_name, tmp_files)
                print(f"LD-pruned GWAS successfully saved to '{out_name}'.")
                # Return a fresh flat scan to prevent downstream AST crashing
                return pl.DataFrame(), pl.scan_csv(out_name, separator="\t")
            else:
                print("No SNPs remained after pruning across any chromosomes.")
                return pl.DataFrame(), None