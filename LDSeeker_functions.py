#!/usr/bin/env python
import os
import re
import gc
import shutil
import polars as pl
import numpy as np
from numba import njit

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
        snp1_col: str = "SNP1",
        snp2_col: str = "SNP2",
        threshold: float | None = None,
        threshold_mode: str = "below"
) -> pl.LazyFrame:
    """
    Massively optimized Polars + Numba LD pruning.
    Now accepts a LazyFrame for ld_pairs_lazy to prevent memory crashes.
    """
    active_lazy = study_lazy
    study_cols = active_lazy.collect_schema().names()

    # 1) Pre-filter based on P-value threshold if provided
    if (threshold is not None) and (p_col is not None) and (p_col in study_cols):
        thresh_expr = pl.lit(threshold, dtype=pl.Float32)
        p_expr = pl.col(p_col).cast(pl.Float32)
        if threshold_mode == "above":
            active_lazy = active_lazy.filter(p_expr > thresh_expr)
        else:
            active_lazy = active_lazy.filter(p_expr < thresh_expr)

    # 2) Extract ordered SNPs map it to Int32 IDs natively
    sort_exprs = []
    if (p_col is not None) and (p_col in study_cols):
        sort_exprs.append(p_col)

    # Use streaming=True to ensure out-of-core chunk processing for massive sorting
    ordered_snps_df = (
        active_lazy.select([snp_col, p_col] if sort_exprs else [snp_col])
        .drop_nulls(subset=[snp_col])
        .unique(subset=[snp_col], maintain_order=False)
        .sort(sort_exprs, descending=False) if sort_exprs else active_lazy.select(snp_col).unique()
    ).collect(streaming=True)

    # Add 0 to N-1 integer IDs mapping to feed Numba natively
    n_snps = len(ordered_snps_df)
    ordered_snps_df = ordered_snps_df.with_columns(
        pl.int_range(0, n_snps, dtype=pl.Int32).alias("_numba_id")
    )

    # 3) Memory-Safe Graph Construction:
    # Drops the massive String columns and R2 values at the Rust layer.
    # We only materialize `u` and `v` Integer Arrays!
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

    # 4) Run extreme speed graph algorithm in C via Numba
    edges_u = ld_mapped["u"].to_numpy()
    edges_v = ld_mapped["v"].to_numpy()

    keep_mask = _numba_greedy_prune(n_snps, edges_u, edges_v)

    # Aggressively return internal array memory to OS before continuing
    del ld_mapped, edges_u, edges_v

    # 5) Filter by the boolean array and return Semi-Join AST node
    kept_snps = ordered_snps_df.filter(pl.Series(keep_mask))[snp_col]
    return _semi_join_filter(active_lazy, snp_col, kept_snps)


# -------------------------------------------------------------------
# TOP-LD (pairwise)
# -------------------------------------------------------------------

def TOP_LD_info_pairwise(rs_series, chrom, population, maf_threshold, R2_threshold, imp_snp_list=None, ld_prune=False):
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

    maf_lazy = _semi_join_filter(maf_lazy, "rsID", snp_series)
    valid_ids1 = _semi_join_filter(maf_lazy, "rsID", rs_series).select("Uniq_ID")
    valid_ids2 = maf_lazy.select("Uniq_ID")

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
                                  ld_prune=False):
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
    maf_df = _semi_join_filter(maf_df, "rsID", snp_series)

    ld_df = ld_raw.filter(pl.col(cols[3]) != pl.col(cols[4]))
    ld_df = _semi_join_filter(ld_df, cols[3], rs_series)
    ld_df = _semi_join_filter(ld_df, cols[4], snp_series)

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
                                        ld_prune=False):
    if R2_threshold is not None and R2_threshold < 0.8:
        print("Pheno Scanner includes data with R2 ≥ 0.8. Setting R2_threshold = 0.8")
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

    maf2_lazy = _semi_join_filter(maf1_lazy, "rsid", snp_series).rename({
        "hg19_coordinates": "pos2(hg19)", maf_pop: "MAF2", "a1": "ALT2", "a2": "REF2"
    }).select(["pos2(hg19)", "rsid", "MAF2", "ALT2", "REF2"])

    maf1_lazy = _semi_join_filter(maf1_lazy, "rsid", rs_series).rename({
        "hg19_coordinates": "ref_hg19_coordinates", "rsid": "ref_rsid", maf_pop: "MAF1", "a1": "ALT1", "a2": "REF1"
    }).select(["ref_hg19_coordinates", "ref_rsid", "MAF1", "ALT1", "REF1"])

    ld_lazy = pl.scan_parquet(ld_file).select(["ref_hg19_coordinates", "ref_rsid", "rsid", "r2", "r", "dprime"]).filter(
        pl.col("ref_rsid") != pl.col("rsid"))
    ld_lazy = _semi_join_filter(ld_lazy, "ref_rsid", rs_series)
    ld_lazy = _semi_join_filter(ld_lazy, "rsid", snp_series)

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

def hg38_1kg_LD_info_pairwise(rs_series, chrom, population, maf_threshold, R2_threshold, imp_snp_list, ld_prune=False):
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
    maf_lazy = _semi_join_filter(maf_lazy, "SNP", snp_series).with_columns(pl.col("MAF").cast(pl.Float32))

    if maf_threshold is not None:
        maf_lazy = maf_lazy.filter(pl.col("MAF") >= pl.lit(maf_threshold, dtype=pl.Float32))

    ld_lazy = ld_lazy.select(["CHR_A", "SNP_A", "CHR_B", "SNP_B", "R"])
    ld_lazy = _semi_join_filter(ld_lazy, "SNP_A", rs_series)
    ld_lazy = _semi_join_filter(ld_lazy, "SNP_B", snp_series)
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


def _generic_vcor_pairwise(rs_series, R2_threshold, maf_threshold, ld_file, imp_snp_list, ld_prune, is_phased=False):
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
    ld_lazy = _semi_join_filter(ld_lazy, "ID_A", rs_series)
    ld_lazy = _semi_join_filter(ld_lazy, "ID_B", snp_series)

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
                                            ld_prune=False):
    print(f"Building query plan: 1KGP high-cov files ({population}) chr{chrom}...")
    ld_file = f'D:/ref/1KGP_high_coverage/LD_{population}_r_unphased/{population}_chr{chrom}_r_unphased.vcor.parquet'
    return _generic_vcor_pairwise(rs_series, R2_threshold, maf_threshold, ld_file, imp_snp_list, ld_prune,
                                  is_phased=False)


def hg38_1kg_process_high_coverage_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list,
                                            ld_prune=False):
    return hg38_1kg_LD_info_high_coverage_pairwise(rs_series, r2threshold, population, maf_input, chromosome,
                                                   imp_snp_list, ld_prune)


def HGDP_LD_info_pairwise(rs_series, R2_threshold, population, maf_threshold, chrom, imp_snp_list, ld_prune=False):
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
                                  is_phased=True)


def HGDP_process_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list, ld_prune=False):
    return HGDP_LD_info_pairwise(rs_series, r2threshold, population, maf_input, chromosome, imp_snp_list, ld_prune)


def UKBB_LD_info_pairwise(rs_series, R2_threshold, population, maf_threshold, chrom, imp_snp_list, ld_prune=False):
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
    ld_lazy = _semi_join_filter(ld_lazy, "snp1", rs_series)
    ld_lazy = _semi_join_filter(ld_lazy, "snp2", snp_series)

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
# Non-Pairwise Functions (Wrappers for backward compatibility)
# -------------------------------------------------------------------
def hg38_1kg_LD_info_high_coverage(rs_series, chrom, pop, maf, R2, imp): return hg38_1kg_LD_info_high_coverage_pairwise(
    rs_series, R2, pop, maf, chrom, imp)


def hg38_1kg_process_high_coverage(rs_series, r2, pop, maf, chrom, imp): return hg38_1kg_LD_info_high_coverage(
    rs_series, chrom, pop, maf, r2, imp)


def HGDP_LD_info(rs_series, chrom, pop, maf, R2, imp): return HGDP_LD_info_pairwise(rs_series, R2, pop, maf, chrom, imp)


def HGDP_process(rs_series, r2, pop, maf, chrom, imp): return HGDP_LD_info(rs_series, chrom, pop, maf, r2, imp)


def UKBB_LD_info(rs_series, chrom, pop, maf, R2, imp): return UKBB_LD_info_pairwise(rs_series, R2, pop, maf, chrom, imp)


def UKBB_process(rs_series, r2, pop, maf, chrom, imp): return UKBB_LD_info(rs_series, chrom, pop, maf, r2, imp)


def hg38_1kg_LD_info(rs_series, chrom, pop, maf, R2, imp): return hg38_1kg_LD_info_pairwise(rs_series, chrom, pop, maf,
                                                                                            R2, imp)


def hg38_1kg_process(rs_series, r2, pop, maf, chrom, imp): return hg38_1kg_LD_info(rs_series, chrom, pop, maf, r2, imp)


def Hap_Map_LD_info_dask(rs_series, chrom, pop, maf, R2, imp): return Hap_Map_LD_info_dask_pairwise(rs_series, chrom,
                                                                                                    pop, maf, R2, imp)


def Hap_Map_process(rs_series, r2, pop, maf, chrom, imp): return Hap_Map_LD_info_dask(rs_series, chrom, pop, maf, r2,
                                                                                      imp)


def pheno_Scanner_LD_info_dask(rs_series, chrom, pop, maf, R2, imp): return pheno_Scanner_LD_info_dask_pairwise(
    rs_series, chrom, pop, maf, R2, imp)


def pheno_Scanner_process(rs_series, r2, pop, maf, chrom, imp): return pheno_Scanner_LD_info_dask(rs_series, chrom, pop,
                                                                                                  maf, r2, imp)


def TOP_LD_info(rs_series, chrom, pop, maf, R2, imp=None): return TOP_LD_info_pairwise(rs_series, chrom, pop, maf, R2,
                                                                                       imp)


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


def _get_study_info_chunked(file_path: str, is_pairwise: bool = False):
    """Loads the GWAS input file in chunks to extract SNPs and Chromosomes without memory spikes."""
    print("Scanning input file in chunks to preserve memory...")

    # Using read_csv_batched to iterate over massive files natively chunk-by-chunk
    reader = pl.read_csv_batched(file_path, separator="\t")

    chroms_set = set()
    snps_chunks = []
    study_cols = None

    while True:
        batches = reader.next_batches(5)  # Process 5 chunks into memory at a time
        if not batches:
            break
        for chunk in batches:
            if study_cols is None:
                study_cols = chunk.columns

            if "CHR" in study_cols:
                chroms_set.update(chunk.get_column("CHR").drop_nulls().unique().to_list())
            if "SNP" in study_cols:
                snps_chunks.append(chunk.get_column("SNP").drop_nulls().unique())

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


def process_data(file_path, r2threshold, population, maf_input, ref_file, imp_snp_list):
    # Determine execution strategy based on file size threshold (500MB)
    file_size_mb = os.path.getsize(file_path) / (1024 * 1024)
    is_large_file = file_size_mb > 500

    with pl.StringCache():
        chroms, rs_series, study_cols = _get_study_info_chunked(file_path, is_pairwise=False)

        processors = {
            "1000G_hg38": hg38_1kg_process, "1000G_hg38_high_cov": hg38_1kg_process_high_coverage,
            "UKBB": UKBB_process, "HGDP": HGDP_process, "TOP_LD": TOP_LD_process,
            "Pheno_Scanner": pheno_Scanner_process, "Hap_Map": Hap_Map_process
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
                if is_large_file and batch_size < len(chroms):
                    tmp_file = f"LD_info_tmp_batch_{i}.txt"
                    print(f"Streaming batch {i // batch_size + 1} (chroms: {batch_chroms}) to disk...")
                    batch_lazy.sink_csv(tmp_file, separator="\t")
                    tmp_files.append(tmp_file)
                else:
                    print("Streaming highly parallelized query across all chromosomes...")
                    batch_lazy.sink_csv("LD_info_chr_all.txt", separator="\t")
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
        ld_prune: bool = False, ld_prune_p_col: str | None = "P", ld_prune_out_prefix: str = "LD_pruned",
        ld_prune_threshold: float | None = None, ld_prune_mode: str = "below"
):
    # Determine execution strategy based on file size threshold (500MB)
    file_size_mb = os.path.getsize(file_path) / (1024 * 1024)
    is_large_file = file_size_mb > 500

    with pl.StringCache():
        chroms, rs_series, study_cols = _get_study_info_chunked(file_path, is_pairwise=True)

        # We still define the lazy frame to execute sinks/filters downstream in an out-of-core manner.
        study_lazy = pl.scan_csv(file_path, separator="\t")

        panel_cfg = {
            "1000G_hg38": {"fn": hg38_1kg_process_pairwise, "snp1_col": "SNP_A", "snp2_col": "SNP_B"},
            "1000G_hg38_high_cov": {"fn": hg38_1kg_process_high_coverage_pairwise, "snp1_col": "rsID1",
                                    "snp2_col": "rsID2"},
            "HGDP": {"fn": HGDP_process_pairwise, "snp1_col": "rsID1", "snp2_col": "rsID2"},
            "UKBB": {"fn": UKBB_process_pairwise, "snp1_col": "rsID1", "snp2_col": "rsID2"},
            "TOP_LD": {"fn": TOP_LD_process_pairwise, "snp1_col": "rsID1", "snp2_col": "rsID2"},
            "Hap_Map": {"fn": Hap_Map_process_pairwise, "snp1_col": "rsID1", "snp2_col": "rsID2"},
            "Pheno_Scanner": {"fn": pheno_Scanner_process_pairwise, "snp1_col": "rsID1", "snp2_col": "rsID2"},
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
                    lazy_data = proc_fn(rs_series, r2threshold, population, maf_input, chrom, imp_snp_list, ld_prune=False)
                    if lazy_data is not None:
                        lazy_results_list.append(lazy_data)
                        print(f"Evaluation graph for chr{chrom} constructed.")

                if lazy_results_list:
                    batch_lazy = pl.concat(lazy_results_list, how="vertical_relaxed")
                    if is_large_file and batch_size < len(chroms):
                        tmp_file = f"LD_info_pairwise_tmp_batch_{i}.txt"
                        print(f"Streaming batch {i // batch_size + 1} (chroms: {batch_chroms}) to disk...")
                        batch_lazy.sink_csv(tmp_file, separator="\t")
                        tmp_files.append(tmp_file)
                    else:
                        print("Streaming highly parallelized pairwise query across all chromosomes...")
                        batch_lazy.sink_csv("LD_info_chr_all_pairwise.txt", separator="\t")
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
                    lazy_data = proc_fn(rs_series, r2threshold, population, maf_input, chrom, imp_snp_list, ld_prune=True)
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
                            print(f"WARNING: 'CHR' column missing in GWAS data. Pruning against chr{chrom} LD will output ALL chromosomes in this file!")
                        chr_study_lazy = study_lazy

                    kept_gwas_lazy = ld_prune_pairs(
                        study_lazy=chr_study_lazy,
                        ld_pairs_lazy=lazy_data,
                        snp_col="SNP", p_col=ld_prune_p_col,
                        snp1_col=snp1_col, snp2_col=snp2_col, threshold=ld_prune_threshold, threshold_mode=ld_prune_mode
                    )
                    kept_lazy_list.append(kept_gwas_lazy)

                if kept_lazy_list:
                    batch_lazy = pl.concat(kept_lazy_list, how="vertical_relaxed")
                    if is_large_file and batch_size < len(chroms):
                        tmp_file = f"{ld_prune_out_prefix}_tmp_batch_{i}.txt"
                        print(f"Streaming pruning batch {i // batch_size + 1} (chroms: {batch_chroms}) to disk...")
                        batch_lazy.sink_csv(tmp_file, separator="\t")
                        tmp_files.append(tmp_file)
                    else:
                        print("Streaming highly parallelized LD pruning across all chromosomes...")
                        batch_lazy.sink_csv(out_name, separator="\t")
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
