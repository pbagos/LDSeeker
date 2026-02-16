#!/usr/bin/env python
import re
import pandas as pd
import gc
import polars as pl
from collections import defaultdict

pd.set_option('display.width', None)
# -------------------------------------------------------------------
# LD pruning helper
# -------------------------------------------------------------------

def ld_prune_pairs(
    study_df: pd.DataFrame,
    ld_pairs,
    snp_col: str = "SNP",
    p_col: str | None = "P",
    snp1_col: str = "SNP1",
    snp2_col: str = "SNP2",
    threshold: float | None = None,
    threshold_mode: str = "below"
):
    """
    LD pruning on a GWAS table given pairwise LD information.

    Parameters
    ----------
    study_df : pandas.DataFrame
        Original GWAS / summary statistics. Must contain `snp_col` and
        optionally `p_col`.

    ld_pairs : polars.DataFrame or pandas.DataFrame
        Pairwise LD table (already filtered at desired R²). Must contain
        `snp1_col` and `snp2_col`.

    snp_col : str
        Name of SNP ID column in `study_df` (e.g. "SNP" or "rsID").

    p_col : str or None
        Column used to rank SNPs (e.g. "P"). Lower values are considered
        more important. If None or missing, the input order of `study_df`
        is used.

    snp1_col, snp2_col : str
        Column names of the two SNP IDs in `ld_pairs`.

    threshold : float or None, optional
        If provided, filter `study_df` based on `p_col` before pruning.

    threshold_mode : str, optional
        "below" to keep SNPs with P < threshold (default).
        "above" to keep SNPs with P > threshold.

    Returns
    -------
    kept_df : pandas.DataFrame
        View of `study_df` restricted to SNPs that REMAIN after LD pruning
        (the LD-independent set).

    pruned_df : pandas.DataFrame
        View of `study_df` restricted to SNPs that were REMOVED by
        LD pruning (in LD with some more important SNP).
    """

    # 0) Pre-filter based on threshold if provided
    # We create 'active_df' which is the subset of study_df we want to process
    active_df = study_df.copy()

    if (threshold is not None) and (p_col is not None) and (p_col in study_df.columns):
        if threshold_mode == "above":
            active_df = active_df[active_df[p_col] > threshold]
        else:
            # Default is "below"
            active_df = active_df[active_df[p_col] < threshold]

    # 1) Build ranking: smallest p first if p_col exists, else original order
    if (p_col is not None) and (p_col in active_df.columns):
        ranking = (
            active_df[[snp_col, p_col]]
            .dropna(subset=[p_col])
            .drop_duplicates(subset=[snp_col])
            .sort_values(p_col, ascending=True)
        )
        ordered_snps = ranking[snp_col].tolist()
    else:
        ranking = active_df[[snp_col]].drop_duplicates(subset=[snp_col])
        ordered_snps = ranking[snp_col].tolist()

    # 2) Build adjacency list (LD graph) from pairwise LD table
    neighbors: dict[str, set[str]] = defaultdict(set)

    if isinstance(ld_pairs, pl.DataFrame):
        s1 = ld_pairs.select(snp1_col).to_series().to_list()
        s2 = ld_pairs.select(snp2_col).to_series().to_list()
    else:
        # assume pandas-like
        s1 = ld_pairs[snp1_col].tolist()
        s2 = ld_pairs[snp2_col].tolist()

    for a, b in zip(s1, s2):
        if pd.isna(a) or pd.isna(b):
            continue
        neighbors[a].add(b)
        neighbors[b].add(a)

    kept_snps: set[str] = set()
    excluded: set[str] = set()   # SNPs that cannot be kept anymore

    # 3) Greedy pruning: walk SNPs in ranked order;
    #    keep if not excluded, then exclude its neighbors.
    for snp in ordered_snps:
        if snp in excluded or snp in kept_snps:
            continue
        kept_snps.add(snp)
        for nb in neighbors.get(snp, ()):
            excluded.add(nb)

    # 4) SNPs with no LD edges (not in neighbors) → keep them
    # Iterate specifically over the filtered 'active_df' unique SNPs
    for snp in active_df[snp_col].drop_duplicates():
        if (snp not in kept_snps) and (snp not in excluded):
            kept_snps.add(snp)

    all_snps = set(active_df[snp_col].unique())
    pruned_snps = all_snps - kept_snps

    # 5) Split the original GWAS into "kept" and "pruned-out"
    # Note: We return subsets of the *filtered* active_df
    kept_df = active_df[active_df[snp_col].isin(kept_snps)].copy()
    pruned_df = active_df[active_df[snp_col].isin(pruned_snps)].copy()

    return kept_df, pruned_df



# -------------------------------------------------------------------
# TOP-LD (pairwise)
# -------------------------------------------------------------------

def TOP_LD_info_pairwise(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list=None):
    """
       Optimized version of TOP_LD_info to minimize memory footprint.
       Uses semi-joins for early filtering of LD pairs.
       """
    print(f"Loading TOP-LD files for {population} chr{chrom}...")

    # 1) Combine rsIDs up-front to define the universe of interesting SNPs
    target_rsids = set(rs_list)
    if imp_snp_list:
        all_rsids = list(target_rsids | set(imp_snp_list))
    else:
        all_rsids = list(target_rsids)

    maf_path = (
        f"D:/ref/TOP_LD/{population}/SNV/"
        f"{population}_chr{chrom}_no_filter_0.2_1000000_info_annotation.parquet"
    )
    ld_path = (
        f"D:/ref/TOP_LD/{population}/SNV/"
        f"{population}_chr{chrom}_no_filter_0.2_1000000_LD.parquet"
    )

    # 2) Scan MAF and filter immediately
    # We create a lazy frame of ONLY the SNPs we care about.
    maf_lazy = (
        pl.scan_parquet(maf_path)
        .filter((pl.col("MAF") >= maf_threshold) & (pl.col("rsID").is_in(all_rsids)))
        .select(["Uniq_ID", "rsID", "MAF", "REF", "ALT"])
    )

    # We need a lightweight version (just IDs) for filtering the LD table efficiently
    valid_ids = maf_lazy.select("Uniq_ID")

    # 3) Optimize LD loading: Use semi-joins to drop irrelevant rows early
    # "semi" join acts like a filter: keep rows in LEFT where key exists in RIGHT.
    # We filter LD rows where *both* Uniq_ID_1 and Uniq_ID_2 exist in our relevant MAF list.
    ld_lazy = (
        pl.scan_parquet(ld_path)
        .filter(pl.col("R2") >= R2_threshold)
        .join(valid_ids, left_on="Uniq_ID_1", right_on="Uniq_ID", how="semi")
        .join(valid_ids, left_on="Uniq_ID_2", right_on="Uniq_ID", how="semi")
    )

    # 4) Prepare two MAF views for the final enrichment
    # Since we've already used semi-joins to filter LD, these inner joins will only
    # attach data to the rows that survived, rather than exploding the search space.
    maf1 = maf_lazy.rename({
        "Uniq_ID": "Uniq_ID_1", "rsID": "rsID1", "MAF": "MAF1",
        "REF": "REF1", "ALT": "ALT1"
    })
    maf2 = maf_lazy.rename({
        "Uniq_ID": "Uniq_ID_2", "rsID": "rsID2", "MAF": "MAF2",
        "REF": "REF2", "ALT": "ALT2"
    })

    # 5) Join LD ↔ MAF1 ↔ MAF2
    joined = (
        ld_lazy
        .join(maf1, on="Uniq_ID_1", how="inner")
        .join(maf2, on="Uniq_ID_2", how="inner")
    )

    # 6) Apply specific directional filtering if imp_snp_list is provided
    # The semi-joins ensured ID1 and ID2 are in the "all_rsids" bucket.
    # This step ensures the specific relationship (rsID1 in rs_list AND rsID2 in imp_snp_list).
    if imp_snp_list:
        joined = joined.filter(
            (pl.col("rsID1").is_in(list(target_rsids))) &
            (pl.col("rsID2").is_in(imp_snp_list))
        )

    # 7) Select + rename final output columns
    final_lazy = joined.select([
        pl.col("Uniq_ID_1").alias("pos1"),
        pl.col("Uniq_ID_2").alias("pos2"),
        "R2", "+/-corr", "Dprime",
        "rsID1", "rsID2",
        "MAF1", "MAF2",
        "REF1", "ALT1", "REF2", "ALT2"
    ])

    # 8) Execute without streaming
    # Note: We removed streaming=True because it is deprecated for Parquet in newer Polars versions.
    # The semi-join optimizations above sufficiently reduce memory usage for the standard engine.
    print("Executing query...")
    result = final_lazy.collect()

    if result.is_empty():
        print("No SNPs found.")

    # 9) FINAL FILTER: ensure both SNP_A and SNP_B are in the original snp_set
    result = result.filter(
        pl.col("rsID1").is_in(rs_list) &
        pl.col("rsID2").is_in(rs_list)
    )

    # 10) Cleanup
    del maf_lazy, maf1, maf2, ld_lazy, joined, final_lazy
    gc.collect()

    return result


def TOP_LD_process_pairwise(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data
    outputData = TOP_LD_info_pairwise(
        list(study_df['SNP']),
        chromosome,
        population,
        maf_input,
        r2threshold,
        imp_snp_list
    )
    return outputData


# -------------------------------------------------------------------
# HapMap (pairwise)
# -------------------------------------------------------------------

def Hap_Map_LD_info_dask_pairwise(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    print(f"Loading Hap Map files ({population}) ...")

    # 1) Validate population
    valid_pops = ['YRI', 'CHB', 'JPT', 'CEU', 'MKK', 'LWK',
                  'CHD', 'GIH', 'TSI', 'MEX', 'ASW']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths
    maf_file = f'D:/ref/Hap_Map/allele_freqs_chr{chrom}_{population}_phase3.2_nr.b36_fwd.parquet'
    ld_file = f'D:/ref/Hap_Map/ld_chr{chrom}_{population}.parquet'

    # 3) Read and preprocess MAF table
    maf_df = (
        pl.read_parquet(maf_file)
        .with_columns(pl.col("otherallele_freq").alias("MAF"))
        .select(["rs#", "refallele", "otherallele", "MAF"])
        .rename({
            "rs#": "rsID",
            "refallele": "REF",
            "otherallele": "ALT"
        })
        .filter(pl.col("MAF") >= float(maf_threshold))
    )

    # 4) Read and preprocess LD table
    ld_raw = pl.read_parquet(ld_file)
    cols = ld_raw.columns

    ld_df = (
        ld_raw
        .filter(
            (pl.col(cols[3]) != pl.col(cols[4])) &
            (pl.col(cols[6]) >= float(R2_threshold))
        )
        .rename({
            cols[0]: "pos1",
            cols[1]: "pos2",
            cols[2]: "pop",
            cols[3]: "rsID1",
            cols[4]: "rsID2",
            cols[5]: "Dprime",
            cols[6]: "R2"
        })
        .select(["pos1", "pos2", "pop", "rsID1", "rsID2", "Dprime", "R2"])
    )

    # 5) Join LD ↔ MAF on rsID1 and rsID2
    tmp = ld_df.join(maf_df, left_on="rsID1", right_on="rsID", how="inner")
    merged = tmp.join(
        maf_df,
        left_on="rsID2",
        right_on="rsID",
        how="inner",
        suffix="_2"
    )

    # 6) Select & reorder columns, and rename duplicates
    result = merged.select([
        pl.col("pos1"),
        pl.col("pos2"),
        pl.col("rsID1"),
        pl.col("rsID2"),
        pl.col("MAF").alias("MAF1"),
        pl.col("MAF_2").alias("MAF2"),
        pl.col("REF").alias("REF1"),
        pl.col("REF_2").alias("REF2"),
        pl.col("ALT").alias("ALT1"),
        pl.col("ALT_2").alias("ALT2"),
        pl.col("R2"),
        pl.col("Dprime")
    ])

    # 7) Apply the user’s SNP-list filter
    if imp_snp_list:
        result = result.filter(
            pl.col("rsID1").is_in(rs_list) &
            pl.col("rsID2").is_in(imp_snp_list)
        )
    else:
        result = result.filter(pl.col("rsID1").is_in(rs_list))

    # 8) FINAL FILTER: ensure both SNP_A and SNP_B are in the original snp_set
    result = result.filter(
        pl.col("rsID1").is_in(rs_list) &
        pl.col("rsID2").is_in(rs_list)
    )

    # 9) Clean up
    del ld_raw, ld_df, maf_df, tmp, merged
    gc.collect()

    return result


def Hap_Map_process_pairwise(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    outputData = Hap_Map_LD_info_dask_pairwise(
        list(study_df['SNP']),
        chromosome,
        population,
        maf_input,
        r2threshold,
        imp_snp_list
    )
    return outputData


# -------------------------------------------------------------------
# PhenoScanner (pairwise)
# -------------------------------------------------------------------

def pheno_Scanner_LD_info_dask_pairwise(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    # 1) Enforce minimum R2 threshold
    if R2_threshold < 0.8:
        print("Pheno Scanner includes data with R2 ≥ 0.8. Setting R2_threshold = 0.8")
        R2_threshold = 0.8

    print("Building lazy plan for Pheno Scanner files...")

    # 2) Paths & pop mapping
    maf_file = "D:/ref/Pheno_Scanner/1000G.parquet"
    ld_file = f"D:/ref/Pheno_Scanner/1000G_{population}/1000G_{population}_chr{chrom}.parquet"
    pop_map = {"EUR": "eur", "EAS": "eas", "AFR": "afr", "AMR": "amr", "SAS": "sas"}
    maf_pop = pop_map.get(population)
    if maf_pop is None:
        raise ValueError(f"Unsupported population: {population}")

    # 3) First MAF scan (for ref side → MAF1)
    maf1_lazy = (
        pl.scan_parquet(maf_file)
        .select(["hg19_coordinates", "chr", "rsid", maf_pop, "a1", "a2"])
        .filter(
            (pl.col("chr") == chrom) &
            (pl.col(maf_pop) != "-")
        )
        .with_columns(pl.col(maf_pop).cast(pl.Float64))
        .filter(pl.col(maf_pop) >= maf_threshold)
        .rename({
            "hg19_coordinates": "ref_hg19_coordinates",
            "rsid": "ref_rsid",
            maf_pop: "MAF1",
            "a1": "ALT1",
            "a2": "REF1",
        })
        .select(["ref_hg19_coordinates", "ref_rsid", "MAF1", "ALT1", "REF1"])
    )

    # 4) LD scan
    ld_lazy = (
        pl.scan_parquet(ld_file)
        .select(["ref_hg19_coordinates", "ref_rsid", "rsid", "r2", "r", "dprime"])
        .filter(
            (pl.col("ref_rsid") != pl.col("rsid")) &
            (pl.col("r2") >= R2_threshold)
        )
    )

    # 5) Second MAF scan (for query side → MAF2), keep `rsid` for join!
    maf2_lazy = (
        pl.scan_parquet(maf_file)
        .select(["hg19_coordinates", "chr", "rsid", maf_pop, "a1", "a2"])
        .filter(
            (pl.col("chr") == chrom) &
            (pl.col(maf_pop) != "-")
        )
        .with_columns(pl.col(maf_pop).cast(pl.Float64))
        .filter(pl.col(maf_pop) >= maf_threshold)
        .rename({
            "hg19_coordinates": "pos2(hg19)",
            maf_pop: "MAF2",
            "a1": "ALT2",
            "a2": "REF2",
        })
        # note: we keep `rsid` here so we can join on it
        .select(["pos2(hg19)", "rsid", "MAF2", "ALT2", "REF2"])
    )

    # 6) Join the three pieces
    joined = (
        ld_lazy
        .join(maf1_lazy, on="ref_rsid", how="inner")
        .join(maf2_lazy, on="rsid", how="inner")
    )

    # 7) Rename LD columns & the `rsid` from maf2 → final names
    renamed = joined.rename({
        "ref_rsid": "rsID1",
        "rsid": "rsID2",
        "ref_hg19_coordinates": "pos1(hg19)",
        "r2": "R2",
    })

    # 8) Filter by user-supplied SNP lists
    filtered = renamed.filter(
        pl.col("rsID2").is_in(rs_list) &
        (pl.col("rsID1").is_in(imp_snp_list) if imp_snp_list else pl.lit(True))
    )

    # 9) Extract numeric coords
    with_pos = filtered.with_columns([
        pl.col("pos1(hg19)").str.extract(r".*:(\d+)$", 1).cast(pl.Int64),
        pl.col("pos2(hg19)").str.extract(r".*:(\d+)$", 1).cast(pl.Int64),
    ])

    # 10) Final projection/order
    final_lazy = with_pos.select([
        "rsID1", "pos1(hg19)", "rsID2", "dprime",
        "pos2(hg19)", "R2", "r",
        "MAF1", "MAF2", "ALT1", "REF1", "ALT2", "REF2",
    ])

    # 11) Execute in streaming, collect to Polars
    final_pl = final_lazy.collect()

    final_pl = final_pl.filter(
        pl.col("rsID1").is_in(rs_list) &
        pl.col("rsID2").is_in(rs_list)
    )

    # 12) Cleanup
    gc.collect()
    if final_pl.is_empty():
        print("No SNPs found")

    final_pd = final_pl.with_columns(
        pl.lit(chrom).alias("CHR")
    )

    return final_pl


def pheno_Scanner_process_pairwise(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    outputData = pheno_Scanner_LD_info_dask_pairwise(
        list(study_df['SNP']),
        chromosome,
        population,
        maf_input,
        r2threshold,
        imp_snp_list
    )
    return outputData


# -------------------------------------------------------------------
# 1000G hg38 (pairwise)
# -------------------------------------------------------------------

def hg38_1kg_LD_info_pairwise(
    rs_list: list[str],
    chrom: int,
    population: str,
    maf_threshold: float | None,
    R2_threshold: float | None,
    imp_snp_list: list[str] | None
) -> pl.DataFrame:
    """
    Fetch pairwise LD (R and R2) for all variant pairs on `chrom` in `population`
    where both SNPs are in `rs_list` (or in `imp_snp_list` if provided),
    apply MAF and R2 thresholds, and return a Polars DataFrame with allele info.
    """
    print(f"Loading 1000 Genomes Project (hg38) files ({population}) ...")

    # 1) Validate population
    valid_pops = ['AMR', 'EUR', 'AFR', 'SAS', 'EAS']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths
    maf_file = f'D:/ref/1000G_hg38/1000G_{population}_0_01.parquet'
    ld_file = f'D:/ref/1000G_hg38/{population}/chr{chrom}_merged.parquet'

    # 3) Load & filter MAF data
    print("Loading and filtering MAF file...")
    tmp = pl.read_parquet(maf_file)
    maf_df = tmp.rename({
        old: new for old, new in zip(tmp.columns[:6],
                                     ['CHR', 'SNP', 'MAF', 'POS', 'REF', 'ALT'])
    })
    if maf_threshold is not None:
        maf_df = maf_df.filter(pl.col("MAF") >= maf_threshold)

    # decide which SNP-set to use
    snp_set = set(imp_snp_list) if imp_snp_list else set(rs_list)
    maf_df = maf_df.filter(pl.col("SNP").is_in(snp_set))

    # 4) Load & filter LD data
    print("Loading and filtering LD data...")
    ld_df = (
        pl.read_parquet(ld_file, columns=["CHR_A", "SNP_A", "CHR_B", "SNP_B", "R"])
        .with_columns((pl.col("R") ** 2).alias("R2"))
        .filter(
            pl.col("SNP_A").is_in(snp_set) &
            pl.col("SNP_B").is_in(snp_set)
        )
    )
    if R2_threshold is not None:
        ld_df = ld_df.filter(pl.col("R2") >= R2_threshold)

    # 5) Annotate allele A
    merged = (
        ld_df.join(
            maf_df,
            left_on=["CHR_A", "SNP_A"],
            right_on=["CHR", "SNP"],
            how="left"
        )
        .rename({
            "CHR_A": "CHR",
            "MAF": "MAF_A",
            "POS": "POS_A",
            "REF": "REF_A",
            "ALT": "ALT_A"
        })
    )

    # 6) Annotate allele B
    merged = (
        merged.join(
            maf_df,
            left_on=["CHR_B", "SNP_B"],
            right_on=["CHR", "SNP"],
            how="left"
        )
        .rename({
            "MAF": "MAF_B",
            "POS": "POS_B",
            "REF": "REF_B",
            "ALT": "ALT_B"
        })
        .drop(["CHR_B"])
    )

    # 7) Drop any pairs where one allele lacked MAF info
    merged = merged.drop_nulls(subset=["MAF_A", "MAF_B"])

    # 8) (Re-)apply MAF threshold on both alleles, if requested
    if maf_threshold is not None:
        merged = merged.filter(
            (pl.col("MAF_A") >= maf_threshold) &
            (pl.col("MAF_B") >= maf_threshold)
        )

    # 9) Reorder columns for readability
    preferred = [
        "SNP_A", "CHR", "POS_A", "SNP_B", "POS_B",
        "REF_A", "ALT_A", "MAF_A", "REF_B", "ALT_B", "MAF_B",
        "R", "R2"
    ]
    cols = merged.columns
    final_cols = [c for c in preferred if c in cols] + [c for c in cols if c not in preferred]

    # 10) FINAL FILTER: ensure both SNP_A and SNP_B are in the original snp_set
    merged = merged.filter(
        pl.col("SNP_A").is_in(snp_set) &
        pl.col("SNP_B").is_in(snp_set)
    )

    return merged.select(final_cols)


def hg38_1kg_process_pairwise(
    study_df,
    r2threshold: float | None,
    population: str,
    maf_input: float | None,
    chromosome: int,
    imp_snp_list: list[str] | None = None
) -> pl.DataFrame:
    """
    Wrapper to fetch pairwise LD for SNPs in `study_df['SNP']` on `chromosome`.
    """
    rs_list = list(study_df['SNP'])
    return hg38_1kg_LD_info_pairwise(
        rs_list=rs_list,
        chrom=chromosome,
        population=population,
        maf_threshold=maf_input,
        R2_threshold=r2threshold,
        imp_snp_list=imp_snp_list
    )
#                                           study_df,r2threshold, population, maf_input,chrom,imp_snp_list
def hg38_1kg_LD_info_high_coverage_pairwise(rs_list, R2_threshold, population, maf_threshold, chrom, imp_snp_list):
    print(f"Loading 1000 Genomes Project (hg38) - High Coverage files ({population}) ...")

    # 1) Validate population
    valid_pops = ['AMR', 'EUR', 'AFR', 'SAS', 'EAS']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")



    snp_set = set(imp_snp_list) if imp_snp_list else set(rs_list)
    # 2) File paths
    ld_file = f'D:/ref/1KGP_high_coverage/LD_{population}_r_unphased/{population}_chr{chrom}_r_unphased.vcor.parquet'

    # 4) Load LD data and filter with maf_threshold
    print("Loading and filtering LD data...")

    #ld_df = pl.read_parquet(ld_file)

    ld_df = pl.scan_parquet(ld_file)

    # ✅ filter EARLY (before loading everything)
    cols = ld_df.collect_schema().names()
    maf1_col = "MAF_1" if "MAF_1" in cols else ("MAF1" if "MAF1" in cols else "NONMAJ_FREQ_A")
    ld_df = ld_df.filter(pl.col(maf1_col) >= maf_threshold)

    # rename (only those that exist)
    rename_map = {
        "#CHROM_A": "CHR",
        "CHROM_B": "CHR2",
        'POS_A': 'POS1',
        'POS_B': 'POS2',
        'MAJ_A': 'MAJ1',
        'MAJ_B': 'MAJ2',
        'NONMAJ_A': 'NONMAJ1',
        'NONMAJ_B': 'NONMAJ2',
        "ID_A": "rsID1",
        "ID_B": "rsID2",
        "NONMAJ_FREQ_A": "MAF1",
        "NONMAJ_FREQ_B": "MAF2",
        "UNPHASED_R": "R",
    }
    #keep only chr1 column
    rename_map = {k: v for k, v in rename_map.items() if k in cols}

    ld_df = ld_df.rename(rename_map)

    # ✅ Filter for the rsID list BEFORE collecting
    ld_df = ld_df.filter(pl.col("rsID1").is_in(rs_list))

    # Convert "R" to a float before squaring it
    ld_df = ld_df.with_columns(
        (pl.col("R").cast(pl.Float64) ** 2).alias("R2")
    )

    if R2_threshold is not None:
        ld_df = ld_df.filter(pl.col("R2") >= R2_threshold)

    # --- NEW LOGIC START ---
    # Keep only CHR1 (renamed from CHROM_A or #CHROM_A)
    # Check if CHR2 exists in the schema before dropping
    current_cols = ld_df.collect_schema().names()
    if "CHR2" in current_cols:
        ld_df = ld_df.drop("CHR2")


    # 10) FINAL FILTER: ensure both SNP_A and SNP_B are in the original snp_set
    ld_df = ld_df.filter(
        pl.col("rsID1").is_in(snp_set) &
        pl.col("rsID2").is_in(snp_set)
    )


    # Now collect
    ld_df = ld_df.collect()

    return ld_df


def hg38_1kg_process_high_coverage_pairwise(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    outputData = hg38_1kg_LD_info_high_coverage_pairwise(
        list(study_df['SNP']),
        chromosome,
        population,
        maf_input,
        r2threshold,
        imp_snp_list
    )
    return outputData

def HGDP_LD_info_pairwise(rs_list, R2_threshold, population, maf_threshold, chrom, imp_snp_list):
    print(f"Loading Human Genome Diversity Project (HGDP) files ({population}) ...")
    snp_set = set(imp_snp_list) if imp_snp_list else set(rs_list)

    # 1) Validate population
    valid_pops = ['AMERICA', 'EUROPE', 'AFRICA', 'CENTRAL_SOUTH_ASIA', 'EAST_ASIA','MIDDLE_EAST','OCEANIA']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths
    ld_file = f'D:/ref/HGDP/LD_{population}_r_phased/{population}_chr{chrom}_r_phased.vcor.parquet'

    # 4) Load LD data and filter with maf_threshold
    print("Loading and filtering LD data...")

    #ld_df = pl.read_parquet(ld_file)

    ld_df = pl.scan_parquet(ld_file)

    # ✅ filter EARLY (before loading everything)
    cols = ld_df.collect_schema().names()
    maf1_col = "MAF_1" if "MAF_1" in cols else ("MAF1" if "MAF1" in cols else "NONMAJ_FREQ_A")
    ld_df = ld_df.filter(pl.col(maf1_col) >= maf_threshold)

    # rename (only those that exist)
    rename_map = {
        "#CHROM_A": "CHR",
        "CHROM_B": "CHR2",
         'POS_A': 'POS1',
        'POS_B': 'POS2',
        'MAJ_A': 'MAJ1',
        'MAJ_B': 'MAJ2',
        'NONMAJ_A': 'NONMAJ1',
        'NONMAJ_B': 'NONMAJ2',
        "ID_A": "rsID1",
        "ID_B": "rsID2",
        "NONMAJ_FREQ_A": "MAF1",
        "NONMAJ_FREQ_B": "MAF2",
        "DPRIME": "Dprime",
        "PHASED_R": "R",
    }
    #keep only chr1 column
    rename_map = {k: v for k, v in rename_map.items() if k in cols}

    ld_df = ld_df.rename(rename_map)

    #   Filter for the rsID list BEFORE collecting
    ld_df = ld_df.filter(pl.col("rsID1").is_in(rs_list))

    # Convert "R" to a float before squaring it
    ld_df = ld_df.with_columns(
        (pl.col("R").cast(pl.Float64) ** 2).alias("R2")
    )

    if R2_threshold is not None:
        ld_df = ld_df.filter(pl.col("R2") >= R2_threshold)


    # Keep only CHR1 (renamed from CHROM_A or #CHROM_A)
    # Check if CHR2 exists in the schema before dropping
    current_cols = ld_df.collect_schema().names()
    if "CHR2" in current_cols:
        ld_df = ld_df.drop("CHR2")

    # 10) FINAL FILTER: ensure both SNP_A and SNP_B are in the original snp_set
    ld_df = ld_df.filter(
        pl.col("rsID1").is_in(snp_set) &
        pl.col("rsID2").is_in(snp_set)
    )
    # Now collect
    ld_df = ld_df.collect()

    return ld_df


def HGDP_process_pairwise(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    outputData = HGDP_LD_info_pairwise(
        list(study_df['SNP']),
        r2threshold,
        population,
        maf_input,
        chromosome,
        imp_snp_list
    )
    return outputData


def UKBB_LD_info_pairwise(rs_list, R2_threshold, population, maf_threshold, chrom, imp_snp_list):
    print(f"Loading UK Biobank (UKBB) files ({population}) ...")
    snp_set = set(imp_snp_list) if imp_snp_list else set(rs_list)

    # 1) Validate population
    valid_pops = ['AMR', 'EUR', 'AFR', 'CSA', 'EAS','MID']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths
    ld_file = f'D:/ref/UKBB/{population}/chr_{chrom}_ld.parquet'

    # 4) Load LD data and filter with maf_threshold
    print("Loading and filtering LD data...")

    #ld_df = pl.read_parquet(ld_file)

    ld_df = pl.scan_parquet(ld_file)

    # ✅ filter EARLY (before loading everything)
    cols = ld_df.collect_schema().names()
    maf1_col = "maf1" if "maf1" in cols else ("MAF1" if "MAF1" in cols else "maf1")
    ld_df = ld_df.filter(pl.col(maf1_col) >= maf_threshold)

    # rename (only those that exist)
    rename_map = {
        "chr": "CHR",

         'bp1': 'POS1',
        'bp2': 'POS2',
        'a1_1': 'A1_1',
        'a1_2': 'A1_2',
        'a2_1': 'A2_1',
        'a2_2': 'A2_2',
        "snp1": "rsID1",
        "snp2": "rsID2",
        "maf1": "MAF1",
        "maf2": "MAF2",
        "r2": "R2",
        "r": "R",
    }
    #keep only chr1 column
    rename_map = {k: v for k, v in rename_map.items() if k in cols}

    ld_df = ld_df.rename(rename_map)

    #   Filter for the rsID list BEFORE collecting
    ld_df = ld_df.filter(pl.col("rsID1").is_in(rs_list))

    # Convert "R" to a float before squaring it
    ld_df = ld_df.with_columns(
        (pl.col("R").cast(pl.Float64) ** 2).alias("R2")
    )

    if R2_threshold is not None:
        ld_df = ld_df.filter(pl.col("R2") >= R2_threshold)

    # 10) FINAL FILTER: ensure both SNP_A and SNP_B are in the original snp_set
    ld_df = ld_df.filter(
        pl.col("rsID1").is_in(snp_set) &
        pl.col("rsID2").is_in(snp_set)
    )
    # Now collect
    ld_df = ld_df.collect()


    return ld_df


def UKBB_process_pairwise(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    outputData = UKBB_LD_info_pairwise(
        list(study_df['SNP']),
        r2threshold,
        population,
        maf_input,
        chromosome,
        imp_snp_list
    )
    return outputData




def hg38_1kg_LD_info_high_coverage(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    print(f"Loading 1000 Genomes Project (hg38) - High Coverage files ({population}) ...")

    # 1) Validate population
    valid_pops = ['AMR', 'EUR', 'AFR', 'SAS', 'EAS']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths
    ld_file = f'D:/ref/1KGP_high_coverage/LD_{population}_r_unphased/{population}_chr{chrom}_r_unphased.vcor.parquet'

    # 4) Load LD data and filter with maf_threshold
    print("Loading and filtering LD data...")

    #ld_df = pl.read_parquet(ld_file)

    ld_df = pl.scan_parquet(ld_file)

    # ✅ filter EARLY (before loading everything)
    cols = ld_df.collect_schema().names()
    maf1_col = "MAF_1" if "MAF_1" in cols else ("MAF1" if "MAF1" in cols else "NONMAJ_FREQ_A")
    ld_df = ld_df.filter(pl.col(maf1_col) >= maf_threshold)

    # rename (only those that exist)
    rename_map = {
        "#CHROM_A": "CHR",
        "CHROM_B": "CHR2",
        'POS_A': 'POS1',
        'POS_B': 'POS2',
        'MAJ_A': 'MAJ1',
        'MAJ_B': 'MAJ2',
        'NONMAJ_A': 'NONMAJ1',
        'NONMAJ_B': 'NONMAJ2',
        "ID_A": "rsID1",
        "ID_B": "rsID2",
        "NONMAJ_FREQ_A": "MAF1",
        "NONMAJ_FREQ_B": "MAF2",
        "UNPHASED_R": "R",
    }
    #keep only chr1 column
    rename_map = {k: v for k, v in rename_map.items() if k in cols}

    ld_df = ld_df.rename(rename_map)

    # ✅ Filter for the rsID list BEFORE collecting
    ld_df = ld_df.filter(pl.col("rsID1").is_in(rs_list))

    # Convert "R" to a float before squaring it
    ld_df = ld_df.with_columns(
        (pl.col("R").cast(pl.Float64) ** 2).alias("R2")
    )

    if R2_threshold is not None:
        ld_df = ld_df.filter(pl.col("R2") >= R2_threshold)

    # --- NEW LOGIC START ---
    # Keep only CHR1 (renamed from CHROM_A or #CHROM_A)
    # Check if CHR2 exists in the schema before dropping
    current_cols = ld_df.collect_schema().names()
    if "CHR2" in current_cols:
        ld_df = ld_df.drop("CHR2")
    # --- NEW LOGIC END ---

    # Now collect
    ld_df = ld_df.collect()

    return ld_df


def hg38_1kg_process_high_coverage(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    outputData = hg38_1kg_LD_info_high_coverage(
        list(study_df['SNP']),
        chromosome,
        population,
        maf_input,
        r2threshold,
        imp_snp_list
    )
    return outputData


def HGDP_LD_info(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    print(f"Loading Human Genome Diversity Project (HGDP) files ({population}) ...")

    # 1) Validate population
    valid_pops = ['AMERICA', 'EUROPE', 'AFRICA', 'CENTRAL_SOUTH_ASIA', 'EAST_ASIA','MIDDLE_EAST','OCEANIA']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths
    ld_file = f'D:/ref/HGDP/LD_{population}_r_phased/{population}_chr{chrom}_r_phased.vcor.parquet'

    # 4) Load LD data and filter with maf_threshold
    print("Loading and filtering LD data...")

    #ld_df = pl.read_parquet(ld_file)

    ld_df = pl.scan_parquet(ld_file)

    # ✅ filter EARLY (before loading everything)
    cols = ld_df.collect_schema().names()
    maf1_col = "MAF_1" if "MAF_1" in cols else ("MAF1" if "MAF1" in cols else "NONMAJ_FREQ_A")
    ld_df = ld_df.filter(pl.col(maf1_col) >= maf_threshold)

    # rename (only those that exist)
    rename_map = {
        "#CHROM_A": "CHR",
        "CHROM_B": "CHR2",
         'POS_A': 'POS1',
        'POS_B': 'POS2',
        'MAJ_A': 'MAJ1',
        'MAJ_B': 'MAJ2',
        'NONMAJ_A': 'NONMAJ1',
        'NONMAJ_B': 'NONMAJ2',
        "ID_A": "rsID1",
        "ID_B": "rsID2",
        "NONMAJ_FREQ_A": "MAF1",
        "NONMAJ_FREQ_B": "MAF2",
        "DPRIME": "Dprime",
        "PHASED_R": "R",
    }
    #keep only chr1 column
    rename_map = {k: v for k, v in rename_map.items() if k in cols}

    ld_df = ld_df.rename(rename_map)

    #   Filter for the rsID list BEFORE collecting
    ld_df = ld_df.filter(pl.col("rsID1").is_in(rs_list))

    # Convert "R" to a float before squaring it
    ld_df = ld_df.with_columns(
        (pl.col("R").cast(pl.Float64) ** 2).alias("R2")
    )

    if R2_threshold is not None:
        ld_df = ld_df.filter(pl.col("R2") >= R2_threshold)


    # Keep only CHR1 (renamed from CHROM_A or #CHROM_A)
    # Check if CHR2 exists in the schema before dropping
    current_cols = ld_df.collect_schema().names()
    if "CHR2" in current_cols:
        ld_df = ld_df.drop("CHR2")



    # Now collect
    ld_df = ld_df.collect()

    return ld_df


def HGDP_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    outputData = HGDP_LD_info(
        list(study_df['SNP']),
        chromosome,
        population,
        maf_input,
        r2threshold,
        imp_snp_list
    )
    return outputData


def UKBB_LD_info(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    print(f"Loading UK Biobank (UKBB) files ({population}) ...")

    # 1) Validate population
    valid_pops = ['AMR', 'EUR', 'AFR', 'CSA', 'EAS','MID']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths
    ld_file = f'D:/ref/UKBB/{population}/chr_{chrom}_ld.parquet'

    # 4) Load LD data and filter with maf_threshold
    print("Loading and filtering LD data...")

    #ld_df = pl.read_parquet(ld_file)

    ld_df = pl.scan_parquet(ld_file)

    # ✅ filter EARLY (before loading everything)
    cols = ld_df.collect_schema().names()
    maf1_col = "maf1" if "maf1" in cols else ("MAF1" if "MAF1" in cols else "maf1")
    ld_df = ld_df.filter(pl.col(maf1_col) >= maf_threshold)

    # rename (only those that exist)
    rename_map = {
        "chr": "CHR",

         'bp1': 'POS1',
        'bp2': 'POS2',
        'a1_1': 'A1_1',
        'a1_2': 'A1_2',
        'a2_1': 'A2_1',
        'a2_2': 'A2_2',
        "snp1": "rsID1",
        "snp2": "rsID2",
        "maf1": "MAF1",
        "maf2": "MAF2",
        "r2": "R2",
        "r": "R",
    }
    #keep only chr1 column
    rename_map = {k: v for k, v in rename_map.items() if k in cols}

    ld_df = ld_df.rename(rename_map)

    #   Filter for the rsID list BEFORE collecting
    ld_df = ld_df.filter(pl.col("rsID1").is_in(rs_list))

    # Convert "R" to a float before squaring it
    ld_df = ld_df.with_columns(
        (pl.col("R").cast(pl.Float64) ** 2).alias("R2")
    )

    if R2_threshold is not None:
        ld_df = ld_df.filter(pl.col("R2") >= R2_threshold)


    # Now collect
    ld_df = ld_df.collect()

    return ld_df


def UKBB_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    outputData = UKBB_LD_info(
        list(study_df['SNP']),
        chromosome,
        population,
        maf_input,
        r2threshold,
        imp_snp_list
    )
    return outputData



# -------------------------------------------------------------------
# 1000G hg38 (non-pairwise)
# -------------------------------------------------------------------

def hg38_1kg_LD_info(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    print(f"Loading 1000 Genomes Project (hg38) files ({population}) ...")

    # 1) Validate population
    valid_pops = ['AMR', 'EUR', 'AFR', 'SAS', 'EAS']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths
    maf_file = f'D:/ref/1000G_hg38/1000G_{population}_0_01.parquet'
    ld_file = f'D:/ref/1000G_hg38/{population}/chr{chrom}_merged.parquet'

    # 3) Load and filter MAF file early
    print("Loading and filtering MAF file...")
    tmp = pl.read_parquet(maf_file)
    maf_df = tmp.rename({
        old: new for old, new in zip(
            tmp.columns[:6],
            ['CHR', 'SNP', 'MAF', 'POS', 'REF', 'ALT']
        )
    })

    if maf_threshold is not None:
        maf_df = maf_df.filter(pl.col("MAF") >= maf_threshold)

    if imp_snp_list:
        maf_df = maf_df.filter(pl.col("SNP").is_in(imp_snp_list))
    elif rs_list:
        maf_df = maf_df.filter(pl.col("SNP").is_in(rs_list))

    maf_snps = maf_df.select("SNP").to_series().to_list()

    # 4) Load LD data
    print("Loading and filtering LD data...")
    ld_df = pl.read_parquet(ld_file, columns=["CHR_A", "SNP_A", "CHR_B", "SNP_B", "R"])

    if imp_snp_list or rs_list:
        filter_snps = set(imp_snp_list or rs_list)
        ld_df = ld_df.filter(
            pl.col("SNP_A").is_in(filter_snps) | pl.col("SNP_B").is_in(filter_snps)
        )

    ld_df = ld_df.with_columns((pl.col("R") ** 2).alias("R2"))

    if R2_threshold is not None:
        ld_df = ld_df.filter(pl.col("R2") >= R2_threshold)

    # 5) Join with allele A
    merged = ld_df.join(
        maf_df,
        left_on=["CHR_A", "SNP_A"],
        right_on=["CHR", "SNP"],
        how="left"
    ).rename({
        "CHR_A": "CHR",
        "MAF": "MAF_A",
        "POS": "POS_A",
        "REF": "REF_A",
        "ALT": "ALT_A"
    })

    # 6) Join with allele B
    merged = merged.join(
        maf_df,
        left_on=["CHR_B", "SNP_B"],
        right_on=["CHR", "SNP"],
        how="left"
    ).rename({
        "MAF": "MAF_B",
        "POS": "POS_B",
        "REF": "REF_B",
        "ALT": "ALT_B"
    }).drop(["CHR_B"])

    # 7) Final MAF filtering
    if maf_threshold is not None:
        merged = merged.filter(
            (pl.col("MAF_A") >= maf_threshold) & (pl.col("MAF_B") >= maf_threshold)
        )

    # 8) Reorder columns
    preferred_order = [
        "SNP_A", "CHR", "POS_A", "SNP_B", "POS_B",
        "REF_A", "ALT_A", "MAF_A", "REF_B", "ALT_B", "MAF_B",
        "R", "R2"
    ]
    existing_columns = merged.columns
    final_order = [col for col in preferred_order if col in existing_columns] + [
        col for col in existing_columns if col not in preferred_order
    ]
    merged = merged.select(final_order)

    return merged


def hg38_1kg_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    outputData = hg38_1kg_LD_info(
        list(study_df['SNP']),
        chromosome,
        population,
        maf_input,
        r2threshold,
        imp_snp_list
    )
    return outputData


# -------------------------------------------------------------------
# HapMap (non-pairwise)
# -------------------------------------------------------------------

def Hap_Map_LD_info_dask(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    print(f"Loading Hap Map files ({population}) ...")

    # 1) Validate population
    valid_pops = ['YRI', 'CHB', 'JPT', 'CEU', 'MKK', 'LWK',
                  'CHD', 'GIH', 'TSI', 'MEX', 'ASW']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths
    maf_file = f'D:/ref/Hap_Map/allele_freqs_chr{chrom}_{population}_phase3.2_nr.b36_fwd.parquet'
    ld_file = f'D:/ref/Hap_Map/ld_chr{chrom}_{population}.parquet'

    # 3) Read and preprocess MAF table
    maf_df = (
        pl.read_parquet(maf_file)
        .with_columns(pl.col("otherallele_freq").alias("MAF"))
        .select(["rs#", "refallele", "otherallele", "MAF"])
        .rename({
            "rs#": "rsID",
            "refallele": "REF",
            "otherallele": "ALT"
        })
        .filter(pl.col("MAF") >= float(maf_threshold))
    )

    # 4) Read and preprocess LD table
    ld_raw = pl.read_parquet(ld_file)
    cols = ld_raw.columns

    ld_df = (
        ld_raw
        .filter(
            (pl.col(cols[3]) != pl.col(cols[4])) &
            (pl.col(cols[6]) >= float(R2_threshold))
        )
        .rename({
            cols[0]: "pos1",
            cols[1]: "pos2",
            cols[2]: "pop",
            cols[3]: "rsID1",
            cols[4]: "rsID2",
            cols[5]: "Dprime",
            cols[6]: "R2"
        })
        .select(["pos1", "pos2", "pop", "rsID1", "rsID2", "Dprime", "R2"])
    )

    # 5) Join LD ↔ MAF on rsID1 and rsID2
    tmp = ld_df.join(maf_df, left_on="rsID1", right_on="rsID", how="inner")
    merged = tmp.join(
        maf_df,
        left_on="rsID2",
        right_on="rsID",
        how="inner",
        suffix="_2"
    )

    result = merged.select([
        pl.col("pos1"),
        pl.col("pos2"),
        pl.col("rsID1"),
        pl.col("rsID2"),
        pl.col("MAF").alias("MAF1"),
        pl.col("MAF_2").alias("MAF2"),
        pl.col("REF").alias("REF1"),
        pl.col("REF_2").alias("REF2"),
        pl.col("ALT").alias("ALT1"),
        pl.col("ALT_2").alias("ALT2"),
        pl.col("R2"),
        pl.col("Dprime")
    ])

    if imp_snp_list:
        result = result.filter(
            pl.col("rsID1").is_in(rs_list) &
            pl.col("rsID2").is_in(imp_snp_list)
        )
    else:
        result = result.filter(pl.col("rsID1").is_in(rs_list))


    result = result.with_columns(
        pl.lit(chrom).alias("CHR")
    )
    del ld_raw, ld_df, maf_df, tmp, merged
    gc.collect()

    return result


def Hap_Map_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    outputData = Hap_Map_LD_info_dask(
        list(study_df['SNP']),
        chromosome,
        population,
        maf_input,
        r2threshold,
        imp_snp_list
    )
    return outputData


# -------------------------------------------------------------------
# PhenoScanner (non-pairwise)
# -------------------------------------------------------------------

def pheno_Scanner_LD_info_dask(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    if R2_threshold < 0.8:
        print("Pheno Scanner includes data with R2 ≥ 0.8. Setting R2_threshold = 0.8")
        R2_threshold = 0.8

    print("Building lazy plan for Pheno Scanner files...")

    maf_file = "D:/ref/Pheno_Scanner/1000G.parquet"
    ld_file = f"D:/ref/Pheno_Scanner/1000G_{population}/1000G_{population}_chr{chrom}.parquet"
    pop_map = {"EUR": "eur", "EAS": "eas", "AFR": "afr", "AMR": "amr", "SAS": "sas"}
    maf_pop = pop_map.get(population)
    if maf_pop is None:
        raise ValueError(f"Unsupported population: {population}")

    maf1_lazy = (
        pl.scan_parquet(maf_file)
        .select(["hg19_coordinates", "chr", "rsid", maf_pop, "a1", "a2"])
        .filter(
            (pl.col("chr") == chrom) &
            (pl.col(maf_pop) != "-")
        )
        .with_columns(pl.col(maf_pop).cast(pl.Float64))
        .filter(pl.col(maf_pop) >= maf_threshold)
        .rename({
            "hg19_coordinates": "ref_hg19_coordinates",
            "rsid": "ref_rsid",
            maf_pop: "MAF1",
            "a1": "ALT1",
            "a2": "REF1",
        })
        .select(["ref_hg19_coordinates", "ref_rsid", "MAF1", "ALT1", "REF1"])
    )

    ld_lazy = (
        pl.scan_parquet(ld_file)
        .select(["ref_hg19_coordinates", "ref_rsid", "rsid", "r2", "r", "dprime"])
        .filter(
            (pl.col("ref_rsid") != pl.col("rsid")) &
            (pl.col("r2") >= R2_threshold)
        )
    )

    maf2_lazy = (
        pl.scan_parquet(maf_file)
        .select(["hg19_coordinates", "chr", "rsid", maf_pop, "a1", "a2"])
        .filter(
            (pl.col("chr") == chrom) &
            (pl.col(maf_pop) != "-")
        )
        .with_columns(pl.col(maf_pop).cast(pl.Float64))
        .filter(pl.col(maf_pop) >= maf_threshold)
        .rename({
            "hg19_coordinates": "pos2(hg19)",
            maf_pop: "MAF2",
            "a1": "ALT2",
            "a2": "REF2",
        })
        .select(["pos2(hg19)", "rsid", "MAF2", "ALT2", "REF2"])
    )

    joined = (
        ld_lazy
        .join(maf1_lazy, on="ref_rsid", how="inner")
        .join(maf2_lazy, on="rsid", how="inner")
    )

    renamed = joined.rename({
        "ref_rsid": "rsID1",
        "rsid": "rsID2",
        "ref_hg19_coordinates": "pos1(hg19)",
        "r2": "R2",
    })

    filtered = renamed.filter(
        pl.col("rsID2").is_in(rs_list) &
        (pl.col("rsID1").is_in(imp_snp_list) if imp_snp_list else pl.lit(True))
    )

    with_pos = filtered.with_columns([
        pl.col("pos1(hg19)").str.extract(r".*:(\d+)$", 1).cast(pl.Int64),
        pl.col("pos2(hg19)").str.extract(r".*:(\d+)$", 1).cast(pl.Int64),
    ])

    final_lazy = with_pos.select([
        "rsID1", "pos1(hg19)", "rsID2", "dprime",
        "pos2(hg19)", "R2", "r",
        "MAF1", "MAF2", "ALT1", "REF1", "ALT2", "REF2",
    ])

    final_pl = final_lazy.collect()
    final_pd = final_pl.to_pandas()

    gc.collect()
    if final_pd.empty:
        print("No SNPs found")

    final_pd["CHR"] = chrom


    return final_pd


def pheno_Scanner_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    outputData = pheno_Scanner_LD_info_dask(
        list(study_df['SNP']),
        chromosome,
        population,
        maf_input,
        r2threshold,
        imp_snp_list
    )
    return outputData


# -------------------------------------------------------------------
# TOP-LD (non-pairwise)
# -------------------------------------------------------------------

def TOP_LD_info(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list=None):
    """
    Optimized version of TOP_LD_info to minimize memory footprint.
    Uses semi-joins for early filtering of LD pairs.
    """
    print(f"Loading TOP-LD files for {population} chr{chrom}...")

    # 1) Combine rsIDs up-front to define the universe of interesting SNPs
    target_rsids = set(rs_list)
    if imp_snp_list:
        all_rsids = list(target_rsids | set(imp_snp_list))
    else:
        all_rsids = list(target_rsids)

    maf_path = (
        f"D:/ref/TOP_LD/{population}/SNV/"
        f"{population}_chr{chrom}_no_filter_0.2_1000000_info_annotation.parquet"
    )
    ld_path = (
        f"D:/ref/TOP_LD/{population}/SNV/"
        f"{population}_chr{chrom}_no_filter_0.2_1000000_LD.parquet"
    )

    # 2) Scan MAF and filter immediately
    # We create a lazy frame of ONLY the SNPs we care about.
    maf_lazy = (
        pl.scan_parquet(maf_path)
        .filter((pl.col("MAF") >= maf_threshold) & (pl.col("rsID").is_in(all_rsids)))
        .select(["Uniq_ID", "rsID", "MAF", "REF", "ALT"])
    )

    # We need a lightweight version (just IDs) for filtering the LD table efficiently
    valid_ids = maf_lazy.select("Uniq_ID")

    # 3) Optimize LD loading: Use semi-joins to drop irrelevant rows early
    # "semi" join acts like a filter: keep rows in LEFT where key exists in RIGHT.
    # We filter LD rows where *both* Uniq_ID_1 and Uniq_ID_2 exist in our relevant MAF list.
    ld_lazy = (
        pl.scan_parquet(ld_path)
        .filter(pl.col("R2") >= R2_threshold)
        .join(valid_ids, left_on="Uniq_ID_1", right_on="Uniq_ID", how="semi")
        .join(valid_ids, left_on="Uniq_ID_2", right_on="Uniq_ID", how="semi")
    )

    # 4) Prepare two MAF views for the final enrichment
    # Since we've already used semi-joins to filter LD, these inner joins will only
    # attach data to the rows that survived, rather than exploding the search space.
    maf1 = maf_lazy.rename({
        "Uniq_ID": "Uniq_ID_1", "rsID": "rsID1", "MAF": "MAF1",
        "REF": "REF1", "ALT": "ALT1"
    })
    maf2 = maf_lazy.rename({
        "Uniq_ID": "Uniq_ID_2", "rsID": "rsID2", "MAF": "MAF2",
        "REF": "REF2", "ALT": "ALT2"
    })

    # 5) Join LD ↔ MAF1 ↔ MAF2
    joined = (
        ld_lazy
        .join(maf1, on="Uniq_ID_1", how="inner")
        .join(maf2, on="Uniq_ID_2", how="inner")
    )

    # 6) Apply specific directional filtering if imp_snp_list is provided
    # The semi-joins ensured ID1 and ID2 are in the "all_rsids" bucket.
    # This step ensures the specific relationship (rsID1 in rs_list AND rsID2 in imp_snp_list).
    if imp_snp_list:
        joined = joined.filter(
            (pl.col("rsID1").is_in(list(target_rsids))) &
            (pl.col("rsID2").is_in(imp_snp_list))
        )

    # 7) Select + rename final output columns
    final_lazy = joined.select([
        pl.col("Uniq_ID_1").alias("pos1"),
        pl.col("Uniq_ID_2").alias("pos2"),
        "R2", "+/-corr", "Dprime",
        "rsID1", "rsID2",
        "MAF1", "MAF2",
        "REF1", "ALT1", "REF2", "ALT2"
    ])

    # 8) Execute without streaming
    # Note: We removed streaming=True because it is deprecated for Parquet in newer Polars versions.
    # The semi-join optimizations above sufficiently reduce memory usage for the standard engine.
    print("Executing query...")
    result = final_lazy.collect()

    if result.is_empty():
        print("No SNPs found.")

    # 9) Cleanup
    del maf_lazy, valid_ids, maf1, maf2, ld_lazy, joined, final_lazy
    gc.collect()

    result = result.with_columns(
        pl.lit(chrom).alias("CHR")
    )
    return result

def TOP_LD_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    outputData = TOP_LD_info(
        list(study_df['SNP']),
        chromosome,
        population,
        maf_input,
        r2threshold,
        imp_snp_list
    )
    return outputData






def _normalize_chr_values(values):
    """
    Convert a list of chromosome-like values to a sorted list of ints in [1..22].
    Accepts: 1, "1", "chr1", "CHR1", etc. Ignores X/Y/MT and anything outside 1..22.
    """
    out = set()
    for v in values:
        if v is None:
            continue
        s = str(v).strip()

        # Extract first integer occurrence (handles "chr1", "CHR:1", etc.)
        m = re.search(r"(\d+)", s)
        if not m:
            continue
        try:
            c = int(m.group(1))
        except ValueError:
            continue
        if 1 <= c <= 22:
            out.add(c)

    return sorted(out)

def _concat_frames(frames):
    """
    Concatenate a list of DataFrames (all Pandas or all Polars).
    Returns a DataFrame of the same type, or None if empty.
    """
    frames = [f for f in frames if f is not None]
    if not frames:
        return None

    first = frames[0]
    if isinstance(first, pl.DataFrame):
        # Ensure all are Polars
        pl_frames = []
        for f in frames:
            if isinstance(f, pl.DataFrame):
                pl_frames.append(f)
            elif isinstance(f, pd.DataFrame):
                pl_frames.append(pl.from_pandas(f))
            else:
                raise TypeError(f"Unsupported frame type: {type(f)}")
        return pl.concat(pl_frames, how="vertical", rechunk=True)

    if isinstance(first, pd.DataFrame):
        # Ensure all are Pandas
        pd_frames = []
        for f in frames:
            if isinstance(f, pd.DataFrame):
                pd_frames.append(f)
            elif isinstance(f, pl.DataFrame):
                pd_frames.append(f.to_pandas())
            else:
                raise TypeError(f"Unsupported frame type: {type(f)}")
        return pd.concat(pd_frames, axis=0, ignore_index=True)

    raise TypeError(f"Unsupported frame type: {type(first)}")

def _write_frame_as_tsv(df, path):
    """
    Write either Pandas or Polars DataFrame to TSV.
    """
    if df is None:
        # Write an empty file with no header (or you can raise)
        with open(path, "w", encoding="utf-8") as f:
            f.write("")
        return

    if isinstance(df, pl.DataFrame):
        df.write_csv(path, separator="\t")
        return

    if isinstance(df, pd.DataFrame):
        df.to_csv(path, sep="\t", index=False)
        return

    raise TypeError(f"Unsupported frame type: {type(df)}")


def process_data(file_path, r2threshold, population, maf_input, ref_file, imp_snp_list):
    """
    Non-pairwise LD extraction across chromosomes.

    Behavior:
    - If CHR column exists: process the chromosomes present (filtered to 1..22).
    - If CHR column does NOT exist: process chromosomes 1..22.
    - Merge all chromosome results at the end into LD_info_chr_all.txt
    """
    final_results_list = []

    # 1) Read with Polars (fast for large TSVs)
    study_pl = pl.read_csv(file_path, separator="\t")
    ##print(study_pl)

    # 2) Decide chromosomes
    if "CHR" in study_pl.columns:
        chroms_raw = study_pl.select(pl.col("CHR").unique()).to_series().to_list()
        chroms = _normalize_chr_values(chroms_raw)

        # If CHR exists but parsing yielded nothing, fall back to 1..22
        if not chroms:
            chroms = list(range(1, 23))
    else:
        chroms = list(range(1, 23))

    # 3) Convert to Pandas once (your downstream functions expect Pandas)
    study_df = study_pl.to_pandas()

    ref_panel = ref_file

    if ref_panel == "1000G_hg38":
        for chrom in chroms:
            final_data = hg38_1kg_process(
                study_df, r2threshold, population, maf_input, chrom, imp_snp_list
            )
            final_results_list.append(final_data)

        final_df = _concat_frames(final_results_list)
        _write_frame_as_tsv(final_df, "LD_info_chr_all.txt")
        print("Check 'LD_info_chr_all.txt' for results")

    elif ref_panel == "1000G_hg38_high_cov":
        for chrom in chroms:
            final_data = hg38_1kg_process_high_coverage(
                study_df, r2threshold, population, maf_input, chrom, imp_snp_list
            )
            final_results_list.append(final_data)

        final_df = _concat_frames(final_results_list)
        _write_frame_as_tsv(final_df, "LD_info_chr_all.txt")
        print("Check 'LD_info_chr_all.txt' for results")
    elif ref_panel == "UKBB":
        for chrom in chroms:
            final_data = UKBB_process(
                study_df, r2threshold, population, maf_input, chrom, imp_snp_list
            )
            final_results_list.append(final_data)

        final_df = _concat_frames(final_results_list)
        _write_frame_as_tsv(final_df, "LD_info_chr_all.txt")
        print("Check 'LD_info_chr_all.txt' for results")
    elif ref_panel == "HGDP":
        for chrom in chroms:
            final_data = HGDP_process(
                study_df, r2threshold, population, maf_input, chrom, imp_snp_list
            )
            final_results_list.append(final_data)

        final_df = _concat_frames(final_results_list)
        _write_frame_as_tsv(final_df, "LD_info_chr_all.txt")
        print("Check 'LD_info_chr_all.txt' for results")

    elif ref_panel == "TOP_LD":
        for chrom in chroms:
            final_data = TOP_LD_process(
                study_df, r2threshold, population, maf_input, chrom, imp_snp_list
            )
            final_results_list.append(final_data)

        final_df = _concat_frames(final_results_list)
        _write_frame_as_tsv(final_df, "LD_info_chr_all.txt")
        print("Check 'LD_info_chr_all.txt' for results")

    elif ref_panel == "Pheno_Scanner":
        for chrom in chroms:
            final_data = pheno_Scanner_process(
                study_df, r2threshold, population, maf_input, chrom, imp_snp_list
            )

            # per-chr file (optional, but you had it before; keep it)
            if isinstance(final_data, pl.DataFrame):
                final_data.write_csv(f"LD_info_chr{chrom}.txt", separator="\t")
            elif isinstance(final_data, pd.DataFrame):
                final_data.to_csv(f"LD_info_chr{chrom}.txt", sep="\t", index=False)
            else:
                raise TypeError(f"Unsupported frame type from pheno_Scanner_process: {type(final_data)}")

            print(f"Check 'LD_info_chr{chrom}.txt' for LD information")
            final_results_list.append(final_data)

        final_df = _concat_frames(final_results_list)
        _write_frame_as_tsv(final_df, "LD_info_chr_all.txt")
        print("Check 'LD_info_chr_all.txt' for merged results")

    elif ref_panel == "Hap_Map":
        for chrom in chroms:
            final_data = Hap_Map_process(
                study_df, r2threshold, population, maf_input, chrom, imp_snp_list
            )
            final_results_list.append(final_data)

        final_df = _concat_frames(final_results_list)
        _write_frame_as_tsv(final_df, "LD_info_chr_all.txt")
        print("Check 'LD_info_chr_all.txt' for results")

    else:
        raise ValueError(f"Unsupported ref_panel: {ref_panel}")

def _normalize_autosomal_chroms(values) -> list[int]:
    """
    Convert a list of mixed chromosome encodings into sorted autosomal ints [1..22].

    Accepts:
      - 1, "1", "01"
      - "chr1", "CHR1"
      - ignores X/Y/MT and anything not 1..22
    """
    out = set()
    for v in values:
        if v is None:
            continue

        # Polars may give numpy types; cast to str safely
        s = str(v).strip()

        # Strip common 'chr' prefix
        s = re.sub(r"^(chr|CHR|Chr)", "", s).strip()

        # Keep only digits if it's like "1.0" or "01"
        # (If you have weird encodings, adapt this part)
        if re.fullmatch(r"\d+(\.0+)?", s):
            s = s.split(".")[0]

        if s.isdigit():
            c = int(s)
            if 1 <= c <= 22:
                out.add(c)

    return sorted(out)


def process_data_pairwise(
    file_path,
    r2threshold,
    population,
    maf_input,
    ref_file,
    imp_snp_list,
    ld_prune: bool = False,
    ld_prune_p_col: str | None = "P",
    ld_prune_out_prefix: str = "LD_pruned",
    ld_prune_threshold: float | None = None,
    ld_prune_mode: str = "below"
):
    """
    Pairwise LD extraction across chromosomes.

    If CHR column exists in the input GWAS, only those chromosomes will be processed.
    Otherwise, defaults to 1..22.

    If ld_prune == True, also perform LD pruning on the input GWAS
    and return (final_ld_df, kept_gwas_df, pruned_gwas_df).
    """
    # -------------------------
    # Read GWAS once
    # -------------------------
    study_pl = pl.read_csv(file_path, separator="\t")
    #print(study_pl)

    # -------------------------
    # Decide which chromosomes to process
    # -------------------------
    if "CHR" in study_pl.columns:
        raw_chroms = (
            study_pl.select(pl.col("CHR").unique())
            .to_series()
            .to_list()
        )
        chroms = _normalize_autosomal_chroms(raw_chroms)

        # Fallback if CHR exists but is empty/badly encoded
        if not chroms:
            chroms = list(range(1, 23))
    else:
        chroms = list(range(1, 23))

    # Convert to pandas if your downstream functions expect pandas
    study_df = study_pl.to_pandas()

    ref_panel = ref_file
    final_results_list = []
    final_df = None
    kept_gwas = None
    pruned_gwas = None

    # -------------------------
    # Small refactor to avoid repetition
    # -------------------------
    panel_cfg = {
        "1000G_hg38": {
            "fn": hg38_1kg_process_pairwise,
            "snp1_col": "SNP_A",
            "snp2_col": "SNP_B",
        },

        "1000G_hg38_high_cov": {
            "fn": hg38_1kg_LD_info_high_coverage_pairwise,
            "snp1_col": "rsID1",
            "snp2_col": "rsID2",
        },

        "HGDP": {
            "fn": HGDP_process_pairwise,
            "snp1_col": "rsID1",
            "snp2_col": "rsID2",
        },

        "UKBB": {
            "fn": UKBB_process_pairwise,
            "snp1_col": "rsID1",
            "snp2_col": "rsID2",
        },

        "TOP_LD": {
            "fn": TOP_LD_process_pairwise,
            "snp1_col": "rsID1",
            "snp2_col": "rsID2",
        },
        "Hap_Map": {
            "fn": Hap_Map_process_pairwise,
            "snp1_col": "rsID1",
            "snp2_col": "rsID2",
        },
        "Pheno_Scanner": {
            "fn": pheno_Scanner_process_pairwise,
            "snp1_col": "rsID1",
            "snp2_col": "rsID2",
        },
    }

    if ref_panel not in panel_cfg:
        raise ValueError(f"Unsupported ref_panel: {ref_panel}")

    proc_fn = panel_cfg[ref_panel]["fn"]
    snp1_col = panel_cfg[ref_panel]["snp1_col"]
    snp2_col = panel_cfg[ref_panel]["snp2_col"]

    # -------------------------
    # Process only the needed chromosomes
    # -------------------------
    for chrom in chroms:
        final_data = proc_fn(study_df,r2threshold, population, maf_input,chrom,imp_snp_list)

        # Be tolerant if your proc functions sometimes return pandas
        if final_data is None:
            continue
        if isinstance(final_data, pl.DataFrame):
            final_results_list.append(final_data)
        else:
            # assume pandas
            final_results_list.append(pl.from_pandas(final_data))

    if not final_results_list:
        # Nothing came back (e.g., no SNPs matched / no LD pairs found)
        final_df = pl.DataFrame()
    else:
        final_df = pl.concat(final_results_list, how="vertical_relaxed")

    final_df.write_csv("LD_info_chr_all_pairwise.txt", separator="\t")
    print("Check 'LD_info_chr_all_pairwise.txt' for results")

    # -------------------------
    # Optional LD pruning
    # -------------------------
    if ld_prune and final_df is not None and final_df.height > 0:
        kept_gwas, pruned_gwas = ld_prune_pairs(
            study_df=study_df,
            ld_pairs=final_df,
            snp_col="SNP",
            p_col=ld_prune_p_col,
            snp1_col=snp1_col,
            snp2_col=snp2_col,
            threshold=ld_prune_threshold,
            threshold_mode=ld_prune_mode
        )

        kept_gwas.to_csv(f"{ld_prune_out_prefix}_kept.txt", sep="\t", index=False)
        pruned_gwas.to_csv(f"{ld_prune_out_prefix}_pruned.txt", sep="\t", index=False)
        print(
            f"LD-pruned GWAS written to '{ld_prune_out_prefix}_kept.txt' (kept) "
            f"and '{ld_prune_out_prefix}_pruned.txt' (pruned)."
        )

    return final_df, kept_gwas, pruned_gwas