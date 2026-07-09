import time
import os
import argparse
import pandas as pd
import LDSeeker_functions

start = time.time()


def main():
    print(r"""
      .---.
     /     \
    |       |   _     ____   ____            _
     \     /   | |   |  _ \ / ___|  ___  ___| | _____ _ __
      `---'__  | |   | | | |\___ \ / _ \/ _ \ |/ / _ \ '__|
           \ \ | |___| |_| | ___) |  __/  __/   <  __/ |
            \ \|_____|____/ |____/ \___|\___|_|\_\___|_|
             `
    """)

    # Parse arguments
    global imp_snp_list
    imp_snp_list = []
    version = '1.0.0'
    print("---------------------------------------------------------------------------------")
    print("LDSeeker : Exploring linkage disequilibrium in GWAS using multiple panels")
    print("Version " + version + "; February 2026")
    print("Copyright (C) 2026 Pantelis Bagos")
    print("Freely distributed under the GNU General Public Licence (GPLv3)")
    print("---------------------------------------------------------------------------------")

    parser = argparse.ArgumentParser(description="Process data in chunks.")
    parser.add_argument('--file-path', type=str, required=True, help='Input file path')
    parser.add_argument('--r2threshold', type=float, required=True, help='R2 threshold')
    parser.add_argument('--pop', type=str, required=True, help='Population')

    parser.add_argument('--maf', type=float, required=True, help='MAF input value')
    parser.add_argument(
        '--ref',
        type=str,
        required=False,
        help='LD Reference files. Available reference panels: Pheno_Scanner, 1000G_hg38 (High Coverage, Standard), TOP_LD, Hap_Map, UKBB (UK Biobank), HGDP (Human Genome Diversity Project) and LASI-DAD',
        default='1000G_hg38'
    )
    parser.add_argument(
        '--pairwise',
        type=str,
        required=False,
        help='LD found pairwise (YES, NO)',
        default='NO'
    )

    parser.add_argument(
        '--imp_list',
        type=str,
        required=False,
        help='A filename to define SNPs to impute (each SNP has a new line, no header)'
    )

    # --- LD pruning related arguments (used only when --pairwise YES) ---
    prune_group = parser.add_argument_group(
        'LD pruning (pairwise mode only)',
        'Greedy pairwise-LD pruning of the input GWAS. Effective only when --pairwise YES.'
    )
    prune_group.add_argument(
        '--ld-prune',
        type=str,
        required=False,
        default='NO',
        help='Apply LD pruning to the input GWAS using pairwise LD (YES, NO). '
             'Effective only when --pairwise YES.'
    )
    prune_group.add_argument(
        '--ld-prune-metric',
        type=str,
        required=False,
        default='P',
        choices=['P', 'Z'],
        help='Which metric to prioritize when pruning. "P" keeps the lowest P-value, '
             '"Z" keeps the highest absolute Z-value. Default is "P".'
    )
    prune_group.add_argument(
        '--ld-prune-prefix',
        type=str,
        required=False,
        default='LD_pruned',
        help='Prefix for LD-pruned GWAS output files '
             '(e.g. LD_pruned_kept.txt).'
    )

    # --- P-Value Arguments (pruning) ---
    prune_group.add_argument(
        '--ld-prune-p-col',
        type=str,
        required=False,
        default='P',
        help='Column name in the input file used for P-values.'
    )
    prune_group.add_argument(
        '--ld-prune-p-threshold',
        type=float,
        required=False,
        default=None,
        help='Threshold value to filter SNPs by P-value before pruning (e.g. 0.05). '
             'Rows not meeting the criterion are ignored.'
    )
    prune_group.add_argument(
        '--ld-prune-p-threshold-mode',
        type=str,
        required=False,
        default='below',
        choices=['below', 'above'],
        help='Filter mode relative to P-value threshold: "below" (keep P < threshold) '
             'or "above" (keep P > threshold).'
    )

    # --- Z-Value Arguments (pruning) ---
    prune_group.add_argument(
        '--ld-prune-z-col',
        type=str,
        required=False,
        default='Z',
        help='Column name in the input file used for Z-values.'
    )
    prune_group.add_argument(
        '--ld-prune-z-threshold',
        type=float,
        required=False,
        default=2.0,
        help='Threshold value to filter SNPs by absolute Z-value before pruning '
             '(e.g. 2.0). Rows with |Z| <= threshold are ignored.'
    )

    # --- Significance filtering of the INPUT variants (applies to both modes) ---
    # Wald z = BETA/SE, two-tailed p = erfc(|z|/sqrt(2)). Keep either the
    # significant or the non-significant variants BEFORE LD retrieval. This is
    # independent of --pairwise and independent of the LD-pruning block above.
    sig_group = parser.add_argument_group(
        'Significance filtering (input variants)',
        'Keep only significant or only non-significant variants of the input file, '
        'judged on the Wald statistic (z = BETA/SE) and its two-tailed p-value. '
        'Applies to both standard and --pairwise YES runs.'
    )
    sig_group.add_argument(
        '--significance',
        type=str,
        required=False,
        default=None,
        choices=['significant', 'nonsignificant'],
        help='Which variants of the input file to keep before LD retrieval. '
             '"significant" keeps associated variants, "nonsignificant" keeps the '
             'complement. Omit to disable significance filtering (default).'
    )
    sig_group.add_argument(
        '--significance-metric',
        type=str,
        required=False,
        default='P',
        choices=['P', 'Z'],
        help='Metric used to define significance. "P" uses the two-tailed p-value '
             'derived from z (z = Z column, else BETA/SE); "Z" filters directly on '
             '|z|. Default is "P".'
    )
    sig_group.add_argument(
        '--significance-threshold',
        type=float,
        required=False,
        default=5e-8,
        help='Cutoff for --significance. With --significance-metric P this is a '
             'p-value (e.g. 5e-8 genome-wide, 0.05 nominal); with metric Z it is an '
             'absolute z cutoff (e.g. 5.4513 ~ p 5e-8, 1.96 ~ p 0.05). Default 5e-8.'
    )
    sig_group.add_argument(
        '--beta-col',
        type=str,
        required=False,
        default='BETA',
        help='Column name for the effect size (BETA), used to derive z = BETA/SE '
             'when no Z column is present. Default "BETA".'
    )
    sig_group.add_argument(
        '--se-col',
        type=str,
        required=False,
        default='SE',
        help='Column name for the standard error (SE), used to derive z = BETA/SE '
             'when no Z column is present. Default "SE".'
    )
    sig_group.add_argument(
        '--sig-z-col',
        type=str,
        required=False,
        default='Z',
        help='Column name for a precomputed z-score, preferred over BETA/SE when '
             'present, for significance filtering. Default "Z".'
    )
    sig_group.add_argument(
        '--sig-p-col',
        type=str,
        required=False,
        default='P',
        help='Column name for a precomputed p-value, used as a fallback for '
             'significance filtering only when neither Z nor BETA/SE is available. '
             'Default "P".'
    )

    args = parser.parse_args()

    # Accessing arguments
    file_path = args.file_path
    r2threshold = args.r2threshold
    population = args.pop
    maf_input = args.maf
    ref_file = args.ref
    imp_snp_list_path = args.imp_list
    pairwise = args.pairwise.upper()

    ld_prune_flag = args.ld_prune.upper()
    ld_prune = (ld_prune_flag == 'YES')

    # Significance filtering (shared by both entry points)
    sig_kwargs = dict(
        sig_keep=args.significance,
        sig_metric=args.significance_metric,
        sig_threshold=args.significance_threshold,
        z_col=args.sig_z_col,
        beta_col=args.beta_col,
        se_col=args.se_col,
        p_col=args.sig_p_col,
    )

    if imp_snp_list_path is not None:
        imp_snp_list = list(pd.read_csv(imp_snp_list_path, header=None)[0])

    if not os.path.exists(file_path):
        print(f"Error: File {file_path} not found.")
        return

    if pairwise == 'NO':
        # Standard case: no pairwise LD matrix, so no LD pruning here.
        # Significance filtering of the input variants still applies.
        LDSeeker_functions.process_data(
            file_path,
            r2threshold,
            population,
            maf_input,
            ref_file,
            imp_snp_list,
            **sig_kwargs
        )
    else:
        # Find pairwise LD among the given variants.
        # If ld_prune == True, LDSeeker_functions.process_data_pairwise()
        # is expected to perform LD pruning and write the pruned GWAS files.
        # Significance filtering of the input variants applies here too and,
        # when pruning is on, feeds the same kept set into the pruner.
        LDSeeker_functions.process_data_pairwise(
            file_path,
            r2threshold,
            population,
            maf_input,
            ref_file,
            imp_snp_list,
            ld_prune=ld_prune,
            ld_prune_p_col=args.ld_prune_p_col,
            ld_prune_p_threshold=args.ld_prune_p_threshold,
            ld_prune_p_threshold_mode=args.ld_prune_p_threshold_mode,
            ld_prune_z_col=args.ld_prune_z_col,
            ld_prune_z_threshold=args.ld_prune_z_threshold,
            ld_prune_metric=args.ld_prune_metric,
            ld_prune_out_prefix=args.ld_prune_prefix,
            **sig_kwargs
        )


if __name__ == "__main__":
    main()

end = time.time()
print(f"Total Time: {end - start} seconds")