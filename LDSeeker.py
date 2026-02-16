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
        help='LD Reference files. Available reference panels: Pheno_Scanner, 1000G_hg38 (High Coverage, Standard), TOP_LD, Hap_Map, UKBB (UK Biobank), HGDP (Human Genome Diversity Project), LASI-DAD and ChinaMAP',
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
    parser.add_argument(
        '--ld-prune',
        type=str,
        required=False,
        default='NO',
        help='Apply LD pruning to the input GWAS using pairwise LD (YES, NO). '
             'Effective only when --pairwise YES.'
    )
    parser.add_argument(
        '--ld-prune-col',
        type=str,
        required=False,
        default='P',
        help='Column name in the input file used to rank SNPs for LD pruning '
             '(e.g. P for p-values). If the column is missing, input order is used.'
    )
    parser.add_argument(
        '--ld-prune-prefix',
        type=str,
        required=False,
        default='LD_pruned',
        help='Prefix for LD-pruned GWAS output files '
             '(e.g. LD_pruned_kept.txt, LD_pruned_pruned.txt).'
    )
    parser.add_argument(
        '--ld-prune-threshold',
        type=float,
        required=False,
        default=None,
        help='Threshold value to filter SNPs before pruning (e.g. 0.05). '
             'Rows not meeting the criterion are ignored.'
    )
    parser.add_argument(
        '--ld-prune-mode',
        type=str,
        required=False,
        default='below',
        choices=['below', 'above'],
        help='Filter mode relative to threshold: "below" (keep P < threshold) '
             'or "above" (keep P > threshold).'
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
    ld_prune_col = args.ld_prune_col
    ld_prune_prefix = args.ld_prune_prefix

    # New arguments for threshold filtering
    ld_prune_threshold = args.ld_prune_threshold
    ld_prune_mode = args.ld_prune_mode

    if imp_snp_list_path is not None:
        imp_snp_list = list(pd.read_csv(imp_snp_list_path, header=None)[0])

    if not os.path.exists(file_path):
        print(f"Error: File {file_path} not found.")
        return

    if pairwise == 'NO':
        # Standard case: no pairwise LD matrix, so no LD pruning here
        LDSeeker_functions.process_data(
            file_path,
            r2threshold,
            population,
            maf_input,
            ref_file,
            imp_snp_list
        )
    else:
        # Find pairwise LD among the given variants
        # If ld_prune == True, LDSeeker_functions.process_data_pairwise()
        # is expected to perform LD pruning and write the pruned GWAS files.
        LDSeeker_functions.process_data_pairwise(
            file_path,
            r2threshold,
            population,
            maf_input,
            ref_file,
            imp_snp_list,
            ld_prune=ld_prune,
            ld_prune_p_col=ld_prune_col,
            ld_prune_out_prefix=ld_prune_prefix,
            ld_prune_threshold=ld_prune_threshold,
            ld_prune_mode=ld_prune_mode
        )


if __name__ == "__main__":
    main()

end = time.time()
print(f"Total Time: {end - start} seconds")