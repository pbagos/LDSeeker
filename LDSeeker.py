#!/usr/bin/env python

import time
import os
import argparse
import pandas as pd
import LDSeeker_functions

start = time.time()


def main():

    # Parse arguments
    global imp_snp_list
    imp_snp_list = []
    version = '1.0.0'
    print("---------------------------------------------------------------------------------")
    print("LDSeeker : Search LD among a set of variants")
    print("Version " + version + "; December 2025")
    print("Copyright (C) 2025 Pantelis Bagos")
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
        help='LD Reference files (Pheno_Scanner, 1000G_hg38, TOP_LD, Hap_Map or all_panels)',
        default='all_panels'
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
            ld_prune_out_prefix=ld_prune_prefix
        )


if __name__ == "__main__":
    main()

end = time.time()
print(f"Total Time: {end - start} seconds")
