import dask.dataframe as dd
import numpy as np
import pandas as pd
import polars as pl
pd.set_option('display.width', None)


def Hap_Map_LD_info_dask(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    print(f"Loading Hap Map files ({population}) ...")

    if population not in ['YRI', 'CHB', 'JPT', 'CEU', 'MKK', 'LWK', 'CHD', 'GIH', "TSI", 'MEX', "ASW"]:
        print("This population is not available in HapMap files. Please select a different population...")
        exit()
    maf_file = f'ref/Hap_Map/allele_freqs_chr{chrom}_{population}_phase3.2_nr.b36_fwd.txt.gz'
    ld_file = f'ref/Hap_Map/ld_chr{chrom}_{population}.txt.gz'
    maf_df = dd.read_csv(maf_file, sep='\s+', blocksize=None)
    # Calculating the Minor Allele Frequency (MAF)
    # maf_df['MAF'] = maf_df[['refallele_freq', 'otherallele_freq']].min(axis=1)

    maf_df['MAF'] = maf_df['otherallele_freq']

    # Renaming columns in the DataFrame
    maf_df = maf_df.rename(columns={'refallele': 'REF', 'otherallele': 'ALT'})

    # Filter MAF DataFrame using Dask
    # maf_df = dd.read_csv(maf_file, sep='\s+',blocksize=None, usecols=['rs#', 'chrom', 'pos', 'otherallele_freq'],  )

    # maf_df['het-freq'] = maf_df['het-freq'].astype(float)
    maf_df = maf_df[maf_df['MAF'] >= float(maf_threshold)]

    # Process LD DataFrame using Dask
    ld_df = dd.read_csv(ld_file, blocksize=None, sep='\s+', header=None)

    ld_df = ld_df[(ld_df[3] != ld_df[4]) & (ld_df[6] >= R2_threshold)]

    maf_df = maf_df.rename(columns={'rs#': 'rsID'})

    # Define the new column names
    new_column_names = {
        0: 'pos1',
        1: 'pos2',
        2: 'pop',
        3: 'rsID1',
        4: 'rsID2',
        5: 'Dprime',
        6: 'R2'
    }

    # Rename the columns
    ld_df = ld_df.rename(columns=new_column_names)

    merged_df = dd.merge(ld_df, maf_df, left_on='rsID1', right_on='rsID')
    merged_df = dd.merge(merged_df, maf_df, left_on='rsID2', right_on='rsID')

    merged_df = merged_df[
        ['pos1', 'pos2', 'rsID1', 'rsID2', 'MAF_x', "MAF_y", 'REF_x', 'REF_y', 'ALT_x', 'ALT_y', "R2", "Dprime"]]
    merged_df = merged_df.rename(
        columns={'MAF_x': 'MAF1', 'MAF_y': 'MAF2', 'REF_x': 'REF1', 'REF_y': 'REF2', 'ALT_x': 'ALT1', 'ALT_y': 'ALT2'})

    if imp_snp_list:
        final_result = merged_df[merged_df['rsID1'].isin(rs_list) & merged_df['rsID2'].isin(imp_snp_list)]

    else:
        final_result = merged_df[merged_df['rsID1'].isin(rs_list)]
    final_result = final_result.compute()  # Important: This triggers the actual computation
    if final_result.empty:
        print("No SNPs found")
        exit()

    final_result.reset_index(inplace=True, drop=True)
    #final_result.to_csv('LD_info_Hap_Map_chr' + str(chrom) + '.txt', sep="\t", index=False)

    final_result.rename(columns={"MAF1": "ALT_AF1", "MAF2": "ALT_AF2"}).to_csv(
        'LD_info_Hap_Map_chr' + str(chrom) + '.txt', sep="\t", index=False
    )

    return final_result


def Hap_Map_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data
    outputData = Hap_Map_LD_info_dask(list(study_df['snp']), chromosome, population, maf_input, r2threshold,
                                      imp_snp_list)


    return outputData


def pheno_Scanner_LD_info_dask(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    if R2_threshold < 0.8:
        print("Pheno Scanner include data with a R2 threshold >= 0.8. The R2 threshold will be set to 0.8")
        R2_threshold = 0.8

    print("Loading Pheno Scanner files...")

    maf_file = 'ref/Pheno_Scanner/1000G.txt'
    ld_file = f'ref/Pheno_Scanner/1000G_{population}/1000G_{population}_chr{chrom}.txt.gz'
    population_map = {'EUR': 'eur', 'EAS': 'eas', 'AFR': 'afr', 'AMR': 'amr', 'SAS': 'sas'}

    maf_pop = population_map.get(population, None)
    if maf_pop is None:
        raise ValueError(f"Unsupported population: {population}")

    # Filter MAF DataFrame using Dask
    maf_df = dd.read_csv(maf_file, sep='\s+', blocksize=None,
                         usecols=['hg19_coordinates', 'chr', 'rsid', maf_pop, 'a1', 'a2'],
                         dtype={maf_pop: 'object'}
                         )
    maf_df = maf_df[(maf_df['chr'] == chrom) & (maf_df[maf_pop] != '-')]
    maf_df[maf_pop] = maf_df[maf_pop].astype(float)

    # # Calculate the frequency of the second allele for each population
    # maf_df[str(maf_pop)+'2'] = 1 -  maf_df[maf_pop]
    #
    # maf_df[maf_pop] = maf_df[[maf_pop, str(maf_pop)+'2']].min(axis=1)
    maf_df = maf_df[maf_df[maf_pop] >= float(maf_threshold)]

    # Process LD DataFrame using Dask
    ld_df = dd.read_csv(ld_file, sep='\s+', blocksize=None,
                        usecols=['ref_hg19_coordinates', 'ref_rsid', 'rsid', 'r2', 'r', 'dprime'],
                        dtype={'r2': 'float64', 'dprime': 'float64', 'r': 'float64'})

    ld_df = ld_df[(ld_df['ref_rsid'] != ld_df['rsid']) & (ld_df['r2'] >= R2_threshold)]
    merged_df = dd.merge(ld_df, maf_df.rename(
        columns={'hg19_coordinates': 'ref_hg19_coordinates', 'rsid': 'ref_rsid', maf_pop: 'MAF1', 'a1': 'ALT1',
                 'a2': 'REF1'}), on='ref_rsid')
    merged_df = dd.merge(merged_df, maf_df.rename(columns={maf_pop: 'MAF2', 'a1': 'ALT2', 'a2': 'REF2'}), on='rsid')

    # Convert to Pandas DataFrame by computing, to finalize and filter based on rs_list
    final_result = merged_df.compute()  # Important: This triggers the actual computation
    # print(final_result.head())
    final_result = final_result.rename(
        columns={'ref_rsid': 'rsID1', 'rsid': 'rsID2', 'ref_hg19_coordinates_x': 'pos1(hg19)',
                 'hg19_coordinates': 'pos2(hg19)', 'r2': 'R2'})
    final_result = final_result[
        ['rsID1', 'pos1(hg19)', 'rsID2', 'dprime', 'pos2(hg19)', 'R2', 'r', 'MAF1', 'MAF2', 'ALT1', 'REF1', 'ALT2',
         'REF2']]

    if imp_snp_list:
        final_result = final_result[final_result['rsID2'].isin(rs_list) & final_result['rsID1'].isin(imp_snp_list)]

    else:
        final_result = final_result[final_result['rsID2'].isin(rs_list)]
    if final_result.empty:
        print("No SNPs found")
        exit()
    # Split the 'location' column at ':' and keep the part after it
    final_result['pos1(hg19)'] = final_result['pos1(hg19)'].str.split(':').str[1]
    final_result['pos2(hg19)'] = final_result['pos2(hg19)'].str.split(':').str[1]
    final_result.reset_index(inplace=True, drop=True)
    # final_result.to_csv('LD_info_chr' + str(chrom) + '.txt', sep="\t", index=False)
   # final_result.to_csv('LD_info_Pheno_Scanner_chr_' + str(chrom) + '.txt', sep="\t", index=False)
    final_result.rename(columns={"MAF1": "ALT_AF1", "MAF2": "ALT_AF2"}).to_csv(
        'LD_info_Pheno_Scanner_chr_' + str(chrom) + '.txt', sep="\t", index=False
    )

    return final_result


def pheno_Scanner_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data
    outputData = pheno_Scanner_LD_info_dask(list(study_df['snp']), chromosome, population, maf_input, r2threshold,
                                            imp_snp_list)

    return outputData


def TOP_LD_info(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    # print("Loading TOP-LD files from Parquet...")
    #
    # # Updated paths to Parquet format
    # maf_file = f'{population}_chr{chrom}_no_filter_0.2_1000000_info_annotation.parquet'
    # ld_file = f'{population}_chr{chrom}_no_filter_0.2_1000000_LD.parquet'
    #
    # # Load MAF DataFrame from Parquet
    # maf_df = dd.read_parquet(maf_file, columns=['Position', 'rsID', 'MAF', 'REF', 'ALT'])
    # maf_df = maf_df[maf_df['MAF'] >= maf_threshold]
    #
    # # Load LD DataFrame from Parquet
    # ld_df = dd.read_parquet(ld_file, columns=['SNP1', 'SNP2', 'R2', '+/-corr', 'Dprime'])
    # ld_df = ld_df[ld_df['R2'] >= R2_threshold]
    #
    # # Rename columns for consistency
    # maf_df = maf_df.rename(columns={'Position': 'SNP'})
    #
    # # Merge LD with MAF data
    # merged_df = dd.merge(ld_df, maf_df.rename(
    #     columns={'SNP': 'SNP1', 'rsID': 'rsID1', 'MAF': 'MAF1', 'REF': 'REF1', 'ALT': 'ALT1'}), on='SNP1')
    # merged_df = dd.merge(merged_df, maf_df.rename(
    #     columns={'SNP': 'SNP2', 'rsID': 'rsID2', 'MAF': 'MAF2', 'REF': 'REF2', 'ALT': 'ALT2'}), on='SNP2')
    #
    # # Final selected and renamed columns
    # final_df = merged_df[
    #     ['SNP1', 'SNP2', 'R2', '+/-corr', 'Dprime', 'rsID1', 'rsID2', 'MAF1', 'MAF2', 'REF1', 'ALT1', 'REF2', 'ALT2']]
    # final_df = final_df.rename(columns={'SNP1': 'pos1', 'SNP2': 'pos2'})
    #
    # # Filter based on SNP list
    # if imp_snp_list:
    #     result = final_df[final_df['rsID1'].isin(rs_list) & final_df['rsID2'].isin(imp_snp_list)].compute()
    # else:
    #     result = final_df[final_df['rsID1'].isin(rs_list)].compute()
    #
    # if result.empty:
    #     print("No SNPs found")
    #     exit()
    #
    # result.reset_index(inplace=True, drop=True)
    # result.rename(columns={"MAF1": "ALT_AF1", "MAF2": "ALT_AF2"}).to_csv(
    #     f'LD_info_TOP_LD_chr{chrom}.txt', sep="\t", index=False
    # )
    #
    # return result
    # Filepaths
    maf_file = f"{population}_chr{chrom}_no_filter_0.2_1000000_info_annotation.parquet"
    ld_file = f"{population}_chr{chrom}_no_filter_0.2_1000000_LD.parquet"

    # 1) Lazy load and filter MAF table
    maf_lazy = (
        pl.scan_parquet(maf_file)
        .select([
            pl.col("Position").alias("SNP"),
            pl.col("rsID"),
            pl.col("MAF"),
            pl.col("REF"),
            pl.col("ALT"),
        ])
        .filter(pl.col("MAF") >= maf_threshold)
    )

    # 2) Lazy load and filter LD table
    ld_lazy = (
        pl.scan_parquet(ld_file)
        .select([
            pl.col("SNP1"),
            pl.col("SNP2"),
            pl.col("R2"),
            pl.col("+/-corr"),
            pl.col("Dprime"),
        ])
        .filter(pl.col("R2") >= R2_threshold)
    )

    # 3) Join LD to MAF on SNP1
    df = ld_lazy.join(
        maf_lazy.rename({
            "SNP": "SNP1",
            "rsID": "rsID1",
            "MAF": "MAF1",
            "REF": "REF1",
            "ALT": "ALT1",
        }),
        on="SNP1",
    )

    # 4) Join the result to MAF on SNP2
    df = df.join(
        maf_lazy.rename({
            "SNP": "SNP2",
            "rsID": "rsID2",
            "MAF": "MAF2",
            "REF": "REF2",
            "ALT": "ALT2",
        }),
        on="SNP2",
    )

    # 5) Filter by your SNP lists
    if imp_snp_list:
        df = df.filter(
            (pl.col("rsID1").is_in(rs_list)) &
            (pl.col("rsID2").is_in(imp_snp_list))
        )
    else:
        df = df.filter(pl.col("rsID1").is_in(rs_list))

    # 6) Collect into memory
    result = df.collect()

    # 7) Check for emptiness
    if result.height == 0:
        print("No SNPs found matching your criteria.", file=sys.stderr)
        sys.exit(1)

    # 8) Rename for output
    result = result.rename({
        "SNP1": "pos1",
        "SNP2": "pos2",
        "MAF1": "ALT_AF1",
        "MAF2": "ALT_AF2",
    })
    print(result)
    # 9) Write to TSV
    out_fname = f"LD_info_TOP_LD_chr{chrom}.txt"
    result.write_csv(out_fname, separator="\t")
    print(f"Wrote LD info to {out_fname}")

    return result

def TOP_LD_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data

    outputData = TOP_LD_info(list(study_df['snp']), chromosome, population, maf_input, r2threshold, imp_snp_list)

    return outputData


def process_data(file_path, r2threshold, population, maf_input, ref_file, imp_snp_list):
    final_results_list = []

    study_df = pd.read_csv(file_path, sep="\t")
    chroms = list(set(study_df['chr']))
    ref_panel = ref_file
    # Check if all required columns are present



    # Depending on the reference panel...
    if ref_panel == 'TOP_LD':
        for chrom in chroms:
            final_data = TOP_LD_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)


            # final_data.to_csv("imputation_results_chr" + str(chrom) + ".txt", sep="\t", index=False)
            print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            # final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)

            final_df.to_csv("LD_info_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")

    if ref_panel == 'Pheno_Scanner':
        for chrom in chroms:
            final_data = pheno_Scanner_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)



            final_data.to_csv("LD_info_chr" + str(chrom) + ".txt", sep="\t", index=False)

            print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")

            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)


            final_df.to_csv("imputation_results_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")

    if ref_panel == 'Hap_Map':
        for chrom in chroms:
            final_data = Hap_Map_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)


            final_data.to_csv("LD_info_chr" + str(chrom) + ".txt", sep="\t", index=False)

            print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)

            final_df.to_csv("LD_info_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")

    if ref_panel == 'all_panels':

        print(f"Checking all LD sources")
        for chrom in chroms:
            # For HapMap, we need to take all the panels and merge them...
            if population == 'EUR':
                pop_hm = "CEU"
                final_data_hm = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

                # pop_hm = "TSI"
                # final_data_hm_TSI = Hap_Map_process(study_df, r2threshold, pop_hm , maf_input, chrom, imp_snp_list)
                #

                # final_data_hm = pd.concat([final_data_hm_CEU, final_data_hm_TSI])
                # Keep the largest R2 value if a snp is common in any of the panels
                final_data_hm = final_data_hm.groupby('snp').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(
                    drop=True)

            if population == 'AFR':
                pop_hm = "YRI"
                final_data_hm_YRI = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

                pop_hm = "MKK"
                final_data_hm_MKK = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

                pop_hm = "LWK"
                final_data_hm_LWK = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

                pop_hm = "ASW"
                final_data_hm_ASW = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

                final_data_hm = pd.concat([final_data_hm_YRI, final_data_hm_MKK, final_data_hm_LWK, final_data_hm_ASW])

            if population == 'EAS':
                pop_hm = "CHB"
                final_data_hm_CHB = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)
                pop_hm = "JPT"
                final_data_hm_JPT = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)
                pop_hm = "CHD"
                final_data_hm_CHD = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

                final_data_hm = pd.concat([final_data_hm_CHB, final_data_hm_JPT, final_data_hm_CHD])

            if population == 'SAS':
                pop_hm = 'GIH'
                final_data_hm = Hap_Map_process(study_df, r2threshold, pop_hm, maf_input, chrom, imp_snp_list)

            # Keep the largest R2 value if a snp is common in any of the panels
            final_data_hm = final_data_hm.groupby('snp').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(drop=True)
            final_data_hm['source'] = 'HapMap'

            # final_data_hm = Hap_Map_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            final_data_ps = pheno_Scanner_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            final_data_tld = TOP_LD_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            final_data_ps['source'] = 'Pheno Scanner'
            final_data_tld['source'] = 'TOP-LD'

            # final_data = pd.concat([final_data_hm, final_data_ps, final_data_tld])
            final_data = pd.concat([final_data_hm, final_data_ps, final_data_tld], ignore_index=True)

            # Keep the largest R2 value if a snp is common in any of the panels
            if imp_snp_list == True:
                final_data = final_data.loc[final_data.groupby('snp')['R2'].idxmax()]
            else:
                final_data = final_data.groupby('snp').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(drop=True)
            print(len(final_data))

            data = pd.read_csv(file_path, sep="\t")
            data['z'] = data['beta'] / data['SE']
            data['imputed'] = 0
            data['source'] = 'GWAS'
            final_data = pd.concat([final_data, data], ignore_index=True)

            print(f"Total Imputed SNPs: {len(final_data[final_data['imputed'] == 1])} SNPs")

            print(f"Total : {len(final_data)} SNPs")

            final_data.to_csv("imputation_results_chr" + str(chrom) + ".txt", sep="\t", index=False)

            print("Check 'imputation_results_chr" + str(chrom) + ".txt' for the results")
            print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)
            # Separate the DataFrame into two based on the 'imputed' column.

            # Remove duplicates in the 'final_df_init' DataFrame based on the 'snp' column.
            final_df = final_df.drop_duplicates(subset="snp")

            # Concatenate the two DataFrames back together. You might consider resetting the index.

            final_df.to_csv("LD_info_chr_all.txt", sep="\t", index=False)
            print("Check 'LD_info_chr_all.txt' for results")