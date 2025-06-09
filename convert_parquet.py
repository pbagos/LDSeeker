
import pandas as pd
import os

# Input and output file paths
input_csv_gz_path = "EUR_chr22_no_filter_0.2_1000000_info_annotation.csv.gz"
output_parquet_path = "EUR_chr22_no_filter_0.2_1000000_info_annotation.parquet"

# Read the CSV.GZ file
print(f"Reading compressed CSV file: {input_csv_gz_path}")
df = pd.read_csv(input_csv_gz_path, compression='gzip')

# Optionally preview the data
print("DataFrame preview:")
print(df.head())

# Save to Parquet
print(f"Writing to Parquet file: {output_parquet_path}")
df.to_parquet(output_parquet_path, engine='pyarrow', index=False)

print("Conversion complete.")
