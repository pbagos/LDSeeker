# **LDSeeker**

  .---.  
 /     \\  
|       |   \_     \_\_\_\_   \_\_\_\_            \_  
 \\     /   | |   |  \_ \\ / \_\_\_|  \_\_\_  \_\_\_| | \_\_\_\_\_ \_ \_\_  
  \`---'\_\_  | |   | | | |\\\_\_\_ \\ / \_ \\/ \_ \\ |/ / \_ \\ '\_\_|  
       \\ \\ | |\_\_\_| |\_| | \_\_\_) |  \_\_/  \_\_/   \<  \_\_/ |  
        \\ \\|\_\_\_\_\_|\_\_\_\_/ |\_\_\_\_/ \\\_\_\_|\\\_\_\_|\_|\\\_\\\_\_\_|\_|  
         \`

**LDSeeker** is a Python-based tool for exploring linkage disequilibrium (LD) in Genome-Wide Association Studies (GWAS) using multiple large-scale reference panels. It leverages **Polars** for high-performance data processing.

**Web Tool:** A web-based version of LDSeeker is available at: [https://compgen.dib.uth.gr/LDSeeker/](https://compgen.dib.uth.gr/LDSeeker/)

LDSeeker can perform two main tasks:

1. **LD Annotation:** Find LD partners for your input SNPs in a reference panel.  
2. **Pairwise LD & Pruning:** Calculate LD between SNPs in your dataset and perform LD pruning (clumping) to identify independent signals.

## **TABLE OF CONTENTS**

* [Features](https://www.google.com/search?q=%23features)  
* [Requirements](https://www.google.com/search?q=%23requirements)  
* [Data Configuration](https://www.google.com/search?q=%23data-configuration)  
* [Usage](https://www.google.com/search?q=%23usage)  
* [Arguments](https://www.google.com/search?q=%23arguments)  
* [Input File Format](https://www.google.com/search?q=%23input-file-format)  
* [Output Files](https://www.google.com/search?q=%23output-files)  
* [License](https://www.google.com/search?q=%23license)

## **FEATURES**

* **Multi-Panel Support:** Seamlessly switch between major reference panels (1000 Genomes, UK Biobank, HGDP, TOP-LD, etc.).  
* **High Performance:** Built on polars for extremely fast processing of large Parquet/TSV datasets.  
* **Pairwise LD Calculation:** Calculate ![][image1] and ![][image2] between provided variants.  
* **LD Pruning (Clumping):** automated "greedy" pruning to filter GWAS results down to independent loci based on P-values and LD thresholds.  
* **Flexible Filtering:** Filter by ![][image1], MAF, and P-value thresholds.

## **REQUIREMENTS**

* Python 3.9+  
* **Polars**  
* **Pandas**

You can install the dependencies via pip:

pip install polars pandas

## **DATA CONFIGURATION**

**⚠️ Note:** The current version of LDSeeker\_functions.py expects reference panel data (Parquet files) to be located in specific directories (e.g., ref/).

Before running the tool, you must either:

1. Move your reference data to ref/ folder.  
2. **OR** (Recommended) Edit LDSeeker\_functions.py and update the file paths in the \_file variable definitions (e.g., maf\_file, ld\_file) within the specific function for the population you are using.

### **Supported Reference Panels**

The tool expects Parquet files for the following panels:

* **1000G\_hg38** (Standard & High Coverage)  
* **UKBB** (UK Biobank)  
* **HGDP** (Human Genome Diversity Project)  
* **TOP\_LD**  
* **Hap\_Map**  
* **Pheno\_Scanner**

## **USAGE**

### **1\. Basic LD Annotation (Non-Pairwise)**

Finds all variants in the reference panel that are in LD with your input SNPs.

python LDSeeker.py \\  
  \--file-path input\_gwas.txt \\  
  \--r2threshold 0.6 \\  
  \--pop EUR \\  
  \--maf 0.01 \\  
  \--ref 1000G\_hg38

### **2\. Pairwise LD Calculation**

Calculates LD specifically *between* the SNPs provided in your input file.

python LDSeeker.py \\  
  \--file-path input\_gwas.txt \\  
  \--r2threshold 0.1 \\  
  \--pop EUR \\  
  \--maf 0.01 \\  
  \--ref 1000G\_hg38 \\  
  \--pairwise YES

### **3\. LD Pruning (Clumping)**

Calculates pairwise LD and then removes SNPs that are in high LD with a more significant SNP (lower P-value).

python LDSeeker.py \\  
  \--file-path input\_gwas.txt \\  
  \--r2threshold 0.2 \\  
  \--pop EUR \\  
  \--maf 0.01 \\  
  \--ref UKBB \\  
  \--pairwise YES \\  
  \--ld-prune YES \\  
  \--ld-prune-col P \\  
  \--ld-prune-threshold 5e-8 \\  
  \--ld-prune-mode below

## **ARGUMENTS**

| Argument | Required | Default | Description |
| :---- | :---- | :---- | :---- |
| \--file-path | Yes | \- | Path to input GWAS summary statistics (Tab-separated). |
| \--r2threshold | Yes | \- | Minimum ![][image1] value (0.0 \- 1.0). |
| \--pop | Yes | \- | Population code (e.g., EUR, AMR, AFR, EAS, SAS). |
| \--maf | Yes | \- | Minor Allele Frequency threshold (e.g., 0.01). |
| \--ref | No | 1000G\_hg38 | Reference panel. Options: 1000G\_hg38, 1000G\_hg38\_high\_cov, UKBB, HGDP, TOP\_LD, Hap\_Map, Pheno\_Scanner. |
| \--pairwise | No | NO | Calculate pairwise LD between input SNPs (YES or NO). |
| \--imp\_list | No | \- | Path to a file containing a specific list of SNPs to impute/filter (no header). |
| \--ld-prune | No | NO | Apply LD pruning? (Requires \--pairwise YES). |
| \--ld-prune-col | No | P | Column name to use for ranking SNPs (usually P-value). |
| \--ld-prune-prefix | No | LD\_pruned | Prefix for pruning output files. |
| \--ld-prune-threshold | No | None | Filter input SNPs by this value before pruning (e.g., 0.05). |
| \--ld-prune-mode | No | below | Keep rows below or above the prune threshold. |

## **INPUT FILE FORMAT**

The input file (--file-path) should be a tab-separated text file.

* **Required Column:** SNP (rsID).  
* **Optional Columns:**  
  * CHR (Chromosome): If present, processing is faster as it targets specific chromosomes.  
  * P (P-value): Required if using \--ld-prune (or specify another column with \--ld-prune-col).

## **OUTPUT FILES**

Depending on the mode, LDSeeker generates the following TSV files:

* **LD\_info\_chr\_all.txt**: Results for standard LD annotation.  
* **LD\_info\_chr\_all\_pairwise.txt**: Matrix/List of pairwise LD values.  
* **LD\_pruned\_kept.txt**: (Pruning mode) The independent SNPs retained after pruning.  
* **LD\_pruned\_pruned.txt**: (Pruning mode) The SNPs removed due to high LD.

## **LICENSE**

Copyright (C) 2026 Pantelis Bagos.

Freely distributed under the **GNU General Public Licence (GPLv3)**.
 
