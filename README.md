# **LDSeeker**

 

**LDSeeker** is a high-performance tool for exploring Linkage Disequilibrium (LD) in Genome-Wide Association Studies (GWAS), using multiple reference panels

## **Web Tool:** 

A web-based version is available at [https://compgen.dib.uth.gr/LDSeeker/](https://compgen.dib.uth.gr/LDSeeker/)

## **Download Reference panels and resources**
 
| Resource | Web address |
| :--- | :--- |
| LDSeeker GitHub repository | https://github.com/pbagos/LDSeeker/ |
| LDSeeker web tool | https://compgen.dib.uth.gr/LDSeeker/ |
| UKBB LD reference panel | https://zenodo.org/records/18847278 |
| 1000 Genomes Project – high coverage LD reference panel | https://zenodo.org/records/18849097 |
| 1000 Genomes Project LD reference panel | https://zenodo.org/records/18861527 |
| TOP-LD LD reference panel | http://195.251.108.185/ref_panels/TOP_LD |
| HapMap LD reference panel | https://zenodo.org/records/20213914 |
| LASI-DAD LD reference panel | https://zenodo.org/records/20256884 |
| HGDP LD reference panel | http://195.251.108.185/HGDP/ |
| SNP to Gene mapping reference | https://zenodo.org/records/20161132 |

**NOTE: The Pheno Scanner LD reference panel is available upon request from the Pheno Scanner database (http://www.phenoscanner.medschl.cam.ac.uk)
## **Key Capabilities**

1. **LD Annotation:** Identify LD partners for input SNPs within a reference panel.  
2. **Pairwise LD & Pruning:** Calculate LD between SNPs in a dataset and perform pruning to identify independent signals.


## **Table of Contents**

* [Features](#features)  
* [Installation Guide](#installation-guide)  
* [Data Configuration](#data-configuration)  
  * [Supported Reference Panels](#supported-reference-panels)  
* [Usage](#usage)  
  * [1. Basic LD Annotation](#1-basic-ld-annotation-non-pairwise)  
  * [2. Pairwise LD Calculation](#2-pairwise-ld-calculation)  
  * [3. LD Pruning](#3-ld-pruning)  
* [Arguments](#arguments)  
* [Input File Format](#input-file-format)  
* [Output Files](#output-files)  
* [License](#license)
  
## **Features**
* **Multiple LD Reference Panels:** Seamlessly switch between major reference panels (1000 Genomes Project, UK Biobank, HGDP, TOP-LD, etc.).  
* **Pairwise LD Calculation:** Calculate ![][image1] and ![][image2] between provided variants.  
* **LD Pruning:** Pruning to filter GWAS results down to independent loci based on LD and  P-value thresholds.  


## **Installation Guide**

LDSeeker is written in Python (ver. 3.10)

1)	Clone or download LDSeeker from: https://github.com/pbagos/LDSeeker
  ```
  git clone  https://github.com/pbagos/LDSeeker
  ```

2)	After downloading the .zip folder of LDSeeker from GitHub, extract it to a working directory. 

3)	Το install the requirements, pip needs to be installed. Download the script for pip, from: https://bootstrap.pypa.io/get-pip.py.

4)	Open a terminal/command prompt, cd to the folder containing the get-pip.py file and run:
    ```
    python get-pip.py
    ```
5)	To install the mentioned requirements with pip, open a terminal/command prompt and run:
    ```
    pip install -r  requirements.txt
    ```

## **Data Configuration**

**⚠️ Important:** LDSeeker\_functions.py expects reference panel data (Parquet files) to be located in specific directories.

Before running the tool, choose one of the following options:

1. **Option A (Standard):** Create a ref/ directory in the project root and move your reference data there.  
2. **Option B (Custom):** Edit LDSeeker\_functions.py and update the \_file variable definitions (e.g., maf\_file, ld\_file) to point to your custom data paths.

### **Supported Reference Panels**

The tool expects Parquet files for the following panels:

### **Supported Reference Panels**

The tool expects Parquet files for the following reference panels and populations:

| Reference Panel | Label | Population Sample | Number of Samples | Number of Variants |
| :--- | :--- | :--- | ---: | ---: |
| **TOP-LD (hg38)** | EUR | European | 13,160 | 69,524,944 |
| | AFR | African | 1,335 | 60,392,677 |
| | SAS | South Asian | 239 | 22,309,649 |
| | EAS | East Asian | 844 | 35,538,656 |
| **PhenoScanner (Phase 3, hg19 & hg38)** | EUR | European | 503 | 11,159,862 |
| | AFR | African | 661 | 11,159,862 |
| | SAS | South Asian | 489 | 11,159,862 |
| | EAS | East Asian | 504 | 11,159,862 |
| | AMR | American | 347 | 11,159,862 |
| **HapMap (hg16)** | ASW | African ancestry in Southwest USA | 90 | 1,561,113 |
| | CEU | Utah residents with Northern and Western European ancestry from the CEPH collection | 180 | 1,412,161 |
| | CHB | Han Chinese in Beijing, China | 90 | 1,328,283 |
| | CHD | Chinese in Metropolitan Denver, Colorado | 100 | 1,305,880 |
| | GIH | Gujarati Indians in Houston, Texas | 100 | 1,407,540 |
| | JPT | Japanese in Tokyo, Japan | 91 | 1,296,969 |
| | LWK | Luhya in Webuye, Kenya | 100 | 1,529,438 |
| | MEX | Mexican ancestry in Los Angeles, California | 90 | 1,409,947 |
| | MKK | Maasai in Kinyawa, Kenya | 180 | 1,419,626 |
| | TSI | Toscans in Italy | 100 | 1,419,920 |
| | YRI | Yoruba in Ibadan, Nigeria | 180 | 1,501,085 |
| **1000 Genomes Project Phase 3 (hg38)** | EUR | European | 632 | 8,193,280 |
| | AFR | African | 893 | 13,876,891 |
| | SAS | South Asian | 601 | 8,579,150 |
| | EAS | East Asian | 585 | 7,245,426 |
| | AMR | Admixed American | 490 | 9,347,814 |
| **1000 Genomes Project Phase 3 (hg38), high coverage** | EUR | European | 632 | 9,263,406 |
| | AFR | African | 893 | 15,952,015 |
| | SAS | South Asian | 601 | 9,439,730 |
| | EAS | East Asian | 585 | 7,989,898 |
| | AMR | Admixed American | 490 | 10,106,451 |
| **UKBB (hg19)** | EUR | European | 362,446 | 1,431,634 |
| | AFR | African | 6,255 | 1,259,175 |
| | CSA | Central and South Asian | 8,284 | 1,379,991 |
| | EAS | East Asian | 2,700 | 1,146,910 |
| | AMR | Admixed American | 987 | 1,286,943 |
| | MID | Middle Eastern | 1,567 | 1,307,269 |
| **LASI-DAD (hg38)** | IND | Indian | 2,680 | 119,948,583 |
| **HGDP** | EUR | European | 155 | 10,144,721 |
| | AFR | African | 104 | 18,062,644 |
| | CSA | Central South Asian | 197 | 11,281,764 |
| | EAS | East Asian | 223 | 9,631,805 |
| | AMR | Admixed American | 61 | 7,592,946 |
| | MID | Middle Eastern | 161 | 11,731,317 |
| | OCN | Oceanian | 28 | 8,147,809 |


## **Usage**

### **1\. Basic LD Annotation (Non-Pairwise)**

Finds all variants in the reference panel that are in LD with your input SNPs.
```
python LDSeeker.py --file-path input_gwas.txt --r2threshold 0.6 --pop EUR --maf 0.01 --ref 1000G_hg38
```
### **2\. Pairwise LD Calculation**

Calculates LD specifically *between* the SNPs provided in your input file.

```
python LDSeeker.py --file-path input_gwas.txt --r2threshold 0.6 --pop EUR --maf 0.01 --ref 1000G_hg38  --pairwise YES
```
### **3\. LD Pruning**

Calculates pairwise LD and removes SNPs that are in high LD.
```
python LDSeeker.py --file-path input_gwas.txt --r2threshold 0.2 --pop EUR --maf 0.01 --ref UKBB --pairwise YES --ld-prune YES 
```
## Arguments

### Core

| Argument | Required | Default | Choices | Description |
| :-- | :-- | :-- | :-- | :-- |
| `--file-path` | **Yes** | — | — | Path to input GWAS summary statistics (tab-separated). |
| `--r2threshold` | **Yes** | — | — | Minimum r² value (0.0–1.0). |
| `--pop` | **Yes** | — | — | Population code (e.g. EUR, AMR, AFR, EAS, SAS). |
| `--maf` | **Yes** | — | — | Minor Allele Frequency threshold (e.g. 0.01). |
| `--ref` | No | `1000G_hg38` | — | Reference panel: `1000G_hg38`, `1000G_hg38_high_cov`, `TOP_LD`, `Hap_Map`, `Pheno_Scanner`, `UKBB`, `HGDP`, `LASI_DAD`. |
| `--pairwise` | No | `NO` | — | Return only pairwise LD among input SNPs (`YES`/`NO`). |
| `--imp_list` | No | — | — | File of SNPs to impute/add as targets (one per line, no header). |

## LD pruning (only active with `--pairwise YES --ld-prune YES`)

| Argument | Required | Default | Choices | Description |
| :-- | :-- | :-- | :-- | :-- |
| `--ld-prune` | No | `NO` | — | Apply pairwise-LD pruning (`YES`/`NO`). |
| `--ld-prune-metric` | No | `P` | `P`, `Z` | Rank SNPs by lowest P (`P`) or highest \|Z\| (`Z`). |
| `--ld-prune-prefix` | No | `LD_pruned` | — | Prefix for pruned output (e.g. `LD_pruned_kept.txt`). |
| `--ld-prune-p-col` | No | `P` | — | Column name for P-values. |
| `--ld-prune-p-threshold` | No | `None` | — | Pre-filter by P before pruning (e.g. 5e-8); rows failing it are dropped. |
| `--ld-prune-p-threshold-mode` | No | `below` | `below`, `above` | Keep P `below` or `above` the threshold. |
| `--ld-prune-z-col` | No | `Z` | — | Column name for Z-values. |
| `--ld-prune-z-threshold` | No | `1.96` | — | Pre-filter by \|Z\| before pruning; rows with \|Z\| ≤ threshold dropped. |

## Significance filtering of input variants (applies to both standard and `--pairwise YES` runs)

| Argument | Required | Default | Choices | Description |
| :-- | :-- | :-- | :-- | :-- |
| `--significance` | No | — (off) | `significant`, `nonsignificant` | Keep only significant or only non-significant variants before LD retrieval. Omit to disable. |
| `--significance-metric` | No | `P` | `P`, `Z` | Define significance by two-tailed p from z (`P`) or directly by \|z\| (`Z`). |
| `--significance-threshold` | No | `5e-08` | — | p-value cutoff (metric `P`) or \|z\| cutoff (metric `Z`). |
| `--beta-col` | No | `BETA` | — | Effect-size column, for deriving z = BETA/SE when no Z column. |
| `--se-col` | No | `SE` | — | Standard-error column, for deriving z = BETA/SE. |
| `--sig-z-col` | No | `Z` | — | Precomputed z column (preferred over BETA/SE) for significance. |
| `--sig-p-col` | No | `P` | — | Precomputed p column, fallback when neither Z nor BETA/SE exists. |

## Notes

- `--pairwise`, `--ld-prune`, and `--imp_list` are not restricted by argparse `choices`. The YES/NO flags are tested as `== 'NO'` / `== 'YES'` after uppercasing, so `--pairwise MAYBE` is treated as pairwise mode (anything ≠ `NO`) and `--ld-prune TRUE` is treated as NO (anything ≠ `YES`). Add `choices=['YES','NO']` to enforce.
- `--ref` is not validated by argparse; an unsupported panel errors at runtime with `Unsupported ref_panel`.
- Metric `P` and metric `Z` are interchangeable at matched thresholds, since two-tailed p < α ⟺ \|z\| > Φ⁻¹(1 − α/2). For example `--significance-metric Z --significance-threshold 5.4513` selects the same variants as `--significance-metric P --significance-threshold 5e-8`.
## **Input File Format**

The input file (--file-path) must be a **tab-separated** text file.

| SNP | CHR | P |
| :---- | :---- | :---- |
| rs12345 | 1 | 0.05 |
| rs67890 | 1 | 1.2e-8 |

* **SNP (Required):** The rsID of the variant.  
* **CHR (Optional):** Chromosome. If present, processing is faster as it targets specific chromosomes.  
* **P (Conditional):** P-value. Required if using \--ld-prune (or specify another column with \--ld-prune-col).

## **Output Files**

Depending on the mode selected, LDSeeker generates the following TSV files:

1. **LD\_info\_chr\_all.txt**  
   * Results for standard LD annotation (LD partners found in reference).  
2. **LD\_info\_chr\_all\_pairwise.txt**  
   * Matrix/List of pairwise LD values between input SNPs.  
3. **LD\_pruned\_kept.txt**  
   * (Pruning mode) The independent SNPs retained after pruning.  
4. **LD\_pruned\_pruned.txt**  
   * (Pruning mode) The SNPs removed due to high LD with a lead variant.

## **License**

Freely distributed under the **GNU General Public Licence (GPLv3)**.

[image1]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAA8AAAAVCAYAAACZm7S3AAABAklEQVR4XmNgoBNglZOT8wXiWdLS0sLokvgAC1BTGZBmBHGA7L9ArI2mBjsQExMTB2mQlZX1A/GB7B1A3ICmDCdgBCouFxIS4gNx5OXllwLxHHRFRAGgxv9AihVdnCAAarwuLi7OjS5OEAgKCvLD/ApiwyXQgx9FEgiAgTUZaGsOlMsM5BeBWUDBB6CgB/kFSE8CxaWMjIwqkH8RiBWBWBIkh4yBaozBmkFBj6wAZDJQ7DWKImwAKOkCVKQJxNFQzdEgcSC9DCh3GF09VgBUeBqo4YCoqCgPuhxBALIVloJIBlD/KaGLEwTAaNKHBhTpAJRugfgKujhRABinnOiJZBgDAEB+OYlaoqoJAAAAAElFTkSuQmCC>

[image2]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABUAAAAZCAYAAADe1WXtAAABPklEQVR4Xu2UvUoDURCFN6gQjFbZdWH/SSeks7DRzsZG7HwE3yCClY0vIIhgkxewtbPSxgeI+ASS0kohFpozMiPryb2rrZgPhixzzszcO2w2CP4daZp28bPA+U+Korj0xCDP831YWlwjQH9HHHE+SJJkuSzLjx/ikesE0TC0z/kvYBhLYHKPNeRuoU3QYNtyuHomJ637ZtATnXJeqKqqrfrEchhwjKYnNdsMi1KUZdkuC4atQp51yHXj1TFxSwqiKFphzag3hX8Hz2/s+QZMV1bgIgzD1XrTX4Gmr4g7zhu45p40hOeFNS9acMZ5QzQ96T1rTvCehlqwzpqhL/kojuM11pzAfNi0K1z9QPUl1ny00HTY1BTaMzw3nPcC84YsH/Hg0C50LeesOcFuOlrgDDR8wh9hM/B9gebM+aNMAdsSZUMPB9lJAAAAAElFTkSuQmCC>
