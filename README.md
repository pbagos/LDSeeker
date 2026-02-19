# **LDSeeker**

 

**LDSeeker** is a high-performance tool for exploring Linkage Disequilibrium (LD) in Genome-Wide Association Studies (GWAS), using multiple reference panels

## **Web Tool:** 

A web-based version is available at [https://compgen.dib.uth.gr/LDSeeker/](https://compgen.dib.uth.gr/LDSeeker/)

## Download Reference panels 
You can download all the LD reference panels in .parquet format from [here](http://195.251.108.185/ref/).

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

| Reference Panel | Label | Population Sample | Number of Samples |
| :--- | :--- | :--- | :--- |
| **TOP-LD (hg38)** | EUR | European | 13,160 |
| | AFR | African | 1,335 |
| | SAS | South Asian | 239 |
| | EAS | East Asian | 844 |
| **PhenoScanner (Phase 3, hg19 & hg38)** | EUR | European | 503 |
| | AFR | African | 661 |
| | SAS | South Asian | 489 |
| | EAS | East Asian | 504 |
| | AMR | American | 347 |
| **HapMap (hg16)** | ASW | African ancestry in Southwest USA | 90 |
| | CEU | Utah residents with N/W European ancestry | 180 |
| | CHB | Han Chinese in Beijing, China | 90 |
| | CHD | Chinese in Metropolitan Denver, Colorado | 100 |
| | GIH | Gujarati Indians in Houston, Texas | 100 |
| | JPT | Japanese in Tokyo, Japan | 91 |
| | LWK | Luhya in Webuye, Kenya | 100 |
| | MEX | Mexican ancestry in Los Angeles, California | 90 |
| | MKK | Maasai in Kinyawa, Kenya | 180 |
| | TSI | Toscans in Italy | 100 |
| | YRI | Yoruba in Ibadan, Nigeria | 180 |
| **RAISS (1000G Phase 3, hg38)** | EUR | European | 632 |
| | AFR | African | 893 |
| | SAS | South Asian | 601 |
| | EAS | East Asian | 585 |
| | AMR | Admixed American | 490 |
| **UKBB (hg19)** | EUR | European | 362,446 |
| | AFR | African | 6,255 |
| | CSA | Central and South Asian | 8,284 |
| | EAS | East Asian | 2,700 |
| | AMR | Admixed American | 987 |
| | MID | Middle Eastern | 1,567 |
| **LASI-DAD** | IND | Indian | 2,680 |
| **HGDP** | EUR | European | 155 |
| | AFR | African | 104 |
| | CSA | Central South Asian | 197 |
| | EAS | East Asian | 223 |
| | AMR | Admixed American | 61 |
| | MID | Middle Eastern | 161 |
| | OCN | Oceanian | 28 |


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

Calculates pairwise LD and removes SNPs that are in high LD with a more significant SNP (lower P-value).
```
python LDSeeker.py --file-path input_gwas.txt --r2threshold 0.2 --pop EUR --maf 0.01 --ref UKBB --pairwise YES --ld-prune YES --ld-prune-col P --ld-prune-threshold 5e-8 --ld-prune-mode below
```
## **Arguments**

| Argument | Required | Default | Description |
| :---- | :---- | :---- | :---- |
| \--file-path | **Yes** | \- | Path to input GWAS summary statistics (Tab-separated). |
| \--r2threshold | **Yes** | \- | Minimum ![][image1] value (0.0 \- 1.0). |
| \--pop | **Yes** | \- | Population code (e.g., EUR, AMR, AFR, EAS, SAS). |
| \--maf | **Yes** | \- | Minor Allele Frequency threshold (e.g., 0.01). |
| \--ref | No | 1000G\_hg38 | Reference panel (e.g., UKBB, HGDP, TOP\_LD, Hap_Map, 1000G\_hg38\_high\_cov,1000G\_hg38). |
| \--pairwise | No | NO | Calculate pairwise LD between input SNPs (YES or NO). |
| \--imp\_list | No | \- | Path to a file containing a specific list of SNPs to impute/filter (no header). |
| \--ld-prune | No | NO | Apply LD pruning? (Requires \--pairwise YES). |
| \--ld-prune-col | No | P | Column name to use for ranking SNPs (usually P-value). |
| \--ld-prune-prefix | No | LD\_pruned | Prefix for pruning output files. |
| \--ld-prune-threshold | No | None | Filter input SNPs by this value before pruning (e.g., 0.05). |
| \--ld-prune-mode | No | below | Keep rows below or above the prune threshold. |

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
