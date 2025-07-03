# BacTax-RNA

### Author: Rodgee Mae Guden  
**Marine Biology Research Group, Ghent University**  
**Contact:** rodgeemae.guden@ugent.be  

---

## Overview

**BacTax-RNA** is a modular Bash workflow for bacterial taxonomic analysis from metatranscriptomic RNA-Seq data. The pipeline focuses on high-throughput processing of RNA reads to extract and classify rRNA sequences, providing taxonomic profiles and OTU tables for bacterial community studies.

This workflow is designed for metatranscriptomic datasets with an emphasis on removing rRNA, classifying bacterial rRNA reads, and generating OTU tables suitable for downstream community composition, diversity, and differential abundance analyses.

---

## Key Features

✔️ Quality Control with FastQC & MultiQC  
✔️ Adapter trimming with Fastp  
✔️ SortMeRNA-based rRNA filtering  
✔️ Kraken2-based taxonomic classification of rRNA reads  
✔️ Bacterial-specific filtering of classification results  
✔️ OTU table generation for genus and species levels  
✔️ Krona visualization  
✔️ Diversity and differential abundance analysis ready  

---

## Dependencies

| Tool         | Purpose                              | Installation link                           |
|--------------|--------------------------------------|----------------------------------------------|
| `fastqc`     | Read quality control                | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| `multiqc`    | Aggregate QC reports                | [MultiQC](https://multiqc.info/)            |
| `fastp`      | Adapter trimming                    | [fastp](https://github.com/OpenGene/fastp)  |
| `sortmerna`  | rRNA filtering                      | [SortMeRNA](https://github.com/biocore/sortmerna) |
| `kraken2`    | Taxonomic classification            | [Kraken2](https://github.com/DerrickWood/kraken2) |
| `krona`      | Taxonomy visualization              | [Krona](https://github.com/marbl/Krona)     |
| `R`          | OTU table processing & diversity analysis | [R Project](https://www.r-project.org/) |

Custom R scripts are included for OTU table correction and diversity analyses.

---

## Usage

### 1. Make the script executable

```bash
chmod +x metatranscriptomics_rrna_pipeline.sh
```

### 2. Run the workflow

Example for running trimming, rRNA removal, and taxonomic analysis:

```bash
bash metatranscriptomics_rrna_pipeline.sh -s qc,trim_fastp,postqc_trim_fastp,download_db,index_db,remove_rrna_unmerged,postqc_sortmerna,summary_sortmerna,classify_rrna,krona_kraken,filter_bacteria_rrna,generate_bacteria_otu_R,autofix_otu_header,run_diversity_analysis,run_deseq2
```

Steps are modular. You can run selected steps by specifying them with commas.

---

## Output

The workflow generates:

- Quality reports before and after trimming  
- Trimmed FASTQ files  
- rRNA-removed reads  
- Kraken2 taxonomic reports for rRNA reads  
- Krona visualizations  
- Bacteria-filtered Kraken2 reports  
- OTU tables at genus and species levels  
- Diversity and differential abundance analysis results  

---

## Notes

- Paths to input/output directories and databases are configured at the top of the script. Modify them to fit your setup.  
- The workflow assumes paired-end FASTQ files ending with `_R1.fastq.gz` and `_R2.fastq.gz`.  
- The SortMeRNA database is automatically downloaded and indexed if requested.  
- Email notifications are implemented but currently disabled by default.  

---

## Citation & License

© 2025 Rodgee Mae Guden  
**Marine Biology Research Group, Ghent University**  
All rights reserved.

For permission to reuse or adapt this workflow, contact the author.

---

## Acknowledgments

This workflow was developed to support RNA-based bacterial taxonomic profiling for the ORCHESTRA project.

---

## Questions & Issues

Feel free to open an [Issue](https://github.com/yourusername/BacTax-RNA/issues) or contact the author directly for support.


