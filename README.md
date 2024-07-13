# Pancancer oncRNAs
Orphan non-coding RNAs (oncRNAs) are a class of cancer-specific, non-coding small RNAs. In this project directory, we include the scripts, code, and analysis notebooks used to systematically annotate a pancancer set of oncRNAs within The Cancer Genome Atlas (TCGA). Details and results can be found in our released [preprint](https://www.biorxiv.org/content/10.1101/2024.03.19.585748v1.full). 

## Package Dependencies
Running our oncRNA discovery pipeline will require installing the packages listed below (here we also include the versions used in our [preprint](https://www.biorxiv.org/content/10.1101/2024.03.19.585748v1.full). We recommend installing the following packages and dependencies using the conda/mamba ecosystem. <br>

```
- mongo
- pymongo
- scipy
- pandas
- numpy
- statsmodels
- bedtools
- pysam
- joblib
```

## Notebooks
Notebooks used for the systematic annotation of oncRNAs in TCGA datasets. **Note:** Access to the controlled TCGA sequencing files used in this project will require obtainig appropriate authorization.

1. `01_preprocess_data.ipynb` – notebook used to preprocess the TCGA sequencing files (bam). 
2. `02_TCGA_oncRNA_analysis.ipynb` – notebook containing our analytical framework for calling oncRNAs.


### License
MIT license

### Citing
[Systematic annotation of orphan RNAs reveals blood-accessible molecular barcodes of cancer identity and cancer-emergent oncogenic drivers
](https://www.biorxiv.org/content/10.1101/2024.03.19.585748v1.full) (preprint).

