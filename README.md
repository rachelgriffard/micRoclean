# micRoclean: Decontamination for Low-Biomass Microbiome Data

<img src="https://github.com/rachelgriffard/micRoclean_development/blob/main/micRoclean.svg" height = "200" align = "right">

micRoclean contains two decontamination pipelines aimed at low-biomass microbiome data.

**For questions on installation or usage, please submit an issue via GitHub or contact Rachel Griffard (rgriffard@kumc.edu).**

Please **download** the [vignette file](https://github.com/rachelgriffard/micRoclean/blob/main/vignette.html) file in this repository for a detailed run through of this package functionality.

## Installation
To install the micRoclean package, users should use the *install_github* function from the **devtools** package. The full command is as follows:
```
devtools::install_github("rachelgriffard/micRoclean")
```

The latest opitma release is available for download from the [repository](https://github.com/rachelgriffard/micRoclean).

## Usage

**Please download the [vignette file](https://github.com/rachelgriffard/micRoclean/blob/main/vignette.html) file in this repository for a detailed run through of this package functionality.**

<img src="https://github.com/rachelgriffard/micRoclean_development/blob/main/FlowChart.png" width = "1000" align = "center">

### micRoclean input

1. *Count matrix* - A samples (n) by features (p) matrix
```
head(counts)
```
| | taxa_1 | taxa_2 | taxa_3 | taxa_4 |
| :-------------: | ------------- |------------- |------------- |------------- |
| Sample_1  |  0 | 5 | 0|20 |
| Sample_2  |  15 | 5 | 0|0 |
| Sample_3  |  0 | 13 | 0| 200 |
| Sample_4  |  4 | 5 | 0| 0 |
3. *Metadata* - A metadata matrix with samples (n) as rows and two required columns *is_control* and *sample_type*. Two optional columns can be included named *batch* and *sample_well*. It is important that the naming scheme of the columns for the metadata match as seen below.
  * *is_control* - Boolean variable where *TRUE* indicates a extraction negative control sample and *FALSE* otherwise.
  * *sample_type* - Sample types named by string, indicating which samples should be read together.
  * (optional) *batch* - String indicating batch.
  * (optional) *sample_well* - String that indicates the well location of the sample. Must be in LETTER-NUMBER format, e.g. A1.
```
head(metadata)
```
| | is_control | sample_type | batch | sample_well |
| :-------------: | ------------- |------------- |------------- |------------- |
| Sample_1  |  FALSE | plasma | A| A2|
| Sample_2  |  FALSE | plasma | A| A4|
| Sample_3  |  TRUE | DNA extraction control | B| B3| 
| Sample_4  |  FALSE | plasma | B| B1|

