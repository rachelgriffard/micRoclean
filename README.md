# **micRoclean**: Decontamination for low-biomass microbiome data

<img src="https://github.com/rachelgriffard/micRoclean_development/blob/main/micRoclean.svg" height = "200" align = "right">

micRoclean contains two  pipelines aimed at decontaminating low-biomass microbiome data.

**For questions on installation or usage, please submit an issue or discussion via GitHub**

Please **download** the [vignette file](https://github.com/rachelgriffard/micRoclean/blob/main/vignette.html) file in this repository for a detailed run through of this package functionality.

## Installation
To install the micRoclean package, users should use the *install_github* function from the **devtools** package. The full command is as follows:
```
devtools::install_github("rachelgriffard/micRoclean")
library(micRoclean)
```

The latest micRoclean release is available for download from the [repository](https://github.com/rachelgriffard/micRoclean).

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
| Sample_5  |  0 | 1 | 6| 0 |
| Sample_6  |  6 | 14 | 21| 2 |
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
| Sample_1  |  FALSE | plasma_1 | A| A2|
| Sample_2  |  FALSE | plasma_1 | B| A4|
| Sample_3  |  TRUE | DNA extraction control | B| B3| 
| Sample_4  |  FALSE | plasma_2 | A| B1|
| Sample_5  |  TRUE | DNA extraction control | B| B4|
| Sample_6  |  FALSE | plasma_2 | B| C12|

### Pipeline 1
This pipeline should be used when the user:
1. Has sample well information available
2. Wants to primarily characterize the original composition of the sample prior to contamination

To run this pipeline, the user can input their data as such:
```
pipeline_1_results = pipeline1(counts = counts,
                               meta = metadata)
```
Once run, the pipeline will return a list object with:
1. *Decontaminated counts matrix (decontaminated_count)* - A samples (n) by features (p) matrix with the decontaminated counts
2. *Filtering loss value (FL)* - A numeric value between 0 and 1

### Pipeline 2
This pipeline should be used when the user:
1. Wants to primarily identify potential biomarkers
2. Does not have sample well information available

To run this pipeline, the user can input their data as such:
```
pipeline_2_results = pipeline2(counts = counts,
                               meta = metadata,
                               blocklist = bl,
                               technical_replicates = tr,
                               remove_if = 1, #optional
                               step2_threshold = 0.5) #optional
```
Where:
* *blocklist* - character string
* *technical_replicates* - Data frame indicating pairs of technical replicates across batches by sample name. See example below where Sample 1 and Sample 6 are technical replicates.

| *batch_1* | *batch_2* |
| ------------- | ------------- |
| Sample_1  |  Sample_6 |
| Sample_2  |  Sample_4 |

Once run, the pipeline will return a list object with:
1. *Decontaminated counts matrix (decontaminated_count)* - A samples (n) by non-contaminant features (p - c) matrix with the decontaminated counts
2. *Filtering loss value (FL)* - A numeric value between 0 and 1
3. *Contaminant ID (contaminant_id)* - Dataframe with features (p) by removal steps and boolean value indicating TRUE if tagged as contaminant in that step, FALSE otherwise.
4. *Removed (removed)* - Character vector of all samples tagged as contaminants and removed from the decontaminated count matrix

Optionally, users can input their results list object from pipeline2 into the *visualize_pipeline* function. If *interactive* is set to true, the resulting visualization is interactive.
```
visualize_pipeline(pipeline_2_results,
                   interactive = FALSE)
```
![VennExample](https://github.com/rachelgriffard/micRoclean_development/assets/95938614/3f26fedf-47b4-4d1d-bd73-23ca6f32d963)


### Filtering loss
