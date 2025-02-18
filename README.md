# *micRoclean*: Decontamination for low-biomass microbiome data

<img src="https://github.com/rachelgriffard/micRoclean/blob/main/images/micRoclean.svg" height = "200" align = "right">

micRoclean contains two pipelines within the *micRoclean()* function for decontaminating low-biomass microbiome data.

**For questions on installation or usage, please submit an [issue](https://github.com/rachelgriffard/micRoclean/issues) or [discussion](https://github.com/rachelgriffard/micRoclean/discussions) via GitHub.**

Please **download** the [vignette file](https://github.com/rachelgriffard/micRoclean/tree/main/vignettes) file in this repository for a detailed run through of this package functionality.

## Installation
To install the micRoclean package, users should use the *install_github* function from the **devtools** package. The full command is as follows:
```
devtools::install_github("rachelgriffard/micRoclean")

library(micRoclean)
```

The latest micRoclean release is available for download from the [repository](https://github.com/rachelgriffard/micRoclean).

## Usage

**Please download the [vignette file](https://github.com/rachelgriffard/micRoclean/tree/main/vignettes) file in this repository for a detailed run through of this package functionality.**

<img src="https://github.com/user-attachments/assets/a78fa445-3c67-42af-a981-0580e636df50" width = "1000" align = "center">

### micRoclean input

1. *Count matrix* - A samples (n) by features (p) matrix generated from 16S-rRNA sequencing data
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
2. *Metadata* - A metadata matrix with samples (n) as rows and two required columns *is_control* and *sample_type*. Two optional columns can be included named *batch* and *sample_well*. It is important that the naming scheme of the columns for the metadata match as seen below.
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
| Sample_2  |  FALSE | plasma | B| A4|
| Sample_3  |  TRUE | DNA extraction control | B| B3| 
| Sample_4  |  FALSE | plasma | A| B1|
| Sample_5  |  TRUE | DNA extraction control | B| B4|
| Sample_6  |  FALSE | plasma | B| C12|

### Original Composition Pipeline
This pipeline should be used when the user:
1. Has sample well information available and/or
2. Wants to primarily characterize the original composition of the sample prior to contamination and/or
3. Has only one batch OR has multiple batches with controls in each batch

Furthermore, users must have control samples present in each batch for this method to be used.

This pipeline implements the [SCRuB method](https://www.nature.com/articles/s41587-023-01696-w) for decontamination (Austin et al., 2023). To run this pipeline, the user can input their data as such:
```
orig.composition_results = micRoclean(counts = counts,
                               meta = metadata,
                               research_goal = 'orig.composition',
                               control_name = 'Control')
```
Once run, the pipeline will return a list object with:
1. *Decontaminated counts matrix (decontaminated_count)* - A samples (n) by features (p) matrix with the decontaminated counts
2. *Filtering loss value (FL)* - A numeric value between 0 and 1

### Biomarker Identification Pipeline
This pipeline should be used when the user:
1. Wants to primarily identify potential biomarkers
2. Does not have sample well information available

The Biomarker Identification Pipeline contains multiple steps to identify potential contaminants, as visualized here:
<img src = "https://github.com/rachelgriffard/micRoclean/assets/95938614/f7c3290e-026a-4d53-9939-e737ea400da1" align = "center" height = 400>

To run this pipeline, the user can input their data as such:
```
biomarkerID_results = micRoclean(counts = counts,
                               meta = metadata,
                               research_goal = 'biomarker',
                               blocklist = bl,
                               technical_replicates = tr,
                               control_name = 'Control',
                               remove_if = 1, #optional
                               step2_threshold = 0.5) #optional
```
Where:
* *blocklist* - character string
* *technical_replicates* - Data frame indicating pairs of technical replicates across batches by sample name. See example below where Sample 1 and Sample 6 are technical replicates.
* *control_name* - String name of control in metadata
* *remove_if* - Number of steps identifying feature as contaminant to remove from final dataframe

| *batch_1* | *batch_2* |
| ------------- | ------------- |
| Sample_1  |  Sample_6 |
| Sample_2  |  Sample_4 |

Once run, the pipeline will return a list object with:
1. *Decontaminated counts matrix (decontaminated_count)* - A samples (n) by non-contaminant features (p - c) matrix with the decontaminated counts
2. *Filtering loss value (FL)* - A numeric value between 0 and 1
3. *Contaminant ID (contaminant_id)* - Dataframe with features (p) by removal steps and boolean value indicating TRUE if tagged as contaminant in that step, FALSE otherwise.
4. *Removed (removed)* - Character vector of all samples tagged as contaminants and removed from the decontaminated count matrix

For convenience, the default blocklist from [Eisenhofer et al. (2019)](https://www.sciencedirect.com/science/article/pii/S0966842X18302531?via%3Dihub) is included below and can be copied, if desired:
```
bl = c('Actinomyces','Corynebacterium','Arthrobacter',
       'Rothia','Propionibacterium','Atopobium',
       'Sediminibacterium','Porphyromonas','Prevotella',
       'Chryseobacterium','Capnocytophaga','Chryseobacterium',
       'Flavobacterium','Pedobacter','UnclassifiedTM7',
       'Bacillus','Geobacillus','Brevibacillus','Paenibacillus',
       'Staphylococcus','Abiotrophia','Granulicatella',
       'Enterococcus','Lactobacillus','Streptococcus',
       'Clostridium','Coprococcus','Anaerococcus','Dialister','Megasphaera',
       'Veillonella','Fusobacterium','Leptotrichia','Brevundimonas','Afipia',
       'Bradyrhizobium','Devosia','Methylobacterium','Mesorhizobium','Phyllobacterium',
       'Rhizobium','Methylobacterium','Phyllobacterium','Roseomonas','Novosphingobium	',
       'Sphingobium','Sphingomonas','Achromobacter','Burkholderia','Acidovorax',
       'Comamonas','Curvibacter','Pelomonas','Cupriavidus','Duganella',
       'Herbaspirillum','Janthinobacterium','Massilia','Oxalobacter','Ralstonia',
       'Leptothrix','kingella','Neisseria','Escherichia','Haemophilus',
       'Acinetobacter','Enhydrobacter','Pseudomonas','Stenotrophomonas','Xanthomonas')
```

Optionally, users can input their results list object from the Biomarker Identification Pipeline into the *visualize_pipeline* function. If *interactive* is set to true, the resulting visualization will be interactive within an HTML pop up.
```
visualize_pipeline(biomarkerID_results,
                   interactive = FALSE)
```

<img src = "https://github.com/rachelgriffard/micRoclean_development/assets/95938614/3f26fedf-47b4-4d1d-bd73-23ca6f32d963" align = "center">

### Well to well contamination (well2well)

Well-to-well leakage is a common contamination where biological samples leak into controls. For batches where users do not have well location information, the well2well function in micRoclean assigns pseudo-locations in a 96-well plate, assuming common order of samples vertically or horizontally. The function then estimates the proportion of each control that originates from a biological sample, indicating well-to-well leakage, through SCRuB package spatial functionality (Austin et al., 2023). If the level of well-to-well contamination is higher than 0.1, the function will return a warning message indicating to the user that they should obtain the well location information for their data and run through pipeline1, which can account for this form of contamination. The user can choose to ignore this message and continue with the pipeline2 analysis, but this is discouraged.

For more information about this, users are encouraged to read the well-to-well functionality of the [SCRuB package (Austin et al., 2023](https://www.nature.com/articles/s41587-023-01696-w).

### Filtering loss (FL)
First introduced for use in a filtering method PERfect by Smirnova, Huzurbazar, and Jafari (2019), the filtering loss (FL) statistic is implemented in the micRoclean package to quantify the impact due to filtering features out in the above pipelines. The filtering loss value is between zero and one, indicating low to high contribution respectively from the removed reads to the total convariance structure. As the filtering loss value gets closer to one, users should be concerned about potential overfiltering.

Filtering loss for removal of reads $J$ is defined as

$$FL(J) = 1 - \frac{\|\|Y^T Y\|\|_F^2}{\|\|X^TX\||\_F^2}$$

where $X$ is the n x p count matrix and $Y$ is the n x q matrix resulting from the partial removal of reads or whole removal of features after applying the decontamination method. As full features may not be removed, $q ≤ p$. 

For more detailed information, users are suggested to read the methods section 2.1 of the [Smirnova, Huzurbazar, and Jafari (2109)](https://doi.org/10.1093/biostatistics/kxy020) publication.

## References
Austin, G. I., Park, H., Meydan, Y., Seeram, D., Sezin, T., Lou, Y. C., Firek, B. A., Morowitz, M. J., Banfield, J. F., Christiano, A. M., Pe'er, I., Uhlemann, A. C., Shenhav, L., & Korem, T. (2023). Contamination source modeling with SCRuB improves cancer phenotype prediction from microbiome data. Nature biotechnology, 41(12), 1820–1828. [https://doi.org/10.1038/s41587-023-01696-w](https://doi.org/10.1038/s41587-023-01696-w) 

Davis, N.M., Proctor, D.M., Holmes, S.P. et al. Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome 6, 226 (2018). [https://doi.org/10.1186/s40168-018-0605-2](https://doi.org/10.1186/s40168-018-0605-2)

Smirnova, E., Huzurbazar, S., & Jafari, F. (2019). PERFect: PERmutation Filtering test for microbiome data. Biostatistics (Oxford, England), 20(4), 615–631. [https://doi.org/10.1093/biostatistics/kxy020](https://doi.org/10.1093/biostatistics/kxy020)

Zozaya-Valdés, E., Wong, S. Q., Raleigh, J., Hatzimihalis, A., Ftouni, S., Papenfuss, A. T., Sandhu, S., Dawson, M. A., & Dawson, S. J. (2021). Detection of cell-free microbial DNA using a contaminant-controlled analysis framework. Genome biology, 22(1), 187. <https://doi.org/10.1186/s13059-021-02401-3>
