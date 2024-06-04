# Function development for microbiome decontamination pipeline (micRoclean)
# 20240305
# Rachel Griffard

# Required libraries
library(phyloseq) # object wrapper
library(SummarizedExperiment) # object wrapper
library(tidyverse)
library(plotly) # for interactive feature
library(SCRuB) # pipeline 1
library(PERfect) # DFL
library(decontam) # pipeline 2 step2
library(microDecon) # pipeline 2
library(ANCOMBC) # pipeline 2
library(ggVennDiagram) # function 3 - pipeline 2 - comparison across removed
library(shiny) # function 3
library(ANCOMBC) # pipeline 2 step1
library(kappa) # pipeline 2 step3


# Function 0: "Sort" function based on order names by strings
## Determine if there is similarity based on this order with
## pop-up window containing information if so
## Provide "warning sign" about potential well to well interaction

#' @name well2well
#' @usage Determine the potential impact of spatial well to well relationships
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta Metadata with column one 'is_control' indicating TRUE if control, FALSE if not and 'sample_type' with sample name
#' @param seed Random seed
#' @return 
#' @exportClass 

well2well = function(counts, meta, seed = 42) {
  # basic horiz/vert sort for now
  set.seed(seed)
  
  # append potential horizontal and vertical well orders together
  vert = data.frame('sample_well' = )
  meta_vert = cbind(meta, vert)
  
  horiz = data.frame('sample_well')
  meta_horiz = cbind(meta, horiz)
  
  # create SCRuB objects
  SCRuB_vert = SCRuB(data = counts,
                     metadata = meta_vert)
    
  SCRuB_horiz = SCRuB(data = counts,
                      metadata = meta_vert)
  
  # output
  return(list('Vertical' = data.frame(SCRuB_vert$decontaminated_samples),
              'Horizontal' = data.frame(SCRuB_horiz$decontaminated_samples)))
}

# Function 1: Pipeline 1 - SCRuB + PERfect
## Useful if data contains well information, significant number of negative controls
## within each batch included
## Filters partial features

#' @name pipeline1
#' @usage
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta Data frame with columns is_control, sample_type, and (optional) sample_well
#' @param control_order Vector ordering run of sample_type controls, default NA
#' @param seed Random seed
#' @return List object with original matrix, decontaminated matrix, and corresponding
#' filtering loss (FL) statistics
#' @exportClass list

pipeline1 = function(counts, meta, control_order = NA, seed = 42) {
  
  set.seed(seed)
  
  # SCRuB
  scr_out = SCRuB(counts, 
                  metadata = meta,
                  control_order = control_order)
  sc_counts = data.frame(scr_out$decontaminated_samples)
  
  sc_FL = FL()
  
  # extract FL values from SCRuB data
  
  
  # Create deliverable
  deliv = list(
    SCRuB_FL = sc_FL,
    SCRuB_res = scr_out
  )
  
  return(deliv)
}


# Function 2: Pipeline 2 - Mahoney
## Useful if data does not contain well information, incomplete negative controls,
## multiple batches
## Filters entire features, not partial


#' @name pipeline2
#' @usage 
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta dataframe with columns is_control, sample_type
#' @param blocklist Vector of known previously identified contaminant features
#' @param remove_if Threshold for number of steps feature must be identified as potential contaminant to be removed from final cleaned count matrix. Default set to 1.
#' @param step2_threshold Threshold value for prevalence method of decontam
#' @param technical_replicates Vector identifying technical replicates across batches
#' @param seed Random seed
#' @return List object with original matrix, decontaminated matrix, and 
#' filtering loss (FL)
#' @exportClass list

pipeline2 = function(counts, meta, blocklist, remove_if = 1, step2_threshold = 0.5,
                     technical_replicates, seed = 42) {
  
  set.seed(seed)
  
  # Step 0: W2W check
  w2w = well2well(counts, meta, seed = seed)
  
  # Step 1: Remove features that showed different abundance in different batches
  ## ancombc comparison across batches
  
  s1_res = step1(counts, meta)
  
  # Step 2: Remove features that are differentially abundant in negative controls
  ## decontam
  
  s2_res = step2(counts, meta, step2_threshold)
  
  # Step 3: Remove if DA in diff batches for technical replicates
  
  s3_res = step3(counts, meta, technical_replicates)
  
  # Step 4: Remove known 'blocklist' of contaminants
  
  s4_res = step4(counts, meta, blocklist)
  
  # Prune failed taxa from final object
  
  res = data.frame('step1' = s1_res,
                   'step2' = s2_res,
                   'step3' = s3_res,
                   'step4' = s4_res)
  rownames(res) = colnames(counts)
  
  # return column with summed cases where feature was true
  res$remove = rowSums(res)
  
  # transpose to same as counts matrix
  res = t(res)
  
  # remove features above specified threshold from original counts frame
  counts = rbind(counts, res['remove',])
  rownames(counts)[nrow(counts)] = 'remove'
  final_counts = counts[counts['remove'<remove_if,]==TRUE,]
  
  removed = counts[counts['remove'>remove_if,]==TRUE,]
  
  # determine filtering loss value
  
  FL = FL(counts, removed)
  
  # Create deliverable
  
  deliv = list('results' = t(res),
               'new_counts' = final_counts,
               'filtering loss' = FL)
  
  return(deliv)
}

# Function 2A: Step 1 Pipeline 2

#' @name step1
#' @usage Run step 1 of pipeline 2 to identify features that are differentially abundant
#' between batches
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta Matrix with columns is_control, sample_type, and batch
#' @return  Vector of features tagged as contaminants
#' @exportClass vector

step1 = function(counts, meta) {
  
  phyloseq = wrap_phyloseq(counts, meta)
  
  # run differential analysis
  s1 = ancombc(phyloseq = phyloseq, assay_name = "counts", 
                     group = "batch", p_adj_method = "BH",  lib_cut = 0,
                     formula = "batch", 
                     struc_zero = TRUE, neg_lb = FALSE,
                     tol = 1e-5, max_iter = 100, conserve = FALSE,
                     alpha = 0.05, global = TRUE) 
  # create results matrix
  s1_res = do.call(cbind, s1$res)
  
  # identify column for diff results
  col = ncol(s1_res)
  
  # return indices for which differentially abundant across batches
  ind = which(s1_res[,col]==TRUE)
  
  # return list of tagged contaminant features
  return(s1_res[ind,1])
}

# Function 2B: Step 2 Pipeline 2

#' @name step2
#' @usage Run step 2 of pipeline 2 to identify features that are expressed higher in negative
#' controls and lower in samples
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta Matrix with columns is_control, sample_type, and batch
#' @param threshold Threshold value for prevalence method of decontam
#' @return Vector of features tagged as contaminants
#' @exportClass vector

step2 = function(counts, meta, threshold) {

  # subset to only batches that contain negative controls
  
  # create phyloseq object
  phyloseq =  wrap_phyloseq(counts, meta)
  
  # run decontam prevalence method
  s2_res = isContaminant(phyloseq, method="prevalence", neg="is_control", threshold=threshold)
  
  # return indices for which features identified as contaminant by decontam prevalence method
  ind = which(s2_res$contaminant)
  
  # return list of tagged contaminant features  
  return(rownames(s2_res[ind,]))

}

# Function 2C: Step 3 Pipeline 2

#' @name step3
#' @usage Run step 3 of pipeline 2 to identify features with different abundance in technical
#' replicates across batches
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta Matrix with columns is_control, sample_type, and batch
#' @param technical_replicates Matrix identifying technical replicates across batches with batch as column and rows matching replicates
#' @return Vector of features tagged as contaminants
#' @exportClass vector

step3 = function(counts, meta, technical_replicates) {
  
  # create phyloseq object
  phyloseq =  wrap_phyloseq(counts, meta)
  
  # wrap dataframes for technical replicates
  for (i in 1:dplyr::n_distinct(meta$batch)) {
    
  }
  
  ################################## start at line 353 from original_pipeline2.R
  # genus_new_sh = prune_samples(rs$Batch1, genus)
  # genus_old_sh = prune_samples(rs$Batch2, genus)
  # 
  # genus_new.df = otu_table(genus_new_sh) # Extract OTU table from phyloseq object
  # genus_new.df.PA = transform_sample_counts(genus_new.df, function(abund) 1*(abund>0)) # Set as present/absence
  # 
  # genus_old.df = otu_table(genus_old_sh) # Extract OTU table from phyloseq object
  # genus_old.df.PA = transform_sample_counts(genus_old.df, function(abund) 1*(abund>0))
  
  # get kappa statistic for all technical replicates shared btwn extraction batches
  kappa_results = data.frame(value = numeric(), statistic = numeric(), p.value = numeric())
  
}

# Function 2D: Step 4 Pipeline 2

#' @name step4
#' @usage Run step 4 of pipeline 2 to identify features previously identified as contaminants
#' based on blocklist
#' 
#' @param blocklist Vector of known previously identified contaminant features
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta Matrix with columns is_control, sample_type, and batch
#' @return Vector of features tagged as contaminants
#' @exportClass vector

step4 = function(counts, meta, blocklist) {
  allTaxa = colnames(counts)
  
  return(allTaxa[(allTaxa %in% blocklist)])
}


# Function 3: Visualize pipeline results
## Venn diagram pipeline2, comparison of PERfect results both

#' @name visualize_pipeline
#' @usage Visualize results from the pipelines within the package
#'
#' @param pipeline_output Output of pipeline object
#' @param interactive TRUE if user wants interactive plot output
#' @return Visualizations relating to pipeline object
#' @exportClass list

visualize_pipeline = function(pipeline_output, interactive = FALSE)  {
  
  if (pipeline1) {

  }
  
  if (pipeline2) {
    # Venn comparison of contaminant taxa removed across steps
    x = list(s1_rem, s2_rem, s3_rem, s4_rem)
    ggVennDiagram(x, stroke.size =1,
                  category.names = c("Step 1", "Step 2", "Step 3", "Step 4"),
                  edge_lty = "solid", set_size = 6,
                  label_alpha = 0.5, label_percent_digit = 1) +
      scale_x_continuous(expand = expansion(mult = .2)) +
      ggplot2::scale_color_grey(start=0, end=0) +
      scale_fill_distiller(direction=1) +
      labs(title="Taxa Removal by Step") +
      theme(legend.position="none", plot.title=element_text(size=25, hjust = 0.5)) +
      scale_fill_distiller(palette = "Spectral")
  }
}

# Function 4: "Unwrap" phyloseq or SummarizedExperiment object

#' @name unwrap_phyloseq
#' @usage 
#'
#' @param phyloseq Phyloseq object to unwrap into required count and meta matrices for pipeline functions
#' @return List containing counts and metadata data frames for input into pipeline functions
#' @exportClass list

unwrap_phyloseq = function(phyloseq) {
  counts = data.frame(t(phyloseq@otu_table)) # requires data frame first, will not coerce to matrix from phyloseq object
  meta = data.frame(phyloseq@sam_data)
  
  return(list(counts = as.matrix(counts), # adjust to matrix for returnable
              meta = as.matrix(meta)))
}

# Function 5: "Wrap" phyloseq or SummarizedExperiment object

#' @name wrap_phyloseq
#' @usage 
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta dataframe with columns is_control, sample_type
#' @return List containing counts and metadata data frames for input into pipeline functions
#' @exportClass phyloseq

wrap_phyloseq = function(counts, meta) {
  counts = t(counts) # transpose to fit with expectation of phyloseq object
  OTU = otu_table(counts, taxa_are_rows = TRUE)
  META = sample_data(meta)
  
  tax_mat = matrix(rownames(counts),nrow=nrow(counts),ncol=1)
  rownames(tax_mat) = rownames(counts)
  TAX = tax_table(tax_mat)
  
  return(phyloseq(OTU, META, TAX))
}

#' @name NP_Order
#' @usage (adjusted from katiasmirn/PERfect/FiltLoss.R)
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @return Data frame with features as row names and associated filtering loss value
#' @exportClass data.frame

NP_Order = function(counts){
  #arrange counts in order of increasing number of samples taxa are present in
  NP = names(sort(apply(counts, 2, Matrix::nnzero)))
  return(NP)
}

# Function 6: Filtering loss statistic 

#' @name FL
#' @usage Determine the filtering loss (FL) for count
#' data based on removed features (adjusted from katiasmirn/PERfect/FiltLoss.R)
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param removed Vector of features to be removed as contaminants
#' @return Data frame with features as row names and associated filtering loss value
#' @exportClass data.frame

FL = function(counts, removed){
  #X - data matrix
  #Ind - only jth taxon removed to calculate FL
  p = dim(counts)[2]#total #of taxa
  Norm_Ratio = rep(1, p)
  
  # Check the format of X
  if(!(class(counts) %in% c("matrix"))){counts = as.matrix(counts)}
  
  #Order columns by importance
  Order.vec = NP_Order(counts)

  counts = counts[,Order.vec] #properly order columns of X
  
  Order_Ind = seq_len(length(Order.vec))
  Netw = t(counts)%*%counts
  
  #Taxa at the top of the list have smallest number of connected nodes
  for (i in seq_len(p)){
    Ind = Order_Ind[-i]
    
    #define matrix X_{-J}'X_{-J} for individual filtering loss
    Netw_R = Netw[Ind, Ind]
    
    Norm_Ratio[i] =  sum(Netw_R*Netw_R)
  }
  
  FL = 1 - Norm_Ratio/sum(Netw*Netw)
  
  names(FL) =  colnames(counts)

  FL_value = sum(FL[names(FL) %in% removed])

  return(FL_value)
  
  # return(FL = data.frame(FL))
}

# ? Function ?: Shiny visualization for output
## Later step once others are finished


