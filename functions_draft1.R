# Function development for microbiome decontamination pipeline (micRoclean)
# 20240305
# Rachel Griffard

# Required libraries
library(phyloseq) # object wrapper
library(SummarizedExperiment) # object wrapper
library(tidyverse)
library(plotly) # for interactive feature
library(SCRuB) # well2well, pipeline 1
library(decontam) # pipeline 2 step2
library(microDecon) # pipeline 2
library(ANCOMBC) # pipeline 2
library(ggVennDiagram) # function 3 - pipeline 2 - comparison across removed
library(shiny) # function 3
library(ANCOMBC) # pipeline 2 step1
library(irr) # pipeline 2 step3


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
  
  # plate wells
  well = data.frame()
  
  for (i in 1:8) { # rows
    row = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
    
    for (j in 1:12) { # columns
      well[i,j] = paste0(row[i], j, sep = '')
    }
  }
  
  # string for well assignments
  vert = unname(unlist(well))
  horiz = unname(unlist(data.frame(t(well))))
  
  # append potential horizontal and vertical well orders together

    # order batches based on naming convention (number in the end of the string)
  meta[order(as.numeric(sub(".*[^0-9](\\d+)$", "\\1", rownames(meta)))),]

    # restart at each batch (different plates)
  num_b = table(meta$batch)
  
  sample_well = c(vert[1:num_b[1]], vert[1:num_b[2]])
  meta_vert = cbind(meta, sample_well)
  meta_vert = subset(meta_vert, select = c(is_control, sample_type, sample_well))
  
  sample_well = c(horiz[1:num_b[1]], horiz[1:num_b[2]])
  meta_horiz = cbind(meta, sample_well)
  meta_horiz = subset(meta_horiz, select = c(is_control, sample_type, sample_well))
  
  # create SCRuB objects
  SCRuB_vert = SCRuB(counts,
                     meta_vert)
    
  SCRuB_horiz = SCRuB(counts,
                      meta_horiz)
  
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
    deliv = list('contaminant_id' = res,
                 'decontaminated_count' = final_counts,
                 'removed' = removed,
                 'filtering_loss' = FL,
                 'pipeline' = 'pipeline1')
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
#' @return List object with contaminant matrix, decontaminated count matrix, character string of all removed contaminants,
#' and filtering loss (FL) statistic
#' @exportClass list

pipeline2 = function(counts, meta, blocklist, technical_replicates, remove_if = 1,
                     step2_threshold = 0.5, seed = 42) {
  
  set.seed(seed)
  
  # Step 0: W2W check
#  w2w = well2well(counts, meta, seed = seed)
  
  # Step 1: Remove features that showed different abundance in different batches
  ## ancombc comparison across batches
  
  s1_res = step1(counts, meta)
  
  # Step 2: Remove features that are differentially abundant in negative controls
  ## decontam
  
  s2_res = step2(counts, meta, step2_threshold)
  
  # Step 3: Remove if DA in diff batches for technical replicates
  
  s3_res = step3(counts, technical_replicates)
  
  # Step 4: Remove known 'blocklist' of contaminants
  
  s4_res = step4(counts, blocklist)
  
  # Create dataframe indicating TRUE if contaminant and FALSE if not tagged
  res = data.frame('feature' = colnames(counts))
  
  res = data.frame('feature' = colnames(counts),
                   'step1' = ifelse(colnames(counts) %in% s1_res, TRUE, FALSE),
                   'step2' = ifelse(colnames(counts) %in% s2_res, TRUE, FALSE),
                #   'step3' = ifelse(colnames(counts) %in% s3_res, TRUE, FALSE),
                   'step4' = ifelse(colnames(counts) %in% s4_res, TRUE, FALSE))
  
  # return column with summed cases where feature was true
  
  # transpose to same as counts matrix
  res2 = res
  res2$remove = rowSums(res2[,-c(1)])
  res2 = t(res2)
  
  # remove features above specified threshold from original counts frame
  counts_rem = rbind(counts, res2['remove',])
  rownames(counts_rem)[nrow(counts_rem)] = 'remove'
  final_counts = counts[,counts_rem['remove',]<remove_if]
  
  removed = setdiff(colnames(counts), colnames(final_counts))
  
  # determine filtering loss value
  
  FL = FL(counts, removed)
  
  # Create deliverable
  
  deliv = list('contaminant_id' = res,
               'decontaminated_count' = final_counts,
               'removed' = removed,
               'filtering_loss' = FL,
               'pipeline' = 'pipeline2')
  
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
  suppressWarnings({
  s1 = ancombc(phyloseq = phyloseq, assay_name = "counts", 
                     group = "batch", p_adj_method = "BH", lib_cut = 0,
                     formula = "batch", 
                     struc_zero = TRUE, neg_lb = FALSE,
                     tol = 1e-5, max_iter = 100, conserve = FALSE,
                     alpha = 0.05, global = TRUE) 
  })
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

step3 = function(counts, technical_replicates) {
  
  # wrap dataframe for technical replicates in each batch ordered by match (line ~336 original_pipeline2.R)
  count_replicate = list()
  PA_replicate = list()
  
  for (i in 1:dplyr::n_distinct(meta$batch)) {
    vals = technical_replicates[,i]
    count_replicate[[i]] = data.frame(t(counts[vals,]))
    
    # create presence absence matrices
    j = as.matrix(count_replicate[[i]])
    j[j!=0] = 1
    PA_replicate[[i]] = j
  }
  
  # create dataframe to contain results from IRR kappa
  kappa_results = data.frame(value = numeric(), statistic = numeric(), p.value = numeric())
  
  # get kappa values using for loop
  batch1.df.PA = PA_replicate[[1]]
  batch2.df.PA = PA_replicate[[2]]
  
  for (i in nrow(PA_replicate[[1]])) {
    k = kappa2(t(rbind(batch1.df.PA[i,], batch2.df.PA[i,])), "unweighted")
    kappa_results[i,"value"] = k$value 
    kappa_results[i,"statistic"] = k$statistic
    kappa_results[i,"p.value"] = k$p.value
  }
  
  row.names(kappa_results) = colnames(counts)
  
  kappa_results.no_NA = subset(kappa_results, !is.na(value) & !is.na(p.value))
  kappa_res_remove = subset(kappa_results.no_NA, p.value >= 0.05 | value <= 0.4)
  
  return(rownames(kappa_res_remove))
}

# Function 2D: Step 4 Pipeline 2

#' @name step4
#' @usage Run step 4 of pipeline 2 to identify features previously identified as contaminants
#' based on blocklist
#' 
#' @param blocklist List of known previously identified contaminant features
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta Matrix with columns is_control, sample_type, and batch
#' @return Vector of features tagged as contaminants
#' @exportClass vector

step4 = function(counts, blocklist) {
  allTaxa = colnames(counts)
  
  if (class(blocklist)!='character') {
    warning('blocklist parameter must be formatted as a character vector')
    break
  }
  
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
  
  if (pipeline_output$pipeline == 'pipeline1') {

  }
  
  if (pipeline_output$pipeline == 'pipeline2') {
    # import data
    s1_rem = pipeline_output$contaminant_id$feature[pipeline_output$contaminant_id$step1==TRUE]
    s2_rem = pipeline_output$contaminant_id$feature[pipeline_output$contaminant_id$step2==TRUE]
    s3_rem = pipeline_output$contaminant_id$feature[pipeline_output$contaminant_id$step3==TRUE]
    s4_rem = pipeline_output$contaminant_id$feature[pipeline_output$contaminant_id$step4==TRUE]
    # Venn comparison of contaminant taxa removed across steps
    x = list(s1_rem, s2_rem, s3_rem, s4_rem)
    p = ggVennDiagram(x, stroke.size =1,
                  category.names = c("Step 1", "Step 2", "Step 3", "Step 4"),
                  edge_lty = "solid", set_size = 6,
                  label_alpha = 0.5, label_percent_digit = 1) +
      scale_x_continuous(expand = expansion(mult = .2)) +
      ggplot2::scale_color_grey(start=0, end=0) +
      scale_fill_distiller(direction=1) +
      labs(title="Taxa Removal by Step") +
      theme(legend.position="none", plot.title=element_text(size=25, hjust = 0.5)) +
      scale_fill_distiller(palette = "Spectral")
    
    if (interactive == FALSE) {
      return(p)
    }
    
    if (interactive == TRUE) {
      return(plotly::ggplotly(p))
    }
    
    else {
      warning('interactive must be set to TRUE or FALSE.')
    }
  }
  
  else {
    warning('Rerun data through pipeline and ensure object in visualize_pipeline is output from pipeline1 or pipeline2.')
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
#' data based on removed features (adjusted from katiasmirn/PERfect/FL_J.R)
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param removed Vector of features to be removed as contaminants
#' @return Data frame with features as row names and associated filtering loss value
#' @exportClass data.frame

FL = function(counts, removed){
  
  # Check the format of J
  if(class(removed) != "character")
    stop('removed argument must be a character vector containing names of taxa to be removed')
  
  Ind <- which(colnames(counts) %in%  removed)
  X_R <- counts[,-Ind]
  #calculate corresponding norm
  Netw <- t(as.matrix(counts))%*%as.matrix(counts)
  Netw_R <- t(as.matrix(X_R))%*%as.matrix(X_R)
  #FL <-  1 - (psych::tr(t(Netw_R)%*%Netw_R)/psych::tr(t(Netw)%*%Netw))
  FL <-  1 - (sum(Netw_R*Netw_R)/sum(Netw*Netw))
  return(FL)
}

# ? Function ?: Shiny visualization for output
## Later step once others are finished


