#' micRoclean
#' Primary function for decontamination low biomass microbiome data
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta Data frame with columns is_control, sample_type, and (optional) sample_well
#' @param research_goal Define primary research goal. Options: 'biomarker' or 'orig.composition'
#' @param control_name Character name of controls within metadata
#' @param control_order Vector ordering run of sample_type controls, default NA
#' @param blocklist Vector of known previously identified contaminant features
#' @param remove_if Threshold for number of steps feature must be identified as potential contaminant to be removed from final cleaned count matrix. Default set to 1.
#' @param step2_threshold Threshold value for prevalence method of decontam
#' @param technical_replicates Vector identifying technical replicates across batches
#' @param seed Random seed
#'
#' @return List object with original matrix, decontaminated matrix, character string of all removed contaminants (if 'biomarker'), and corresponding
#' filtering loss (FL) statistics
#' @export

micRoclean = function(counts, meta, research_goal, control_name, control_order = NA, blocklist, remove_if = 1, step2_threshold = 0.5, technical_replicates, seed = 42) {
  if (research_goal == 'orig.composition') {
     res = pipeline1(counts, meta, control_name, control_order, seed)
     return(res)
  }

  if (research_goal == 'biomarker') {
     res = pipeline2(counts, meta, blocklist, control_name, technical_replicates, remove_if,
                     step2_threshold, seed)
     return(res)
  }
  else {
    stop(paste0('research_goal must be: biomarker or orig.composition'))
  }
}

#' unwrap_phyloseq
#'
#' Take input phyloseq object and "unwrap" data into matrices expected as input for functions in
#' micRoclean package.
#'
#' @family management
#'
#' @import phyloseq
#' @param phyloseq Phyloseq object to unwrap into required count and meta matrices for pipeline functions
#' @return List containing counts and metadata data frames for input into pipeline functions
#' @export

unwrap_phyloseq = function(phyloseq) {
  counts = data.frame(t(phyloseq@otu_table)) # requires data frame first, will not coerce to matrix from phyloseq object
  meta = data.frame(phyloseq@sam_data)

  return(list(counts = as.matrix(counts), # adjust to matrix for returnable
              meta = as.matrix(meta)))
}

#' wrap_phyloseq
#'
#' Wrap count matrix and metadata into a phyloseq object
#'
#' @family management
#'
#' @import phyloseq
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta dataframe with columns is_control, sample_type
#' @return List containing counts and metadata data frames for input into pipeline functions
#' @export

wrap_phyloseq = function(counts, meta) {
  counts = t(counts) # transpose to fit with expectation of phyloseq object
  OTU = otu_table(counts, taxa_are_rows = TRUE)
  META = sample_data(meta)

  tax_mat = matrix(rownames(counts),nrow=nrow(counts),ncol=1)
  rownames(tax_mat) = rownames(counts)
  TAX = tax_table(tax_mat)

  return(phyloseq(OTU, META, TAX))
}


#' FL
#' Determine the filtering loss (FL) for count
#' data based on removed features (adjusted from katiasmirn/PERfect/FL_J.R)
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param new_counts Count matrix with samples as rows and features as counts after being partially filtered by SCRuB method
#' @param removed Vector of features to be removed as contaminants
#' @return Data frame with features as row names and associated filtering loss value
#' @export

FL = function(counts, new_counts = NULL){

  X_R = new_counts

  #calculate corresponding norm
  Netw = t(as.matrix(counts))%*%as.matrix(counts)
  Netw_R = t(as.matrix(X_R))%*%as.matrix(X_R)

  FL =  1-(norm(Netw_R), "F")^2/norm(Netw), "F")^2)
  return(FL)
}
