#' @name unwrap_phyloseq
#' @usage Take input phyloseq object and "unwrap" data into matrices expected as input for functions in 
#' micRoclean package.
#' 
#' @family management
#'
#' @import phyloseq
#' @param phyloseq Phyloseq object to unwrap into required count and meta matrices for pipeline functions
#' @return List containing counts and metadata data frames for input into pipeline functions
#' @exportClass list

unwrap_phyloseq = function(phyloseq) {
  counts = data.frame(t(phyloseq@otu_table)) # requires data frame first, will not coerce to matrix from phyloseq object
  meta = data.frame(phyloseq@sam_data)
  
  return(list(counts = as.matrix(counts), # adjust to matrix for returnable
              meta = as.matrix(meta)))
}

#' @name wrap_phyloseq
#' @usage 
#' 
#' @family management
#'
#' @import phyloseq
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
#' @usage Order taxa by counts, required inside FL function (adjusted from katiasmirn/PERfect/FiltLoss.R)
#' 
#' @family management
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @return Data frame with features as row names and associated filtering loss value
#' @exportClass data.frame

NP_Order = function(counts){
  #arrange counts in order of increasing number of samples taxa are present in
  NP = names(sort(apply(counts, 2, Matrix::nnzero)))
  return(NP)
}

#' @name FL
#' @usage Determine the filtering loss (FL) for count
#' data based on removed features (adjusted from katiasmirn/PERfect/FL_J.R)
#' 
#' @family management
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param removed Vector of features to be removed as contaminants
#' @return Data frame with features as row names and associated filtering loss value
#' @exportClass data.frame

FL = function(counts, removed){
  
  # Check the format of J
  if(class(removed) != "character")
    stop('removed argument must be a character vector containing names of taxa to be removed')
  
  Ind =which(colnames(counts) %in%  removed)
  X_R =counts[,-Ind]
  #calculate corresponding norm
  Netw =t(as.matrix(counts))%*%as.matrix(counts)
  Netw_R =t(as.matrix(X_R))%*%as.matrix(X_R)
  
  FL = 1 - (sum(Netw_R*Netw_R)/sum(Netw*Netw))
  return(FL)
}