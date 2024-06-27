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
#' @exportClass list

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


#' FL
#' Determine the filtering loss (FL) for count
#' data based on removed features (adjusted from katiasmirn/PERfect/FL_J.R)
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param new_counts Count matrix with samples as rows and features as counts after being partially filtered by SCRuB method
#' @param removed Vector of features to be removed as contaminants
#' @return Data frame with features as row names and associated filtering loss value
#' @exportClass data.frame

FL = function(counts, new_counts = NULL, removed = NULL){

  # for pipeline 1
  if (is.null(new_counts) == FALSE) { #pipeline 1
    X_R = new_counts
  }

  # for pipeline 2
  if (is.null(removed) == FALSE) {

    # Check the format of removed vector
    if(class(removed) != "character")
      stop('removed argument must be a character vector containing names of taxa to be removed')

    Ind = which(colnames(counts) %in%  removed)
    X_R = counts[,-Ind]
  }

  else {
    warning('Ensure the correct new_counts or removed values are input for this method')
  }

  #calculate corresponding norm
  Netw = t(as.matrix(counts))%*%as.matrix(counts)
  Netw_R = t(as.matrix(X_R))%*%as.matrix(X_R)

  FL =  1 - (sum(Netw_R*Netw_R)/sum(Netw*Netw))
  return(FL)
}
