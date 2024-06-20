#' @name pipeline1
#' @usage Using the SCRuB method, this pipeline is best used when the user wants
#' to closely characterize the composition of the data prior to contamination.
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
  scr_out = SCRuB::SCRuB(counts, 
                         metadata = meta[,-3],
                         control_order = control_order)
  sc_counts = data.frame(scr_out$decontaminated_samples)
  
  sc_FL = FL(counts, new_counts = sc_counts)
  
  # extract FL values from SCRuB data
  
  
  # Create deliverable
  return(list('decontaminated_count' = sc_counts,
              'filtering_loss' = sc_FL,
              'pipeline' = 'pipeline1')
  )
}
