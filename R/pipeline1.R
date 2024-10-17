#' pipeline1
#'
#' Using the SCRuB method, this pipeline is best used when the user wants
#' to closely characterize the composition of the data prior to contamination.
#'
#' @family pipeline1
#'
#' @import tidyverse
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta Data frame with columns is_control, sample_type, and (optional) sample_well
#' @param control_name Character name of controls within metadata
#' @param control_order Vector ordering run of sample_type controls, default NA
#' @param seed Random seed
#' @return List object with original matrix, decontaminated matrix, and corresponding
#' filtering loss (FL) statistics
#' @export

pipeline1 = function(counts, meta, control_name, control_order = NA, seed = 42) {

  if (sum(colnames(meta) %in% 'batch')==0) {
    meta$batch = rep('a')
  }

  # check to ensure each batch contains some controls
  for(b in meta$batch) {
    index = meta %>% filter(batch == b) %>% row.names() # select within batch

    if (sum(meta[index, 'is_control']==TRUE) == 0) {
      stop(paste0('To use pipeline1, all batches must contain controls. Batch named  ',
                  b, ' does not contain samples specified as controls.', sep = '')) # break loop if missing
    }
  }

  if (sum((colnames(meta) %in% 'sample_well'))==0) {
    well2well(counts, meta, control_name = control_name, seed = seed)
  }

  set.seed(seed)

  batch = unique(meta$batch)

  # SCRuB
  sc_outs = list() # get count matrices from batches scrubbed separately
  for(b in batch) {
    index = meta %>% filter(batch == b) %>% row.names() # select only of one batch

    sc_outs[[b]] = SCRuB::SCRuB(counts[index,],
                                meta[index,] %>%
                                  select(any_of(c('is_control', 'sample_type', 'sample_well'))))$decontaminated_samples
  }

  sc_counts = do.call(rbind, sc_outs) # append batches back together

  # Identify filtering loss
    # remove row of blank sample from counts for FL comparison
  bl = rownames(meta)[meta$is_control==TRUE]
  counts_nobl = counts[!rownames(counts) %in% bl,]

  sc_FL = FL(counts_nobl, new_counts = sc_counts)


  # Create deliverable
  return(list('decontaminated_count' = sc_counts,
              'filtering_loss' = sc_FL,
              'pipeline' = 'pipeline1')
  )
}
