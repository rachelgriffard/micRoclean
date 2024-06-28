#' well2well
#'
#' @name well2well
#' @usage Determine the potential impact of spatial well to well relationships
#'
#' @import dplyr
#' @importFrom SCRuB SCRuB
#' @importFrom vegan vegdist mantel
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta Metadata with column one 'is_control' indicating TRUE if control, FALSE if not and 'sample_type' with sample name
#' @param seed Random seed
#' @export

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
  vert = unname(unlist(well)) # vertical alignment
  horiz = unname(unlist(data.frame(t(well)))) # horizontal alignment

  # order samples by name convention
  meta = meta %>%
    arrange(batch, as.numeric(str_extract(rownames(meta), "\\d+$")))

  # append potential horizontal and vertical well orders together

  # order batches based on naming convention (number in the end of the string)
  meta[order(as.numeric(sub(".*[^0-9](\\d+)$", "\\1", rownames(meta)))),]

  # restart at each batch (different plates)
  num_b = table(meta$batch)

  sample_well = c(vert[1:num_b[1]], vert[1:num_b[2]])
  meta_vert = cbind(meta, sample_well)

  sample_well = c(horiz[1:num_b[1]], horiz[1:num_b[2]])
  meta_horiz = cbind(meta, sample_well)

  # order counts by name convention for SCRuB function
  counts = as.data.frame(counts) %>%
    add_column(meta$batch) %>%
    arrange(`meta$batch`, as.numeric(str_extract(rownames(counts), "\\d+$"))) %>%
    mutate(`meta$batch` = NULL)

  # create SCRuB objects by batch
  sc_outs_vert = list()
  for(b in unique(meta_vert$batch)) {
    index = meta_vert %>% filter(batch == b) %>% row.names()

    try((sc_outs_vert[[b]] = SCRuB(counts[index,],
                                   meta_vert[index,] %>%
                                     select(is_control, sample_type, sample_well))$p), silent = TRUE)
  }

  sc_outs_horiz = list()
  for(b in unique(meta_horiz$batch)) {
    index = meta_horiz %>% filter(batch == b) %>% row.names()

    try((sc_outs_horiz[[b]] = SCRuB(counts[index,],
                                    meta_horiz[index,] %>%
                                      select(is_control, sample_type, sample_well))$p), silent = TRUE)
  }

  sc_outs = list()
  for(b in unique(meta$batch)) {
    index = meta %>% filter(batch == b) %>% row.names()

    try((sc_outs[[b]] = SCRuB(counts[index,],
                              meta[index,] %>%
                                select(is_control, sample_type))$p), silent = TRUE)
  }

  # append batches back together for full, decontaminated dataframe
  SCRuB_vert = unlist(sc_outs_vert) # append batches back together

  SCRuB_horiz = unlist(sc_outs_horiz) # append batches back together

  SCRuB = unlist(sc_outs) # append batches back together


  # determine if vert/horiz significantly different from without spatial
  res_vert = irr::icc(as.matrix(cbind(SCRuB_vert, SCRuB)))

  ifelse(res_vert$p<0.01, warning('Strong evidence of well to well contamination. User is encouraged to rerun pipeline1 with well location information.'), NA)

  res_horiz = irr::icc(as.matrix(cbind(SCRuB_horiz, SCRuB)))

  ifelse(res_horiz$p<0.01, warning('Strong evidence of well to well contamination. User is encouraged to rerun pipeline1 with well location information.'), return(paste0('No strong evidence of well to well contamination.')))
}
