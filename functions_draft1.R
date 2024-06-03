# Function development for microbiome decontamination pipeline
# 20240305
# Rachel Griffard

# Required libraries throughout package
library(phyloseq) # object wrapper
library(SummarizedExperiment) # object wrapper
library(tidyverse)
library(plotly) # for interactive feature
library(SCRuB) # pipeline 1
library(PERfect) # DFL
library(decontam) # pipeline 2
library(microDecon) # pipeline 2
library(ANCOMBC) # pipeline 2
library(ggVennDiagram) # function 3 - pipeline 2 - comparison across removed
library(shiny) # function 3


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
#' @return List object with original matrix, decontaminated matrix, and 
#' difference in filtering loss (DFL) statistics for both
#' @exportClass list

pipeline1 = function(counts, meta, control_order = NA, seed = 42) {
  
  set.seed(seed)
  
  # PERfect implementation before SCRuBbing
  pre_sim = PERFect_sim(X = counts)
  
  # SCRuB
  scr_out = SCRuB(counts, 
                  metadata = meta,
                  control_order = control_order)
  sc_counts = data.frame(scr_out$decontaminated_samples)
  
  # Extract DFL values ordered by contribution
  pre_DFL = DFL(counts)
  post_DFL = DFL(sc_counts) # post-scrub
  
  
  # Create deliverable
  deliv = list(
    preSCRuB_DFL = DFL,
    postSCRuB_DFL = sc_DFL,
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
#' @param blocklist Vector of known contaminant features to remove
#' @param remove_if 
#' @param step2 Method for step 2 ('decontam' or 'microDecon')
#' @param step2_threshold 
#' @param seed Random seed
#' @return List object with original matrix, decontaminated matrix, and 
#' difference in filtering loss (DFL) statistics for both
#' @exportClass list

pipeline2 = function(counts, meta, blocklist, remove_if = 1, step2 = 'decontam', step2_threshold = 0.5, seed = 42) {
  
  set.seed(seed)
  
  # Step 0: W2W check
  ## vert and horiz output
  w2w = well2well(counts, meta, seed = seed)
  
  sim_0 = FL(counts)
  
  
  # Step 1: Remove features that showed different abundance in different batches
  ## ancombc comparison across batches
  
  s1_res = step1(counts, meta)
  
  # Step 2: Remove features that are differentially abund in negative controls
  ## decontam OR microDecon
  
  s2_res = step2(counts, meta, step2, step2_threshold)
  
  # Step 3: Remove if DA in diff batches for technical replicates
  
  s3_res = step3(counts, meta)
  
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
#' @usage 
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta dataframe with columns is_control, sample_type
#' @return List object with original matrix, decontaminated matrix, and 
#' difference in filtering loss (DFL) statistics for both
#' @exportClass data.frame

step1 = function(counts, meta) {
  
  micro = wrap_phyloseq(counts, meta)
  
  micro_s1 = ancombc(phyloseq = micro, assay_name = "counts", 
                     group = "batch", p_adj_method = "BH",  lib_cut = 0,
                     formula = "batch", 
                     struc_zero = TRUE, neg_lb = FALSE,
                     tol = 1e-5, max_iter = 100, conserve = FALSE,
                     alpha = 0.05, global = TRUE) 
  
  micro_s1_res = do.call(cbind, micro_s1$res)
  
  micro_s1_res[ncol(micro_s1_res)-2==TRUE,] # check this / add eval statement to determine if still true
  
  s1_rem = rownames(micro_s1_res[micro_s1_res$`diff_abn.batch2. New`==TRUE,])
  
  index3 = grep(paste(s1_rem,collapse="$|"), colnames(counts))
  counts_s1 = counts[,-index3]
}

# Function 2B: Step 2 Pipeline 2

#' @name step2
#' @usage 
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta dataframe with columns is_control, sample_type
#' @param method Method for step 2 ('decontam' or 'microDecon')
#' @return List object with original matrix, decontaminated matrix, and 
#' difference in filtering loss (DFL) statistics for both
#' @exportClass data.frame

step2 = function(counts, meta, method, threshold = 0.5) {
  ## Option 1: decontam
  if (step2 == 'decontam') {
    # subset to only batches that contain negative controls
    
    decontam_output = isContaminant(phyloseq, method="prevalence", neg="is.neg", threshold=threshold)
    index4 = which(decontam_output$contaminant)
    
    micro_s2 = rownames(decontam_output[decontam_output$contaminant==TRUE,])
  }
  
  ## Option 2: microDecon
  if (step2 == 'microDecon') {
    
  }
  
  ## Option 3: Invalid
  else {
    paste('Must select method for step two: microDecon OR decontam')
    break
  }
}

# Function 2C: Step 3 Pipeline 2

step3 = function(counts, meta) {
  
}

# Function 2D: Step 4 Pipeline 2

step4 = function(counts, meta, blocklist) {
  
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
    DFL = rownames_to_column(DFL, "feature")
    sc_DFL = rownames_to_column(sc_DFL, "feature")
    
    DFL_com = DFL %>%
      left_join(sc_DFL, by = c('feature' = 'feature')) # join all DFL values into data.frame
    
    DFL_pl = DFL_com %>%
      pivot_longer(!feature, names_to = 'Group', values_to = 'DFL') %>%
      mutate(Group = ifelse(Group == 'DFL.x', 'Pre-SCRuB', 'Post-SCRuB')) %>%
      mutate(DFL = ifelse(is.na(DFL)==TRUE, 0, DFL)) %>%
      mutate(order(DFL)) %>%
      ggplot(aes(fct_reorder(feature, DFL), DFL, fill = Group)) +
      geom_bar(stat = 'identity', position = 'dodge2') +
      xlab('Feature') + ggtitle('DFL Comparison Pre- and Post-SCRuB Method') +
      coord_flip()
    
    ifelse(interactive == TRUE, ggplotly(DFL_pl), DFL_pl)
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
#' @param phyloseq Phyloseq object to unwrap into required data frames for pipeline functions
#' @return List containing counts and metadata data frames for input into pipeline functions
#' @exportClass list

unwrap_phyloseq = function(phyloseq) {
  
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
  OTU = otu_table(as.matrix(counts), taxa_are_rows = TRUE)
  META = sample_data(meta)
  
  tax_mat = matrix(rownames(counts),nrow=nrow(counts),ncol=1)
  rownames(tax_mat) = rownames(counts)
  TAX = tax_table(tax_mat)
  
  return(phyloseq(OTU, META, TAX))
}

#' @name NP_Order
#' @usage (adjusted from katiasmirn/PERfect/FiltLoss.R)
#'
#' @param Counts Count matrix with samples as rows and features as counts
#' @return Data frame with features as row names and associated filtering loss value
#' @exportClass data.frame

NP_Order = function(Counts){
  #arrange counts in order of increasing number of samples taxa are present in
  NP = names(sort(apply(Counts, 2, Matrix::nnzero)))
  return(NP)
}

# Function 6: Filtering loss statistic 

#' @name FL
#' @usage Determine the filtering loss (FL) for count
#' data based on removed features (adjusted from katiasmirn/PERfect/FiltLoss.R)
#'
#' @param X Count matrix with samples as rows and features as counts
#' @param removed Vector of features to be removed as contaminants
#' @return Data frame with features as row names and associated filtering loss value
#' @exportClass data.frame

FL = function(X, removed){
  #X - data matrix
  #Ind - only jth taxon removed to calculate FL
  p = dim(X)[2]#total #of taxa
  Norm_Ratio = rep(1, p)
  
  # Check the format of X
  if(!(class(X) %in% c("matrix"))){X = as.matrix(X)}
  
  #Order columns by importance
  Order.vec = NP_Order(X)

  X = X[,Order.vec] #properly order columns of X
  
  Order_Ind = seq_len(length(Order.vec))
  Netw = t(X)%*%X
  
  #Taxa at the top of the list have smallest number of connected nodes
  for (i in seq_len(p)){
    Ind = Order_Ind[-i]
    
    #define matrix X_{-J}'X_{-J} for individual filtering loss
    Netw_R = Netw[Ind, Ind]
    
    Norm_Ratio[i] =  sum(Netw_R*Netw_R)
  }
  
  FL = 1 - Norm_Ratio/sum(Netw*Netw)
  
  names(FL) =  colnames(X)

  FL_value = sum(FL[rownames(FL) %in% removed,])

  return(FL_value)
  
  # return(FL = data.frame(FL))
}

# ? Function ?: Shiny visualization for output
## Later step once others are finished


