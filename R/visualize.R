#' visualize_pipeline

#' Visualize results from the pipelines within the package
#'
#' @param pipeline_output Output of pipeline object
#' @param interactive TRUE if user wants interactive plot output
#' @return Visualizations relating to pipeline object
#' @export

visualize_pipeline = function(pipeline_output, interactive = FALSE)  {

  if (pipeline_output$pipeline == 'pipeline1') {
    warning('micRoclean does not currently have functionality to visualize SCRuB decontamination')
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
