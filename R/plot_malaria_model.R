#' @importFrom utils globalVariables
utils::globalVariables(c("time", "population", "count"))

#' Plot Malaria Model Results
#'
#' This function takes the output of the malaria model and creates a plot showing the dynamics of each population over time.
#'
#' @param results The output from the `run` method of the malaria model.
#' @param colors A vector of colors for the population lines. If NULL, default colors will be used.
#' @param line_width The width of the population lines.
#' @param legend_position The position of the legend ("top", "bottom", "left", "right", or "none").
#' @export
#' @import ggplot2
#' @import tidyr
#' @import scales
#' @import dplyr
plot_malaria_model <- function(results, colors = NULL, line_width = 1, legend_position = "right") {
  results_df <- as.data.frame(results)
  colnames(results_df)[1] <- "time"  # Rename the first column to "time"
  results_long <- tidyr::pivot_longer(results_df, -time, names_to = "population", values_to = "count")

  # Ensure variables are defined for ggplot
  if (is.null(colors)) {
    colors <- scales::hue_pal()(length(unique(results_long$population)))
  }

  legend_position <- match.arg(legend_position, c("top", "bottom", "left", "right", "none"))

  # Create plot
  p <- ggplot(results_long, aes(x = time, y = count, color = population)) +
    geom_line(linewidth = line_width) +
    scale_color_manual(values = colors) +
    labs(x = "Time", y = "Proportion", title = "Malaria Model Dynamics") +
    theme_minimal()

  if (legend_position != "none") {
    p <- p + theme(legend.position = legend_position)
  } else {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}
