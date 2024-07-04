#' @import ggplot2
#' @importFrom data.table as.data.table :=
NULL

#' Plot Range Model Results
#'
#' @param model_results A data frame containing the model results from `range_model`.
#' @param x_var A string specifying the variable to be used for the x-axis.
#' @param y_var A string specifying the variable to be used for the y-axis.
#' @param color_var A string specifying the variable to be used for color classification.
#' @param output_vars A character vector specifying the output variables to plot. If NULL, all output variables will be plotted.
#' @param colors A named character vector specifying the colors for each variable. If NULL, default colors will be used.
#' @param line_types A named character vector specifying the line types for each variable. If NULL, default line types will be used.
#' @param line_widths A named numeric vector specifying the line widths for each variable. If NULL, default line widths will be used.
#' @param title The title of the plot.
#' @param subtitle The subtitle of the plot.
#'
#' @return A ggplot object.
#' @export
plot_range <- function(model_results, x_var, y_var, color_var, output_vars = NULL, colors = NULL, line_types = NULL,
                       line_widths = NULL, title = "Range Model Results", subtitle = NULL) {

  dt <- data.table::as.data.table(model_results)

  if (is.null(output_vars)) {
    output_vars <- setdiff(names(dt), c("t", names(dt)[1:2]))
  }

  # Check if required columns exist
  required_cols <- c(x_var, y_var, color_var)
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Ensure the required variables are numeric
  dt[[x_var]] <- as.numeric(dt[[x_var]])
  dt[[y_var]] <- as.numeric(dt[[y_var]])

  # Convert color_var to factor
  dt[[color_var]] <- as.factor(dt[[color_var]])

  # Calculate mean values for each combination of x_var and color_var
  dt_summary <- dt[, list(mean_y = mean(get(y_var))), by = c(x_var, color_var)]

  p <- ggplot(dt_summary, aes_string(x = x_var, y = "mean_y", color = color_var, group = color_var)) +
    geom_line() +
    theme_minimal() +
    labs(title = title, subtitle = subtitle, x = x_var, y = y_var) +
    scale_color_discrete(name = color_var)

  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }

  if (!is.null(line_types)) {
    p <- p + scale_linetype_manual(values = line_types)
  }

  if (!is.null(line_widths)) {
    p <- p + scale_size_manual(values = line_widths)
  }

  return(p)
}
