#' @import ggplot2
#' @importFrom stats time as.formula setNames
#' @importFrom ggplot2 aes facet_wrap
utils::globalVariables(c("value", "variable"))

#' Plot Malaria Model Results
#'
#' This function plots the results of a malaria model simulation.
#'
#' @param model_results A data frame containing the model results.
#' @param res_time The time at which resistance was introduced.
#' @param ton The time at which treatment was turned on.
#' @param toff The time at which treatment was turned off.
#' @param output_vars A character vector specifying the output variables to plot. If NULL, all output variables will be plotted.
#' @param colors A named character vector specifying the colors for each variable. If NULL, default colors will be used.
#' @param line_types A named character vector specifying the line types for each variable. If NULL, default line types will be used.
#' @param line_widths A named numeric vector specifying the line widths for each variable. If NULL, default line widths will be used.
#' @param title The title of the plot.
#' @param subtitle The subtitle of the plot. If NULL, a default subtitle will be generated.
#' @param show_lines A character vector specifying which vertical lines to show. Can include "res_time", "ton", and "toff".
#' @param vline_colors A named character vector specifying the colors for each vertical line. Default colors are provided.
#' @param vline_types A named character vector specifying the line types for each vertical line. Default line types are provided.
#' @param vline_widths A named numeric vector specifying the line widths for each vertical line. Default widths are provided.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \donttest{
#' library(AMRSpreadModel)
#'
#' params <- list(
#'   S0 = 0.9, Ds0 = 0.05, As0 = 0.02, Ts0 = 0.03, DR0 = 0, AR0 = 0,
#'   TR0 = 0, Sv0 = 0.9, Ev_s0 = 0.08, Iv_s0 = 0.02, Ev_r0 = 0, Iv_r0 = 0,
#'   m = 10, a = 0.3, b = 0.5876259, phi = 0.7, fT = 0.1, rD = 0.2,
#'   rA = 0.01, rTs = 0.2, rTR_true = 0.01, e = 0.132, mu = 0.132,
#'   n = 10, c_A = 0.05, c_D = 0.06, c_T = 0.02, ton = 8000,
#'   toff = 99999999, res_time = 3000, init_res = 0.2
#' )
#'
#' model <- malaria_model(params = params)
#' if (!is.null(model)) {
#'   results <- model$run(0:500)
#'   plot_model(results, res_time = params$res_time, ton = params$ton, toff = params$toff)
#' } else {
#'   print("Failed to create model. Check parameters and try again.")
#' }
#' }
plot_model <- function(model_results, res_time, ton, toff = NULL, output_vars = NULL,
                       colors = NULL, line_types = NULL, line_widths = NULL, title = "Malaria Model Results",
                       subtitle = NULL, show_lines = c("res_time", "ton", "toff"),
                       vline_colors = NULL,
                       vline_types = NULL,
                       vline_widths = NULL) {

  col_names <- colnames(model_results)

  if (is.null(output_vars)) {
    output_vars <- col_names[col_names != "t"]
  } else {
    if (!all(output_vars %in% col_names)) {
      stop("Invalid output variables. Please make sure the variables are valid column names in the model results.")
    }
  }

  results_df <- data.frame(time = model_results[, "t"], model_results[, output_vars, drop = FALSE])
  results_long <- tidyr::pivot_longer(results_df, cols = -time, names_to = "variable", values_to = "value")

  # Create a data frame for vertical lines
  vertical_lines <- data.frame(
    time = c(res_time, ton, toff),
    variable = c("res_time", "ton", "toff"),
    value = max(results_long$value, na.rm = TRUE)
  )

  # Filter out NA times
  vertical_lines <- vertical_lines[!is.na(vertical_lines$time), ]

  p <- ggplot() +
    geom_line(data = results_long, aes(x = time, y = value, color = variable, linetype = variable, size = variable), na.rm = TRUE) +
    geom_segment(data = vertical_lines, aes(x = time, xend = time, y = 0, yend = value, color = variable, linetype = variable, size = variable), show.legend = TRUE)

  if (is.null(subtitle)) {
    subtitle_parts <- c()
    if ("res_time" %in% show_lines && !is.null(res_time)) {
      subtitle_parts <- c(subtitle_parts, paste("Resistance introduced at time", res_time))
    }
    if ("ton" %in% show_lines && !is.null(ton)) {
      subtitle_parts <- c(subtitle_parts, paste("Treatment turned on at time", ton))
    }
    if ("toff" %in% show_lines && !is.null(toff)) {
      subtitle_parts <- c(subtitle_parts, paste("Treatment turned off at time", toff))
    }
    subtitle <- paste(subtitle_parts, collapse = ", ")
  }

  p <- p + labs(x = "Time", y = "Value", color = "Variable", linetype = "Variable", size = "Variable",
                title = title, subtitle = subtitle)

  all_variables <- unique(c(output_vars, vertical_lines$variable))

  # Default colors
  default_colors <- scale_color_hue()$palette(length(all_variables))
  names(default_colors) <- all_variables

  # Default line types
  default_line_types <- setNames(rep("solid", length(output_vars)), output_vars)
  default_vline_types <- setNames(rep("dotted", length(vertical_lines$variable)), vertical_lines$variable)

  # Default line widths
  default_line_widths <- setNames(rep(0.5, length(output_vars)), output_vars)
  default_vline_widths <- setNames(rep(1, length(vertical_lines$variable)), vertical_lines$variable)

  if (!is.null(colors)) {
    color_values <- c(setNames(colors, output_vars), setNames(vline_colors, vertical_lines$variable))
    color_values <- color_values[all_variables]
  } else {
    color_values <- default_colors
  }

  if (!is.null(line_types)) {
    linetype_values <- c(setNames(line_types, output_vars), setNames(vline_types, vertical_lines$variable))
    linetype_values <- linetype_values[all_variables]
  } else {
    linetype_values <- c(default_line_types, default_vline_types)
  }

  if (!is.null(line_widths)) {
    linewidth_values <- c(setNames(line_widths, output_vars), setNames(vline_widths, vertical_lines$variable))
    linewidth_values <- linewidth_values[all_variables]
  } else {
    linewidth_values <- c(default_line_widths, default_vline_widths)
  }

  p <- p + scale_color_manual(values = color_values)
  p <- p + scale_linetype_manual(values = linetype_values)
  p <- p + scale_size_manual(values = linewidth_values)
  p <- p + theme_minimal()

  return(p)
}
