#' @import parallel
#' @importFrom data.table rbindlist
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export range_model
NULL

#' Summarize Model Results
#'
#' @param results A data frame containing the model results.
#' @param time_points A vector of specific time points to summarize results.
#' @return A summarized data frame.
summarize_model_results <- function(results, time_points) {
  summarized_results <- results[results$t %in% time_points, ]
  return(summarized_results)
}

#' Run Malaria Model for Ranges of Parameters
#'
#' @param param_ranges A named list of parameter ranges. Each element should be a vector of values to run the model over.
#' @param params A list of other parameters to pass to the `malaria_model` function
#' @param output_vars A character vector specifying the output variables to return. If NULL, all output variables will be returned.
#' @param use_eir_ft Logical. If TRUE, uses EIR and ft to generate parameters. If FALSE, uses provided params.
#' @param run_time Numeric. The total time to run the model. Default is 5000.
#' @param time_points Numeric vector. Specific time points to output results. If NULL, all time points will be returned.
#' @param time_scale Character. Either "day" or "year". Default is "day". This affects both input and output time units.
#' @param ton Time at which treatment is turned on
#' @param toff Time at which treatment is turned off
#' @param init_res Initial resistance level
#' @param res_time Time at which resistance is introduced
#' @param rTR_true True treatment rate for resistant parasites
#' @param verbose Logical. If TRUE, prints detailed logs. Default is FALSE.
#' @return A data frame of results for each combination of parameter values
#' @export
range_model <- function(param_ranges, params = NULL, output_vars = NULL, use_eir_ft = FALSE,
                        run_time = 5000, time_points = NULL, time_scale = "day",
                        ton = 5000, toff = 50000, init_res = 0.01, res_time = 3000, rTR_true = 0.1,
                        verbose = FALSE) {
  if (verbose) {
    cat("Generating parameter grid...\n")
  }
  param_grid <- tryCatch({
    expand.grid(param_ranges)
  }, error = function(e) {
    stop("Error in expand.grid with param_ranges: ", e$message)
  })

  n_combinations <- nrow(param_grid)
  if (verbose) {
    cat(paste("Number of parameter combinations:", n_combinations), "\n")
    cat("Parameter grid:\n")
    print(param_grid)
  }

  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = n_combinations, style = 3)

  # Convert run_time and time_points to days if time_scale is "year"
  if (time_scale == "year") {
    run_time <- run_time * 365
    if (!is.null(time_points)) {
      time_points <- time_points * 365
    }
  }

  if (is.null(time_points)) {
    time_points <- 0:run_time
  }

  if (verbose) {
    cat(paste("Run time:", run_time), "\n")
    cat("Time points:\n")
    print(time_points)
  }

  results_list <- list()

  for (i in 1:n_combinations) {
    row <- param_grid[i, ]
    tryCatch({
      if (use_eir_ft) {
        if (!all(c("EIR", "ft") %in% names(row))) {
          stop("When use_eir_ft is TRUE, param_ranges must include 'EIR' and 'ft'")
        }
        if (verbose) {
          cat(paste("Generating parameters for EIR =", row[["EIR"]], "and ft =", row[["ft"]]), "\n")
        }
        current_params <- phi_eir_rel(row[["EIR"]], row[["ft"]], ton, toff, init_res, res_time, row[["rTR_true"]])
        if (is.null(current_params)) {
          stop("phi_eir_rel returned NULL")
        }
        current_params$rTR_true <- row[["rTR_true"]]
      } else {
        if (is.null(params)) {
          stop("When use_eir_ft is FALSE, params must be provided")
        }
        current_params <- modifyList(params, as.list(row))
      }

      if (verbose) {
        cat("Creating malaria model...\n")
      }
      model <- malaria_model(params = current_params, ton = ton, toff = toff, init_res = init_res, res_time = res_time, rTR_true = current_params$rTR_true, verbose = verbose)
      if (is.null(model)) {
        stop("malaria_model returned NULL")
      }
      if (verbose) {
        cat("Running model...\n")
      }
      results <- as.data.frame(model$run(0:run_time))
      if (verbose) {
        cat("Summarizing results...\n")
      }
      summarized_results <- summarize_model_results(results, time_points)

      # Add all varied parameters to the results
      for (param_name in names(row)) {
        summarized_results[[param_name]] <- row[[param_name]]
      }

      summarized_results$prevalence_res <- (summarized_results$AR + summarized_results$DR + summarized_results$TR) /
        (summarized_results$As + summarized_results$Ds + summarized_results$Ts + summarized_results$AR + summarized_results$DR + summarized_results$TR)

      if (verbose) {
        cat("Results generated successfully\n")
      }
      results_list[[i]] <- summarized_results
    }, error = function(e) {
      message("Error in processing row: ", paste(row, collapse = ", "))
      message("Error message: ", e$message)
      return(NULL)
    })

    # Update progress bar
    setTxtProgressBar(pb, i)
  }

  # Close progress bar
  close(pb)

  # Remove NULL results (if any)
  results_list <- results_list[!sapply(results_list, is.null)]
  if (verbose) {
    cat(paste("Number of successful runs:", length(results_list)), "\n")
  }

  if (length(results_list) == 0) {
    warning("No successful model runs. Check your parameters and model.")
    return(data.table::data.table())
  }

  results_df <- data.table::rbindlist(results_list, fill = TRUE)

  # Convert time scale to years if necessary
  if (time_scale == "year") {
    results_df$t <- results_df$t / 365
  }

  if (verbose) {
    cat(paste("Number of rows in final results:", nrow(results_df)), "\n")
    cat(paste("Columns in final results:", paste(names(results_df), collapse = ", ")), "\n")
  }

  return(results_df)
}
