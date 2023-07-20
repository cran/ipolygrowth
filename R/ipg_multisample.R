#' @title Growth parameter estimates for multiple samples
#'
#' @description
#' This function estimates growth parameters for multiple (biological) samples. Technical replicates (multiple time series) are allowed.
#'
#' @param data Input data frame containing the time and dependent variable (y) of multiple biological samples.
#' Data needs to be in long format (i.e. one row per time point per sample).
#' Remove rows with missing in the dependent variable (y).
#' @param id Unique identifier to indicate each biological sample.
#' @param time.name Name of the time variable. The variable needs to be numeric.
#' @param y.name Name of the dependent variable (y). The variable needs to be numeric.
#' @param epsilon Convergence threshold for max y time calculation.
#' The input represents the fraction of the range of the observed dependent variable (y).
#' It needs to be between 0 and 1, and a small value is recommended.
#' The input can be either a single value or a vector of different values (in the same length and order as id) for multiple samples.
#' Default is 0.2% for all samples.
#'
#' @details
#' The function uses the same approach to estimate growth parameters as in `ipg_singlesample()`.
#'
#' @return A list that contains a table of estimates, polynomial models, a table of beta coefficients, and a table of fitted values, all by sample ID.
#' Growth parameters include peak growth rate, peak growth time, doubling time (at the peak growth), lag time, max y, and max y time.
#'
#' @importFrom rlang .data
#'
#' @examples
#' library(dplyr)
#' data <- growthrates::bactgrowth
#' data <- data %>% mutate(id = paste(strain, conc, sep = "-"))
#' out.multisample <- ipg_multisample(data, "id", "time", "value", 0.2/100)
#' @export
#'
#'
ipg_multisample <- function(data, id, time.name, y.name, epsilon = 0.2/100) {
  n <- . <- NULL
  X <- 0


  ls.df <- split(data, data[, id], sep = "^")                                                       # split datasets; name is separated by ^

  df.factor <- data.frame(order = names(ls.df)) %>%
    tidyr::separate_wider_delim("order", delim = "^", names = id)                                   # unique combination of stratified id by row;

  if (length(epsilon == 1)) {
    epsilon <- rep(epsilon, length(ls.df))                                                          # assign epsilon for all samples
  } else {
    epsilon <- epsilon
  }

  ls.model <- vector(mode = 'list', length = nrow(df.factor))                                       # create lists to store model and fitted values
  names(ls.model) <- do.call(paste, df.factor[names(df.factor)])                                    # match component name with stratified id

  output.estimates <- data.frame(
    peak.growth.rate = NA_real_, peak.growth.time = NA_real_,doubling.time = NA_real_,
    start.peak.growth = NA_real_, fit.max.y = NA_real_, fit.max.y.time = NA_real_)                  # create a dataset to store param estimates together

  output.betas <- data.frame(
    B0 = NA_real_, B1 = NA_real_, B2 = NA_real_, B3 = NA_real_, B4 = NA_real_)                      # create a dataset to store beta coefficient

  output.fitted <- data.frame(time = NULL, fit = NULL)                                              # create a dataset to store fitted values

  duptimes <- data[, c(id, time.name)] %>%                                                          # There will be issue when missing in y
    dplyr::distinct() %>%
    dplyr::group_by(dplyr::across(tidyselect::all_of(id))) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last")
    duptimes <- duptimes$n                                                                          # get duplicated times for each data subset
    df.factor2 <- df.factor %>% dplyr::slice(rep(1:nrow(.), duptimes))                              # duplicate id for fitted value

  for (i in 1:length(ls.df)) {
    tryCatch(
      ipg_singlesample(ls.df[[i]], time.name, y.name, epsilon[i]),
      message = function(cnd) {
        X <<- X+1
      },
      finally = {
        temp <- suppressMessages(ipg_singlesample(ls.df[[i]], time.name, y.name, epsilon[i]))  # requires a vector of threshold
        ls.model[[i]] <- temp$model
        output.estimates[i, ] <- temp$estimates
        output.betas[i, ] <- temp$betas
        output.fitted <- rbind(output.fitted, temp$fitted)
      }
    )
  }

  if (X>0) {
    message(paste0("max y time is equal to the largest value of \"", time.name, "\" in ", X, " samples."))
  }


  ## clean output ##----
  colnames(output.estimates) <- c("peak growth rate", "peak growth time", "doubling time", "lag time", "max y", "max y time")
  output.estimates <- cbind(df.factor, output.estimates)
  output.betas <- cbind(df.factor, output.betas)
  output.fitted <- cbind(df.factor2, output.fitted)

  ls.output <- list(
    model = ls.model,
    estimates = output.estimates,
    betas = output.betas,
    fitted = output.fitted
  )
  return(ls.output)
}
