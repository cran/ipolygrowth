#' @title Growth parameter estimates for single sample
#'
#' @description
#' This function estimates growth parameters for a single (biological) sample. Technical replicates (multiple time series) are allowed.
#'
#' @param data Input data frame containing the time and dependent variable (y) from a single biological sample.
#' Data needs to be in long format (i.e. one row per time point).
#' Remove rows with missing in the dependent variable (y).
#' @param time.name Name of the time variable. The variable needs to be numeric.
#' @param y.name Name of the dependent variable (y). The variable needs to be numeric.
#' @param epsilon Convergence threshold for max y time calculation.
#' The input represents the fraction of the range of the observed dependent variable (y).
#' The input needs to be between 0 and 1, and a small value is recommended. Default is 0.2%.
#'
#' @details
#' This function calculates growth curve parameters for a single sample.
#' A 4th degree polynomial is fit to the input data using ordinary least squares estimation.
#' Peak growth time is identified using the second derivative of the estimated polynomial function.
#' Peak growth rate is calculated using the first derivative at peak growth time.
#' Doubling time at peak growth is calculated using the equation: ln(2)/peak growth rate .
#' Lag time is determined using linear interpolation of the peak growth rate to identify the start of the exponential growth phase.
#' Max y time is identified by convergence of the dependent variable where the growth curve reaches an asymptote, with convergence threshold epsilon.
#' Max y is the value of the fitted polynomial function at max y time.
#'
#' @return A list that contains a table of estimates, the polynomial model, a table of beta coefficients, and a table of fitted values.
#' Growth parameters include peak growth rate, peak growth time, doubling time (at the peak growth), lag time, max y, and max y time.
#'
#' @importFrom rlang .data
#'
#' @examples
#' library(dplyr)
#' data <- growthrates::bactgrowth
#' df.singlesample <- data %>% dplyr::filter(strain == "D", conc == 0)
#' out.singlesample <- ipg_singlesample(data = df.singlesample, time.name = "time", y.name = "value")
#' @export
#'
#'
ipg_singlesample <- function(data, time.name, y.name, epsilon = 0.2/100) {
  stopifnot("Y is not found in the data" = !is.null(data[[y.name]]))
  stopifnot("Y must be a numeric variable" = is.numeric(data[[y.name]]))
  stopifnot("Time is not found in the data" = !is.null(data[[time.name]]))
  stopifnot("Time must be a numeric variable" = is.numeric(data[[time.name]]))
  tryCatch(stopifnot(epsilon>0, epsilon<1),
           error = function(e) {
             e$message <- paste("epsilon is outside the recommended range (0-1)")
             print(e$class)
             stop(e)
             })

  time <- fit <- NULL                                                                               # to avoid notes from CRAN check

      gcdata <- data %>%
        dplyr::mutate(time = .data[[time.name]],
                      time_sq = time^2,
                      time_cb = time^3,
                      time_qt = time^4) %>%
        dplyr::arrange(time)

      ## polynomial model ##----
      curve <- stats::lm(stats::as.formula(paste(y.name, "~ time + time_sq + time_cb + time_qt")),
                         data = gcdata)                                                             # 4th order polynomial regression (assumes independent obs)

      gcdata2 <- gcdata %>%
        dplyr::select(time, tidyselect::any_of(y.name)) %>%
        stats::na.omit() %>%                                                                        # to avoid mismatch in length if any missing in y
        dplyr::mutate(fit = stats::predict(curve)) %>%
        dplyr::select(time, fit) %>%
        dplyr::distinct() %>%                                                                       # create a data frame for predicted log OD and time
        dplyr::mutate(diffOD = c(0, diff(fit, lag = 1)))                                            # difference in fitted log OD

      tb.beta <- data.frame(B0 = NA_real_, B1 = NA_real_, B2 = NA_real_, B3 = NA_real_, B4 = NA_real_)
      tb.beta$B0 <- A0 <- as.numeric(stats::coef(curve)[1])                                         # extract betas
      tb.beta$B1 <- A1 <- as.numeric(stats::coef(curve)[2])
      tb.beta$B2 <- A2 <- as.numeric(stats::coef(curve)[3])
      tb.beta$B3 <- A3 <- as.numeric(stats::coef(curve)[4])
      tb.beta$B4 <- A4 <- as.numeric(stats::coef(curve)[5])

      tb.result <- data.frame(
        peak.growth.rate = NA_real_, peak.growth.time = NA_real_, doubling.time = NA_real_,
        start.peak.growth = NA_real_, fit.max.y = NA_real_, fit.max.y.time = NA_real_)


      ## calculate peak growth rate ##----
        # first derivative of 4th degree polynomial is A1 + A2 * (2 * x) + A3 * (3 * x^2) + A4 * (4 * x^3)
        # second derivative of 4th degree polynomial is A2 * 2 + A3 * (3 * (2 * x)) + A4 * (4 * (3 * x^2))
      x.roots <- polyroot(c(A2 * 2, A3 * (3 * (2)), A4 * (4 * (3))))                                # find roots of second derivative; peak growth times
      x <- c(x.roots,0)                                                                             # add 0 to calculate the slope at 0

      pgr <- Re(eval(A1 + A2 * (2 * x) + A3 * (3 * x^2) + A4 * (4 * x^3)))                          # calculate growth rate (slope) at peak growth times (real number only)
      max <- max(gcdata$time, na.rm = TRUE)
      pgr.index <- which(pgr == max(pgr[pgr>0 & Re(x) >= 0 & Re(x) <= max]))[1]                     # choose growth time with maximum growth rate and peak growth time between the start and end point

      tb.result$peak.growth.time <- pgt <- Re(x[pgr.index])                                         # peak growth time
      tb.result$peak.growth.rate <- pgr[pgr.index]                                                  # peak growth rate at peak growth time
      tb.result$doubling.time <-
        Re(log(2)/eval(A1 + A2*(2*x[pgr.index]) + A3*(3*x[pgr.index]^2) + A4*(4*x[pgr.index]^3)))   # doubling time at peak growth = ln(2)/growth rate


      ## calculate lag time (start of peak growth) ##----
        #  the time when the line of intercept from the polynomial model intersects with the linear line of peak growth at peak growth time
      if (tb.result$peak.growth.time == 0) {
        tb.result$start.peak.growth <- 0
      } else {
        y0 <- gcdata2$fit[1]                                                                        # intercept of polynomial curve
        y1 <- eval(expression(A0 + A1 * pgt + A2 * pgt^2 + A3 * pgt^3 + A4 * pgt^4))                # calculate predicted y at peak growth time
        tb.result$start.peak.growth <- solve(pgr[pgr.index], pgr[pgr.index]*pgt-(y1-y0))            # solve for x0 when (y1-y0)/(pgt-x0) = pgr
      }


      ## calculate max y ##----
      r <- gcdata %>% dplyr::select(tidyselect::any_of(y.name)) %>% range()                         # get range of y
      threshold <- (r[2] - r[1]) * epsilon                                                          # calculate cutoff based on range

      if (any((gcdata2$diffOD < threshold) & (gcdata2$diffOD > 0), na.rm=TRUE)) {                   # look for difference in fitted y b/w 0 and threshold

        if (any(gcdata2$time[(gcdata2$diffOD < threshold) & (gcdata2$diffOD > 0)]
                > tb.result$peak.growth.time, na.rm=TRUE)) {                                        # look for time of that larger than peak growth time

          temp.time <- gcdata2$time[(gcdata2$diffOD < threshold) & (gcdata2$diffOD > 0)]            # find all times when fitted y is b/w 0 and threshold
          temp <- temp.time[temp.time > tb.result$peak.growth.time][1]                              # select the first one that is later than the peak growth time
          tb.result$fit.max.y <- gcdata2$fit[gcdata2$time == temp]                                  # return the corresponding fitted OD value
        } else {
          tb.result$fit.max.y <- max(gcdata2$fit)
        }
      } else {
        tb.result$fit.max.y <- max(gcdata2$fit)
      }

      tb.result$fit.max.y.time <- gcdata2$time[gcdata2$fit == tb.result$fit.max.y][1]               # record time of max y for all scenario

      if (tb.result$fit.max.y.time == max(gcdata2$time)){
        message(paste0("max y time is equal to the largest value of \"", time.name, "\""))
      }


      ## clean output ##----
      colnames(tb.result) <- c("peak growth rate", "peak growth time", "doubling time", "lag time", "max y", "max y time")

      output <- list(estimates = tb.result,
                     model = curve,
                     betas = tb.beta,
                     fitted = gcdata2 %>% dplyr::select(time, fit))

      return(output)
}
