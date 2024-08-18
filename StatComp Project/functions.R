# Chenyao Yu (s2156882)

# Place your function definitions that may be needed in the report.Rmd, including function documentation.
# You can also include any needed library() calls here

#' Construct Bayesian CI for Monte Carlo p-values
#'
#' @param x How many times we observe a randomised test statistic as extreme as or more extreme than the observed test statistic
#' @param N Total number of replications
#' @param level The confidence level required
#'
#' @return A list contains upper bar and lower bar of the CI

p_value_CI <- function(x, N, level = 0.95) {
  return(list(lower = qbeta((1 - level) / 2, 
                            shape1 = 1 + x, 
                            shape2 = 1 + N - x),
              upper = qbeta(1 - (1 - level) / 2, 
                            shape1 = 1 + x, 
                            shape2 = 1 + N - x)))
}

#' Estimate model for the square root of monthly average precipitation
#'
#' @param k Frequency for adding covariates to capture seasonal variability
#' @param subset Specification of the rows to be used: defaults to all rows
#'
#' @return An object of class "lm"

precipitation_estimate <- function(k = 0, subset = NULL) {
  formula0 <- "Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear"
  
  m0 <- lm(formula0,
           data = ghcnd_month, subset = subset)
  
  formula_k <- paste(
    "+ I(cos(2 * pi *",
    1:k,
    "* DecYear)) + I(sin(2 * pi *",
    1:k,
    "* DecYear))", collapse = " ")
  
  mk <- lm(paste(formula0,
                 formula_k),
           data = ghcnd_month, subset = subset)
  
  if (k == 0) {
    return(m0)
  } else {
    return(mk)
  }
}


