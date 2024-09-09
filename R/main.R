#---------------------------------------------------------------
# PSPA: post-prediction adaptive inference
# Jiacheng Miao, Xinran Miao, Yixuan Wu, Jiwei Zhao, and Qiongshi Lu
# Available from https://arxiv.org/abs/2311.14220
#---------------------------------------------------------------

#' PSPA M-Estimation for ML-predicted labels
#'
#' \code{pspa_y} function conducts post-prediction M-Estimation.
#' @param X_lab Array or DataFrame containing observed covariates in labeled data.
#' @param X_unlab Array or DataFrame containing observed or predicted covariates in unlabeled data.
#' @param Y_lab Array or DataFrame of observed outcomes in labeled data.
#' @param Yhat_lab Array or DataFrame of predicted outcomes in labeled data.
#' @param Yhat_unlab Array or DataFrame of predicted outcomes in unlabeled data.
#' @param alpha Specifies the confidence level as 1 - alpha for confidence intervals.
#' @param weights weights vector PSPA linear regression (d-dimensional, where d equals the number of covariates).
#' @param quant quantile for quantile estimation
#' @param intercept Boolean indicating if the input covariates' data contains the intercept (TRUE if the input data contains)
#' @param method indicates the method to be used for M-estimation. Options include "mean", "quantile", "ols", "logistic", and "poisson".
#' @return  A summary table presenting point estimates, standard error, confidence intervals (1 - alpha), P-values, and weights.
#' @examples
#' data <- sim_data_y()
#' X_lab <- data$X_lab
#' X_unlab <- data$X_unlab
#' Y_lab <- data$Y_lab
#' Yhat_lab <- data$Yhat_lab
#' Yhat_unlab <- data$Yhat_unlab
#' pspa_y(Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
#'       alpha = 0.05, method = "mean")
#' pspa_y(Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
#'       alpha = 0.05, quant = 0.75, method = "quantile")
#' pspa_y(X_lab = X_lab, X_unlab = X_unlab,
#'       Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
#'       alpha = 0.05, method = "ols")
#' @export

# A general function for M-estimation with ML-predicted labels
pspa_y <- function(X_lab = NA, X_unlab = NA, Y_lab, Yhat_lab, Yhat_unlab,
                  alpha = 0.05, weights = NA, quant = NA, intercept = FALSE,
                  method) {
  # Common values
  if (method %in% c("ols", "logistic", "poisson")) {
    if (intercept) {
      X_lab <- as.matrix(X_lab)
      X_unlab <- as.matrix(X_unlab)
    } else {
      X_lab <- cbind(1, as.matrix(X_lab))
      X_unlab <- cbind(1, as.matrix(X_unlab))
    }
  }

  Y_lab <- as.matrix(as.numeric(unlist(Y_lab)))
  Yhat_lab <- as.matrix(as.numeric(unlist(Yhat_lab)))
  Yhat_unlab <- as.matrix(as.numeric(unlist(Yhat_unlab)))
  n <- nrow(Y_lab)
  N <- nrow(Yhat_unlab)
  if (method %in% c("mean", "quantile")) {
    q <- 1
  } else {
    q <- ncol(X_lab)
  }

  # Initial values
  est <- est_ini(X_lab, Y_lab, quant, method)
  if (is.na(sum(weights))) {
    current_w <- rep(0, q)
    optimized_weight <- optim_weights(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab, current_w, est, quant, method)
    current_w <- optimized_weight
  } else{
    current_w <- weights

  }

  # Calculate the final coefficients and standard errors
  est <- optim_est(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab, current_w, est, quant, method)
  final_sigma <- Sigma_cal(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab, current_w, est, quant, method)
  standard_errors <- sqrt(diag(final_sigma) / n)
  p_values <- 2 * pnorm(abs(est / standard_errors), lower.tail = FALSE)

  # Create output table
  lower_ci <- est - qnorm(1 - alpha / 2) * standard_errors
  upper_ci <- est + qnorm(1 - alpha / 2) * standard_errors
  output_table <- data.frame(
    Estimate = est, Std.Error = standard_errors,
    Lower.CI = lower_ci, Upper.CI = upper_ci, P.value = p_values, Weight = current_w
  )

  colnames(output_table) <- c("Estimate", "Std.Error", "Lower.CI", "Upper.CI", "P.value", "Weight")

  return(output_table)
}

#' PSPA M-Estimation for ML-predicted covariates
#'
#' \code{pspa_x} function conducts post-prediction M-Estimation.
#' @param X_lab Array or DataFrame containing observed covariates in labeled data.
#' @param Xhat_lab Array or DataFrame containing predicted covariates in labeled data.
#' @param Xhat_unlab Array or DataFrame containing predicted covariates in unlabeled data.
#' @param Y_lab Array or DataFrame of observed outcomes in labeled data.
#' @param Y_unlab Array or DataFrame of observed outcomes in unlabeled data.
#' @param alpha Specifies the confidence level as 1 - alpha for confidence intervals.
#' @param weights weights vector PSPA linear regression (d-dimensional, where d equals the number of covariates).
#' @param quant quantile for quantile estimation
#' @param intercept Boolean indicating if the input covariates' data contains the intercept (TRUE if the input data contains)
#' @param method indicates the method to be used for M-estimation. Options include "mean", "quantile", "ols", "logistic", and "poisson".
#' @return  A summary table presenting point estimates, standard error, confidence intervals (1 - alpha), P-values, and weights.
#' @examples
#' data <- sim_data_x()
#' X_lab <- data$X_lab
#' Xhat_lab <- data$Xhat_lab
#' Xhat_unlab <- data$Xhat_unlab
#' Y_lab <- data$Y_lab
#' Y_unlab <- data$Y_unlab
#' pspa_x(X_lab = X_lab, Xhat_lab = Xhat_lab,
#'       Xhat_unlab = Xhat_unlab, Y_lab = Y_lab, Y_unlab = Y_unlab,
#'       alpha = 0.05, method = "ols")
#' @export
           
pspa_x <- function(X_lab, Xhat_lab, Xhat_unlab, Y_lab, Y_unlab,
                  alpha = 0.05, weights = NA, quant = NA, intercept = FALSE,
                  method) {
  # Common values
  if (method %in% c("ols", "logistic", "poisson")) {
    if (intercept) {
      X_lab <- as.matrix(X_lab)
      Xhat_lab <- as.matrix(Xhat_lab)
      Xhat_unlab <- as.matrix(Xhat_unlab)
    } else {
      X_lab <- cbind(1, as.matrix(X_lab))
      Xhat_lab <- cbind(1, as.matrix(Xhat_lab))
      Xhat_unlab <- cbind(1, as.matrix(Xhat_unlab))
    }
  }

  Y_lab <- as.matrix(as.numeric(unlist(Y_lab)))
  Y_unlab <- as.matrix(as.numeric(unlist(Y_unlab)))
  n <- nrow(Y_lab)
  N <- nrow(Y_unlab)
  if (method %in% c("mean", "quantile")) {
    q <- 1
  } else {
    q <- ncol(X_lab)
  }

  # Initial values
  est <- est_ini(X_lab, Y_lab, quant, method)
  current_w <- if (is.na(sum(weights))) rep(0, q) else weights

  current_w <- optim_weights_x(X_lab, Xhat_lab, Xhat_unlab, Y_lab, Y_unlab, current_w, est, quant, method)

  # Calculate the final coefficients and standard errors
  est <- optim_est_x(X_lab, Xhat_lab, Xhat_unlab, Y_lab, Y_unlab, current_w, est, quant, method)
  A_lab <- A(X_lab, Y_lab, quant, est, method)
  Ahat_lab <- A(Xhat_lab, Y_lab, quant, est, method)
  Ahat_unlab <- A(Xhat_unlab, Yhat_unlab, quant, est, method)

  final_sigma <- Sigma_cal_x(X_lab, Xhat_lab, Xhat_unlab, Y_lab, Y_unlab, current_w, est, quant, A_lab, Ahat_lab, Ahat_unlab, method)
  standard_errors <- sqrt(diag(final_sigma) / n)
  p_values <- 2 * pnorm(abs(est / standard_errors), lower.tail = FALSE)

  # Create output table
  lower_ci <- est - qnorm(1 - alpha / 2) * standard_errors
  upper_ci <- est + qnorm(1 - alpha / 2) * standard_errors
  output_table <- data.frame(
    Estimate = est, Std.Error = standard_errors,
    Lower.CI = lower_ci, Upper.CI = upper_ci, P.value = p_values, Weight = current_w
  )

  colnames(output_table) <- c("Estimate", "Std.Error", "Lower.CI", "Upper.CI", "P.value", "Weight")

  return(output_table)
}
