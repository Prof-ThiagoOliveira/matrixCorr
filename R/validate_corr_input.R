#' Validation and normalization for correlation
#'
#' @description
#' Validates and normalizes input for correlation computations. Accepts either a
#' numeric matrix or a data frame, filters numeric columns, checks dimensions and
#' missing values, and returns a numeric (double) matrix with preserved
#' column names.
#'
#' @details
#' Rules enforced:
#' \itemize{
#'   \item Input must be a matrix or data.frame.
#'   \item Only numeric (integer or double) columns are retained (data.frame path).
#'   \item At least two numeric columns are required.
#'   \item All columns must have the same length and \eqn{\ge} 2 observations.
#'   \item Missing values are not allowed.
#'   \item Returns a \code{double} matrix; integer input is converted once.
#' }
#'
#' @param data A matrix or data frame. Non-numeric columns are dropped
#'   (data.frame path). For matrix input, storage mode must be integer or double.
#'
#' @return A numeric matrix (type \code{double}) with column names preserved.
#'
#' @seealso [pearson_corr()], [spearman_rho()], [kendall_tau()]
#'
#' @author Thiago de Paula Oliveira
#' @useDynLib matrixCorr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
validate_corr_input <- function(data) {
  validate_corr_input_cpp(data)
}
