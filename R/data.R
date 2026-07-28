#' Default value labels that are generally considered as missing value.
#' 
#' @source Hand picked labels that appear in the SÜP, HAP and POP study data.
#' @format ## `Default missing  valuelabels`
#' Character vector of value labels.
"missing_labels_default"

#' Default numeric values generally considered as missing value.
#'
#' @source numeric value(s) that are standardly used for missing values.
#' @format ## `Default numercial missing values`
"missing_values_default"

#' Variables feeding the Multidimensional Prognostic Index (MPI).
#'
#' @source Domain definition of the NAPKON SUEP MPI, see `data-raw/mpi_config.R`.
#' @format ## `mpi_domains`
#' Data frame with one row per domain variable and the columns:
#' \describe{
#'   \item{domain}{Domain the variable contributes to.}
#'   \item{table}{Table holding the variable, file name without ".csv".}
#'   \item{variable}{Variable name.}
#'   \item{start}{Variable holding the start of the validity interval, or NA to
#'     fall back to "<variable>_start".}
#'   \item{end}{Variable holding the end of the validity interval, or NA to fall
#'     back to "<variable>_end".}
#'   \item{companions}{Comma separated further variables that share the interval
#'     of `variable` and are needed for scoring, or NA.}
#'   \item{aggregate}{How to combine records valid at the same visit: "max" for
#'     the highest code, "count" for the number of distinct entries.}
#' }
"mpi_domains"

#' Scoring rules of the Multidimensional Prognostic Index (MPI).
#'
#' @source Scoring definition of the NAPKON SUEP MPI, see
#'   `data-raw/mpi_config.R`.
#' @format ## `mpi_scoring`
#' Data frame with one row per domain and the columns:
#' \describe{
#'   \item{domain}{Domain name. Also fixes the order of the domain score
#'     columns.}
#'   \item{type}{Scoring rule: "sum", "bmi", "living_situation" or "none".}
#'   \item{cut_0}{Cut-off for domain score 0, no problem.}
#'   \item{cut_0.5}{Cut-off for domain score 0.5, minor problem.}
#'   \item{cut_1}{Cut-off for domain score 1, major problem.}
#' }
#' For "sum" domains a cut-off is the inclusive lower bound of its score. For
#' "bmi" domains the cut-offs are read literally as "lo-hi" and "<x,>y" bands.
"mpi_scoring"