# ------------------------------------------------------------------------------
# MULTIDIMENSIONAL PROGNOSTIC INDEX (MPI) CONFIGURATION
# ------------------------------------------------------------------------------
#
# Two objects define the MPI:
#
#   mpi_domains  which variables feed each domain, and where they live
#   mpi_scoring  how each domain is turned into a 0 / 0.5 / 1 domain score
#
# They are defined here as source rather than shipped as data/*.rda, so that a
# change to a threshold reviews as a readable diff and needs no regeneration
# step, and so that the numbers sit next to the rules that consume them in
# primary-coding-mpi.R.
#
# Both are arguments of build_mpi_df(), so a study that codes its variables
# differently or adapts the cut-offs to its own data supplies its own tables
# instead of editing this file. mpi_validate_config() checks such a table, and
# build_mpi_df() runs it on whatever it is handed.
#
# These replace the FileGrouping.rds / Scoring.rds files the original analysis
# script read from disk. Copies of those are kept in data-raw/ for
# data-raw/compare-mpi-implementations.R.
#
# ------------------------------------------------------------------------------


#' Stack the domain blocks, filling in the columns a block leaves out
#'
#' Lets each domain below list only the columns it actually uses.
#'
#' @param ... Data frames, one per domain.
#' @return Data frame with the full set of domain columns.
#' @noRd
mpi_bind <- function(...) {

  columns <- c("domain", "table", "variable",
               "start", "end", "companions", "aggregate", "multi", "floor")

  blocks <- lapply(list(...), function(block) {
    for (column in columns) {
      if (!column %in% names(block)) {
        block[[column]] <- switch(column,
                                  aggregate = "max",
                                  multi = "na",
                                  floor = 0,
                                  NA_character_)
      }
      if (!column %in% c("aggregate", "floor")) {
        block[[column]] <- as.character(block[[column]])
      }
    }
    block[["floor"]] <- as.numeric(block[["floor"]])
    block[, columns]
  })

  out <- do.call(rbind, blocks)
  rownames(out) <- NULL

  return(out)
}


#' Variables feeding the Multidimensional Prognostic Index (MPI)
#'
#' @description
#' One row per domain variable, naming the table it lives in and how its
#' records are matched to a visit. Passed to [build_mpi_df()] as `domains`.
#'
#' @details
#' A variable carrying a validity interval contributes while
#' `start <= visit date < end`. When `start` and `end` are NA the companions
#' `"<variable>_start"` and `"<variable>_end"` are used if the table has them,
#' and otherwise the variable is not treated as interval data: it contributes
#' to its own visit when the table carries a visit id, and to every visit of
#' the patient when it does not.
#'
#' @format Data frame with the columns:
#' \describe{
#'   \item{domain}{Domain the variable contributes to. Must match a domain of
#'     [mpi_scoring].}
#'   \item{table}{Table holding the variable, file name without ".csv".}
#'   \item{variable}{Variable name.}
#'   \item{start}{Variable holding the start of the validity interval, or NA.}
#'   \item{end}{Variable holding the end of the validity interval, or NA.}
#'   \item{companions}{Comma separated further variables that share the
#'     interval of `variable` and are needed for scoring, or NA.}
#'   \item{aggregate}{How to combine several records valid at the same visit:
#'     "max" for the highest code, "count" for the number of distinct
#'     entries.}
#'   \item{multi}{How to read a cell that holds several codes separated by
#'     "|": "max" reduces it to the highest code, "na" reads it as missing.
#'     Only the social variables combine, so everything else is "na" and an
#'     ambiguous cell is reported rather than quietly collapsed.}
#'   \item{floor}{The code that means "no problem". A "sum" domain adds
#'     `max(0, value - floor)`, so a graded scale starting at 1 contributes
#'     nothing when it is answered at its best. 0 for a yes/no flag, which
#'     leaves it unchanged.}
#' }
#'
#' @seealso [mpi_scoring] for the cut-offs, [mpi_validate_config()] to check a
#'   table of your own, [mpi_variable_report()] to see how it resolves against
#'   an export.
#' @export
mpi_domains <- mpi_bind(

  ## Comorbidity ---------------------------------------------------------------
  data.frame(
    domain = "comorbidity",
    table = c("fv2_1",
              "fv2_2", "fv2_2", "fv2_2", "fv2_2", "fv2_2",
              "fv2_3", "fv2_3", "fv2_3", "fv2_3", "fv2_3",
              "fv2_4", "fv2_4", "fv2_4",
              "fv2_5",
              "fv2_6"),
    variable = c("gec_cv",
                 "gec_pd", "gec_ckd", "aki", "gec_ulc", "pankr",
                 "gec_on_solid", "gec_on_hema", "gec_hiv", "gec_rh", "gec_dm",
                 "gec_ne", "gec_ne_dem", "de",
                 "gec_tx",
                 "diag_oth"),
    stringsAsFactors = FALSE
  ),

  ## Cognition -----------------------------------------------------------------
  data.frame(
    domain = "cognition",
    table = c("fv3_3", "fv3_3", "prom", "prom", "prom", "promext"),
    variable = c("sy_ne_cogn", "gec_sy_ne_conf",
                 "pro_pc2r", "pro_pc35r", "pro_pc36r",
                 "cfs_seid_crit4"),
    stringsAsFactors = FALSE
  ),

  ## Nutrition -----------------------------------------------------------------
  # Scored from BMI, so height and weight are the inputs rather than a sum.
  data.frame(
    domain = "nutrition",
    table = c("bv1", "fuv1"),
    variable = c("gec_height", "gec_weight"),
    stringsAsFactors = FALSE
  ),

  ## Activities of daily living ------------------------------------------------
  data.frame(
    domain = "adl",
    table = c("fv3_3", "prom", "promext", "promext"),
    variable = c("sy_ne_atax", "eq5d_mob", "pro_29_pfa21", "pro_29_pfa23"),
    stringsAsFactors = FALSE
  ),

  ## Mental health -------------------------------------------------------------
  # The only domain that mixes yes/no flags with graded scales. EQ-5D and
  # PROMIS-29 items start at 1 for "no problem", so summing their raw codes
  # would make the domain total depend on how many items a patient answered
  # rather than on how the patient is: answering all of them at their best
  # already reaches the 0.5 cut-off. `floor` is the code that means "no
  # problem", and each item contributes max(0, value - floor).
  #
  # pss_* are left at 0 because the coding of the PSS items in this export is
  # not established. Set them to 1 if they turn out to run 1-5 rather than
  # 0-4, and lower the mental cut-offs by four to match.
  data.frame(
    domain = "mental",
    table = c("fv2_4", "fv2_4", "prom",
              "promext", "promext", "promext", "promext", "promext",
              "promext", "promext", "promext", "promext", "promext",
              "promext", "promext"),
    variable = c("gec_ne_depr", "gec_ne_anx", "eq5d_anx",
                 "pro_29_edanx01", "pro_29_edanx40", "pro_29_edanx41",
                 "pro_29_edanx53", "pro_29_eddep04", "pro_29_eddep06",
                 "pro_29_eddep29", "pro_29_eddep41",
                 "pss_1", "pss_2", "pss_3", "pss_4"),
    floor = c(0, 0, 1,
              1, 1, 1, 1, 1,
              1, 1, 1, 0, 0,
              0, 0),
    stringsAsFactors = FALSE
  ),

  ## Social / living situation -------------------------------------------------
  # home_adults and home_childr share the validity interval of home and are
  # needed by the living_situation rules.
  # The only domain whose cells may hold several codes at once. `multi = "max"`
  # keeps the highest of them, which for these ordered categories is the least
  # independent living situation. Everywhere else such a cell reads as missing.
  data.frame(
    domain = "social",
    table = c("ehome", "eresid"),
    variable = c("home", "resid"),
    companions = c("home_adults,home_childr", NA),
    multi = "max",
    stringsAsFactors = FALSE
  ),

  ## Medications ---------------------------------------------------------------
  # Interval data: a drug counts towards a visit when the visit date falls in
  # [start, end). Several drugs are typically active at once, hence "count".
  data.frame(
    domain = "medications",
    table = c("emed_i", "emed_ih", "emed_oth", "emed_oth"),
    variable = c("gec_med_i_drug", "med_ih_drug",
                 "gec_med_drug", "gec_med_drug_oth"),
    start = c("gec_med_i_start", "gec_med_ih_start",
              "gec_med_start", "gec_med_start"),
    end = c("med_i_end", "med_ih_end", "gec_med_end", "gec_med_end"),
    aggregate = "count",
    stringsAsFactors = FALSE
  ),

  ## Biomarkers ----------------------------------------------------------------
  # Collected but not scored: no cut-offs are defined for this domain, so it
  # contributes NA to every patient. Kept so that cut-offs can be added later
  # without changing the variable list.
  data.frame(
    domain = "biomarker",
    table = c(rep("stv2", 5), rep("fuv2", 3), rep("stv2", 17)),
    variable = c("gec_albumin", "gec_tropt", "gec_tropt_est",
                 "gec_tropi", "gec_tropi_est",
                 "tsh", "tsh_est", "bc_hba1c",
                 "gec_crp", "gec_crp_est", "gec_pct", "gec_pct_est",
                 "gec_ferritin", "gec_il6", "gec_il6_est",
                 "gec_coag_ddimer", "gec_coag_fibrinog", "gec_bilitotal",
                 "gec_crea", "urea", "gec_ast", "gec_alt",
                 "bloodcount", "gec_bc_hb", "gec_bc_wbc"),
    stringsAsFactors = FALSE
  )
)


#' Scoring rules of the Multidimensional Prognostic Index (MPI)
#'
#' @description
#' One row per domain, fixing both the rule that turns the domain variables
#' into a 0 / 0.5 / 1 domain score and the order of the domain score columns.
#' Passed to [build_mpi_df()] as `scoring`.
#'
#' @details
#' The cut-offs below are the ones adapted to the variables a NAPKON SUEP
#' export actually carries. They are not the cut-offs of the original MPI,
#' which scores instruments this export does not collect. A study with other
#' data is expected to adapt them again and pass its own table.
#'
#' `type` selects the rule:
#' \describe{
#'   \item{sum}{Sum the domain variables, then compare against the cut-offs. A
#'     cut-off is the *inclusive lower bound* of its score, so the highest
#'     cut-off that is still `<=` the sum decides. This matches the bands of
#'     the original MPI, where a comorbidity count of 0 scores 0 rather than
#'     going unscored. Each variable contributes `max(0, value - floor)`, see
#'     the `floor` column of [mpi_domains].}
#'   \item{bmi}{Derive BMI from `gec_weight` and `gec_height` and read the
#'     cut-offs literally as `"lo-hi"` and `"<x,>y"` bands. Only `cut_0.5` and
#'     `cut_1` are evaluated; `cut_0` records the normal band for the reader.}
#'   \item{living_situation}{Apply the `home` / `resid` decision rules. The
#'     cut-off columns are unused.}
#'   \item{none}{No cut-offs defined, the domain scores NA throughout.}
#' }
#'
#' @format Data frame with the columns:
#' \describe{
#'   \item{domain}{Domain name. Also fixes the order of the domain score
#'     columns.}
#'   \item{type}{Scoring rule: "sum", "bmi", "living_situation" or "none".}
#'   \item{cut_0}{Cut-off for domain score 0, no problem.}
#'   \item{cut_0.5}{Cut-off for domain score 0.5, minor problem.}
#'   \item{cut_1}{Cut-off for domain score 1, major problem.}
#' }
#'
#' @seealso [mpi_domains] for the variables, [mpi_validate_config()] to check a
#'   table of your own.
#' @export
mpi_scoring <- data.frame(
  domain  = c("comorbidity", "cognition", "nutrition", "adl",
              "mental", "social", "medications", "biomarker"),
  type    = c("sum", "sum", "bmi", "sum",
              "sum", "living_situation", "sum", "none"),
  cut_0   = c("0", "0", "23-33", "0",
              "0", NA, "0", NA),
  cut_0.5 = c("1", "4", "<23,>33", "4",
              "3", NA, "2", NA),
  cut_1   = c("3", "7", "<20,>38", "7",
              "6", NA, "4", NA),
  stringsAsFactors = FALSE
)


#' Check an MPI configuration
#'
#' @description
#' Checks that a domain definition and a scoring definition fit together and
#' can be scored, and stops naming what is wrong when they cannot.
#'
#' @details
#' [build_mpi_df()] runs this on whatever it is given, so a study supplying its
#' own tables finds a mistake immediately rather than through a domain that
#' quietly scores NA for every patient.
#'
#' Cut-offs are only checked where they are read: numeric and non-decreasing
#' for a "sum" domain, parseable as a band for the `cut_0.5` and `cut_1` of a
#' "bmi" domain. The columns of a "living_situation" or "none" domain are not
#' looked at.
#'
#' @param domains Domain definition, see [mpi_domains].
#' @param scoring Domain scoring definition, see [mpi_scoring].
#'
#' @return `TRUE`, invisibly. Called for the error it raises.
#'
#' @examples
#' mpi_validate_config()
#' @export
mpi_validate_config <- function(domains = epicodr::mpi_domains,
                                scoring = epicodr::mpi_scoring) {

  fail <- function(...) stop(paste0("mpi_validate_config(): ", ...),
                             call. = FALSE)
  listed <- function(x) paste(unique(x), collapse = ", ")

  ## The columns the collectors and the scoring rules read ---------------------
  missing_columns <- setdiff(
    c("domain", "table", "variable", "start", "end", "companions", "aggregate"),
    names(domains)
  )
  if (length(missing_columns) > 0) {
    fail("`domains` has no column(s) ", listed(missing_columns), ".")
  }

  missing_columns <- setdiff(
    c("domain", "type", "cut_0", "cut_0.5", "cut_1"),
    names(scoring)
  )
  if (length(missing_columns) > 0) {
    fail("`scoring` has no column(s) ", listed(missing_columns), ".")
  }

  ## Values the switches in build_mpi_df() can act on --------------------------
  unknown <- setdiff(domains$aggregate, c("max", "count"))
  if (length(unknown) > 0) {
    fail("unknown aggregate ", listed(unknown), ". Use \"max\" or \"count\".")
  }

  unknown <- setdiff(scoring$type, c("sum", "bmi", "living_situation", "none"))
  if (length(unknown) > 0) {
    fail("unknown scoring type ", listed(unknown), ".")
  }

  if ("multi" %in% names(domains)) {
    unknown <- setdiff(stats::na.omit(domains$multi), c("max", "na"))
    if (length(unknown) > 0) {
      fail("unknown multi ", listed(unknown), ". Use \"max\" to reduce a cell ",
           "holding several codes to the highest one, or \"na\" to read it ",
           "as missing.")
    }
  }

  # `floor` is optional, so a table written before it existed still validates.
  if ("floor" %in% names(domains)) {
    floors <- suppressWarnings(as.numeric(domains$floor))
    bad <- is.na(floors) | floors < 0
    if (any(bad)) {
      fail("the floor of ", listed(domains$variable[bad]),
           " must be a number of at least 0: it is the code that means ",
           "\"no problem\", so it can be neither missing nor negative.")
    }
  }

  ## The two tables have to line up --------------------------------------------
  duplicated_domains <- scoring$domain[duplicated(scoring$domain)]
  if (length(duplicated_domains) > 0) {
    fail("`scoring` has more than one rule for domain(s) ",
         listed(duplicated_domains), ".")
  }

  unknown <- setdiff(domains$domain, scoring$domain)
  if (length(unknown) > 0) {
    fail("no scoring rule for domain(s) ", listed(unknown), ".")
  }

  # A scored domain without variables can only ever score NA, which is almost
  # always a typo in a domain name rather than an intention.
  unfed <- setdiff(scoring$domain[scoring$type != "none"], domains$domain)
  if (length(unfed) > 0) {
    fail("no variables for scored domain(s) ", listed(unfed),
         ". Give them variables in `domains`, or set their type to \"none\".")
  }

  ## Cut-offs, where they are read ---------------------------------------------
  for (i in which(scoring$type == "sum")) {
    cuts <- suppressWarnings(as.numeric(
      c(scoring$cut_0[i], scoring$cut_0.5[i], scoring$cut_1[i])
    ))
    if (anyNA(cuts)) {
      fail("the cut-offs of the sum domain '", scoring$domain[i],
           "' must all be numbers, got ",
           paste(c(scoring$cut_0[i], scoring$cut_0.5[i], scoring$cut_1[i]),
                 collapse = " / "), ".")
    }
    if (is.unsorted(cuts)) {
      fail("the cut-offs of the sum domain '", scoring$domain[i],
           "' decrease from cut_0 to cut_1 (", paste(cuts, collapse = " / "),
           "). A cut-off is the inclusive lower bound of its score, so they ",
           "must not go down.")
    }
  }

  for (i in which(scoring$type == "bmi")) {
    for (column in c("cut_0.5", "cut_1")) {
      band <- mpi_parse_band(scoring[[column]][i])
      if (is.na(band$lower) && is.na(band$upper)) {
        fail("the ", column, " of the bmi domain '", scoring$domain[i],
             "' holds no band. Expected something like \"<23,>33\".")
      }
    }
  }

  invisible(TRUE)
}
