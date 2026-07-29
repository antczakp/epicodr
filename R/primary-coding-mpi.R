# ------------------------------------------------------------------------------
# MULTIDIMENSIONAL PROGNOSTIC INDEX (MPI)
# ------------------------------------------------------------------------------
#
# The MPI aggregates eight domains into a single 0-1 index. Every domain is
# scored 0 (no problem), 0.5 (minor problem) or 1 (major problem) and the MPI is
# the mean of the available domain scores.
#
# Which variables feed a domain is defined by `mpi_domains`, how a domain is
# scored by `mpi_scoring`. Both ship as package data and can be replaced by the
# caller.
#
# Domain scores are evaluated at visit dates. A variable contributes the value
# that is valid at the visit:
#
#   interval data      records carrying a start/end date contribute while
#                      start <= visit date < end
#   longitudinal data  records carrying a visit id contribute to that visit
#   cross sectional    records without a visit id contribute to every visit
#
# ------------------------------------------------------------------------------


# Low level helpers ------------------------------------------------------------

#' Coerce an export variable to numeric codes
#'
#' Handles the representations a Transferstelle export variable can have after
#' `read_tsExport()`: plain numeric, labelled numeric, factor and character.
#' Cells holding several codes separated by "|" collapse to their maximum.
#'
#' @param x Vector to coerce.
#' @return Double vector.
#' @noRd
mpi_as_numeric <- function(x) {

  if (is.numeric(x) || is.logical(x)) {
    return(as.numeric(x))
  }

  if (is.factor(x)) {
    x <- as.character(x)
  }

  x <- as.character(x)

  # Multi select cells such as "1|3" keep the highest code, as in the original
  # MPI script.
  if (any(grepl("|", x, fixed = TRUE), na.rm = TRUE)) {
    return(vapply(
      strsplit(x, "|", fixed = TRUE),
      function(codes) {
        codes <- suppressWarnings(as.numeric(codes))
        if (all(is.na(codes))) NA_real_ else max(codes, na.rm = TRUE)
      },
      numeric(1)
    ))
  }

  suppressWarnings(as.numeric(x))
}


#' Coerce an export variable to Date
#'
#' Accepts the raw "%Y%m%d" representation as well as variables that have
#' already been converted by `dates_tsExport()`.
#'
#' @param x Vector to coerce.
#' @return Date vector.
#' @noRd
mpi_as_date <- function(x) {

  if (inherits(x, "Date")) {
    return(x)
  }

  if (inherits(x, "POSIXt")) {
    return(as.Date(x))
  }

  if (is.factor(x)) {
    x <- as.character(x)
  }

  x <- as.character(x)
  out <- as.Date(x, format = "%Y%m%d")

  # Fall back to ISO dates for exports that were converted to character.
  retry <- is.na(out) & !is.na(x)
  if (any(retry)) {
    out[retry] <- as.Date(x[retry], format = "%Y-%m-%d")
  }

  return(out)
}


#' Replace generic missing codes by NA
#'
#' @param x Numeric vector.
#' @param missing_values Numeric vector of codes that stand for "missing".
#' @return Numeric vector.
#' @noRd
mpi_mask_missing <- function(x, missing_values) {

  if (length(missing_values) == 0) {
    return(x)
  }

  x[x %in% missing_values] <- NA_real_
  return(x)
}


#' Split a comma separated configuration entry
#'
#' @param x Length one character vector, possibly NA.
#' @return Character vector, empty when `x` is NA or blank.
#' @noRd
mpi_split_config <- function(x) {

  if (length(x) == 0 || is.na(x) || !nzchar(trimws(x))) {
    return(character(0))
  }

  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}


#' Resolve the variable holding an interval bound
#'
#' Uses the configured name when given. Otherwise the table is searched, first
#' for the exact "<variable>_start" / "<variable>_end" companion and then for
#' any column starting with "<variable>_st" / "<variable>_en".
#'
#' The prefix search reproduces what the original MPI script did. Interval
#' bounds are not spelled consistently across the export, which is why they
#' were searched for rather than listed: `emed_i` pairs `gec_med_i_start` with
#' `med_i_end`, and only some of the bounds ever made it into the domain
#' definition.
#'
#' @param configured Configured variable name or NA.
#' @param variable Name of the variable the bound belongs to.
#' @param suffix Either "start" or "end".
#' @param table_names Column names of the table.
#' @return Length one character vector, NA when no bound is available.
#' @noRd
mpi_resolve_bound <- function(configured, variable, suffix, table_names) {

  if (!is.na(configured) && nzchar(configured)) {
    if (configured %in% table_names) {
      return(configured)
    }
    return(NA_character_)
  }

  exact <- paste(variable, suffix, sep = "_")
  if (exact %in% table_names) {
    return(exact)
  }

  # "start" -> "<variable>_st", "end" -> "<variable>_en".
  prefix <- paste0(variable, "_", substr(suffix, 1, 2))
  matches <- table_names[startsWith(table_names, prefix)]

  if (length(matches) == 0) {
    return(NA_character_)
  }

  if (length(matches) > 1) {
    warning(paste0("MPI: several candidates for the ", suffix, " of '",
                   variable, "': ", paste(matches, collapse = ", "),
                   ". Using '", matches[1],
                   "'. Set the bound explicitly in `domains` to choose."))
  }

  return(matches[1])
}


#' Aggregate the records that are valid at a visit
#'
#' @param values Values of the records, one element per record.
#' @param keys Visit key of each record.
#' @param key_levels Visit keys to return values for, in order.
#' @param aggregate Either "max" or "count".
#' @param missing_values Codes that stand for "missing".
#' @return Double vector, one element per element of `key_levels`.
#' @noRd
mpi_aggregate_records <- function(values, keys, key_levels, aggregate,
                                  missing_values) {

  groups <- split(seq_along(keys), factor(keys, levels = key_levels))

  if (aggregate == "count") {
    # Number of distinct entries that are active at the visit. Cells may hold
    # several entries separated by "|". A visit without any active record is a
    # count of zero, not a missing value.
    values <- as.character(values)
    counts <- vapply(groups, function(idx) {
      entries <- unlist(strsplit(values[idx], "|", fixed = TRUE))
      entries <- entries[!is.na(entries)]
      entries <- entries[nzchar(trimws(entries))]
      length(unique(entries))
    }, numeric(1))
    return(unname(counts))
  }

  values <- mpi_mask_missing(mpi_as_numeric(values), missing_values)
  maxima <- vapply(groups, function(idx) {
    group_values <- values[idx]
    if (all(is.na(group_values))) NA_real_ else max(group_values, na.rm = TRUE)
  }, numeric(1))

  return(unname(maxima))
}


#' Carry values forward and backward within a patient
#'
#' Used for the nutrition domain, where height and weight are typically
#' recorded once but the resulting BMI applies to every visit of the patient.
#'
#' @param x Numeric vector, ordered by visit date within patient.
#' @param patient Patient identifier of each element.
#' @return Numeric vector of the same length.
#' @noRd
mpi_fill_within_patient <- function(x, patient) {

  if (length(x) != length(patient)) {
    x <- rep(x, length.out = length(patient))
  }

  for (idx in split(seq_along(x), patient)) {
    values <- x[idx]
    observed <- which(!is.na(values))
    if (length(observed) == 0 || length(observed) == length(values)) {
      next
    }
    # Nearest observed value: last one carried forward, first one carried back.
    nearest <- vapply(seq_along(values), function(i) {
      before <- observed[observed <= i]
      if (length(before) > 0) before[length(before)] else observed[1]
    }, numeric(1))
    x[idx] <- values[nearest]
  }

  return(x)
}


#' Resolve the patient and visit ID column names
#'
#' Falls back to `set_id_names()` only when the caller did not name both IDs.
#' Exports that ship without the secuTrial metadata tables cannot be read by
#' `read_tsExport()` and therefore have no `export_options`, but they can still
#' be scored by naming the two ID columns directly.
#'
#' @param trial_data A secuTrial data object or a named list of data frames.
#' @param pid Column name of the patient ID, or NULL.
#' @param visitid Column name of the visit ID, or NULL.
#' @return List with the entries `pid` and `visitid`.
#' @noRd
mpi_resolve_ids <- function(trial_data, pid, visitid) {

  if (!is.null(pid) && !is.null(visitid)) {
    return(list(pid = pid, visitid = visitid))
  }

  # set_id_names() derives the ID names from export_options, so it can only run
  # on an object that read_tsExport() produced.
  if (!("id_names" %in% names(trial_data$export_options))) {
    if ("data_names" %in% names(trial_data$export_options)) {
      trial_data <- set_id_names(trial_data)
    }
  }

  if (!("id_names" %in% names(trial_data$export_options))) {
    stop(paste0(
      "Could not determine the ID columns. This object has no ",
      "\"export_options\", which happens when an export ships without the ",
      "secuTrial metadata tables and was read outside read_tsExport(). ",
      "Pass the columns directly, e.g. ",
      "pid = \"export_psn\", visitid = \"mnpcvpid\"."
    ))
  }

  list(
    pid = if (is.null(pid)) trial_data$export_options$id_names$pid else pid,
    visitid = if (is.null(visitid)) {
      trial_data$export_options$id_names$visitid
    } else {
      visitid
    }
  )
}


# Visit table ------------------------------------------------------------------

#' Default visit forms of a SUEP export
#'
#' @return Named list of visit form definitions.
#' @noRd
mpi_default_visit_forms <- function() {
  list(
    scv = list(date = "gec_pr_incl_date", mode = NA_character_),
    m   = list(date = "pr_docudate",      mode = "pr_visit_mode"),
    m2  = list(date = "gec_pr_docudate_1", mode = "pr_visit_mode_1")
  )
}


#' Build the visit table underlying the MPI
#'
#' @description
#' Collects patient, visit and documentation date from the visit forms of an
#' export into a single table. The MPI and its outcomes are evaluated at the
#' visit dates returned here.
#'
#' @param trial_data A secuTrial data object.
#' @param pid Column name of the patient ID in trial_data.
#' @param visitid Column name of the visit ID in trial_data.
#' @param visit_forms Named list describing the visit forms. Each element is a
#'   list with the entries `date` and `mode` naming the documentation date and
#'   the visit mode variable of that form.
#'
#' @return Data frame with one row per visit and the columns patient ID, visit
#'   ID, `mnpvisid`, `mnpvislabel`, `mnpvisno`, `visit_mode` and `visit_date`.
#' @export
build_mpi_visit_df <- function(trial_data, pid = NULL, visitid = NULL,
                               visit_forms = mpi_default_visit_forms()) {

  ids <- mpi_resolve_ids(trial_data, pid, visitid)
  pid <- ids$pid
  visitid <- ids$visitid

  carry <- c("mnpvisid", "mnpvislabel", "mnpvisno")

  collected <- list()

  for (form_name in names(visit_forms)) {

    form <- trial_data[[form_name]]
    spec <- visit_forms[[form_name]]

    if (is.null(form) || !all(c(pid, visitid) %in% names(form))) {
      next
    }
    if (is.na(spec$date) || !(spec$date %in% names(form))) {
      warning(paste0("build_mpi_visit_df(): no date variable '", spec$date,
                     "' in form '", form_name, "'. Form skipped."))
      next
    }

    rows <- data.frame(
      .patient = as.character(form[[pid]]),
      .visit_key = as.character(form[[visitid]]),
      stringsAsFactors = FALSE
    )

    for (column in carry) {
      rows[[column]] <- if (column %in% names(form)) {
        as.character(form[[column]])
      } else {
        NA_character_
      }
    }

    rows$visit_mode <- if (!is.na(spec$mode) && spec$mode %in% names(form)) {
      as.character(form[[spec$mode]])
    } else {
      NA_character_
    }

    rows$visit_date <- mpi_as_date(form[[spec$date]])
    rows$visit_form <- form_name

    collected[[form_name]] <- rows
  }

  if (length(collected) == 0) {
    stop("build_mpi_visit_df(): none of the visit forms were found in the export.")
  }

  visits <- do.call(rbind, c(collected, list(make.row.names = FALSE)))

  # A visit may show up in more than one form. Prefer the record that carries a
  # documentation date, then drop the remaining duplicates.
  visits <- visits[order(visits$.visit_key, is.na(visits$visit_date)), ]
  visits <- visits[!duplicated(visits$.visit_key), ]
  visits <- visits[order(visits$.patient, visits$visit_date, visits$.visit_key), ]
  rownames(visits) <- NULL

  names(visits)[names(visits) == ".patient"] <- pid
  names(visits)[names(visits) == ".visit_key"] <- visitid

  return(visits)
}


# Domain values ----------------------------------------------------------------

#' Collect the values of one configured variable at every visit
#'
#' @param trial_data A secuTrial data object.
#' @param config One row of the domain configuration.
#' @param visits Visit table with the helper columns `.patient`, `.key`, `.date`.
#' @param pid Column name of the patient ID in trial_data.
#' @param visitid Column name of the visit ID in trial_data.
#' @param missing_values Codes that stand for "missing".
#' @return Named list of double vectors, one per collected variable, or NULL.
#' @noRd
mpi_collect_variable <- function(trial_data, config, visits, pid, visitid,
                                 missing_values) {

  form <- trial_data[[config$table]]

  if (is.null(form)) {
    warning(paste0("MPI: table '", config$table,
                   "' is not in the export. Variable '", config$variable,
                   "' is treated as missing."))
    return(NULL)
  }

  form_names <- names(form)
  variables <- c(config$variable, mpi_split_config(config$companions))
  available <- variables[variables %in% form_names]

  if (length(available) == 0 || !(pid %in% form_names)) {
    warning(paste0("MPI: variable '", config$variable, "' is not in table '",
                   config$table, "'. It is treated as missing."))
    return(NULL)
  }

  start_var <- mpi_resolve_bound(config$start, config$variable, "start",
                                 form_names)
  end_var <- mpi_resolve_bound(config$end, config$variable, "end", form_names)
  is_interval <- !is.na(start_var) || !is.na(end_var)

  key_levels <- visits$.key

  if (is_interval) {

    records <- data.frame(
      .patient = as.character(form[[pid]]),
      stringsAsFactors = FALSE
    )
    records$.start <- if (is.na(start_var)) {
      as.Date(rep(NA, nrow(form)))
    } else {
      mpi_as_date(form[[start_var]])
    }
    records$.end <- if (is.na(end_var)) {
      as.Date(rep(NA, nrow(form)))
    } else {
      mpi_as_date(form[[end_var]])
    }
    for (variable in available) {
      records[[variable]] <- form[[variable]]
    }

    matched <- merge(
      visits[, c(".patient", ".key", ".date")],
      records,
      by = ".patient"
    )

    # A record is valid at a visit when start <= visit date < end. An open
    # bound does not restrict the interval.
    valid <-
      (is.na(matched$.start) | matched$.start <= matched$.date) &
      (is.na(matched$.end) | matched$.date < matched$.end)
    matched <- matched[valid, , drop = FALSE]

  } else if (visitid %in% form_names) {

    # Longitudinal form: a record belongs to the visit it was documented in.
    matched <- data.frame(
      .key = as.character(form[[visitid]]),
      stringsAsFactors = FALSE
    )
    for (variable in available) {
      matched[[variable]] <- form[[variable]]
    }
    matched <- matched[matched$.key %in% key_levels, , drop = FALSE]

  } else {

    # Cross sectional form: a record applies to every visit of the patient.
    records <- data.frame(
      .patient = as.character(form[[pid]]),
      stringsAsFactors = FALSE
    )
    for (variable in available) {
      records[[variable]] <- form[[variable]]
    }

    matched <- merge(visits[, c(".patient", ".key")], records, by = ".patient")
  }

  values <- lapply(available, function(variable) {
    mpi_aggregate_records(
      values = matched[[variable]],
      keys = matched$.key,
      key_levels = key_levels,
      aggregate = config$aggregate,
      missing_values = missing_values
    )
  })
  names(values) <- available

  return(values)
}


#' Collect the values of every variable of one domain
#'
#' @inheritParams mpi_collect_variable
#' @param domain_config Domain configuration rows of a single domain.
#' @return Data frame with one row per visit and one column per variable.
#' @noRd
mpi_collect_domain <- function(trial_data, domain_config, visits, pid, visitid,
                               missing_values) {

  collected <- list()

  for (i in seq_len(nrow(domain_config))) {
    values <- mpi_collect_variable(
      trial_data = trial_data,
      config = domain_config[i, ],
      visits = visits,
      pid = pid,
      visitid = visitid,
      missing_values = missing_values
    )
    if (is.null(values)) {
      next
    }
    for (variable in names(values)) {
      collected[[variable]] <- values[[variable]]
    }
  }

  if (length(collected) == 0) {
    return(NULL)
  }

  return(as.data.frame(collected, stringsAsFactors = FALSE))
}


# Domain scoring ---------------------------------------------------------------

#' Parse a "<x,>y" cut-off band
#'
#' @param x Length one character vector.
#' @return List with the entries `lower` and `upper`.
#' @noRd
mpi_parse_band <- function(x) {

  empty <- list(lower = NA_real_, upper = NA_real_)

  if (length(x) == 0 || is.na(x) || !nzchar(trimws(x))) {
    return(empty)
  }

  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  lower <- suppressWarnings(as.numeric(sub("^<", "", parts[grepl("^<", parts)])))
  upper <- suppressWarnings(as.numeric(sub("^>", "", parts[grepl("^>", parts)])))

  list(
    lower = if (length(lower) > 0) lower[1] else NA_real_,
    upper = if (length(upper) > 0) upper[1] else NA_real_
  )
}


#' Score a domain from the sum of its variables
#'
#' @description
#' Sums the domain variables and compares the sum against the cut-offs. A
#' cut-off is the inclusive lower bound of its score, so the highest cut-off
#' that is still smaller than or equal to the sum decides the score. A visit
#' without a single observed variable scores NA rather than 0.
#'
#' @param values Data frame of the domain variables, one row per visit.
#' @param cuts Named numeric vector of cut-offs, names are the domain scores.
#' @return Double vector of domain scores.
#' @noRd
mpi_score_sum <- function(values, cuts) {

  values <- as.matrix(values)
  observed <- rowSums(!is.na(values))
  total <- rowSums(values, na.rm = TRUE)

  cuts <- cuts[!is.na(cuts)]
  cuts <- cuts[order(cuts)]

  score <- rep(NA_real_, nrow(values))
  for (i in seq_along(cuts)) {
    score[total >= cuts[i]] <- as.numeric(names(cuts)[i])
  }

  score[observed == 0] <- NA_real_

  return(score)
}


#' Score the nutrition domain from BMI
#'
#' @description
#' Reads the cut-offs literally: a BMI below `cut_1`'s lower bound or above its
#' upper bound scores 1, a BMI below `cut_0.5`'s lower bound or above its upper
#' bound scores 0.5, anything else scores 0.
#'
#' @param bmi Numeric vector of body mass index values.
#' @param cuts Named character vector of cut-offs.
#' @return Double vector of domain scores.
#' @noRd
mpi_score_bmi <- function(bmi, cuts) {

  minor <- mpi_parse_band(cuts[["0.5"]])
  major <- mpi_parse_band(cuts[["1"]])

  score <- rep(NA_real_, length(bmi))
  observed <- !is.na(bmi)
  score[observed] <- 0

  if (!is.na(minor$lower)) score[observed & bmi < minor$lower] <- 0.5
  if (!is.na(minor$upper)) score[observed & bmi > minor$upper] <- 0.5
  if (!is.na(major$lower)) score[observed & bmi < major$lower] <- 1
  if (!is.na(major$upper)) score[observed & bmi > major$upper] <- 1

  return(score)
}


#' Score the social domain from the living situation
#'
#' @description
#' Applies the decision rules of the original MPI script. `home` distinguishes
#' living at home (2) from institutional care (1), `resid` records the type of
#' residence.
#'
#' `resid` in (2, 3) cannot occur together with a known `home`: a patient who
#' has a home is not in one of those residence types, so the combination is
#' absent from the data rather than unhandled. It scores NA, as does any other
#' combination the rules do not cover.
#'
#' @param home,home_adults,home_childr,resid Numeric vectors of equal length.
#' @return Double vector of domain scores.
#' @noRd
mpi_score_living_situation <- function(home, home_adults, home_childr, resid) {

  # Comparison that yields FALSE instead of NA, so it is safe for indexing.
  is_code <- function(x, code) !is.na(x) & x %in% code

  has_home <- !is.na(home)
  has_resid <- !is.na(resid)
  adults_known <- !is.na(home_adults)

  # An unknown number of children is read as "none reported", which keeps the
  # rules decidable whenever the number of adults is known.
  children <- ifelse(is.na(home_childr), 0, home_childr)

  score <- rep(NA_real_, length(home))

  ## Living situation known, type of residence unknown
  home_only <- !has_resid & has_home
  score[home_only & is_code(home, 2) & adults_known &
          (home_adults > 0 | children > 0)] <- 0
  score[home_only & is_code(home, 2) & adults_known &
          home_adults == 0 & children == 0] <- 0.5
  score[home_only & is_code(home, 2) & !adults_known] <- 0.5
  score[home_only & is_code(home, 1)] <- 1

  ## Type of residence known, living situation unknown
  resid_only <- has_resid & !has_home
  score[resid_only & is_code(resid, 1)] <- 1
  score[resid_only & is_code(resid, c(2, 3))] <- 0.5
  score[resid_only & is_code(resid, 4)] <- 1

  ## Both known
  both <- has_resid & has_home
  score[both & is_code(resid, 4)] <- 1
  score[both & is_code(resid, 1) & is_code(home, 2)] <- 0
  score[both & is_code(resid, 1) & is_code(home, 1)] <- 1

  return(score)
}


#' Categorize the Multidimensional Prognostic Index
#'
#' @description
#' Assigns the three conventional MPI risk groups: MPI-1 (low risk, up to
#' 0.33), MPI-2 (moderate risk, up to 0.66) and MPI-3 (high risk, above 0.66).
#'
#' @param mpi Numeric vector of MPI values between 0 and 1.
#' @return Factor of MPI risk groups.
#' @export
categorize_mpi_ecu <- function(mpi) {

  factor(
    ifelse(is.na(mpi), NA_character_,
           ifelse(mpi <= 0.33, "MPI-1",
                  ifelse(mpi <= 0.66, "MPI-2", "MPI-3"))),
    levels = c("MPI-1", "MPI-2", "MPI-3")
  )
}


# Main entry points ------------------------------------------------------------

#' Build the Multidimensional Prognostic Index table
#'
#' @description
#' Computes the MPI and its domain scores for every visit of every patient.
#'
#' Each domain is scored 0 (no problem), 0.5 (minor problem) or 1 (major
#' problem) and the MPI is the mean of the available domain scores, so it
#' ranges from 0 to 1. Domains without observed data are skipped, which is why
#' `ecu_mpi_n_domains` reports how many domains a value rests on.
#'
#' Averaging over too few domains stops being meaningful, so a visit needs at
#' least `min_domains` scored domains to get an MPI at all. Note that the
#' biomarker domain carries no cut-offs and therefore always scores NA, which
#' leaves seven scorable domains rather than eight.
#'
#' @details
#' Domain variables are read at visit dates. Records carrying a start and end
#' date contribute while `start <= visit date < end`, records carrying a visit
#' id contribute to that visit, and records without a visit id contribute to
#' every visit of the patient.
#'
#' Codes listed in `missing_values` are read as missing before a domain is
#' summed. They are not applied to the nutrition domain, where 88 or 99 are
#' plausible measurements. Pass `numeric(0)` to switch this off.
#'
#' @param trial_data A secuTrial data object.
#' @param pid Column name of the patient ID in trial_data.
#' @param visitid Column name of the visit ID in trial_data.
#' @param domains Domain definition, see [mpi_domains].
#' @param scoring Domain scoring definition, see [mpi_scoring].
#' @param visit_forms Named list describing the visit forms, see
#'   [build_mpi_visit_df()].
#' @param missing_values Numeric codes that stand for "missing".
#' @param min_domains Minimum number of scored domains required for an MPI.
#'   Visits below it get NA for `ecu_mpi` and `ecu_mpi_cat`, while their
#'   individual domain scores and `ecu_mpi_n_domains` are still reported.
#'   Defaults to 5, below which the mean rests on too little of the
#'   assessment to be informative. Pass 1 to score every visit that has any
#'   data at all.
#'
#' @return Data frame with one row per visit, the visit identification columns,
#'   one `ecu_mpi_<domain>` column per domain, `ecu_mpi_n_domains`, `ecu_mpi`
#'   and `ecu_mpi_cat`.
#' @export
build_mpi_df <- function(trial_data, pid = NULL, visitid = NULL,
                         domains = epicodr::mpi_domains,
                         scoring = epicodr::mpi_scoring,
                         visit_forms = mpi_default_visit_forms(),
                         missing_values = epicodr::missing_values_default,
                         min_domains = 5) {

  ids <- mpi_resolve_ids(trial_data, pid, visitid)
  pid <- ids$pid
  visitid <- ids$visitid

  unknown <- setdiff(domains$domain, scoring$domain)
  if (length(unknown) > 0) {
    stop(paste0("build_mpi_df(): no scoring rule for domain(s) ",
                paste(unknown, collapse = ", "), "."))
  }

  visits <- build_mpi_visit_df(trial_data, pid, visitid, visit_forms)

  # Helper columns keep the collectors independent of the actual ID names.
  visits$.patient <- as.character(visits[[pid]])
  visits$.key <- as.character(visits[[visitid]])
  visits$.date <- visits$visit_date

  domain_scores <- list()

  for (i in seq_len(nrow(scoring))) {

    domain <- scoring$domain[i]
    type <- scoring$type[i]
    cuts <- c(
      "0" = scoring$cut_0[i],
      "0.5" = scoring$cut_0.5[i],
      "1" = scoring$cut_1[i]
    )

    if (type == "none") {
      domain_scores[[domain]] <- rep(NA_real_, nrow(visits))
      next
    }

    domain_config <- domains[domains$domain == domain, , drop = FALSE]

    # Missing codes are category codes, so they are not applied to the
    # continuous measurements the BMI rule works on.
    domain_missing <- if (type == "bmi") numeric(0) else missing_values

    values <- mpi_collect_domain(
      trial_data = trial_data,
      domain_config = domain_config,
      visits = visits,
      pid = pid,
      visitid = visitid,
      missing_values = domain_missing
    )

    if (is.null(values)) {
      warning(paste0("build_mpi_df(): no data for domain '", domain,
                     "'. It scores NA for every visit."))
      domain_scores[[domain]] <- rep(NA_real_, nrow(visits))
      next
    }

    domain_scores[[domain]] <- switch(
      type,

      # The names of `cuts` carry the domain scores, so keep them.
      sum = mpi_score_sum(
        values,
        stats::setNames(suppressWarnings(as.numeric(cuts)), names(cuts))
      ),

      bmi = {
        # Height and weight are typically recorded once, so carry them to the
        # remaining visits of the patient before deriving the BMI.
        height <- mpi_fill_within_patient(
          if ("gec_height" %in% names(values)) values$gec_height else NA_real_,
          visits$.patient
        )
        weight <- mpi_fill_within_patient(
          if ("gec_weight" %in% names(values)) values$gec_weight else NA_real_,
          visits$.patient
        )
        mpi_score_bmi(calculate_bmi(weight, height), cuts)
      },

      living_situation = {
        pick <- function(name) {
          if (name %in% names(values)) {
            values[[name]]
          } else {
            rep(NA_real_, nrow(visits))
          }
        }
        mpi_score_living_situation(
          home = pick("home"),
          home_adults = pick("home_adults"),
          home_childr = pick("home_childr"),
          resid = pick("resid")
        )
      },

      stop(paste0("build_mpi_df(): unknown scoring type '", type, "'."))
    )
  }

  # Keep the domain order of the scoring table.
  scores <- as.data.frame(domain_scores, stringsAsFactors = FALSE)
  scores <- scores[, scoring$domain, drop = FALSE]
  score_matrix <- as.matrix(scores)

  n_domains <- rowSums(!is.na(score_matrix))
  mpi <- rowMeans(score_matrix, na.rm = TRUE)
  mpi[is.nan(mpi)] <- NA_real_
  mpi[n_domains < min_domains] <- NA_real_

  names(scores) <- paste0("ecu_mpi_", names(scores))

  result <- cbind(
    visits[, c(pid, visitid, "mnpvisid", "mnpvislabel", "mnpvisno",
               "visit_mode", "visit_date")],
    scores,
    data.frame(
      ecu_mpi_n_domains = n_domains,
      ecu_mpi = mpi,
      ecu_mpi_cat = categorize_mpi_ecu(mpi),
      stringsAsFactors = FALSE
    )
  )
  rownames(result) <- NULL

  return(result)
}


#' Build the MPI outcome table
#'
#' @description
#' Derives days spent on a ward and vital status between the baseline visit and
#' each follow-up visit.
#'
#' @details
#' A day is attributed to exactly one ward. Where ward stays overlap, the later
#' record wins, so `ecu_rehosp` equals the sum of the `ecu_rehosp_<level>`
#' columns.
#'
#' A patient who was not on a ward between the two visits gets a count of zero.
#' Ward stays with an incomplete start or end date cannot be placed on the
#' timeline, so they are skipped and the remaining stays of that patient are
#' still counted. Because such a count is a lower bound rather than a complete
#' one, the patient is marked in `ecu_rehosp_incomplete`. Counts are NA only
#' when the ward table itself is missing from the export.
#'
#' The vital status is derived from the date of death alone, so it stays usable
#' whatever the state of the ward records.
#'
#' @param trial_data A secuTrial data object.
#' @param pid Column name of the patient ID in trial_data.
#' @param visitid Column name of the visit ID in trial_data.
#' @param visits Visit table as returned by [build_mpi_visit_df()]. Built from
#'   `trial_data` when NULL.
#' @param visit_forms Named list describing the visit forms, see
#'   [build_mpi_visit_df()].
#' @param ward_table Name of the table holding the ward stays.
#' @param ward_vars Named character vector naming the ward, start and end
#'   variable in `ward_table`.
#' @param death_table Name of the table holding the date of death.
#' @param death_var Name of the date of death variable in `death_table`.
#' @param baseline_label Visit label of the baseline visit.
#' @param follow_up_labels Visit labels the outcomes are derived for.
#' @param ward_levels Ward codes that get their own day count.
#'
#' @return Data frame with one row per follow-up visit, the visit
#'   identification columns, `ecu_rehosp`, one `ecu_rehosp_<level>` column per
#'   ward code, `ecu_rehosp_incomplete`, `ecu_death` and `ecu_death_date`.
#' @export
build_mpi_outcome_df <- function(trial_data, pid = NULL, visitid = NULL,
                                 visits = NULL,
                                 visit_forms = mpi_default_visit_forms(),
                                 ward_table = "eward",
                                 ward_vars = c(ward = "gec_ward",
                                               start = "ward_start",
                                               end = "ward_end"),
                                 death_table = "fv2_6",
                                 death_var = "death_date",
                                 baseline_label = "Baseline",
                                 follow_up_labels = c("3M Follow-Up",
                                                      "12M Follow-Up"),
                                 ward_levels = 1:7) {

  ids <- mpi_resolve_ids(trial_data, pid, visitid)
  pid <- ids$pid
  visitid <- ids$visitid

  if (is.null(visits)) {
    visits <- build_mpi_visit_df(trial_data, pid, visitid, visit_forms)
  }

  visits$.patient <- as.character(visits[[pid]])

  ## Baseline date of every patient -------------------------------------------
  baseline <- visits[visits$mnpvislabel %in% baseline_label &
                       !is.na(visits$visit_date), ]
  baseline <- baseline[!duplicated(baseline$.patient), ]
  baseline_date <- stats::setNames(baseline$visit_date, baseline$.patient)

  ## Follow-up visits the outcomes are reported for ---------------------------
  target <- visits[visits$mnpvislabel %in% follow_up_labels, ]
  target$.baseline <- baseline_date[target$.patient]

  # Without a baseline or a follow-up date there is no interval to look at.
  target <- target[!is.na(target$.baseline) & !is.na(target$visit_date), ]
  rownames(target) <- NULL

  outcome <- data.frame(
    ecu_rehosp = rep(NA_real_, nrow(target)),
    stringsAsFactors = FALSE
  )
  for (level in ward_levels) {
    outcome[[paste0("ecu_rehosp_", level)]] <- rep(NA_real_, nrow(target))
  }

  ## Ward days ----------------------------------------------------------------
  wards <- trial_data[[ward_table]]
  incomplete <- character(0)

  if (nrow(target) > 0 && !is.null(wards) &&
      all(c(pid, ward_vars) %in% names(wards))) {

    ward_records <- data.frame(
      .patient = as.character(wards[[pid]]),
      .ward = mpi_as_numeric(wards[[ward_vars[["ward"]]]]),
      .start = mpi_as_date(wards[[ward_vars[["start"]]]]),
      .end = mpi_as_date(wards[[ward_vars[["end"]]]]),
      stringsAsFactors = FALSE
    )

    # Stays that cannot be placed on the timeline are skipped, but the
    # patients they belong to are flagged so the counts stay interpretable.
    broken <- is.na(ward_records$.start) | is.na(ward_records$.end)
    incomplete <- unique(ward_records$.patient[broken])
    ward_records <- ward_records[!broken, , drop = FALSE]

    # A patient without a countable stay has not been on a ward: zero days.
    outcome[, ] <- 0

    if (nrow(ward_records) > 0) {

      # Expand the stays to one row per day. The later record wins where stays
      # overlap, so a day is counted once.
      day_counts <- as.numeric(ward_records$.end - ward_records$.start) + 1
      day_counts[day_counts < 0] <- 0

      ward_days <- data.frame(
        .patient = rep(ward_records$.patient, day_counts),
        .day = do.call(c, c(
          list(as.Date(character(0))),
          Map(
            function(from, n) {
              if (n > 0) from + seq_len(n) - 1 else as.Date(character(0))
            },
            ward_records$.start, day_counts
          )
        )),
        .ward = rep(ward_records$.ward, day_counts),
        stringsAsFactors = FALSE
      )
      ward_days <- ward_days[!duplicated(
        ward_days[, c(".patient", ".day")], fromLast = TRUE
      ), ]

      matched <- merge(
        target[, c(".patient", visitid, ".baseline", "visit_date")],
        ward_days,
        by = ".patient"
      )
      matched <- matched[matched$.day >= matched$.baseline &
                           matched$.day <= matched$visit_date, , drop = FALSE]

      key_levels <- as.character(target[[visitid]])
      matched_key <- factor(as.character(matched[[visitid]]),
                            levels = key_levels)

      on_ward <- !is.na(matched$.ward) & matched$.ward != 0
      outcome$ecu_rehosp <- as.numeric(table(matched_key[on_ward]))

      for (level in ward_levels) {
        selected <- !is.na(matched$.ward) & matched$.ward == level
        outcome[[paste0("ecu_rehosp_", level)]] <-
          as.numeric(table(matched_key[selected]))
      }
    }

  } else if (nrow(target) > 0) {
    warning(paste0("build_mpi_outcome_df(): table '", ward_table,
                   "' or its variables are missing. Ward days are NA."))
  }

  # Day counts of a patient with an unusable ward date are lower bounds.
  outcome$ecu_rehosp_incomplete <- target$.patient %in% incomplete

  ## Vital status -------------------------------------------------------------
  death_date <- rep(as.Date(NA), nrow(target))
  deaths <- trial_data[[death_table]]

  if (!is.null(deaths) && all(c(pid, death_var) %in% names(deaths))) {
    death_records <- data.frame(
      .patient = as.character(deaths[[pid]]),
      .date = mpi_as_date(deaths[[death_var]]),
      stringsAsFactors = FALSE
    )
    death_records <- death_records[!is.na(death_records$.date), ]
    death_records <- death_records[!duplicated(death_records$.patient), ]
    death_date <- death_records$.date[match(target$.patient,
                                            death_records$.patient)]
  } else {
    warning(paste0("build_mpi_outcome_df(): table '", death_table,
                   "' or variable '", death_var,
                   "' is missing. Vital status is NA."))
  }

  # A patient counts as deceased at a visit only when the date of death is not
  # after that visit.
  died <- !is.na(death_date) & death_date <= target$visit_date

  result <- cbind(
    target[, c(pid, visitid, "mnpvisid", "mnpvislabel", "mnpvisno",
               "visit_date")],
    outcome,
    data.frame(
      ecu_death = died,
      ecu_death_date = death_date,
      stringsAsFactors = FALSE
    )
  )
  rownames(result) <- NULL

  return(result)
}


#' Report how the MPI variables resolve against an export
#'
#' @description
#' Shows, for every variable of the domain definition, whether it was found in
#' the export, how its records are matched to visits, and which variables were
#' picked up as its validity interval.
#'
#' @details
#' Interval bounds are searched for rather than listed, because the export does
#' not spell them consistently. Use this function to confirm that the intended
#' bounds were found before trusting the domain scores, in particular for the
#' medication and living situation domains.
#'
#' @param trial_data A secuTrial data object.
#' @param pid Column name of the patient ID in trial_data.
#' @param visitid Column name of the visit ID in trial_data.
#' @param domains Domain definition, see [mpi_domains].
#'
#' @return Data frame with one row per domain variable and the columns
#'   `domain`, `table`, `variable`, `found`, `matching`, `start`, `end` and
#'   `companions`. `matching` is one of "interval", "longitudinal",
#'   "cross-sectional", "no table" or "no variable".
#' @export
mpi_variable_report <- function(trial_data, pid = NULL, visitid = NULL,
                                domains = epicodr::mpi_domains) {

  ids <- mpi_resolve_ids(trial_data, pid, visitid)
  pid <- ids$pid
  visitid <- ids$visitid

  rows <- lapply(seq_len(nrow(domains)), function(i) {

    config <- domains[i, ]
    form <- trial_data[[config$table]]

    empty <- data.frame(
      domain = config$domain, table = config$table, variable = config$variable,
      found = FALSE, matching = "no table",
      start = NA_character_, end = NA_character_, companions = NA_character_,
      stringsAsFactors = FALSE
    )

    if (is.null(form)) {
      return(empty)
    }

    form_names <- names(form)

    if (!(config$variable %in% form_names)) {
      empty$matching <- "no variable"
      return(empty)
    }

    start_var <- suppressWarnings(
      mpi_resolve_bound(config$start, config$variable, "start", form_names)
    )
    end_var <- suppressWarnings(
      mpi_resolve_bound(config$end, config$variable, "end", form_names)
    )

    matching <- if (!is.na(start_var) || !is.na(end_var)) {
      "interval"
    } else if (visitid %in% form_names) {
      "longitudinal"
    } else {
      "cross-sectional"
    }

    companions <- mpi_split_config(config$companions)
    companions <- companions[companions %in% form_names]

    data.frame(
      domain = config$domain, table = config$table, variable = config$variable,
      found = TRUE, matching = matching,
      start = start_var, end = end_var,
      companions = if (length(companions) > 0) {
        paste(companions, collapse = ",")
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL

  return(result)
}


#' Primary coding Multidimensional Prognostic Index (MPI)
#'
#' @description
#' Adds the following tables to trial_data:
#' `ecu_mpi` with the MPI and its domain scores per visit, and
#' `ecu_mpi_outcome` with ward days and vital status per follow-up visit.
#'
#' @param trial_data A secuTrial data object.
#' @param outcomes Should the rehospitalisation and vital status table be
#'   added?
#' @param pid Column name of the patient ID in trial_data.
#' @param visitid Column name of the visit ID in trial_data.
#' @param ... Further arguments passed on to [build_mpi_df()].
#'
#' @return A secuTrial data object with the MPI tables added.
#' @export
primary_coding_suep_mpi <- function(trial_data, outcomes = TRUE, pid = NULL,
                                    visitid = NULL, ...) {

  ids <- mpi_resolve_ids(trial_data, pid, visitid)
  pid <- ids$pid
  visitid <- ids$visitid

  trial_data[["ecu_mpi"]] <- build_mpi_df(trial_data, pid, visitid, ...) %>%
    labelled::remove_var_label()

  if (outcomes) {
    trial_data[["ecu_mpi_outcome"]] <-
      build_mpi_outcome_df(trial_data, pid, visitid) %>%
      labelled::remove_var_label()
  }

  return(trial_data)
}
