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
#' A cell holding several codes separated by "|" is resolved by `multi`:
#' "max" keeps the highest code, "na" reads the cell as missing because
#' combining the codes of that variable is not defined.
#'
#' @param x Vector to coerce.
#' @param multi How to resolve a cell holding several codes, "max" or "na".
#' @return Double vector.
#' @noRd
mpi_as_numeric <- function(x, multi = c("max", "na")) {

  multi <- match.arg(multi)

  if (is.logical(x)) {
    return(as.numeric(x))
  }

  # Labelled variables keep their codes underneath the labels.
  if (inherits(x, "haven_labelled")) {
    return(as.numeric(unclass(x)))
  }

  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  if (is.factor(x)) {
    # A level is either the code itself ("1", "2") or the label that belongs
    # to it. Use the codes where they parse. Otherwise fall back to the level
    # order: the categorical scales of the export are ordered factors, so the
    # position of a level is its code. as.character() would discard that and
    # turn every labelled variable into NA.
    level_codes <- suppressWarnings(as.numeric(levels(x)))
    if (!any(grepl("|", levels(x), fixed = TRUE))) {
      if (!any(is.na(level_codes))) {
        return(level_codes[as.integer(x)])
      }
      return(as.numeric(as.integer(x)))
    }
  }

  x <- as.character(x)

  # Multi select cells such as "1|3".
  if (any(grepl("|", x, fixed = TRUE), na.rm = TRUE)) {

    if (multi == "na") {
      # Reducing several codes to one is only meaningful where the variable
      # says so. Everywhere else an ambiguous cell is missing rather than
      # quietly collapsed.
      x[grepl("|", x, fixed = TRUE)] <- NA_character_
      return(suppressWarnings(as.numeric(x)))
    }

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
                                  missing_values, multi = "max") {

  groups <- split(seq_along(keys), factor(keys, levels = key_levels))

  if (aggregate == "count") {
    # Number of distinct entries that are active at the visit. Cells may hold
    # several entries separated by "|". A visit without any active record is a
    # count of zero, not a missing value.
    values <- as.character(values)
    counts <- vapply(groups, function(idx) {
      entries <- trimws(unlist(strsplit(values[idx], "|", fixed = TRUE)))
      entries <- entries[!is.na(entries) & nzchar(entries)]
      # An unused slot of a multi select cell is written out as a token rather
      # than left off, so "MDA343|NA" is one drug and not two.
      entries <- entries[!grepl("^(NA|N/A|NaN|NULL)$", entries,
                                ignore.case = TRUE)]
      length(unique(entries))
    }, numeric(1))
    return(unname(counts))
  }

  values <- mpi_mask_missing(mpi_as_numeric(values, multi), missing_values)
  maxima <- vapply(groups, function(idx) {
    group_values <- values[idx]
    if (all(is.na(group_values))) NA_real_ else max(group_values, na.rm = TRUE)
  }, numeric(1))

  return(unname(maxima))
}


#' Carry the last observed value forward within a patient
#'
#' A visit keeps the most recent earlier measurement until a new one is
#' recorded, from which point that value applies. Visits before the first
#' measurement stay missing: a later value says nothing about an earlier
#' visit.
#'
#' Reserved for quantities that are stable or change slowly, which in the MPI
#' means height and weight. It must not be used for the assessment scales,
#' where carrying a value across visits would invent an observation.
#'
#' @param x Numeric vector, ordered by visit date within patient.
#' @param patient Patient identifier of each element.
#' @return Numeric vector of the same length.
#' @noRd
mpi_carry_forward <- function(x, patient) {

  if (length(x) != length(patient)) {
    x <- rep(x, length.out = length(patient))
  }

  for (idx in split(seq_along(x), patient)) {
    values <- x[idx]
    observed <- which(!is.na(values))
    if (length(observed) == 0 || length(observed) == length(values)) {
      next
    }
    # Number of measurements taken at or before each visit. A visit with none
    # keeps its NA, the rest take the most recent one.
    preceding <- findInterval(seq_along(values), observed)
    filled <- rep(NA_real_, length(values))
    known <- preceding > 0
    filled[known] <- values[observed[preceding[known]]]
    x[idx] <- filled
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

  # A visit key has to identify one visit of one patient. An id shared across
  # patients names the visit of the study rather than the visit of the
  # patient, and the deduplication below would quietly collapse the table to
  # one row per visit type instead of saying so.
  patients_per_key <- tapply(visits$.patient, visits$.visit_key,
                             function(x) length(unique(x)))
  shared <- names(patients_per_key)[!is.na(patients_per_key) &
                                      patients_per_key > 1]
  if (length(shared) > 0) {
    stop(paste0(
      "build_mpi_visit_df(): the visit id '", visitid, "' is shared by ",
      "several patients (", length(shared), " of ",
      length(patients_per_key), " values, e.g. '", shared[1],
      "'), so it cannot identify a visit. Pass the per-patient visit id ",
      "instead, e.g. visitid = \"mnpcvpid\"."
    ))
  }

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

  # How a cell holding several codes is read. Defaults to "na", so only a
  # variable that asks for it has its codes reduced to the highest one.
  multi <- if ("multi" %in% names(config) && !is.na(config$multi)) {
    config$multi
  } else {
    "na"
  }

  if (multi == "na" && config$aggregate != "count") {
    for (variable in available) {
      ambiguous <- sum(grepl("|", as.character(matched[[variable]]),
                             fixed = TRUE), na.rm = TRUE)
      if (ambiguous > 0) {
        warning(paste0("MPI: '", variable, "' in table '", config$table,
                       "' holds ", ambiguous, " record(s) with several codes ",
                       "in one cell. Reducing them to one value is not ",
                       "defined for this variable, so they are read as ",
                       "missing. Set multi = \"max\" in `domains` to keep ",
                       "the highest code instead."))
      }
    }
  }

  values <- lapply(available, function(variable) {
    mpi_aggregate_records(
      values = matched[[variable]],
      keys = matched$.key,
      key_levels = key_levels,
      aggregate = config$aggregate,
      missing_values = missing_values,
      multi = multi
    )
  })
  names(values) <- available

  # A variable that holds data but reads as missing throughout is almost
  # always a coding problem rather than an empty column, so say so instead of
  # letting the domain quietly lose a variable.
  for (variable in available) {
    if (length(values[[variable]]) > 0 &&
        all(is.na(values[[variable]])) &&
        any(!is.na(matched[[variable]]))) {
      warning(paste0("MPI: '", variable, "' in table '", config$table,
                     "' holds data but could not be read as a numeric code. ",
                     "It is treated as missing."))
    }
  }

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

  # Read every variable of the domain as a numeric code, the way the original
  # script applied as.numeric() across the whole domain before summing. Doing
  # it here rather than in each collector means a variable cannot reach the
  # scoring rules still carrying labels.
  collected <- lapply(collected, mpi_as_numeric)

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
#' Each variable contributes `max(0, value - floor)`, so a graded scale whose
#' "no problem" code is 1 adds nothing when it is answered at its best.
#' Without this the total would depend on how many items a patient answered
#' rather than on how the patient is. Clamping at zero stops a value below its
#' floor from subtracting from the domain burden.
#'
#' @param values Data frame of the domain variables, one row per visit.
#' @param cuts Named numeric vector of cut-offs, names are the domain scores.
#' @param floors Named numeric vector of "no problem" codes, names are
#'   variables. Variables it does not name keep a floor of zero.
#' @return Double vector of domain scores.
#' @noRd
mpi_score_sum <- function(values, cuts, floors = NULL) {

  # Guard the sum against a column that still carries labels: as.matrix() on a
  # frame with one non-numeric column turns every column into character.
  values <- as.matrix(as.data.frame(lapply(values, mpi_as_numeric)))

  if (!is.null(floors) && nrow(values) > 0) {
    shift <- suppressWarnings(as.numeric(floors[colnames(values)]))
    shift[is.na(shift)] <- 0
    # pmax() leaves NA as NA, so an unobserved variable stays unobserved.
    values <- pmax(values - rep(shift, each = nrow(values)), 0)
  }

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

#' Age at every visit
#'
#' @param trial_data A secuTrial data object.
#' @param visits Visit table carrying the helper columns `.patient` and `.date`.
#' @param pid Column name of the patient ID in trial_data.
#' @param birthdate_table Name of the table holding the date of birth.
#' @param birthdate_var Name of the date of birth variable. The `.date`
#'   companion created by `dates_tsExport()` is used where it exists.
#' @return Double vector of decimal years, one per visit.
#' @noRd
mpi_age_at_visit <- function(trial_data, visits, pid, birthdate_table,
                             birthdate_var) {

  age <- rep(NA_real_, nrow(visits))
  births <- trial_data[[birthdate_table]]
  column <- intersect(c(paste0(birthdate_var, ".date"), birthdate_var),
                      names(births))

  if (is.null(births) || !(pid %in% names(births)) || length(column) == 0) {
    warning(paste0("build_mpi_df(): table '", birthdate_table,
                   "' or variable '", birthdate_var, "' is missing. ",
                   "Age is NA."))
    return(age)
  }

  raw <- births[[column[1]]]
  birth_records <- data.frame(
    .patient = as.character(births[[pid]]),
    .date = mpi_as_date(raw),
    stringsAsFactors = FALSE
  )

  # Name what could not be read rather than letting it vanish as NA.
  unread <- !is.na(raw) & trimws(as.character(raw)) != "" &
    is.na(birth_records$.date)
  if (any(unread)) {
    warning(paste0("build_mpi_df(): ", sum(unread), " value(s) of '",
                   column[1], "' are not a complete date, e.g. '",
                   as.character(raw[unread][1]), "'. Age is NA for them."))
  }

  birth_records <- birth_records[!is.na(birth_records$.date), ]
  birth_records <- birth_records[!duplicated(birth_records$.patient), ]
  birth_date <- birth_records$.date[match(visits$.patient,
                                          birth_records$.patient)]

  # Whole years plus the elapsed share of the current one, so leap years count.
  return(lubridate::time_length(lubridate::interval(birth_date, visits$.date),
                                "years"))
}


#' A value recorded once per patient, at every visit
#'
#' @param trial_data A secuTrial data object.
#' @param visits Visit table carrying the helper column `.patient`.
#' @param pid Column name of the patient ID in trial_data.
#' @param table Name of the table holding the variable.
#' @param variable Name of the variable. The `.factor` companion created by
#'   `factorize_tsExport()` is used where it exists.
#' @param label What the variable is, for the warning when it is missing.
#' @return Vector with one value per visit, a factor where the companion
#'   exists.
#' @noRd
mpi_patient_value <- function(trial_data, visits, pid, table, variable,
                              label) {

  records <- trial_data[[table]]
  column <- intersect(c(paste0(variable, ".factor"), variable), names(records))

  if (is.null(records) || !(pid %in% names(records)) || length(column) == 0) {
    warning(paste0("build_mpi_df(): table '", table, "' or variable '",
                   variable, "' is missing. ", label, " is NA."))
    return(rep(NA, nrow(visits)))
  }

  value <- records[[column[1]]]
  patient <- as.character(records[[pid]])

  # The value is recorded once per patient, so the first documented one stands
  # for every visit.
  documented <- !is.na(value) & trimws(as.character(value)) != ""

  # A column without a single documented value usually means the variable is
  # recorded in another table, so say so rather than return NA silently.
  if (!any(documented)) {
    warning(paste0("build_mpi_df(): variable '", column[1], "' in table '",
                   table, "' holds no documented value. ", label, " is NA."))
  }

  value <- value[documented]
  patient <- patient[documented]
  first <- !duplicated(patient)

  return(value[first][match(visits$.patient, patient[first])])
}


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
#' `ecu_age` is the age in decimal years on the visit date, derived from the
#' date of birth for every visit, so it advances between visits rather than on
#' birthdays only. Use `floor(ecu_age)` for completed years. The age at
#' baseline is the value of the baseline row.
#'
#' `ecu_sex`, `ecu_ethnicity`, `ecu_origin` and `ecu_palliative` are recorded
#' once per patient and carried on every visit. They are the factors
#' `gec_gender.factor`, `gec_ethnicity.factor`, `origin.factor` and
#' `palliative.factor` where the export was factorized and the raw codes
#' otherwise.
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
#' @param birthdate_table Name of the table holding the date of birth.
#' @param birthdate_var Name of the date of birth variable in
#'   `birthdate_table`. Its `.date` companion is used where `dates_tsExport()`
#'   created one.
#' @param sex_table Name of the table holding the sex.
#' @param sex_var Name of the sex variable in `sex_table`. Its `.factor`
#'   companion is used where `factorize_tsExport()` created one.
#' @param ethnicity_table Name of the table holding the ethnicity.
#' @param ethnicity_var Name of the ethnicity variable in `ethnicity_table`.
#'   Its `.factor` companion is used where `factorize_tsExport()` created one.
#' @param origin_table Name of the table holding the origin.
#' @param origin_var Name of the origin variable in `origin_table`. Its
#'   `.factor` companion is used where `factorize_tsExport()` created one.
#' @param palliative_table Name of the table holding the palliative status.
#' @param palliative_var Name of the palliative status variable in
#'   `palliative_table`. Its `.factor` companion is used where
#'   `factorize_tsExport()` created one.
#'
#' @return Data frame with one row per visit, the visit identification columns,
#'   `ecu_age`, `ecu_sex`, `ecu_ethnicity`, `ecu_origin`, `ecu_palliative`, one
#'   `ecu_mpi_<domain>` column per domain, `ecu_mpi_n_domains`, `ecu_mpi` and
#'   `ecu_mpi_cat`.
#' @export
build_mpi_df <- function(trial_data, pid = NULL, visitid = NULL,
                         domains = epicodr::mpi_domains,
                         scoring = epicodr::mpi_scoring,
                         visit_forms = mpi_default_visit_forms(),
                         missing_values = epicodr::missing_values_default,
                         min_domains = 5,
                         birthdate_table = "bv1",
                         birthdate_var = "gec_birthdate",
                         sex_table = "bv1",
                         sex_var = "gec_gender",
                         ethnicity_table = "bv1",
                         ethnicity_var = "gec_ethnicity",
                         origin_table = "bv1",
                         origin_var = "origin",
                         palliative_table = "bv1",
                         palliative_var = "palliative") {

  ids <- mpi_resolve_ids(trial_data, pid, visitid)
  pid <- ids$pid
  visitid <- ids$visitid

  mpi_validate_config(domains, scoring)

  visits <- build_mpi_visit_df(trial_data, pid, visitid, visit_forms)

  # Helper columns keep the collectors independent of the actual ID names.
  visits$.patient <- as.character(visits[[pid]])
  visits$.key <- as.character(visits[[visitid]])
  visits$.date <- visits$visit_date

  age <- mpi_age_at_visit(trial_data, visits, pid, birthdate_table,
                          birthdate_var)
  sex <- mpi_patient_value(trial_data, visits, pid, sex_table, sex_var,
                           "Sex")
  ethnicity <- mpi_patient_value(trial_data, visits, pid, ethnicity_table,
                                 ethnicity_var, "Ethnicity")
  origin <- mpi_patient_value(trial_data, visits, pid, origin_table,
                              origin_var, "Origin")
  palliative <- mpi_patient_value(trial_data, visits, pid, palliative_table,
                                  palliative_var, "Palliative status")

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
        stats::setNames(suppressWarnings(as.numeric(cuts)), names(cuts)),
        floors = if ("floor" %in% names(domain_config)) {
          stats::setNames(suppressWarnings(as.numeric(domain_config$floor)),
                          domain_config$variable)
        } else {
          NULL
        }
      ),

      bmi = {
        # Height and weight are recorded on separate forms, so they rarely
        # share a visit. Carry each of them forward, which lets the BMI be
        # derived from the most recent measurement of either. Height does not
        # change, and a weight applies until the next one is recorded.
        #
        # Forward only: a visit before the first measurement stays missing,
        # because a later weight says nothing about an earlier visit. This is
        # confined to height and weight, the only quantities in the MPI stable
        # enough to justify it.
        height <- mpi_carry_forward(
          if ("gec_height" %in% names(values)) values$gec_height else NA_real_,
          visits$.patient
        )
        weight <- mpi_carry_forward(
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
    data.frame(ecu_age = age, ecu_sex = sex, ecu_ethnicity = ethnicity,
               ecu_origin = origin, ecu_palliative = palliative),
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
#' Derives days spent on a ward between the baseline visit and each follow-up
#' visit. The vital status is derived per patient by [build_mpi_vital_df()],
#' the ward days per patient and ward group by [build_mpi_rehosp_df()].
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
#' @param baseline_label Visit label of the baseline visit.
#' @param follow_up_labels Visit labels the outcomes are derived for.
#' @param ward_levels Ward codes that get their own day count.
#'
#' @return Data frame with one row per follow-up visit, the visit
#'   identification columns, `ecu_rehosp`, one `ecu_rehosp_<level>` column per
#'   ward code and `ecu_rehosp_incomplete`.
#' @export
build_mpi_outcome_df <- function(trial_data, pid = NULL, visitid = NULL,
                                 visits = NULL,
                                 visit_forms = mpi_default_visit_forms(),
                                 ward_table = "eward",
                                 ward_vars = c(ward = "gec_ward",
                                               start = "ward_start",
                                               end = "ward_end"),
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

  result <- cbind(
    target[, c(pid, visitid, "mnpvisid", "mnpvislabel", "mnpvisno",
               "visit_date")],
    outcome
  )
  rownames(result) <- NULL

  return(result)
}


#' Build the patient level rehospitalisation table
#'
#' @description
#' Sums the ward days of [build_mpi_outcome_df()] over groups of wards, once
#' per follow-up visit, and reports them with one row per patient.
#'
#' @details
#' Every count covers the days from the baseline visit up to the follow-up
#' visit, so the 12 month count includes the days of the 3 month count. The
#' columns are named `ecu_rehosp_<group>_<follow-up>`, with the defaults
#' `ecu_rehosp_1_3_3M`, `ecu_rehosp_4_3M`, `ecu_rehosp_1_3_12M` and
#' `ecu_rehosp_4_12M`.
#'
#' A count is NA when the patient has no dated baseline visit or no dated
#' follow-up visit of that kind, or when the ward table is missing from the
#' export. Counts of a patient flagged in `ecu_rehosp_incomplete` of the
#' outcome table are lower bounds.
#'
#' @param trial_data A secuTrial data object.
#' @param pid Column name of the patient ID in trial_data.
#' @param visitid Column name of the visit ID in trial_data.
#' @param visits Visit table as returned by [build_mpi_visit_df()]. Built from
#'   `trial_data` when NULL.
#' @param visit_forms Named list describing the visit forms, see
#'   [build_mpi_visit_df()].
#' @param outcome Outcome table as returned by [build_mpi_outcome_df()]. Built
#'   from `trial_data` when NULL.
#' @param ward_groups Named list of the ward codes whose days are summed. The
#'   names become part of the column names.
#' @param follow_ups Named character vector of follow-up visit labels. The
#'   names become the suffix of the column names.
#'
#' @return Data frame with one row per patient of the visit table, the patient
#'   ID and one count per follow-up visit and ward group.
#' @export
build_mpi_rehosp_df <- function(trial_data, pid = NULL, visitid = NULL,
                                visits = NULL,
                                visit_forms = mpi_default_visit_forms(),
                                outcome = NULL,
                                ward_groups = list("1_3" = 1:3, "4" = 4),
                                follow_ups = c("3M" = "3M Follow-Up",
                                               "12M" = "12M Follow-Up")) {

  ids <- mpi_resolve_ids(trial_data, pid, visitid)
  pid <- ids$pid
  visitid <- ids$visitid

  if (is.null(visits)) {
    visits <- build_mpi_visit_df(trial_data, pid, visitid, visit_forms)
  }

  if (is.null(outcome)) {
    outcome <- build_mpi_outcome_df(
      trial_data, pid, visitid, visits = visits,
      follow_up_labels = unname(follow_ups),
      ward_levels = sort(unique(unlist(ward_groups)))
    )
  }

  ward_columns <- lapply(ward_groups, function(levels) {
    paste0("ecu_rehosp_", levels)
  })
  absent <- setdiff(unlist(ward_columns), names(outcome))
  if (length(absent) > 0) {
    stop(paste0("build_mpi_rehosp_df(): the outcome table has no column ",
                paste(absent, collapse = ", "), "."))
  }

  patients <- unique(as.character(visits[[pid]]))
  result <- data.frame(.patient = patients, stringsAsFactors = FALSE)

  for (follow_up in names(follow_ups)) {

    # Patients without this follow-up visit get a row of NA.
    rows <- outcome[outcome$mnpvislabel %in% follow_ups[[follow_up]], ,
                    drop = FALSE]
    rows <- rows[match(patients, as.character(rows[[pid]])), , drop = FALSE]

    for (group in names(ward_groups)) {
      result[[paste0("ecu_rehosp_", group, "_", follow_up)]] <-
        unname(rowSums(rows[, ward_columns[[group]], drop = FALSE]))
    }
  }

  names(result)[names(result) == ".patient"] <- pid

  return(result)
}


#' Build the MPI vital status table
#'
#' @description
#' Derives the vital status of every patient from the date of death, once at
#' the end of the study and once at each horizon after baseline, by default 3
#' and 12 months.
#'
#' @details
#' Vital status belongs to the patient rather than to a visit, so it does not
#' depend on a follow-up visit having been documented. Patients who die before
#' their first follow-up rarely have one, and they are the events the MPI is
#' meant to predict.
#'
#' `ecu_alive` is FALSE for every patient with a date of death.
#' `ecu_alive_<months>M`, by default `ecu_alive_3M` and `ecu_alive_12M`,
#' censors at that many calendar months after the baseline visit: a death on or
#' before that day makes it FALSE, a later death leaves it TRUE. It is NA when
#' the baseline visit has no date.
#'
#' A death documented without a usable date cannot be told apart from survival
#' and counts as alive.
#'
#' @param trial_data A secuTrial data object.
#' @param pid Column name of the patient ID in trial_data.
#' @param visitid Column name of the visit ID in trial_data.
#' @param visits Visit table as returned by [build_mpi_visit_df()]. Built from
#'   `trial_data` when NULL.
#' @param visit_forms Named list describing the visit forms, see
#'   [build_mpi_visit_df()].
#' @param death_table Name of the table holding the date of death.
#' @param death_var Name of the date of death variable in `death_table`.
#' @param baseline_label Visit label of the baseline visit.
#' @param horizons Months after baseline at which the vital status is censored.
#'   Each gets an `ecu_alive_<months>M` column.
#'
#' @return Data frame with one row per patient of the visit table, the patient
#'   ID, `ecu_death_date`, `ecu_alive` and one `ecu_alive_<months>M` column per
#'   horizon.
#' @export
build_mpi_vital_df <- function(trial_data, pid = NULL, visitid = NULL,
                               visits = NULL,
                               visit_forms = mpi_default_visit_forms(),
                               death_table = "fv2_6",
                               death_var = "death_date",
                               baseline_label = "Baseline",
                               horizons = c(3, 12)) {

  ids <- mpi_resolve_ids(trial_data, pid, visitid)
  pid <- ids$pid
  visitid <- ids$visitid

  if (is.null(visits)) {
    visits <- build_mpi_visit_df(trial_data, pid, visitid, visit_forms)
  }

  visits$.patient <- as.character(visits[[pid]])
  patients <- unique(visits$.patient)

  ## Baseline date of every patient -------------------------------------------
  baseline <- visits[visits$mnpvislabel %in% baseline_label &
                       !is.na(visits$visit_date), ]
  baseline <- baseline[!duplicated(baseline$.patient), ]
  baseline_date <- baseline$visit_date[match(patients, baseline$.patient)]

  ## Date of death ------------------------------------------------------------
  death_date <- rep(as.Date(NA), length(patients))
  deaths <- trial_data[[death_table]]
  known <- !is.null(deaths) && all(c(pid, death_var) %in% names(deaths))

  if (known) {
    death_records <- data.frame(
      .patient = as.character(deaths[[pid]]),
      .date = mpi_as_date(deaths[[death_var]]),
      stringsAsFactors = FALSE
    )
    death_records <- death_records[!is.na(death_records$.date), ]
    death_records <- death_records[!duplicated(death_records$.patient), ]
    death_date <- death_records$.date[match(patients, death_records$.patient)]
  } else {
    warning(paste0("build_mpi_vital_df(): table '", death_table,
                   "' or variable '", death_var,
                   "' is missing. Vital status is NA."))
  }

  ## Vital status -------------------------------------------------------------
  alive <- is.na(death_date)
  if (!known) {
    alive[] <- NA
  }

  result <- data.frame(
    .patient = patients,
    ecu_death_date = death_date,
    ecu_alive = alive,
    stringsAsFactors = FALSE
  )

  for (months in horizons) {
    # Calendar months, rolled back to the last day of a shorter month.
    horizon <- lubridate::add_with_rollback(baseline_date,
                                            lubridate::period(months, "months"))

    alive_at <- is.na(death_date) | death_date > horizon
    alive_at[is.na(horizon) | !known] <- NA

    result[[paste0("ecu_alive_", months, "M")]] <- alive_at
  }

  names(result)[names(result) == ".patient"] <- pid

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
#' `ecu_mpi_outcome` with ward days per follow-up visit.
#'
#' With `outcomes`, every row of `ecu_mpi` also carries the vital status of its
#' patient from [build_mpi_vital_df()] and the ward days of its patient per
#' ward group and follow-up visit from [build_mpi_rehosp_df()], so the baseline
#' MPI and the outcomes it predicts sit on one row.
#'
#' @param trial_data A secuTrial data object.
#' @param outcomes Should the vital status and the rehospitalisation table be
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

  mpi <- build_mpi_df(trial_data, pid, visitid, ...)

  if (outcomes) {
    outcome <- build_mpi_outcome_df(trial_data, pid, visitid)

    patient_tables <- list(
      build_mpi_vital_df(trial_data, pid, visitid),
      build_mpi_rehosp_df(trial_data, pid, visitid, outcome = outcome)
    )
    for (patient_table in patient_tables) {
      columns <- setdiff(names(patient_table), pid)
      mpi[columns] <- patient_table[match(as.character(mpi[[pid]]),
                                          as.character(patient_table[[pid]])),
                                    columns]
    }
  }

  trial_data[["ecu_mpi"]] <- mpi %>%
    labelled::remove_var_label()

  if (outcomes) {
    trial_data[["ecu_mpi_outcome"]] <- outcome %>%
      labelled::remove_var_label()
  }

  return(trial_data)
}
