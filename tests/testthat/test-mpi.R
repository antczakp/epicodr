# Tests for the Multidimensional Prognostic Index.
#
# The fixture is a minimal synthetic export: two patients, three visits each,
# and just enough tables to exercise cross sectional, longitudinal and interval
# data as well as every scoring rule.

make_mpi_fixture <- function() {

  visits <- data.frame(
    export_psn = c("P1", "P1", "P1", "P2", "P2", "P2"),
    mnpcvpid = c("11", "12", "13", "21", "22", "23"),
    mnpvisid = c("1", "2", "3", "4", "5", "6"),
    mnpvislabel = rep(c("Baseline", "3M Follow-Up", "12M Follow-Up"), 2),
    mnpvisno = rep(1:3, 2),
    pr_docudate = c("20210101", "20210401", "20220101",
                    "20210201", "20210501", "20220201"),
    pr_visit_mode = 1,
    stringsAsFactors = FALSE
  )

  trial_data <- list(

    export_options = list(
      id_names = data.frame(pid = "export_psn", visitid = "mnpcvpid",
                            stringsAsFactors = FALSE)
    ),

    # Visit form.
    m = visits,

    ## Comorbidity: cross sectional, one row per patient ----------------------
    # P1 sums to 2 -> score 0.5 (cut-offs 0 / 1 / 3)
    # P2 sums to 0 -> score 0
    fv2_1 = data.frame(export_psn = c("P1", "P2"), gec_cv = c(1, 0),
                       stringsAsFactors = FALSE),
    fv2_2 = data.frame(export_psn = c("P1", "P2"), gec_pd = c(1, 0),
                       stringsAsFactors = FALSE),
    # 99 is a generic missing code and must not inflate the sum.
    fv2_6 = data.frame(export_psn = c("P1", "P2"), diag_oth = c(99, 0),
                       death_date = c(NA, "20210601"),
                       stringsAsFactors = FALSE),

    ## Nutrition: height at baseline only, weight per visit -------------------
    bv1 = data.frame(export_psn = c("P1", "P2"), mnpcvpid = c("11", "21"),
                     gec_height = c(180, 170), stringsAsFactors = FALSE),
    fuv1 = data.frame(export_psn = c("P1", "P2"), mnpcvpid = c("11", "21"),
                      gec_weight = c(85, 40), stringsAsFactors = FALSE),

    ## Social: interval data with companion variables -------------------------
    # P1 lives at home with another adult -> 0. P2 lives alone -> 0.5.
    ehome = data.frame(
      export_psn = c("P1", "P2"),
      home = c(2, 2),
      home_adults = c(1, 0),
      home_childr = c(0, 0),
      home_start = c("20200101", "20200101"),
      home_end = c(NA, NA),
      stringsAsFactors = FALSE
    ),

    ## Medications: interval data, counted ------------------------------------
    # P1 takes two drugs at baseline, one of which stops before the 3M visit.
    emed_oth = data.frame(
      export_psn = c("P1", "P1", "P2"),
      gec_med_drug = c("A", "B", "C"),
      gec_med_start = c("20200101", "20200101", "20200101"),
      gec_med_end = c(NA, "20210301", NA),
      stringsAsFactors = FALSE
    ),

    ## Patient reported outcomes: longitudinal, one row per visit ------------
    # Feeds the cognition, adl and mental domains, so the fixture reaches the
    # five scored domains an MPI needs.
    prom = data.frame(
      export_psn = visits$export_psn,
      mnpcvpid = visits$mnpcvpid,
      pro_pc2r = rep(c(1, 3), each = 3),
      pro_pc35r = rep(c(1, 3), each = 3),
      pro_pc36r = rep(c(1, 3), each = 3),
      eq5d_mob = rep(c(1, 5), each = 3),
      eq5d_anx = rep(c(1, 5), each = 3),
      stringsAsFactors = FALSE
    ),

    ## Ward stays -------------------------------------------------------------
    # P1 spends 5 days on ward 1 between baseline and the 3M visit.
    eward = data.frame(
      export_psn = "P1",
      gec_ward = 1,
      ward_start = "20210201",
      ward_end = "20210205",
      stringsAsFactors = FALSE
    )
  )

  class(trial_data) <- "tsExportdata"
  return(trial_data)
}


# Visit table ------------------------------------------------------------------

test_that("build_mpi_visit_df() returns one row per visit with parsed dates", {

  visits <- build_mpi_visit_df(make_mpi_fixture())

  expect_equal(nrow(visits), 6)
  expect_true(inherits(visits$visit_date, "Date"))
  expect_equal(
    visits$visit_date[visits$mnpcvpid == "11"],
    as.Date("2021-01-01")
  )
  expect_setequal(visits$mnpvislabel,
                  c("Baseline", "3M Follow-Up", "12M Follow-Up"))
})


# Domain scoring ---------------------------------------------------------------

test_that("sum cut-offs are inclusive lower bounds", {

  cuts <- c("0" = 0, "0.5" = 1, "1" = 3)
  values <- data.frame(a = c(0, 1, 2, 3, 4), b = 0)

  expect_equal(mpi_score_sum(values, cuts), c(0, 0.5, 0.5, 1, 1))
})


test_that("the shipped cut-offs map sums to scores as specified", {

  # A cut-off is the inclusive lower bound of its score: below the 0.5
  # cut-off scores 0, reaching it scores 0.5, reaching the 1 cut-off scores 1.
  # Locked against the shipped table so the rule cannot drift.
  score_for <- function(domain, sums) {
    row <- epicodr::mpi_scoring[epicodr::mpi_scoring$domain == domain, ]
    cuts <- stats::setNames(
      as.numeric(c(row$cut_0, row$cut_0.5, row$cut_1)),
      c("0", "0.5", "1")
    )
    mpi_score_sum(data.frame(a = sums), cuts)
  }

  # cut-offs 0 / 1 / 3
  expect_equal(score_for("comorbidity", 0:4), c(0, 0.5, 0.5, 1, 1))
  # cut-offs 0 / 4 / 7
  expect_equal(score_for("cognition", c(0, 3, 4, 6, 7, 9)),
               c(0, 0, 0.5, 0.5, 1, 1))
  expect_equal(score_for("adl", c(0, 3, 4, 6, 7, 9)),
               c(0, 0, 0.5, 0.5, 1, 1))
  # cut-offs 0 / 12 / 15
  expect_equal(score_for("mental", c(0, 11, 12, 14, 15, 20)),
               c(0, 0, 0.5, 0.5, 1, 1))
  # cut-offs 0 / 2 / 4
  expect_equal(score_for("medications", c(0, 1, 2, 3, 4, 9)),
               c(0, 0, 0.5, 0.5, 1, 1))
})


test_that("a domain without a single observed variable scores NA, not 0", {

  cuts <- c("0" = 0, "0.5" = 1, "1" = 3)
  values <- data.frame(a = c(NA_real_, 0), b = c(NA_real_, NA_real_))

  expect_equal(mpi_score_sum(values, cuts), c(NA, 0))
})


test_that("BMI bands are read literally from the cut-off strings", {

  cuts <- c("0" = "23-33", "0.5" = "<23,>33", "1" = "<20,>38")
  bmi <- c(19, 20, 22, 23, 30, 33, 35, 38, 40, NA)

  expect_equal(
    mpi_score_bmi(bmi, cuts),
    c(1, 0.5, 0.5, 0, 0, 0, 0.5, 0.5, 1, NA)
  )
})


test_that("living situation rules cover the documented combinations", {

  # home only
  expect_equal(mpi_score_living_situation(2, 1, 0, NA), 0)
  expect_equal(mpi_score_living_situation(2, 0, 0, NA), 0.5)
  expect_equal(mpi_score_living_situation(2, NA, NA, NA), 0.5)
  expect_equal(mpi_score_living_situation(1, NA, NA, NA), 1)

  # resid only
  expect_equal(mpi_score_living_situation(NA, NA, NA, 1), 1)
  expect_equal(mpi_score_living_situation(NA, NA, NA, 2), 0.5)
  expect_equal(mpi_score_living_situation(NA, NA, NA, 4), 1)

  # both known
  expect_equal(mpi_score_living_situation(2, NA, NA, 1), 0)
  expect_equal(mpi_score_living_situation(1, NA, NA, 1), 1)
  expect_equal(mpi_score_living_situation(2, NA, NA, 4), 1)

  # nothing known
  expect_equal(mpi_score_living_situation(NA, NA, NA, NA), NA_real_)
})


test_that("an unknown number of children does not abort the rules", {

  # home_adults is known, home_childr is not: the original script raised
  # "missing value where TRUE/FALSE needed" here.
  expect_equal(mpi_score_living_situation(2, 0, NA, NA), 0.5)
  expect_equal(mpi_score_living_situation(2, 2, NA, NA), 0)
})


test_that("categorize_mpi_ecu() assigns the three MPI risk groups", {

  expect_equal(
    as.character(categorize_mpi_ecu(c(0, 0.33, 0.5, 0.66, 0.9, NA))),
    c("MPI-1", "MPI-1", "MPI-2", "MPI-2", "MPI-3", NA)
  )
})


# Helpers ----------------------------------------------------------------------

test_that("dates are accepted as %Y%m%d, ISO and Date", {

  expected <- as.Date("2021-03-04")

  expect_equal(mpi_as_date("20210304"), expected)
  expect_equal(mpi_as_date(20210304), expected)
  expect_equal(mpi_as_date("2021-03-04"), expected)
  expect_equal(mpi_as_date(expected), expected)
  expect_true(is.na(mpi_as_date(NA)))
})


test_that("multi select cells collapse to their highest code", {

  expect_equal(mpi_as_numeric(c("1|3", "2", NA)), c(3, 2, NA))
})


test_that("generic missing codes are masked", {

  expect_equal(mpi_mask_missing(c(0, 1, -1, 99), c(-1, 88, 98, 99)),
               c(0, 1, NA, NA))
  expect_equal(mpi_mask_missing(c(0, 99), numeric(0)), c(0, 99))
})


test_that("values are carried forward and backward within a patient", {

  x <- c(NA, 5, NA, NA, 7, NA)
  patient <- c("a", "a", "a", "b", "b", "b")

  expect_equal(mpi_fill_within_patient(x, patient), c(5, 5, 5, 7, 7, 7))
})


# End to end -------------------------------------------------------------------

test_that("build_mpi_df() scores every visit of every patient", {

  result <- suppressWarnings(build_mpi_df(make_mpi_fixture()))

  expect_equal(nrow(result), 6)
  expect_true(all(c("ecu_mpi", "ecu_mpi_cat", "ecu_mpi_n_domains") %in%
                    names(result)))
  expect_true(all(paste0("ecu_mpi_", epicodr::mpi_scoring$domain) %in%
                    names(result)))

  # MPI is a mean of 0 / 0.5 / 1 domain scores.
  expect_true(all(result$ecu_mpi >= 0 & result$ecu_mpi <= 1, na.rm = TRUE))
})


test_that("cross sectional comorbidity applies to every visit and masks 99", {

  result <- suppressWarnings(build_mpi_df(make_mpi_fixture()))

  # P1: gec_cv = 1 plus gec_pd = 1 is 2, diag_oth = 99 is masked -> 0.5
  expect_equal(unique(result$ecu_mpi_comorbidity[result$export_psn == "P1"]),
               0.5)
  # P2: everything 0 -> score 0 rather than NA
  expect_equal(unique(result$ecu_mpi_comorbidity[result$export_psn == "P2"]), 0)
})


test_that("nutrition carries a single measurement to every visit", {

  result <- suppressWarnings(build_mpi_df(make_mpi_fixture()))

  # P1: 85 kg at 1.80 m is a BMI of 26.2 -> inside 23-33 -> 0
  expect_equal(unique(result$ecu_mpi_nutrition[result$export_psn == "P1"]), 0)
  # P2: 40 kg at 1.70 m is a BMI of 13.8 -> below 20 -> 1
  expect_equal(unique(result$ecu_mpi_nutrition[result$export_psn == "P2"]), 1)
})


test_that("interval data expire between visits", {

  result <- suppressWarnings(build_mpi_df(make_mpi_fixture()))
  p1 <- result[result$export_psn == "P1", ]

  # Two drugs at baseline, one of which ends 2021-03-01, so one remains at the
  # later visits. Cut-offs are 0 / 2 / 4, so 2 drugs score 0.5 and 1 scores 0.
  expect_equal(p1$ecu_mpi_medications[p1$mnpvislabel == "Baseline"], 0.5)
  expect_equal(p1$ecu_mpi_medications[p1$mnpvislabel == "3M Follow-Up"], 0)
})


test_that("social companions decide the living situation score", {

  result <- suppressWarnings(build_mpi_df(make_mpi_fixture()))

  # P1 lives with another adult, P2 lives alone.
  expect_equal(unique(result$ecu_mpi_social[result$export_psn == "P1"]), 0)
  expect_equal(unique(result$ecu_mpi_social[result$export_psn == "P2"]), 0.5)
})


test_that("the unscored biomarker domain stays NA and is left out of the mean", {

  result <- suppressWarnings(build_mpi_df(make_mpi_fixture()))

  expect_true(all(is.na(result$ecu_mpi_biomarker)))
  expect_true(all(result$ecu_mpi_n_domains <= nrow(epicodr::mpi_scoring) - 1))
})


test_that("min_domains suppresses an MPI that rests on too little data", {

  result <- suppressWarnings(build_mpi_df(make_mpi_fixture(), min_domains = 8))

  expect_true(all(is.na(result$ecu_mpi)))
  expect_true(all(is.na(result$ecu_mpi_cat)))
})


test_that("an MPI needs five scored domains by default", {

  expect_equal(formals(build_mpi_df)$min_domains, 5)

  # The fixture scores seven domains, biomarker being the unscorable one.
  full <- suppressWarnings(build_mpi_df(make_mpi_fixture()))
  expect_true(all(full$ecu_mpi_n_domains == 7))
  expect_false(any(is.na(full$ecu_mpi)))

  # Drop the patient reported outcomes and only four domains remain.
  sparse_data <- make_mpi_fixture()
  sparse_data$prom <- NULL
  sparse <- suppressWarnings(build_mpi_df(sparse_data))

  expect_true(all(sparse$ecu_mpi_n_domains == 4))
  expect_true(all(is.na(sparse$ecu_mpi)))
  expect_true(all(is.na(sparse$ecu_mpi_cat)))

  # The domain scores themselves are still reported, so it stays visible why
  # the MPI was withheld.
  expect_false(any(is.na(sparse$ecu_mpi_comorbidity)))
  expect_false(any(is.na(sparse$ecu_mpi_nutrition)))

  # Lowering the threshold brings the same visits back.
  relaxed <- suppressWarnings(build_mpi_df(sparse_data, min_domains = 4))
  expect_false(any(is.na(relaxed$ecu_mpi)))
})


test_that("build_mpi_outcome_df() counts ward days inside the visit window", {

  result <- suppressWarnings(build_mpi_outcome_df(make_mpi_fixture()))

  p1_3m <- result[result$export_psn == "P1" &
                    result$mnpvislabel == "3M Follow-Up", ]

  # 2021-02-01 to 2021-02-05 inclusive is 5 days on ward 1.
  expect_equal(p1_3m$ecu_rehosp, 5)
  expect_equal(p1_3m$ecu_rehosp_1, 5)
  expect_equal(p1_3m$ecu_rehosp_2, 0)

  # P2 never went to a ward.
  expect_true(all(result$ecu_rehosp[result$export_psn == "P2"] == 0))
})


test_that("death is only reported from the visit it happened before", {

  result <- suppressWarnings(build_mpi_outcome_df(make_mpi_fixture()))
  p2 <- result[result$export_psn == "P2", ]

  # P2 died 2021-06-01: after the 3M visit (2021-05-01), before the 12M visit.
  expect_false(p2$ecu_death[p2$mnpvislabel == "3M Follow-Up"])
  expect_true(p2$ecu_death[p2$mnpvislabel == "12M Follow-Up"])
})


test_that("no rehospitalisation is zero, not NA", {

  result <- suppressWarnings(build_mpi_outcome_df(make_mpi_fixture()))

  # P2 never went to a ward.
  expect_equal(unique(result$ecu_rehosp[result$export_psn == "P2"]), 0)
  expect_false(any(is.na(result$ecu_rehosp)))
  expect_false(any(result$ecu_rehosp_incomplete))
})


test_that("an incomplete ward stay is skipped but the patient is flagged", {

  trial_data <- make_mpi_fixture()
  # A second stay whose end date was never documented.
  trial_data$eward <- rbind(
    trial_data$eward,
    data.frame(export_psn = "P1", gec_ward = 4,
               ward_start = "20210310", ward_end = NA,
               stringsAsFactors = FALSE)
  )

  result <- suppressWarnings(build_mpi_outcome_df(trial_data))
  p1_3m <- result[result$export_psn == "P1" &
                    result$mnpvislabel == "3M Follow-Up", ]

  # The countable stay still contributes its 5 days.
  expect_equal(p1_3m$ecu_rehosp, 5)
  expect_true(p1_3m$ecu_rehosp_incomplete)

  # P2 is unaffected.
  expect_false(any(result$ecu_rehosp_incomplete[result$export_psn == "P2"]))

  # The vital status does not depend on the ward timeline.
  expect_false(any(is.na(result$ecu_death)))
})


test_that("interval bounds are found by prefix when they are spelled apart", {

  # The export pairs a "gec_..._start" with a bare "..._end", which is why the
  # original script searched for bounds instead of listing them.
  expect_equal(mpi_resolve_bound(NA, "home", "start",
                                 c("home", "home_startdate")),
               "home_startdate")
  expect_equal(mpi_resolve_bound(NA, "home", "end",
                                 c("home", "home_enddate")),
               "home_enddate")

  # An exact companion still wins over a prefix match.
  expect_equal(mpi_resolve_bound(NA, "home", "start",
                                 c("home_startdate", "home_start")),
               "home_start")

  # Nothing to find.
  expect_equal(mpi_resolve_bound(NA, "home", "start", c("home", "resid")),
               NA_character_)

  # A configured name is used as given.
  expect_equal(
    mpi_resolve_bound("med_i_end", "med_ih_drug", "end",
                      c("med_ih_drug", "med_i_end")),
    "med_i_end"
  )
})


test_that("primary_coding_suep_mpi() attaches the MPI tables to the object", {

  # The intended workflow: read_tsExport() -> primary_coding_suep() ->
  # primary_coding_suep_mpi(), attaching tables to the object the same way
  # primary_coding_suep_pcs_score() attaches ecu_pcs_score.
  trial_data <- suppressWarnings(primary_coding_suep_mpi(make_mpi_fixture()))

  expect_s3_class(trial_data, "tsExportdata")
  expect_true(all(c("ecu_mpi", "ecu_mpi_outcome") %in% names(trial_data)))

  expect_equal(nrow(trial_data$ecu_mpi), 6)
  expect_true("ecu_mpi" %in% names(trial_data$ecu_mpi))
  expect_true("ecu_rehosp" %in% names(trial_data$ecu_mpi_outcome))

  # The tables that were already in the object are left alone.
  expect_true(all(c("m", "fv2_1", "ehome", "emed_oth", "eward") %in%
                    names(trial_data)))
  expect_equal(trial_data$m, make_mpi_fixture()$m)

  # IDs are still taken from export_options when they are not named.
  expect_true("export_psn" %in% names(trial_data$ecu_mpi))
})


test_that("primary_coding_suep_mpi() can skip the outcome table", {

  trial_data <- suppressWarnings(
    primary_coding_suep_mpi(make_mpi_fixture(), outcomes = FALSE)
  )

  expect_true("ecu_mpi" %in% names(trial_data))
  expect_false("ecu_mpi_outcome" %in% names(trial_data))
})


test_that("primary_coding_suep_mpi() passes options on to build_mpi_df()", {

  trial_data <- suppressWarnings(
    primary_coding_suep_mpi(make_mpi_fixture(), min_domains = 8)
  )

  expect_true(all(is.na(trial_data$ecu_mpi$ecu_mpi)))
})


test_that("a plain list of tables works when both IDs are named", {

  # Exports without the secuTrial metadata tables cannot go through
  # read_tsExport(), so they arrive as a bare list of data frames read with
  # read.csv2() and carry no export_options at all.
  raw <- make_mpi_fixture()
  raw$export_options <- NULL
  raw <- unclass(raw)

  expect_error(build_mpi_df(raw), "pid")

  result <- suppressWarnings(
    build_mpi_df(raw, pid = "export_psn", visitid = "mnpcvpid")
  )
  reference <- suppressWarnings(build_mpi_df(make_mpi_fixture()))

  expect_equal(result, reference)

  outcome <- suppressWarnings(
    build_mpi_outcome_df(raw, pid = "export_psn", visitid = "mnpcvpid")
  )
  expect_equal(outcome,
               suppressWarnings(build_mpi_outcome_df(make_mpi_fixture())))

  report <- suppressWarnings(
    mpi_variable_report(raw, pid = "export_psn", visitid = "mnpcvpid")
  )
  expect_equal(nrow(report), nrow(epicodr::mpi_domains))
})


test_that("mpi_variable_report() shows how each variable is matched", {

  report <- suppressWarnings(mpi_variable_report(make_mpi_fixture()))

  expect_equal(nrow(report), nrow(epicodr::mpi_domains))

  home <- report[report$variable == "home", ]
  expect_true(home$found)
  expect_equal(home$matching, "interval")
  expect_equal(home$start, "home_start")
  expect_equal(home$companions, "home_adults,home_childr")

  # Cross sectional comorbidity form.
  expect_equal(report$matching[report$variable == "gec_cv"], "cross-sectional")
  # Longitudinal weight and patient reported outcome forms.
  expect_equal(report$matching[report$variable == "gec_weight"],
               "longitudinal")
  expect_equal(report$matching[report$variable == "pro_pc2r"], "longitudinal")
  # Tables the fixture does not have.
  expect_equal(report$matching[report$variable == "cfs_seid_crit4"],
               "no table")
  expect_false(any(report$found[report$table == "stv2"]))
})
