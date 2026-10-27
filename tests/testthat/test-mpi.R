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
                     gec_height = c(180, 170),
                     gec_birthdate = c("19540101", "19480815"),
                     gec_gender = c(1, 2),
                     gec_gender.factor = factor(
                       c("Männlich", "Weiblich"),
                       levels = c("Männlich", "Weiblich", "Divers")
                     ),
                     gec_ethnicity = c(2, 1),
                     gec_ethnicity.factor = factor(c("B", "A"),
                                                   levels = c("A", "B", "C")),
                     origin = c(1, 2),
                     origin.factor = factor(c("Home", "Nursing home"),
                                            levels = c("Home", "Nursing home")),
                     palliative = c(2, 1),
                     palliative.factor = factor(c("No", "Yes"),
                                                levels = c("Yes", "No")),
                     stringsAsFactors = FALSE),
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
  # cut-offs 0 / 3 / 6, applied to sums the floors have already shifted
  expect_equal(score_for("mental", c(0, 2, 3, 5, 6, 9)),
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


test_that("labelled factors are read as their level codes", {

  # An ordered scale whose levels are labels: the position of a level is the
  # code. as.character() would turn all of these into NA.
  eq5d_mob <- factor(
    c("Keine Probleme", "Grosse Probleme", NA, "Leichte Probleme"),
    levels = c("Keine Probleme", "Leichte Probleme", "Maessige Probleme",
               "Grosse Probleme", "Extreme Probleme"),
    ordered = TRUE
  )
  expect_equal(mpi_as_numeric(eq5d_mob), c(1, 4, NA, 2))

  # Levels that are the codes themselves keep their value, not their position.
  coded <- factor(c("0", "3", "7"), levels = c("0", "3", "7"))
  expect_equal(mpi_as_numeric(coded), c(0, 3, 7))

  # An unordered factor is treated the same way, since the export orders its
  # value labels by code.
  expect_equal(mpi_as_numeric(factor(c("nein", "ja"),
                                     levels = c("nein", "ja"))),
               c(1, 2))
})


test_that("labelled vectors are read as their underlying codes", {

  skip_if_not_installed("haven")

  x <- haven::labelled(c(1, 5, NA), labels = c(none = 1, extreme = 5))
  expect_equal(mpi_as_numeric(x), c(1, 5, NA))
})


test_that("a domain of label-coded factors is scored, not dropped", {

  # The reported failure: the variables of a domain arrive as labelled
  # factors, so the domain scored NA where it should have carried a score.
  trial_data <- make_mpi_fixture()
  as_label <- function(x) {
    levels <- c("keine", "leicht", "maessig", "gross", "extrem")
    factor(levels[x], levels = levels, ordered = TRUE)
  }
  trial_data$prom$eq5d_anx <- as_label(trial_data$prom$eq5d_anx)
  trial_data$prom$eq5d_mob <- as_label(trial_data$prom$eq5d_mob)

  result <- suppressWarnings(build_mpi_df(trial_data))
  reference <- suppressWarnings(build_mpi_df(make_mpi_fixture()))

  expect_false(any(is.na(result$ecu_mpi_mental)))
  expect_false(any(is.na(result$ecu_mpi_adl)))

  # Reading the levels recovers exactly the codes the fixture started from.
  expect_equal(result$ecu_mpi_mental, reference$ecu_mpi_mental)
  expect_equal(result$ecu_mpi_adl, reference$ecu_mpi_adl)
  expect_equal(result$ecu_mpi, reference$ecu_mpi)
})


test_that("a column that cannot be read as a code is reported", {

  trial_data <- make_mpi_fixture()
  trial_data$prom$eq5d_anx <- "keine Angabe"

  expect_warning(build_mpi_df(trial_data), "eq5d_anx")
})


test_that("generic missing codes are masked", {

  expect_equal(mpi_mask_missing(c(0, 1, -1, 99), c(-1, 88, 98, 99)),
               c(0, 1, NA, NA))
  expect_equal(mpi_mask_missing(c(0, 99), numeric(0)), c(0, 99))
})


test_that("values are carried forward but never backward", {

  # A visit before the first measurement stays missing.
  x <- c(NA, 5, NA, NA, 7, NA)
  patient <- c("a", "a", "a", "b", "b", "b")

  expect_equal(mpi_carry_forward(x, patient), c(NA, 5, 5, NA, 7, 7))

  # A new measurement takes over from the visit it was recorded at.
  expect_equal(
    mpi_carry_forward(c(80, NA, 60, NA), rep("a", 4)),
    c(80, 80, 60, 60)
  )

  # Patients do not bleed into each other.
  expect_equal(
    mpi_carry_forward(c(80, NA, NA, 60), c("a", "a", "b", "b")),
    c(80, 80, NA, 60)
  )
})


test_that("a new weight takes over the nutrition score from its visit on", {

  trial_data <- make_mpi_fixture()
  # P1 is 1.80 m and weighs 85 kg at baseline, then 40 kg at the 12M visit.
  trial_data$fuv1 <- rbind(
    trial_data$fuv1,
    data.frame(export_psn = "P1", mnpcvpid = "13", gec_weight = 40,
               stringsAsFactors = FALSE)
  )

  result <- suppressWarnings(build_mpi_df(trial_data))
  p1 <- result[result$export_psn == "P1", ]
  p1 <- p1[order(p1$visit_date), ]

  # BMI 26.2 at baseline, carried to the 3M visit, then 12.3 from 12M on.
  expect_equal(p1$ecu_mpi_nutrition, c(0, 0, 1))
})


test_that("a visit before the first measurement has no BMI", {

  trial_data <- make_mpi_fixture()
  # Height and weight documented only at the 12M visit.
  trial_data$bv1$mnpcvpid <- c("13", "23")
  trial_data$fuv1$mnpcvpid <- c("13", "23")

  result <- suppressWarnings(build_mpi_df(trial_data))
  p1 <- result[result$export_psn == "P1", ]
  p1 <- p1[order(p1$visit_date), ]

  expect_true(all(is.na(p1$ecu_mpi_nutrition[1:2])))
  expect_false(is.na(p1$ecu_mpi_nutrition[3]))
})

test_that("a missing height does not abort the nutrition domain", {

  # calculate_bmi() judges the unit from the mean height, which is NaN when
  # every height is missing. That used to error instead of leaving the BMI
  # missing, taking the whole MPI down with it.
  expect_true(all(is.na(calculate_bmi(c(70, 80), c(NA_real_, NA_real_)))))

  trial_data <- make_mpi_fixture()
  trial_data$bv1$gec_height <- NA_real_

  result <- suppressWarnings(build_mpi_df(trial_data))

  expect_true(all(is.na(result$ecu_mpi_nutrition)))
  expect_false(any(is.na(result$ecu_mpi)))
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


# Age --------------------------------------------------------------------------

test_that("ecu_age is the age in decimal years at every visit", {

  result <- suppressWarnings(build_mpi_df(make_mpi_fixture()))

  # P1 was born 1954-01-01, so baseline and 12M fall on birthdays and the 3M
  # visit is 90 days into the 365 day year.
  p1 <- result[result$export_psn == "P1", ]
  expect_equal(p1$ecu_age, c(67, 67 + 90 / 365, 68))

  # P2 was born 1948-08-15. The age moves between visits within one year of
  # life, which completed years would not show.
  p2 <- result[result$export_psn == "P2", ]
  expect_equal(p2$ecu_age, c(72 + 170 / 365, 72 + 259 / 365, 73 + 170 / 365))

  # Completed years are the floor of it.
  expect_equal(floor(p2$ecu_age), c(72, 72, 73))

  # The baseline age is the value of the baseline row.
  expect_equal(result$ecu_age[result$mnpvislabel == "Baseline"],
               c(67, 72 + 170 / 365))
})


test_that("the converted .date companion of the date of birth is used", {

  trial_data <- make_mpi_fixture()
  trial_data$bv1$gec_birthdate.date <- as.Date(c("1954-01-01", "1948-08-15"))
  # The raw variable would give a different age, so this shows which is read.
  trial_data$bv1$gec_birthdate <- c("19600101", "19600101")

  result <- suppressWarnings(build_mpi_df(trial_data))

  expect_equal(result$ecu_age[result$export_psn == "P1"],
               c(67, 67 + 90 / 365, 68))
})


test_that("a date of birth that is not a complete date is reported", {

  trial_data <- make_mpi_fixture()
  trial_data$bv1$gec_birthdate <- c("19540101", "19480000")

  expect_warning(result <- build_mpi_df(trial_data), "gec_birthdate")
  expect_equal(result$ecu_age[result$export_psn == "P1"],
               c(67, 67 + 90 / 365, 68))
  expect_true(all(is.na(result$ecu_age[result$export_psn == "P2"])))
})


test_that("a missing date of birth makes the age NA, not the MPI", {

  trial_data <- make_mpi_fixture()
  trial_data$bv1$gec_birthdate <- NULL

  expect_warning(result <- build_mpi_df(trial_data), "gec_birthdate")
  expect_true(all(is.na(result$ecu_age)))
  expect_equal(result$ecu_mpi,
               suppressWarnings(build_mpi_df(make_mpi_fixture()))$ecu_mpi)
})


test_that("the sex is the .factor companion, carried on every visit", {

  result <- suppressWarnings(build_mpi_df(make_mpi_fixture()))

  expect_s3_class(result$ecu_sex, "factor")
  expect_equal(levels(result$ecu_sex), c("Männlich", "Weiblich", "Divers"))
  expect_equal(as.character(result$ecu_sex[result$export_psn == "P1"]),
               rep("Männlich", 3))
  expect_equal(as.character(result$ecu_sex[result$export_psn == "P2"]),
               rep("Weiblich", 3))
})


test_that("the raw sex code is used without a .factor companion", {

  trial_data <- make_mpi_fixture()
  trial_data$bv1$gec_gender.factor <- NULL

  result <- suppressWarnings(build_mpi_df(trial_data))

  expect_equal(result$ecu_sex[result$export_psn == "P2"], rep(2, 3))
})


test_that("the ethnicity is the .factor companion, carried on every visit", {

  result <- suppressWarnings(build_mpi_df(make_mpi_fixture()))

  expect_s3_class(result$ecu_ethnicity, "factor")
  expect_equal(levels(result$ecu_ethnicity), c("A", "B", "C"))
  expect_equal(as.character(result$ecu_ethnicity[result$export_psn == "P1"]),
               rep("B", 3))
  expect_equal(as.character(result$ecu_ethnicity[result$export_psn == "P2"]),
               rep("A", 3))
})


test_that("a missing ethnicity makes the ethnicity NA, not the MPI", {

  trial_data <- make_mpi_fixture()
  trial_data$bv1$gec_ethnicity <- NULL
  trial_data$bv1$gec_ethnicity.factor <- NULL

  expect_warning(result <- build_mpi_df(trial_data), "gec_ethnicity")
  expect_true(all(is.na(result$ecu_ethnicity)))
  # The sex read from the same table is unaffected.
  expect_false(any(is.na(result$ecu_sex)))
  expect_equal(result$ecu_mpi,
               suppressWarnings(build_mpi_df(make_mpi_fixture()))$ecu_mpi)
})


test_that("an empty ethnicity column warns instead of turning NA silently", {

  trial_data <- make_mpi_fixture()
  trial_data$bv1$gec_ethnicity.factor[] <- NA

  expect_warning(result <- build_mpi_df(trial_data),
                 "'gec_ethnicity.factor' in table 'bv1' holds no documented")
  expect_true(all(is.na(result$ecu_ethnicity)))
})


test_that("origin and palliative status are carried on every visit", {

  result <- suppressWarnings(build_mpi_df(make_mpi_fixture()))

  expect_s3_class(result$ecu_origin, "factor")
  expect_equal(as.character(result$ecu_origin[result$export_psn == "P1"]),
               rep("Home", 3))
  expect_equal(as.character(result$ecu_origin[result$export_psn == "P2"]),
               rep("Nursing home", 3))

  expect_s3_class(result$ecu_palliative, "factor")
  expect_equal(levels(result$ecu_palliative), c("Yes", "No"))
  expect_equal(as.character(result$ecu_palliative[result$export_psn == "P1"]),
               rep("No", 3))
  expect_equal(as.character(result$ecu_palliative[result$export_psn == "P2"]),
               rep("Yes", 3))
})


test_that("missing origin or palliative status makes only that column NA", {

  trial_data <- make_mpi_fixture()
  trial_data$bv1$origin <- NULL
  trial_data$bv1$origin.factor <- NULL
  trial_data$bv1$palliative <- NULL
  trial_data$bv1$palliative.factor <- NULL

  warnings <- character(0)
  result <- withCallingHandlers(
    build_mpi_df(trial_data),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("'origin'", warnings)))
  expect_true(any(grepl("'palliative'", warnings)))
  expect_true(all(is.na(result$ecu_origin)))
  expect_true(all(is.na(result$ecu_palliative)))
  expect_false(any(is.na(result$ecu_ethnicity)))
  expect_equal(result$ecu_mpi,
               suppressWarnings(build_mpi_df(make_mpi_fixture()))$ecu_mpi)
})


test_that("a missing sex makes the sex NA, not the MPI", {

  trial_data <- make_mpi_fixture()
  trial_data$bv1$gec_gender <- NULL
  trial_data$bv1$gec_gender.factor <- NULL

  expect_warning(result <- build_mpi_df(trial_data), "gec_gender")
  expect_true(all(is.na(result$ecu_sex)))
  expect_equal(result$ecu_mpi,
               suppressWarnings(build_mpi_df(make_mpi_fixture()))$ecu_mpi)
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


# Vital status -----------------------------------------------------------------

# Adds patients whose follow-up is cut short by death, one per rule:
#   P3 died two weeks after baseline and never had a follow-up visit.
#   P4 died after the 12 month horizon.
#   P5 died on the last day of the 12 month horizon.
#   P6 has no dated baseline visit.
make_vital_fixture <- function() {

  trial_data <- make_mpi_fixture()

  trial_data$m <- rbind(
    trial_data$m,
    data.frame(
      export_psn = c("P3", "P4", "P4", "P5", "P6", "P6"),
      mnpcvpid = c("31", "41", "42", "51", "61", "62"),
      mnpvisid = c("7", "8", "9", "10", "11", "12"),
      mnpvislabel = c("Baseline", "Baseline", "3M Follow-Up", "Baseline",
                      "Baseline", "3M Follow-Up"),
      mnpvisno = c(1, 1, 2, 1, 1, 2),
      pr_docudate = c("20210101", "20210101", "20210401", "20210101",
                      NA, "20210401"),
      pr_visit_mode = 1,
      stringsAsFactors = FALSE
    )
  )

  trial_data$fv2_6 <- rbind(
    trial_data$fv2_6,
    data.frame(
      export_psn = c("P3", "P4", "P5", "P6"),
      diag_oth = 0,
      death_date = c("20210115", "20220102", "20220101", "20210601"),
      stringsAsFactors = FALSE
    )
  )

  return(trial_data)
}


test_that("build_mpi_vital_df() reports every patient once", {

  result <- build_mpi_vital_df(make_vital_fixture())

  expect_equal(sort(result$export_psn), paste0("P", 1:6))
  expect_true(inherits(result$ecu_death_date, "Date"))
  expect_equal(result$ecu_death_date[result$export_psn == "P2"],
               as.Date("2021-06-01"))
  expect_true(is.na(result$ecu_death_date[result$export_psn == "P1"]))
})


test_that("ecu_alive reports survival to the end of the study", {

  result <- build_mpi_vital_df(make_vital_fixture())
  alive <- stats::setNames(result$ecu_alive, result$export_psn)

  expect_equal(unname(alive[paste0("P", 1:6)]),
               c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
})


test_that("ecu_alive_12M censors deaths after 12 months from baseline", {

  result <- build_mpi_vital_df(make_vital_fixture())
  alive <- stats::setNames(result$ecu_alive_12M, result$export_psn)

  # P2 died 4 months after baseline, P3 before any follow-up visit.
  expect_false(alive[["P2"]])
  expect_false(alive[["P3"]])
  # P4 died one day past the horizon, P5 on its last day.
  expect_true(alive[["P4"]])
  expect_false(alive[["P5"]])
  expect_true(alive[["P1"]])
  # Without a baseline date there is no horizon to censor at.
  expect_true(is.na(alive[["P6"]]))
  expect_false(result$ecu_alive[result$export_psn == "P6"])
})


test_that("ecu_alive_3M censors deaths after 3 months from baseline", {

  result <- build_mpi_vital_df(make_vital_fixture())
  alive <- stats::setNames(result$ecu_alive_3M, result$export_psn)

  expect_equal(names(result),
               c("export_psn", "ecu_death_date", "ecu_alive", "ecu_alive_3M",
                 "ecu_alive_12M"))

  # P3 died two weeks after baseline.
  expect_false(alive[["P3"]])
  # P2 died 4 months after baseline, P4 and P5 after a year.
  expect_true(alive[["P2"]])
  expect_true(alive[["P4"]])
  expect_true(alive[["P5"]])
  expect_true(alive[["P1"]])
  expect_true(is.na(alive[["P6"]]))
})


test_that("a death on the last day of the 3 month horizon counts", {

  trial_data <- make_mpi_fixture()

  # P2's baseline is 2021-02-01, so the horizon ends on 2021-05-01.
  trial_data$fv2_6$death_date[trial_data$fv2_6$export_psn == "P2"] <-
    "20210501"
  result <- build_mpi_vital_df(trial_data)
  expect_false(result$ecu_alive_3M[result$export_psn == "P2"])

  trial_data$fv2_6$death_date[trial_data$fv2_6$export_psn == "P2"] <-
    "20210502"
  result <- build_mpi_vital_df(trial_data)
  expect_true(result$ecu_alive_3M[result$export_psn == "P2"])
})


test_that("the 3 month horizon rolls back to the end of a shorter month", {

  trial_data <- make_mpi_fixture()
  # 30 November plus 3 months has no 30 February to land on.
  trial_data$m$pr_docudate[trial_data$m$mnpcvpid == "21"] <- "20201130"
  trial_data$fv2_6$death_date[trial_data$fv2_6$export_psn == "P2"] <-
    "20210228"

  result <- build_mpi_vital_df(trial_data)

  expect_false(result$ecu_alive_3M[result$export_psn == "P2"])
})


test_that("the 12 month horizon rolls back to the end of a shorter month", {

  trial_data <- make_mpi_fixture()
  trial_data$m$pr_docudate[trial_data$m$mnpcvpid == "21"] <- "20200229"
  trial_data$fv2_6$death_date[trial_data$fv2_6$export_psn == "P2"] <-
    "20210228"

  result <- build_mpi_vital_df(trial_data)

  expect_false(result$ecu_alive_12M[result$export_psn == "P2"])
})


test_that("a missing death table makes the vital status NA, not alive", {

  trial_data <- make_mpi_fixture()
  trial_data$fv2_6 <- NULL

  expect_warning(result <- build_mpi_vital_df(trial_data), "fv2_6")
  expect_true(all(is.na(result$ecu_alive)))
  expect_true(all(is.na(result$ecu_alive_3M)))
  expect_true(all(is.na(result$ecu_alive_12M)))
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
})


# Patient level rehospitalisation ----------------------------------------------

# Adds stays of P1 that exercise the ward groups, on top of the 5 days on
# ward 1 in February 2021:
#   ward 3 for 2 days before the 3M visit,
#   ward 4 for 10 days between the 3M and the 12M visit,
#   ward 5 for 3 days, which belongs to neither group.
make_rehosp_fixture <- function() {

  trial_data <- make_mpi_fixture()

  trial_data$eward <- rbind(
    trial_data$eward,
    data.frame(export_psn = "P1", gec_ward = c(3, 4, 5),
               ward_start = c("20210301", "20210601", "20210701"),
               ward_end = c("20210302", "20210610", "20210703"),
               stringsAsFactors = FALSE)
  )

  return(trial_data)
}


test_that("build_mpi_rehosp_df() sums ward days per group and follow-up", {

  result <- suppressWarnings(build_mpi_rehosp_df(make_rehosp_fixture()))

  expect_equal(sort(result$export_psn), c("P1", "P2"))
  expect_equal(names(result),
               c("export_psn", "ecu_rehosp_1_3_3M", "ecu_rehosp_4_3M",
                 "ecu_rehosp_1_3_12M", "ecu_rehosp_4_12M"))

  p1 <- result[result$export_psn == "P1", ]
  expect_equal(p1$ecu_rehosp_1_3_3M, 7)
  expect_equal(p1$ecu_rehosp_4_3M, 0)
  # The 12 month window starts at baseline, so it includes the 3 month days.
  expect_equal(p1$ecu_rehosp_1_3_12M, 7)
  expect_equal(p1$ecu_rehosp_4_12M, 10)

  # P2 never went to a ward.
  p2 <- result[result$export_psn == "P2", ]
  expect_equal(unlist(p2[, -1], use.names = FALSE), rep(0, 4))
})


test_that("build_mpi_rehosp_df() is NA without a follow-up visit", {

  result <- suppressWarnings(build_mpi_rehosp_df(make_vital_fixture()))
  rehosp <- result[match(paste0("P", 1:6), result$export_psn), ]

  # P3 and P5 only have a baseline visit.
  expect_true(all(is.na(rehosp[rehosp$export_psn %in% c("P3", "P5"), -1])))
  # P4 has a 3M but no 12M visit.
  expect_equal(rehosp$ecu_rehosp_4_3M[rehosp$export_psn == "P4"], 0)
  expect_true(is.na(rehosp$ecu_rehosp_4_12M[rehosp$export_psn == "P4"]))
  # P6 has no dated baseline visit, so there is no window.
  expect_true(all(is.na(rehosp[rehosp$export_psn == "P6", -1])))
})


test_that("build_mpi_rehosp_df() is NA when the ward table is missing", {

  trial_data <- make_mpi_fixture()
  trial_data$eward <- NULL

  expect_warning(result <- build_mpi_rehosp_df(trial_data), "eward")
  expect_true(all(is.na(result[, -1])))
})


test_that("build_mpi_rehosp_df() names a ward group the outcome lacks", {

  outcome <- suppressWarnings(
    build_mpi_outcome_df(make_mpi_fixture(), ward_levels = 1:3)
  )

  expect_error(
    build_mpi_rehosp_df(make_mpi_fixture(), outcome = outcome),
    "ecu_rehosp_4"
  )
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
  expect_false("ecu_alive" %in% names(trial_data$ecu_mpi_outcome))

  # Vital status is carried on every visit of the patient.
  p2 <- trial_data$ecu_mpi[trial_data$ecu_mpi$export_psn == "P2", ]
  expect_equal(p2$ecu_alive, rep(FALSE, 3))
  expect_equal(p2$ecu_alive_3M, rep(TRUE, 3))
  expect_equal(p2$ecu_alive_12M, rep(FALSE, 3))
  expect_equal(p2$ecu_death_date, rep(as.Date("2021-06-01"), 3))
  expect_equal(p2$ecu_age, c(72 + 170 / 365, 72 + 259 / 365, 73 + 170 / 365))

  # So are the ward days per group, while the detail stays in the outcome.
  p1 <- trial_data$ecu_mpi[trial_data$ecu_mpi$export_psn == "P1", ]
  expect_equal(p1$ecu_rehosp_1_3_3M, rep(5, 3))
  expect_equal(p1$ecu_rehosp_4_12M, rep(0, 3))
  expect_false("ecu_rehosp_1" %in% names(trial_data$ecu_mpi))
  expect_true("ecu_rehosp_1" %in% names(trial_data$ecu_mpi_outcome))

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
  expect_false("ecu_alive" %in% names(trial_data$ecu_mpi))
  expect_false("ecu_rehosp_1_3_3M" %in% names(trial_data$ecu_mpi))
  # Age and sex describe the patient, not an outcome, so they stay.
  expect_true(all(c("ecu_age", "ecu_sex", "ecu_ethnicity", "ecu_origin",
                    "ecu_palliative") %in% names(trial_data$ecu_mpi)))
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

  vital <- build_mpi_vital_df(raw, pid = "export_psn", visitid = "mnpcvpid")
  expect_equal(vital, build_mpi_vital_df(make_mpi_fixture()))

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


# Configuration ----------------------------------------------------------------

test_that("the shipped configuration is valid", {

  expect_true(mpi_validate_config())
  expect_true(is.data.frame(epicodr::mpi_domains))
  expect_true(is.data.frame(epicodr::mpi_scoring))
  expect_setequal(epicodr::mpi_domains$domain, epicodr::mpi_scoring$domain)
})


test_that("mpi_validate_config() names what is wrong", {

  domains <- epicodr::mpi_domains
  scoring <- epicodr::mpi_scoring

  # Cut-offs that go down make the higher score unreachable.
  broken <- scoring
  broken$cut_1[broken$domain == "comorbidity"] <- "0"
  expect_error(mpi_validate_config(domains, broken), "decrease")

  broken <- scoring
  broken$cut_0.5[broken$domain == "comorbidity"] <- "many"
  expect_error(mpi_validate_config(domains, broken), "must all be numbers")

  broken <- scoring
  broken$cut_1[broken$domain == "nutrition"] <- "somewhat high"
  expect_error(mpi_validate_config(domains, broken), "holds no band")

  broken <- scoring
  broken$type[broken$domain == "mental"] <- "median"
  expect_error(mpi_validate_config(domains, broken), "unknown scoring type")

  expect_error(
    mpi_validate_config(domains, scoring[scoring$domain != "mental", ]),
    "no scoring rule"
  )

  # A scored domain nothing feeds can only ever score NA.
  expect_error(
    mpi_validate_config(domains[domains$domain != "mental", ], scoring),
    "no variables"
  )

  broken <- domains
  broken$aggregate[1] <- "mean"
  expect_error(mpi_validate_config(broken, scoring), "unknown aggregate")
})


test_that("build_mpi_df() rejects a configuration it cannot score", {

  broken <- epicodr::mpi_scoring
  broken$type[broken$domain == "adl"] <- "median"

  expect_error(build_mpi_df(make_mpi_fixture(), scoring = broken),
               "unknown scoring type")
})


test_that("a visit id shared across patients is rejected", {

  # mnpvisid names the visit of the study, so every patient carries the same
  # values. Deduplicating on it used to collapse the table to one row per
  # visit type, and mpi_collect_domain() then failed on duplicated factor
  # levels rather than on the wrong id.
  trial_data <- make_mpi_fixture()
  trial_data$m$mnpvisid <- rep(c("1", "2", "3"), 2)

  expect_error(
    build_mpi_visit_df(trial_data, pid = "export_psn", visitid = "mnpvisid"),
    "shared by several patients"
  )
  expect_error(
    build_mpi_df(trial_data, pid = "export_psn", visitid = "mnpvisid"),
    "shared by several patients"
  )

  # The per-patient id still works.
  expect_equal(
    nrow(build_mpi_visit_df(trial_data, "export_psn", "mnpcvpid")), 6
  )
})


# Floors -----------------------------------------------------------------------

test_that("a floor is subtracted from a variable before the sum", {

  cuts <- c("0" = 0, "0.5" = 1, "1" = 3)

  # Two graded items answered at their best contribute nothing.
  values <- data.frame(a = c(1, 2, 5), b = c(1, 1, 3))
  floors <- c(a = 1, b = 1)

  # shifted sums: 0, 1, 6
  expect_equal(mpi_score_sum(values, cuts, floors), c(0, 0.5, 1))

  # Without floors the same data reads as 2, 3, 8.
  expect_equal(mpi_score_sum(values, cuts), c(0.5, 1, 1))
})


test_that("a value below its floor does not subtract from the domain", {

  cuts <- c("0" = 0, "0.5" = 1, "1" = 3)
  values <- data.frame(a = c(0, 4), b = c(4, 4))

  # a = 0 is below its floor of 1, so it contributes 0 rather than -1.
  expect_equal(mpi_score_sum(values, cuts, c(a = 1, b = 1)), c(1, 1))
})


test_that("floors leave unobserved variables unobserved", {

  cuts <- c("0" = 0, "0.5" = 1, "1" = 3)
  values <- data.frame(a = c(NA_real_, 1), b = c(NA_real_, NA_real_))

  expect_equal(mpi_score_sum(values, cuts, c(a = 1, b = 1)), c(NA, 0))
})


test_that("a variable the floors do not name keeps a floor of zero", {

  cuts <- c("0" = 0, "0.5" = 1, "1" = 3)
  values <- data.frame(flag = c(0, 1), item = c(1, 3))

  # flag is untouched, item loses its floor of 1: sums are 0 and 3.
  expect_equal(mpi_score_sum(values, cuts, c(item = 1)), c(0, 1))
})


test_that("only the mental domain carries non-zero floors", {

  floors <- epicodr::mpi_domains
  expect_true("floor" %in% names(floors))
  expect_setequal(unique(floors$domain[floors$floor > 0]), "mental")

  # The nine graded items: eq5d_anx and the eight PROMIS-29 items.
  expect_equal(sum(floors$floor), 9)

  # The cut-offs were lowered by that shift, so a patient who answers every
  # graded item lands where the previous 12 / 15 cut-offs put them.
  mental <- epicodr::mpi_scoring[epicodr::mpi_scoring$domain == "mental", ]
  expect_equal(as.numeric(c(mental$cut_0, mental$cut_0.5, mental$cut_1)),
               c(0, 3, 6))
})


test_that("answering every graded mental item at its best scores 0", {

  # The reported problem: EQ-5D and PROMIS items start at 1, so a fully
  # answered mental domain used to reach the 0.5 cut-off on a healthy patient.
  trial_data <- make_mpi_fixture()
  trial_data$prom$eq5d_anx <- 1

  result <- suppressWarnings(build_mpi_df(trial_data))
  expect_true(all(result$ecu_mpi_mental == 0))

  # And a patient at the top of the scale still reaches a problem score.
  trial_data$prom$eq5d_anx <- 5
  result <- suppressWarnings(build_mpi_df(trial_data))
  expect_true(all(result$ecu_mpi_mental == 0.5))
})


test_that("mpi_validate_config() rejects an unusable floor", {

  broken <- epicodr::mpi_domains
  broken$floor[broken$variable == "eq5d_anx"] <- -1
  expect_error(mpi_validate_config(broken), "floor")

  broken$floor[broken$variable == "eq5d_anx"] <- NA_real_
  expect_error(mpi_validate_config(broken), "floor")

  # A domain table written before floors existed still validates.
  legacy <- epicodr::mpi_domains
  legacy$floor <- NULL
  expect_true(mpi_validate_config(legacy))
})


# Multi code cells -------------------------------------------------------------

test_that("only the social variables reduce a multi code cell", {

  d <- epicodr::mpi_domains
  expect_true("multi" %in% names(d))
  expect_setequal(unique(d$domain[d$multi == "max"]), "social")
  expect_true(all(d$multi[d$domain != "social"] == "na"))
})


test_that("mpi_as_numeric() resolves a multi code cell as told", {

  expect_equal(mpi_as_numeric(c("1|3", "2", NA), multi = "max"), c(3, 2, NA))
  expect_equal(mpi_as_numeric(c("1|3", "2", NA), multi = "na"), c(NA, 2, NA))

  # A vector without any multi code cell reads the same either way.
  expect_equal(mpi_as_numeric(c("1", "2"), multi = "na"), c(1, 2))
})


test_that("a multi code cell outside social is reported, not collapsed", {

  trial_data <- make_mpi_fixture()
  trial_data$fv2_1$gec_cv <- c("1|3", "0")

  expect_warning(build_mpi_df(trial_data), "several codes")

  result <- suppressWarnings(build_mpi_df(trial_data))
  # P1 loses gec_cv, so its comorbidity rests on gec_pd alone: sum 1 -> 0.5.
  expect_equal(unique(result$ecu_mpi_comorbidity[result$export_psn == "P1"]),
               0.5)
})


test_that("a social cell holding several codes keeps the highest", {

  trial_data <- make_mpi_fixture()
  # P2 is recorded both as living at home and in an institution.
  trial_data$ehome$home <- c("2", "1|2")

  result <- suppressWarnings(build_mpi_df(trial_data))
  # max(1, 2) is 2, so P2 is still read as living at home, alone -> 0.5.
  expect_equal(unique(result$ecu_mpi_social[result$export_psn == "P2"]), 0.5)
})


test_that("an empty slot of a drug cell is not counted as a drug", {

  # The export writes an unused slot of a multi select cell out as a token,
  # so "MDA343|NA" is one drug and not two.
  trial_data <- make_mpi_fixture()
  trial_data$emed_oth$gec_med_drug <- c("MDA343|NA", "B", "C|NA")

  result <- suppressWarnings(build_mpi_df(trial_data))
  p1 <- result[result$export_psn == "P1", ]

  # Two drugs at baseline, one after the second ends. Cut-offs 0 / 2 / 4.
  expect_equal(p1$ecu_mpi_medications[p1$mnpvislabel == "Baseline"], 0.5)
  expect_equal(p1$ecu_mpi_medications[p1$mnpvislabel == "3M Follow-Up"], 0)
  # P2 has a single drug beside an empty slot.
  expect_equal(unique(result$ecu_mpi_medications[result$export_psn == "P2"]), 0)
})


test_that("mpi_validate_config() rejects an unknown multi", {

  broken <- epicodr::mpi_domains
  broken$multi[1] <- "mean"
  expect_error(mpi_validate_config(broken), "unknown multi")

  # A domain table written before the column existed still validates.
  legacy <- epicodr::mpi_domains
  legacy$multi <- NULL
  expect_true(mpi_validate_config(legacy))
})
