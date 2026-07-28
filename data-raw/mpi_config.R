# Multidimensional Prognostic Index (MPI) configuration ------------------------
#
# The MPI is aggregated from eight domains. Two objects define it:
#
#   mpi_domains  which variables feed each domain, and where they live
#   mpi_scoring  how each domain is turned into a 0 / 0.5 / 1 domain score
#
# These replace the FileGrouping.rds / Scoring.rds files that the original
# analysis script read from disk. Keeping them as package data makes the
# definition reviewable in version control and lets callers override single
# domains without editing code.
#
# Both objects are passed as arguments to build_mpi_df(), so a study that
# codes its variables differently can supply its own tables.


# Domain definition ------------------------------------------------------------
#
# domain      Domain the variable contributes to. Must match mpi_scoring$domain.
# table       tsExport table holding the variable (file name without ".csv").
# variable    Variable name.
# start       Variable holding the start of the validity interval, or NA.
# end         Variable holding the end of the validity interval, or NA.
#             When both are NA, "<variable>_start" / "<variable>_end" are used
#             if they exist in the table, otherwise the variable is not
#             treated as interval data.
# companions  Comma separated further variables that share the interval of
#             `variable` and are needed for scoring, or NA.
# aggregate   How to combine several records that are valid at the same visit:
#             "max"   highest code wins (severity scales, living situation)
#             "count" number of distinct non-missing entries (drug counts)

# Stack the domain blocks, filling in the optional columns a block leaves out.
mpi_bind <- function(...) {
  columns <- c("domain", "table", "variable",
               "start", "end", "companions", "aggregate")
  blocks <- lapply(list(...), function(block) {
    for (column in columns) {
      if (!column %in% names(block)) {
        block[[column]] <- if (column == "aggregate") "max" else NA_character_
      }
      if (column != "aggregate") {
        block[[column]] <- as.character(block[[column]])
      }
    }
    block[, columns]
  })
  do.call(rbind, blocks)
}

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
    stringsAsFactors = FALSE
  ),

  ## Social / living situation -------------------------------------------------
  # home_adults and home_childr share the validity interval of home and are
  # needed by the living_situation rules.
  data.frame(
    domain = "social",
    table = c("ehome", "eresid"),
    variable = c("home", "resid"),
    companions = c("home_adults,home_childr", NA),
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

rownames(mpi_domains) <- NULL


# Domain scoring ---------------------------------------------------------------
#
# domain   Domain name, also fixes the order of the domain score columns.
# type     Scoring rule to apply:
#          "sum"              sum the domain variables, then compare against
#                             the cut-offs. A cut-off is the *inclusive lower
#                             bound* of its score, so the highest cut-off that
#                             is still <= the sum decides the score.
#          "bmi"              derive BMI from gec_weight / gec_height and read
#                             the cut-offs as literal "<x" / ">x" / "lo-hi"
#                             bands.
#          "living_situation" apply the home / resid decision rules.
#          "none"             no cut-offs defined, domain scores NA.
# cut_0    Cut-off for domain score 0   (no problem).
# cut_0.5  Cut-off for domain score 0.5 (minor problem).
# cut_1    Cut-off for domain score 1   (major problem).

mpi_scoring <- data.frame(
  domain = c("comorbidity", "cognition", "nutrition", "adl",
             "mental", "social", "medications", "biomarker"),
  type = c("sum", "sum", "bmi", "sum",
           "sum", "living_situation", "sum", "none"),
  cut_0 = c("0", "0", "23-33", "0",
            "0", NA, "0", NA),
  cut_0.5 = c("1", "4", "<23,>33", "4",
              "12", NA, "2", NA),
  cut_1 = c("3", "7", "<20,>38", "7",
            "15", NA, "4", NA),
  stringsAsFactors = FALSE
)


# Consistency checks -----------------------------------------------------------

stopifnot(
  setequal(mpi_domains$domain, mpi_scoring$domain),
  !anyDuplicated(mpi_scoring$domain),
  all(mpi_domains$aggregate %in% c("max", "count")),
  all(mpi_scoring$type %in% c("sum", "bmi", "living_situation", "none"))
)

usethis::use_data(mpi_domains, overwrite = TRUE)
usethis::use_data(mpi_scoring, overwrite = TRUE)
