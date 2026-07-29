# Compare the packaged MPI against the original analysis script -----------------
#
# The rewrite deliberately changed several things, so the two implementations
# are not expected to agree out of the box. This script localises the
# disagreement to a single domain and a single rule, which is far easier to
# reason about than a difference in the final index.
#
# Run the original script first, so that `ascores` and `ma` are in the session,
# then run the package on the same export:
#
#   trial_data <- primary_coding_suep_mpi(trial_data)
#
# and source this file.


# 1. Per-visit domain scores of the original ----------------------------------
#
# ascores[[patient]] is a matrix of day (rows) by domain (columns), so the
# score of a visit is the row of that visit's documentation date.

original_long <- do.call(rbind, lapply(names(ascores), function(patient) {

  scores <- ascores[[patient]]
  visits <- ma[ma$export_psn == patient, ]
  visits <- visits[visits$datenr %in% rownames(scores), ]

  if (nrow(visits) == 0) {
    return(NULL)
  }

  out <- data.frame(
    export_psn = patient,
    datenr = as.numeric(visits$datenr),
    scores[visits$datenr, , drop = FALSE],
    stringsAsFactors = FALSE
  )
  out$mpi_original <- rowMeans(scores[visits$datenr, , drop = FALSE],
                               na.rm = TRUE)
  out
}))

original_long <- original_long[!duplicated(
  original_long[, c("export_psn", "datenr")]
), ]


# 2. Per-visit domain scores of the package -----------------------------------

packaged <- trial_data$ecu_mpi
packaged$datenr <- as.numeric(packaged$visit_date)

domains <- epicodr::mpi_scoring$domain


# 3. Join and compare ---------------------------------------------------------

both <- merge(
  packaged[, c("export_psn", "datenr", paste0("ecu_mpi_", domains),
               "ecu_mpi_n_domains", "ecu_mpi")],
  original_long,
  by = c("export_psn", "datenr"),
  suffixes = c("", "_orig")
)

cat("visits compared:", nrow(both), "\n")
cat("visits only in the package output:",
    nrow(packaged) - nrow(both), "\n\n")


## Domain by domain -----------------------------------------------------------
#
# `changed` counts visits where the two disagree, treating NA as a value so a
# score that appeared or disappeared is counted too.

differs <- function(a, b) (is.na(a) != is.na(b)) | (!is.na(a) & !is.na(b) & a != b)

summary_rows <- lapply(domains, function(domain) {

  packaged_score <- both[[paste0("ecu_mpi_", domain)]]
  original_score <- both[[domain]]

  if (is.null(original_score)) {
    return(NULL)
  }

  changed <- differs(packaged_score, original_score)

  data.frame(
    domain = domain,
    n_changed = sum(changed),
    pct_changed = round(100 * mean(changed), 1),
    was_na_now_scored = sum(is.na(original_score) & !is.na(packaged_score)),
    was_scored_now_na = sum(!is.na(original_score) & is.na(packaged_score)),
    package_higher = sum(!is.na(packaged_score) & !is.na(original_score) &
                           packaged_score > original_score),
    package_lower = sum(!is.na(packaged_score) & !is.na(original_score) &
                          packaged_score < original_score),
    stringsAsFactors = FALSE
  )
})

cat("--- domain level differences ---\n")
print(do.call(rbind, summary_rows), row.names = FALSE)


## Cross tabulation of the worst offender -------------------------------------
#
# A shifted cut-off shows up as a clean diagonal offset, a masked missing code
# as a scattered drop, and a lost domain as a column of NA.

worst <- do.call(rbind, summary_rows)
worst <- worst$domain[which.max(worst$n_changed)]

cat("\n--- ", worst, ": original (rows) against package (columns) ---\n",
    sep = "")
print(table(
  original = both[[worst]],
  package = both[[paste0("ecu_mpi_", worst)]],
  useNA = "ifany"
))


## Index level ----------------------------------------------------------------

cat("\n--- MPI ---\n")
cat("original: "); print(summary(both$mpi_original))
cat("package:  "); print(summary(both$ecu_mpi))

cat("\nvisits where the package withheld an MPI:",
    sum(is.na(both$ecu_mpi) & !is.na(both$mpi_original)), "\n")
cat("of those, number of scored domains:\n")
print(table(both$ecu_mpi_n_domains[is.na(both$ecu_mpi) &
                                     !is.na(both$mpi_original)]))


# 4. Switch the changes off one at a time -------------------------------------
#
# Reproduces the original more and more closely. If the numbers still differ
# once everything below is switched off, the cause is not one of the documented
# changes and worth reporting.

legacy <- build_mpi_df(
  trial_data,
  missing_values = numeric(0),  # do not mask -1 / 88 / 98 / 99
  min_domains = 1               # score a visit with a single domain
)
legacy$datenr <- as.numeric(legacy$visit_date)

check <- merge(
  legacy[, c("export_psn", "datenr", "ecu_mpi")],
  original_long[, c("export_psn", "datenr", "mpi_original")],
  by = c("export_psn", "datenr")
)

cat("\n--- with masking and the domain minimum switched off ---\n")
cat("visits still differing:",
    sum(differs(check$ecu_mpi, check$mpi_original)), "of", nrow(check), "\n")
cat("(what remains is the inclusive cut-off, the medication counting fix,\n")
cat(" and the all-missing-is-NA rule, none of which have a switch)\n")
