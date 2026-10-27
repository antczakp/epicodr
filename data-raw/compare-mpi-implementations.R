# Compare the packaged MPI against the original analysis script -----------------
#
# The rewrite deliberately changed several things, so the two implementations
# are not expected to agree out of the box. This script localises the
# disagreement to a single domain and a single rule, which is far easier to
# reason about than a difference in the final index.
#
# Run the original script first, so that `ascores` and `ma` are in the session.
# Then score the same export with the package. Two ways, and the first is the
# better one for a comparison because both sides then read byte-identical
# input through read.csv2(), leaving only the algorithms to differ:
#
#   ## a) reuse the tables the original script already read
#   mpi_pid <- "export_psn"
#   mpi_visitid <- "mnpcvpid"
#   trial_data <- data.all
#   trial_data$ecu_mpi <- build_mpi_df(trial_data, mpi_pid, mpi_visitid)
#
#   ## b) go through the package's own reader, which also exercises labelling
#   trial_data <- read_tsExport(zip_file_path, separator = ";", decimal = ",")
#   trial_data <- primary_coding_suep_mpi(trial_data)
#
# and source this file.


# A bare list of tables read with read.csv2() carries no export_options, so the
# ID columns cannot be derived and have to be named. Left NULL for path (b).
if (!exists("mpi_pid")) mpi_pid <- NULL
if (!exists("mpi_visitid")) mpi_visitid <- NULL


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
#
# Both sides are keyed by patient and day, because the original indexes its
# grid by date. Visits that share a date therefore cross-multiply here, and
# the original keeps only the last of them: its grid has one row per day, so a
# second visit on the same date overwrites the first. Counts below are a
# little inflated for those visits. Check how many there are with
#
#   key <- paste(visits[[mpi_pid]], visits$visit_date); table(table(key))
#

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


# 4. The index, on the visits both implementations score ----------------------
#
# The package withholds an MPI below `min_domains`, so the two summaries above
# are computed over different sets of visits and are not comparable. This is
# the like for like figure.

ok <- !is.na(both$ecu_mpi) & !is.na(both$mpi_original)

cat("\n--- MPI on the ", sum(ok), " visits both implementations score ---\n",
    sep = "")
cat("original: "); print(summary(both$mpi_original[ok]))
cat("package:  "); print(summary(both$ecu_mpi[ok]))
cat("visits differing:",
    sum(differs(both$ecu_mpi[ok], both$mpi_original[ok])), "\n")


# 5. What the domain minimum withholds ----------------------------------------
#
# `min_domains` is a deliberate choice rather than a difference to reconcile:
# an MPI averaged over three domains is not on the same scale as one averaged
# over seven, so pooling them biases whatever is computed downstream. This
# shows what the default costs and how the withheld visits would have scored.
#
# `relaxed` and `packaged` come from the same visits in the same order, so the
# withheld rows are taken directly. Joining them on a date would cross-multiply
# the visits that share one.

relaxed <- build_mpi_df(
  trial_data,
  pid = mpi_pid,
  visitid = mpi_visitid,
  min_domains = 1               # score a visit with a single domain
)

withheld <- relaxed[is.na(packaged$ecu_mpi), ]

cat("\n--- visits the domain minimum withholds ---\n")
cat("withheld:", nrow(withheld), "of", nrow(relaxed), "\n")
cat("scored domains:\n")
print(table(withheld$ecu_mpi_n_domains))
cat("MPI they would have been given:\n")
print(summary(withheld$ecu_mpi))
cat("MPI of the visits that are scored:\n")
print(summary(packaged$ecu_mpi))
