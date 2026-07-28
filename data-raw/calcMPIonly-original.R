fs <- dir()
data.all <- lapply(fs, read.csv2)
names(data.all) <- gsub(".csv$","",fs)

# # Specifiy location where to find data zip file
# zip_file_path <- "FileSenderDownload_17-6-2026__10-59-59/tsExport_tsExport_Polidori_NAPKON_Frailty_Index_SUEP_Nachreich.zip"
#
# # Import the zip file data
# data.all <- zip_file_path %>%
#   read_tsExport(separator = ";", decimal = ",")

nms <- sapply(data.all, names)
nms.exp <- names(which(sapply(nms, function(x)
  any(x == "export_psn"))))
nms.vis <- names(which(sapply(nms, function(x)
  any(x == "mnpvisid"))))
nms.start <- names(which(sapply(nms, function(x)
  length(grep(
    "_start$", x
  )) > 0)))

nms.long <- intersect(nms.exp, nms.vis)
nms.cros <- setdiff(nms.exp, nms.vis)
nms.long2 <- setdiff(nms.cros, nms.start)
nms.attr <- setdiff(nms.vis, nms.exp)

grps <- readRDS("../../tsExport_Polidori_NAPKON_Frailty_Index_SUEP/FileGrouping.rds")
scores <- readRDS("../../tsExport_Polidori_NAPKON_Frailty_Index_SUEP/Scoring.rds")

grps <- lapply(grps, function(x) {
  x[, 1] <- gsub(".csv", "", x[, 1])
  return(x)
})

data.res <- lapply(grps, function(x) {
  res <- apply(x, 1, function(z) {
    tmp <- data.all[[z[1]]]
    if (length(z) == 2) {
      return(tmp[, c(
        grep("export_psn", colnames(tmp), value = T),
        grep("^mnp", colnames(tmp), value = T),
        z[2],
        grep(paste0("^", z[2], "_[seac][tndh]"), colnames(tmp), value =
               T)
      )])
    } else{
      tmp2 <- tmp[, c(grep("export_psn", colnames(tmp), value = T), z[2:length(z)])]
      colnames(tmp2) <- c(colnames(tmp2)[1:2], paste(z[2], names(z)[3:length(z)], sep =
                                                       "_"))
      return(tmp2)
    }
  })
  return(res)
})

m.tmp <- data.all$m
m2.tmp <- data.all$m2
b.tmp <- data.all$scv
b.tmp$pr_visit_mode <- NA
b.tmp <- b.tmp[, c(
  "export_psn",
  "mnpvisid",
  "mnpcvpid",
  "mnpvislabel",
  "mnpvisno",
  "pr_visit_mode",
  "gec_pr_incl_date"
)]
m.tmp <- m.tmp[, c(
  "export_psn",
  "mnpvisid",
  "mnpcvpid",
  "mnpvislabel",
  "mnpvisno",
  "pr_visit_mode",
  "pr_docudate"
)]
m2.tmp <- m2.tmp[, c(
  "export_psn",
  "mnpvisid",
  "mnpcvpid",
  "mnpvislabel",
  "mnpvisno",
  "pr_visit_mode_1",
  "gec_pr_docudate_1"
)]
colnames(b.tmp) <- colnames(m2.tmp) <- colnames(m.tmp)
ma <- rbind(b.tmp, m.tmp, m2.tmp)
ma$datenr <- as.character(as.numeric(as.Date(as.character(ma$pr_docudate), format =
                                               "%Y%m%d")))

pat_id <- unique(unlist(sapply(data.res, function(x)
  sapply(x, function(x)
    x$export_psn))))
fullres <- lapply(pat_id, function(id) {
  m <- subset(ma, export_psn == id)
  m <- m[order(m$pr_docudate), ]
  m$pr_docudate <- as.Date(as.character(m$pr_docudate), format = "%Y%m%d")
  rownames(m) <- m$mnpcvpid
  dater <- seq(min(m$pr_docudate), max(m$pr_docudate), by = 1)
  .id <<- id
  
  dg <- lapply(data.res, function(dr) {
    .dr <<- dr
    dr <- lapply(dr, data.frame)
    dats <- unique(unlist(sapply(dr, function(x)
      setdiff(
        names(x),
        c(
          "export_psn",
          "mnpvisid",
          "mnpcvpid",
          "mnpvislabel",
          "mnpvisno",
          "mnpdocid"
        )
      ))))
    dats.s <- grep("_start.*", dats, invert = T, value = T)
    dats.s <- grep("_end.*", dats.s, invert = T, value = T)
    
    mm <- matrix(NA, nrow = length(dater), ncol = length(dats.s))
    rownames(mm) <- dater
    colnames(mm) <- dats.s
    mm <- data.frame(mm)
    
    for (i in 1:length(dr)) {
      .i <<- i
      tmp <- subset(dr[[i]], export_psn == id)
      if (nrow(tmp) == 0)
        next
      is.long <- any(names(tmp) == "mnpcvpid")
      tmp.var <- setdiff(
        names(tmp),
        c(
          "export_psn",
          "mnpvisid",
          "mnpcvpid",
          "mnpvislabel",
          "mnpvisno",
          "mnpdocid"
        )
      )
      tmp.var <- tapply(tmp.var, sapply(strsplit(tmp.var, "_"), function(x)
        x[1]), function(x)
          grep("_start_|_end_", x, value = T, invert = T), simplify = F)
      if (is.long) {
        for (k in 1:length(tmp.var)) {
          mm[as.character(as.numeric(m[as.character(tmp$mnpcvpid), "pr_docudate"])), tmp.var[[k]]] <- tmp[[tmp.var[[k]]]]
        }
      } else{
        for (j in 1:length(tmp.var)) {
          tmp.var.tmp <- tmp.var[[j]]
          is.otherlong <- length(grep("_start|_end", tmp.var.tmp)) > 0
          is.adds <- length(grep("_uk$", tmp.var.tmp)) > 0
          if (is.adds) {
            uks <- grep("_uk$", tmp.var.tmp, value = T)
            uks.var <- gsub("_uk$", "", uks)
            for (k in 1:length(uks)) {
              sel <- which(tmp[, uks[k]] == -1)
              if (length(sel) > 0)
                tmp[sel, uks.var[k]] <- NA
            }
            tmp <- tmp[, -match(uks, colnames(tmp))]
            tmp.var.tmp <- tmp.var.tmp[-match(uks, tmp.var.tmp)]
            mm <- mm[, -match(uks, colnames(mm))]
          }
          if (is.otherlong) {
            rownames.list <- function(x)
              attr(x, "row.names")
            class(mm) <- "list"
            starts <- grep("_start", tmp.var.tmp, value = T)
            ends <- grep("_end", tmp.var.tmp, value = T)
            act <- grep("_start|_end",
                        tmp.var.tmp,
                        value = T,
                        invert = T)
            if (length(starts) == 1 & length(ends) == 1) {
              for (k in 1:nrow(tmp)) {
                if (!is.na(tmp[[starts]][k]) & !is.na(tmp[[ends]][k])) {
                  sel <- which(
                    as.numeric(rownames.list(mm)) >= as.Date(as.character(tmp[[starts]][k]), format =
                                                               "%Y%m%d") &
                      as.numeric(rownames.list(mm)) < as.Date(as.character(tmp[[ends]][k]), format =
                                                                "%Y%m%d")
                  )
                } else if (!is.na(tmp[[starts]][k]) &
                           is.na(tmp[[ends]][k])) {
                  sel <- which(
                    as.numeric(rownames.list(mm)) >= as.Date(as.character(tmp[[starts]][k]), format =
                                                               "%Y%m%d")
                  )
                } else if (is.na(tmp[[starts]][k]) &
                           !is.na(tmp[[ends]][k])) {
                  sel <- which(
                    as.numeric(rownames.list(mm)) < as.Date(as.character(tmp[[ends]][k]), format =
                                                              "%Y%m%d")
                  )
                } else{
                  sel <- 1:length(rownames.list(mm))
                }
                for (l in 1:length(act))
                  mm[[act[l]]][sel] <- lapply(mm[[act[1]]][sel], function(x)
                    unique(c(na.omit(
                      unlist(x)
                    ), tmp[[act[l]]][k])))
              }
            } else if (length(starts) == 1 & length(ends) == 0) {
              for (k in 1:nrow(tmp)) {
                if (!is.na(tmp[[starts]][k])) {
                  sel <- which(
                    as.numeric(rownames.list(mm)) >= as.Date(as.character(tmp[[starts]][k]), format =
                                                               "%Y%m%d")
                  )
                } else{
                  sel <- 1:length(rownames.list(mm))
                }
                for (l in 1:length(act))
                  mm[[act[l]]][sel] <- lapply(mm[[act[1]]][sel], function(x)
                    unique(c(na.omit(
                      unlist(x)
                    ), tmp[[act[l]]][k])))
                
                #mm[sel,act[l]] <- paste(unique(c(strsplit(mm[sel,act[l]],"\\|"),tmp[[act[l]]][k])),collapse="|")
              }
            } else if (length(starts) == 0 & length(ends) == 1) {
              for (k in 1:nrow(tmp)) {
                if (!is.na(tmp[[ends]][k])) {
                  sel <- which(
                    as.numeric(rownames.list(mm)) >= as.Date(as.character(tmp[[ends]][k]), format =
                                                               "%Y%m%d")
                  )
                } else{
                  sel <- 1:length(rownames.list(mm))
                }
                for (l in 1:length(act))
                  mm[[act[l]]][sel] <- lapply(mm[[act[1]]][sel], function(x)
                    unique(c(na.omit(
                      unlist(x)
                    ), tmp[[act[l]]][k])))
                #mm[sel,act[l]] <- paste(unique(c(strsplit(mm[sel,act[l]],"\\|"),tmp[[act[l]]][k])),collapse="|")
              }
            } else{
              .tmp <<- tmp
              stop("What?! How is that possible!")
            }
            if (is.list(mm)) {
              xsel <- which(sapply(mm, class) == "list")
              if (length(xsel) == 0)
                class(mm) <- "data.frame"
              for (k in xsel)
                mm[[k]] <- type.convert(sapply(mm[[k]], function(z)
                  paste(unique(z[z != ""]), collapse = "|")))
              class(mm) <- "data.frame"
            }
          } else{
            for (l in 1:length(tmp.var.tmp))
              mm[, tmp.var.tmp[l]] <- tmp[[tmp.var.tmp[l]]]
          }
        }
      }
    }
    return(mm)
  })
  return(dg)
})

ascores <- lapply(fullres, function(pat) {
  fscs <- matrix(NA, ncol = length(pat), nrow = nrow(pat[[1]]))
  colnames(fscs) <- rownames(scores)
  rownames(fscs) <- rownames(pat[[1]])
  .pat <<- pat
  for (i in 1:length(pat)) {
    .i <<- i
    tmp <- pat[[i]]
    #tmp <- tmp[,intersect(vars.touse,colnames(tmp))]
    if (any(sapply(tmp, class) == "character")) {
      sel <- which(sapply(tmp, class) == "character")
      for (j in 1:length(sel)) {
        if (colnames(tmp)[sel[j]] == "eq5d_mob") {
          tmp[, sel[j]] <- as.numeric(tmp[, sel[j]])
        }
        if (colnames(tmp)[sel[j]] == "home_childr" |
            colnames(tmp)[sel[j]] == "resid_admission" |
            colnames(tmp)[sel[j]] == "home_adults" |
            colnames(tmp)[sel[j]] == "resid" |
            colnames(tmp)[sel[j]] == "home") {
          tmp[, sel[j]] <- sapply(strsplit(tmp[, sel[j]], "\\|"), function(x)
            max(as.numeric(x)))
        }
        if (colnames(tmp)[sel[j]] == "gec_med_i_drug" |
            colnames(tmp)[sel[j]] == "gec_med_drug" |
            colnames(tmp)[sel[j]] == "med_ih_drug" |
            colnames(tmp)[sel[j]] == "gec_med_drug_oth") {
          tmp[, sel[j]] <- lengths(strsplit(tmp[, sel[j]], "\\|"))
        }
      }
    }
    tmp <- data.frame(lapply(tmp, as.numeric))
    if (rownames(scores)[i] == "nutrition") {
      scs <- rep(NA, nrow(tmp))
      bmi <- tmp$gec_weight / ((tmp$gec_height / 100)^2)
      sel <- which(!is.na(bmi))
      if (length(sel) == 1) {
        bmi <- rep(bmi[sel], nrow(tmp))
      } else if (length(sel) == 0) {
        bmi <- rep(NA, nrow(tmp))
      } else{
        stop("Implement")
      }
      if (all(is.na(bmi))) {
        scs <- rep(NA, nrow(tmp))
      } else{
        normal <- sapply(strsplit(scores[i, 1], "-")[[1]], as.numeric)
        p0.5 <- strsplit(scores[i, 2], ",")[[1]]
        p0.5.lower <- as.numeric(gsub("<", "", grep("<", p0.5, value = T)))
        p0.5.higher <- as.numeric(gsub(">", "", grep(">", p0.5, value =
                                                       T)))
        p1 <- strsplit(scores[i, 3], ",")[[1]]
        p1.lower <- as.numeric(gsub("<", "", grep("<", p1, value = T)))
        p1.higher <- as.numeric(gsub(">", "", grep(">", p1, value = T)))
        bins <- c(0, unique(
          c(p1.lower, p0.5.lower, normal, p0.5.higher, p1.higher)
        ), 1000)
        bins.score <- as.numeric(c(
          colnames(scores)[3],
          colnames(scores)[2],
          colnames(scores)[1],
          colnames(scores)[2],
          colnames(scores)[3]
        ))
        scs <- bins.score[as.numeric(cut(
          bmi,
          bins,
          include.lowest = T,
          ordered_result = T
        ))]
      }
    } else if (rownames(scores)[i] == "social") {
      scs <- apply(tmp, 1, function(z) {
        .z <<- z
        if (is.na(z["resid"]) & !is.na(z["home"])) {
          if (!is.na(z["home_adults"])) {
            if (z["home"] == 2 & (z["home_adults"] > 0 | z["home_childr"] > 0))
              return(0)
            if (z["home"] == 2 &
                z["home_adults"] == 0 & z["home_childr"] == 0)
              return(0.5)
          }
          if (z["home"] == 2)
            return(0.5)
          if (z["home"] == 1)
            return(1)
        } else if (!is.na(z["resid"]) & is.na(z["home"])) {
          if (z["resid"] == 1)
            return(1)
          if (z["resid"] == 2 | z["resid"] == 3)
            return(0.5)
          if (z["resid"] == 4)
            return(1)
        } else if (!is.na(z["resid"]) & !is.na(z["home"])) {
          if (z["resid"] == 4)
            return(1)
          if (z["resid"] == 1 & z["home"] == 2)
            return(0)
          if (z["resid"] == 1 & z["home"] == 1)
            return(1)
        } else{
          return(NA)
        }
        return(NA)
      })
    } else if (rownames(scores)[i] == "biomarker") {
      scs <- rep(NA, nrow(tmp))
    } else{
      dsum <- rowSums(tmp, na.rm = T)
      scs <- rep(NA, length(dsum))
      for (j in 1:ncol(scores)) {
        scs[which(dsum > as.numeric(scores[i, j]))] <- as.numeric(colnames(scores)[j])
      }
    }
    fscs[, rownames(scores)[i]] <- scs
  }
  return(fscs)
})

names(ascores) <- pat_id

death_date <- data.frame(
  export_psn = data.all$fv2_6$export_psn,
  date = as.Date(as.character(data.all$fv2_6$death_date), format = "%Y%m%d")
)
wards <- data.all$eward[, c("export_psn", "gec_ward", "ward_start", "ward_end")]
wards$ward_start <- as.Date(as.character(wards$ward_start), format = "%Y%m%d")
wards$ward_end <- as.Date(as.character(wards$ward_end), format = "%Y%m%d")

mpis <- lapply(1:length(ascores), function(x) {
  .n <<- names(ascores)[x]
  .x <<- x
  mpi <- rowMeans(ascores[[x]], na.rm = T)
  mm <- data.frame(MPI = mpi[subset(ma, export_psn == names(ascores)[x])$datenr], subset(ma, export_psn ==
                                                                                           names(ascores)[x]))
  mm$diff <- as.numeric(mm$datenr) - min(as.numeric(mm$datenr))
  mm$death <- NA
  mm$death_date <- NA
  mm$rehosp <- NA
  mm$rehosp_1 <- NA
  mm$rehosp_2 <- NA
  mm$rehosp_3 <- NA
  mm$rehosp_4 <- NA
  mm$rehosp_5 <- NA
  mm$rehosp_6 <- NA
  mm$rehosp_7 <- NA
  selb <- which(mm$mnpvislabel == "Baseline")
  sel3 <- which(mm$mnpvislabel == "3M Follow-Up")
  sel12 <- which(mm$mnpvislabel == "12M Follow-Up")
  if (length(selb) == 0)
    return(NULL)
  if (length(sel3) == 0 & length(sel12) == 0)
    return(NULL)
  sels <- c(sel3, sel12)
  for (s in sels) {
    tln <- as.numeric(mm[selb, "datenr"]):as.numeric(mm[s, "datenr"])
    tl <- rep(0, length(tln))
    names(tl) <- tln
    wtmp <- subset(wards, export_psn == names(ascores)[x])
    dd <- subset(death_date, export_psn == names(ascores)[x])
    if (nrow(wtmp) > 0) {
      for (l in 1:nrow(wtmp)) {
        if (!is.na(wtmp[l, "ward_start"]) & !is.na(wtmp[l, "ward_end"]))
          tl[as.character(as.numeric(wtmp[l, "ward_start"]):as.numeric(wtmp[l, "ward_end"]))] <- wtmp[l, "gec_ward"]
        if (is.na(wtmp[l, "ward_start"]) |
            is.na(wtmp[l, "ward_end"]))
          tl[1:length(tl)] <- NA
      }
    }
    for (j in 1:7) {
      lx <- length(which(tl == j))
      mm[s, paste0("rehosp_", j)] <- lx
    }
    mm[s, "rehosp"] <- length(which(tl != 0))
    if (nrow(dd) > 1)
      stop("WTF")
    if (!is.na(dd[1, 2])) {
      tl[as.character(as.numeric(dd[1, 2]))] <- -1
      mm$death_date <- dd[1, 2]
    }
    mm[s, "death"] <- F
    if (all(!is.na(tl)))
      if (any(tl == -1))
        mm[s, "death"] <- T
    
  }
  
  #tm <- as.numeric(rownames(mm))
  #rownames(mm) <- tm-min(tm)
  return(mm)
})
names(mpis) <- names(ascores)


mpir <- sapply(mpis, function(x) {
  if (!is.null(x)) {
    selb <- which(x$mnpvislabel == "Baseline")
    sel3 <- which(x$mnpvislabel == "3M Follow-Up")
    sel12 <- which(x$mnpvislabel == "12M Follow-Up")
    baseline_mpi <- x[selb, "MPI"]
    if (length(sel3) != 0) {
      m3data <- c(
        MPI = x[sel3, "MPI"],
        rehosp = x[sel3, "rehosp"],
        rehosp_1_3 = sum(x[sel3, c("rehosp_1", "rehosp_2", "rehosp_3")], na.rm =
                           T),
        rehosp_4 = x[sel3, "rehosp_4"],
        death = x[sel3, "death"]
      )
    } else{
      m3data <- c(
        MPI = NA,
        rehosp = NA,
        rehosp_1_3 = NA,
        rehosp_4 = NA,
        death = NA
      )
    }
    names(m3data) <- paste(names(m3data), "_3", sep = "")
    if (length(sel12) != 0) {
      m12data <- c(
        MPI = x[sel12, "MPI"],
        rehosp = x[sel12, "rehosp"],
        rehosp_1_3 = sum(x[sel12, c("rehosp_1", "rehosp_2", "rehosp_3")], na.rm =
                           T),
        rehosp_4 = x[sel12, "rehosp_4"],
        death = x[sel12, "death"]
      )
    } else{
      m12data <- c(
        MPI = NA,
        rehosp = NA,
        rehosp_1_3 = NA,
        rehosp_4 = NA,
        death = NA
      )
    }
    names(m12data) <- paste(names(m12data), "_12", sep = "")
    return(data.frame(c(MPI = baseline_mpi, m3data, m12data)))
  }
})
mpirs <- t(do.call(cbind, mpir[lengths(mpir) != 0]))
rownames(mpirs) <- names(mpis)[lengths(mpir) != 0]