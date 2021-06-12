#' Experimental function: Produce a summary statistics table
#'
#' @param obj a genind object
#'
#' @return a data.frame with summary statistics
#' @export
#'
#' @author Nikolaos Tourvas
#' @description This function allows the user to produce to summary table
#' similar to the one produced by the software GenAlEx for diploid codominant
#' data. WARNING: This function is still experimental.
#'
#' @import poppr
#' @importFrom hierfstat basic.stats
#' @import reshape2
#' @import ShannonGen
summary_by_pop <- function(obj) {
  ## N
  N_by_locus <- basic.stats(obj)[["n.ind.samp"]]
  obj_list <- seppop(obj)
  N <- list()
  for (i in 1:length(obj_list)) {
    N[[i]] <- length(obj_list[[i]]@pop)
  }
  N <- melt(N)
  N <- c(N[, 1], mean(N[, 1]))

  ## na
  na_by_locus <- poppr2hierfstat_out(obj, "allele")
  na <- table_out(obj, na_by_locus, "na")

  ## ne
  ne_by_locus_Hs <- 1 / (1 - (basic.stats(obj)[["Hs"]]))
  ne_Hs <- table_out(obj, ne_by_locus_Hs, "ne")

  ## ## ne
  ## ne_by_locus_He <- 1 / (1 - (basic.stats(obj)[["Hs"]]))
  ## ne_Hs <- table_out(obj, ne_by_locus, "ne")

  ## uHe
  uHe_by_locus <- poppr2hierfstat_out(obj, "Hexp")
  uHe <- table_out(obj, uHe_by_locus, "uHe")

  ## Ho
  Ho_by_locus <- basic.stats(obj)[["Ho"]]
  Ho <- table_out(obj, Ho_by_locus, "Ho")

  ## Shannon's I
  I_by_locus <- ShannonGen(obj, estimator = "sh")
  I <- table_out(obj, as.matrix(I_by_locus$Shannon_1949), "I")

  ## Shannon's I Zahl
  I_z_by_locus <- ShannonGen(obj, estimator = "z")
  I_z <- table_out(obj, as.matrix(I_z_by_locus$Zahl_1977), "I_z")

  ## Fis
  Fis_by_locus <- basic.stats(obj)[["Fis"]]
  Fis <- table_out(obj, Fis_by_locus, "Fis")

  # save only the values - not their SEs
  summary_df <- cbind(N, na, ne_Hs, I, I_z, Ho, uHe, Fis)
  rownames(summary_df) <- c(popNames(obj), "Mean")
  colnames(summary_df) <- c("N", "na", "na_SE", "ne", "ne_SE",
                            "I", "I_SE", "I_z", "I_z_SE",
                            "Ho", "Ho_SE", "uHe", "uHe_SE", "Fis", "Fis_SE")
  summary_df <- round(as.data.frame(summary_df), digits = 3)

  return(summary_df)
}



#' Experimental function: Produce summary statistics by locus
#'
#' @param obj
#'
#' @return
#' @export
#'
summary_by_loc <- function(obj) {
  obj_list <- seppop(obj)

  ## N
  basic_res <- basic.stats(obj)
  N_by_locus <- basic_res$n.ind.samp

  ## na
  na_by_locus <- poppr2hierfstat_out(obj, "allele")
  # na <- table_out(obj, na_by_locus, "na")

  ## ne
  ne_by_locus_Hs <- 1 / (1 - (basic_res$perloc$Hs))
  # ne_Hs <- table_out(obj, ne_by_locus_Hs, "ne")

  ## ## ne
  ## ne_by_locus_He <- 1 / (1 - (basic.stats(obj)[["Hs"]]))
  ## ne_Hs <- table_out(obj, ne_by_locus, "ne")

  ## uHe
  uHe_by_locus <- poppr2hierfstat_out(obj, "Hexp")
  # uHe <- table_out(obj, uHe_by_locus, "uHe")

  ## Ho
  Ho_by_locus <- basic_res$perloc$Ho
  # Ho <- table_out(obj, Ho_by_locus, "Ho")

  ## Shannon's I
  I_by_locus <- ShannonGen(obj, estimator = "sh")

  ## Shannon's I Zahl
  I_z_by_locus <- ShannonGen(obj, estimator = "z")

  ## Fis
  Fis_by_locus <- basic_res$perloc$Fis
  # Fis <- table_out(obj, Fis_by_locus, "Fis") ## better use boot.ppfis

  ## PIC
  PIC <- pic_calc(obj)

  ## PID
  pid <- pid_calc(obj)
  pidsibs <- pid$pidsibs_by_loc
  pid <- pid$pid_by_loc


  summary_by_locus <- data.frame(
    cbind(rowMeans(N_by_locus, na.rm = T),
          rowMeans(na_by_locus, na.rm = T),
          ne_by_locus_Hs,
          rowMeans(I_by_locus$Shannon_1949, na.rm = T),
          rowMeans(I_z_by_locus$Zahl_1977, na.rm = T),
          Ho_by_locus,
          rowMeans(uHe_by_locus, na.rm = T),
          Fis_by_locus,
          PIC[-nrow(PIC), ncol(PIC)],
          pid,
          pidsibs
    )
  )

  colnames(summary_by_locus) <- c(
    "N", "na", "ne", "I", "I_z", "Ho", "uHe", "Fis", "PIC", "PID", "PIDsibs"
  )

  summary_by_locus <- rbind(summary_by_locus, colMeans(summary_by_locus))
  rownames(summary_by_locus) <-
    c(rownames(summary_by_locus)[-nrow(summary_by_locus)], "Mean")

  return(summary_by_locus)
}
