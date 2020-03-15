#' Produce a summary statistics table
#'
#' @param obj a genind object
#'
#' @return a data.frame with summary statistics
#' @export
#'
#' @author Nikolaos Tourvas
#' @describeIn This function allows the user to produce to summary table
#' similar to the one produced by the software GenAlEx.
#'
#' @examples
#' @import poppr
#' @import hierfstat
#' @import reshape2
summary_by_pop <- function(obj){
  ## N
  N_by_locus <- basic.stats(obj)[["n.ind.samp"]]
  obj_list <- seppop(obj)
  N <- list()
  for(i in 1:length(obj_list)){
    N[[i]] <- length(obj_list[[i]]@pop)
  }
  N <- melt(N)
  N <- c(N[,1], sum(N[,1]))

  ## na
  na_by_locus <- poppr2hierfstat_out(obj, "allele")
  na <- table_out(obj, na_by_locus, "na")


  ## uHe
  uHe_by_locus <- poppr2hierfstat_out(obj, "Hexp")
  uHe <- table_out(obj, uHe_by_locus, "uHe")

  ## Ho
  Ho_by_locus <- basic.stats(obj)[["Ho"]]
  Ho <- table_out(obj, Ho_by_locus, "Ho")

  ## ne
  ne_by_locus_Hs <- 1 / (1 - (basic.stats(obj)[["Hs"]]))
  ne_Hs <- table_out(obj, ne_by_locus_Hs, "ne")

  ## ## ne
  ## ne_by_locus_He <- 1 / (1 - (basic.stats(obj)[["Hs"]]))
  ## ne_Hs <- table_out(obj, ne_by_locus, "ne")

  ## Fis
  Fis_by_locus <- basic.stats(obj)[["Fis"]]
  Fis <- table_out(obj, Fis_by_locus, "Fis") ## better use boot.ppfis

  # save only the values - not their SEs
  summary_df <- cbind(N, na[,1], ne_Hs[,1], Ho[,1], uHe[,1], Fis[,1])
  rownames(summary_df) <- c(popNames(obj), "Total")
  colnames(summary_df) <- c("N", "na", "ne", "Ho", "uHe", "Fis")
  summary_df <- round(as.data.frame(summary_df), digits = 3)

  return(summary_df)

}

# knitr::kable(summary_df, "html", digits = 3) %>%
#         kable_styling(bootstrap_options = "striped", full_width = F)




#' Produce summary statistics by locus
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
summary_by_loc <- function(obj){
  obj_list <- seppop(obj)
  #create workaround for hierfstat
  # if only one pop is present, stats are calculated as if there were
  # two identical pops
  if (length(obj_list) == 1){
    obj2 <- obj

    pop_name <- obj@pop
    pop_name <- as.factor(str_replace(pop_name, "^", "x") )
    obj2@pop <- pop_name

    data_set <- repool(obj2, obj)
  } else {
    data_set <- obj
  }

  ## N
  N_by_locus <- basic.stats(data_set)[["n.ind.samp"]]
  data_set_list <- seppop(data_set)
  N <- list()
  for(i in 1:length(data_set_list)){
    N[[i]] <- length(data_set_list[[i]]@pop)
  }
  N <- melt(N)
  N <- c(N[,1], sum(N[,1]))

  ## na
  na_by_locus <- poppr2hierfstat_out(data_set, "allele")
  na <- table_out(data_set, na_by_locus, "na")


  ## uHe
  uHe_by_locus <- poppr2hierfstat_out(data_set, "Hexp")
  uHe <- table_out(data_set, uHe_by_locus, "uHe")

  ## Ho
  Ho_by_locus <- basic.stats(data_set)[["Ho"]]
  Ho <- table_out(data_set, Ho_by_locus, "Ho")

  ## ne
  ne_by_locus_Hs <- 1 / (1 - (basic.stats(data_set)[["Hs"]]))
  ne_Hs <- table_out(data_set, ne_by_locus_Hs, "ne")

  ## ## ne
  ## ne_by_locus_He <- 1 / (1 - (basic.stats(data_set)[["Hs"]]))
  ## ne_Hs <- table_out(data_set, ne_by_locus, "ne")

  ## Fis
  Fis_by_locus <- basic.stats(data_set)[["Fis"]]
  Fis <- table_out(data_set, Fis_by_locus, "Fis") ## better use boot.ppfis

  ## PIC
  PIC <- pic_calc(data_set)

  ## PID
  pid <- pid_calc(obj)
  pidsibs <- pid$pidsibs_by_loc
  pid <- pid$pid_by_loc


  summary_by_locus <- data.frame(
    cbind(rowMeans(N_by_locus, na.rm=T),
          rowMeans(na_by_locus, na.rm=T),
          rowMeans(ne_by_locus_Hs, na.rm=T),
          rowMeans(Ho_by_locus, na.rm=T),
          rowMeans(uHe_by_locus, na.rm=T),
          rowMeans(Fis_by_locus, na.rm=T),
          PIC[-nrow(PIC), ],
          pid,
          pidsibs)
  )

  colnames(summary_by_locus) <- c("N", "na", "ne", "Ho", "uHe", "Fis", "PIC",
                                  "PI", "PIsibs")

  summary_by_locus <- rbind(summary_by_locus, colMeans(summary_by_locus))
  rownames(summary_by_locus) <-
    c(rownames(summary_by_locus)[-nrow(summary_by_locus)], "Mean")

  return(summary_by_locus)
}
