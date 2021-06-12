#' Polymorphic Information Content calculation for genind objects
#'
#' @param obj a codominant genind object
#' @param ploidy The ploidy of the data
#'
#' @return a data.frame with PIC values (per locus & mean) for each population
#' and overall
#' @export
#'
#' @author Nikolaos Tourvas
#' @description This function is a simple wrapper of the PIC function from
#' polysat package. It allows to quickly calculate Polymorphic Information
#' Content or PIC from a genind object. Check PIC help page under polysat
#' package for references.
#' @import poppr
#' @import pegas
#' @import polysat
pic_calc <- function(obj, ploidy = 2){
  # create genambig objects
  # without the need to import poppr
  obj_df <- as.loci(obj)
  popinfo <- as.numeric(obj_df$population)
  obj_df <- obj_df[ ,-1]
  obj_df$Sample.Name <- rownames(obj_df)

  # make it tidy for polysat
  obj_df_long <- suppressWarnings(
    gather(obj_df, Marker, Alleles, -Sample.Name) )

  obj_df_long <- separate(obj_df_long, col=Alleles, sep = "/",
                          into = c("Allele.1", "Allele.2"))

  obj_df_long$Allele.1[is.na(obj_df_long$Allele.1)] <- -9

  write.table(obj_df_long, file = "polysat_input.txt",
              row.names = F, sep = "\t", na = "")

  obj_genambig <- read.GeneMapper("polysat_input.txt")
  invisible(file.remove("polysat_input.txt"))

  # Create genambig object
  PopInfo(obj_genambig) <- popinfo
  Ploidies(obj_genambig) <- rep(ploidy, nInd(obj))

  # as.genambig is imported from poppr
  # obj_genambig <- as.genambig(obj)
  # obj_genambig@Ploidies@pld[is.na(obj_genambig@Ploidies@pld)] <- ploidy

  # calculate pic
  freq_obj <- simpleFreq(obj_genambig)
  pic_obj <- PIC(freq_obj, overall = TRUE)

  pic_obj <- as.data.frame(pic_obj)
  pic_obj_mean <- rowMeans(pic_obj)

  pic_obj$mean <- pic_obj_mean

  pic_obj <- as.data.frame(t(as.matrix(pic_obj)))

  return(pic_obj)
}
