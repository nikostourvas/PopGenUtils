#' Convert a genind object to a genambig object
#'
#' @param obj a genind object
#'
#' @return a genambig object
#' @export
#'
#' @author Nikolaos Tourvas
#' @describeIn This function allows the user to convert a genind object to a
#' genambig object for analyses using the package polysat.
#'
#' @import tidyr
#' @import polysat
#' @import pegas
genind2genambig <- function(obj){
  obj_df <- as.loci(obj)
  obj_df <- obj_df[,-1]
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
  PopInfo(obj_genambig) <- rep(1, nInd(obj))
  Ploidies(obj_genambig) <- rep(2, nInd(obj))

  return(obj_genambig)
}
