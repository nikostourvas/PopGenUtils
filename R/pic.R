pic_calc <- function(obj){
  # GENIND 2 POLYSAT
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

  test <- read.GeneMapper("polysat_input.txt")
  invisible(file.remove("polysat_input.txt"))

  # Create genambig object
  PopInfo(test) <- rep(1, nInd(obj))
  Ploidies(test) <- rep(2, nInd(obj))

  freq_obj <- simpleFreq(test)
  pic_obj <- PIC(freq_obj, overall = F)

  pic_obj <- as.data.frame(pic_obj)
  pic_obj_mean <- pic_obj %>%
    rowMeans()

  pic_obj$mean <- pic_obj_mean

  pic_obj <- as.data.frame(t(as.matrix(pic_obj)))

  return(pic_obj)
}
