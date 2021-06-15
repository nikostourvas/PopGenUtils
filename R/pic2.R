pic_calc2 <- function(obj, ploidy = 2){

  obj_hier <- genind2hierfstat(obj)
  # hierfstat adds a fake pop and ind when there is only one pop
  # in the data set. Remove this extra line
  if (rownames(obj_hier)[length(rownames(obj_hier))] == "dumind") {
    obj_hier <- obj_hier[-length(rownames(obj_hier)), ]
  }

  write.struct(obj_hier,
               fname="str.txt")

  str <- read.table("str.txt")

  # Change the pop column to integer as is required by structure
  str$V1 <- rep(as.numeric(obj@pop), each = ploidy)

  write.table(str, file="str.txt",
              row.names = FALSE,
              col.names = FALSE)

  obj_genambig <- read.Structure(infile = "str.txt",
                                 ploidy = ploidy,
                                 sep = " ",
                                 markernames = FALSE,
                                 labels = FALSE,
                                 extrarows = 0,
                                 extracols = 1,
                                 popinfocol = 1,
                                 getexcols = FALSE,
                                 ploidyoutput = "one")

  invisible(file.remove("str.txt"))

  # calculate pic
  freq_obj <- simpleFreq(obj_genambig)
  pic_obj <- PIC(freq_obj, overall = TRUE)

  pic_obj <- as.data.frame(pic_obj)
  pic_obj_mean <- rowMeans(pic_obj)

  pic_obj$mean <- pic_obj_mean

  pic_obj <- as.data.frame(t(as.matrix(pic_obj)))

  return(pic_obj)
}
