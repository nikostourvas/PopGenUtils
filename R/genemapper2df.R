genemapper2df <- function(folder_name,
                          export_name = "genemapper_export.csv",
                          export_option = "csv") {

# library(tidyverse)

# function
# if Allele.2 spot is empty but Allele.1 is not,
# then copy Allele.1 to Allele.2
replace.empty <- function(x) {
  x %>%
    mutate(Allele.2.edit =
             case_when(
               is.na(Allele.2) ~ Allele.1,
               is.numeric(Allele.2) ~ Allele.2
             )
    )
}

temp <- list.files(path = folder_name,
                   pattern = "*.txt",
                   full.names = TRUE)
df <- lapply(temp, read.csv, sep = "\t",
               na.strings = "")

# Check if all alleles are coded as numbers
for (i in 1:length(df)){
  if(!is.numeric(df[[i]]$Allele.1)) {
    stop(call. = FALSE,
         "An allele is not coded as a number in the txt file of marker ",
         unique(df[[i]]$Marker))
  } else if (!is.numeric(df[[i]]$Allele.2)) {
    stop(call. = FALSE,
         "An allele is not coded as a number in the txt file of marker ",
         unique(df[[i]]$Marker))
  }
}

## TODO: make check that marker names do not have spaces

# sort data.frame based on Sample.File column
df <- lapply(df, arrange, Sample.File)
df <- lapply(df, select,
               Sample.File, Marker, Allele.1, Allele.2)

df_edit <- lapply(df, replace.empty)
df_edit <- lapply(df_edit, select,
               Allele.1, Allele.2.edit)

# vector of marker names
markers <- list()
for (i in 1:length(df)){
  markers[[i]] <- unique(df[[i]]$Marker)
}
markers <- unlist(markers)
markers <- rep(markers, each=2)

# merge multiple data.frames inside the list
# into one data.frame
df_final <- do.call(cbind, df_edit)
colnames(df_final) <- markers # add marker names

# add ind names column
ind_names <- df[[1]]$Sample.File
ind_names <- gsub(pattern = ".fsa",
                  replacement = "",
                  x = ind_names)

# df_final <- df_final %>%
#   add_column(Genotype = ind_names,
#              .before = 1)
rownames(df_final) <- ind_names
df_final[is.na(df_final)] <- 0

if (export_option == "csv") {
  write.csv(x = df_final, file = export_name)
} else if (export_option == "R") {
  return(df_final)
} else {
  stop(call. = FALSE,
  'The export_option parameter should take one of the following values
        "csv" or "R" ')
}
}
