genemapper2polysat <- function(folder_name){

library(tidyverse)
  
temp <- list.files(path = folder_name,
                   pattern = "*.txt",
                   full.names = TRUE)
df <- lapply(temp, read.csv, sep = "\t",
             na.strings = "")

# sort data.frame based on Sample.File column
df <- lapply(df, arrange, Sample.File)
df <- lapply(df, select, 
             Sample.File, Marker, Allele.1, Allele.2)

# merge data.frames
df_final <- do.call(rbind, df)

# remove ".fsa" from ind names
df_final$Sample.File <- gsub(pattern = ".fsa",
                         replacement = "",
                         x = df_final$Sample.File)

colnames(df_final) <- c("Sample.Name",
                        "Marker",
                        "Allele.1",
                        "Allele.2")

df_final$Allele.1[is.na(df_final$Allele.1)] <- -9
# df_final <- subset(df_final, Allele.1 != -9)

df_final$Allele.2[is.na(df_final$Allele.2)] <- ""

return(df_final)
}