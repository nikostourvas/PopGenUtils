#' Generate microsatellite input file for TreeMix
#'
#' @param x a genind object
#' @param export_name Name of the exported input file. Note that it should include the path to the directory you would like to export it to.
#'
#' @return the input file for TreeMix as is required for microsatellite data
#' @export
#'
#' @author Nikolaos Tourvas
#' @import hierfstat
#' @import adegenet
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr %>%
genind2treemix <- function(x, export_name = "."){
  # calculate allele freqs
  counts <- allele.count(genind2hierfstat(x))

  # make dataframe
  counts <- lapply(counts, function(x){
                  as.data.frame(x)
  })

  # convert allele names from factor to integer
  for(i in 1:length(locNames(x))){
      counts[[i]]$x <- as.numeric(as.character(counts[[i]]$x))
  }

  # Mean length - Col 1
  meanlength <- function(x){
    x %>%
      group_by(Var2) %>%
      mutate(mult = x * Freq / sum(Freq)) %>%
      summarise(avg = sum(mult))
  }
  col1 <- lapply(counts, meanlength) %>%
          lapply(as.data.frame)

  col1 <- bind_rows(col1)
  col1$Locus <- rep(locNames(x), each = length(popNames(x)))
  col1 <- arrange(col1, Var2, Locus)
  col1$avg <- round(col1$avg, 1)

  # Variance in length - Col 2
  varinlength <- function(x){
    x %>%
      group_by(Var2) %>%
      summarise(varinlength = var(x))
  }

    # remove alleles with 0 occurences in each pop
  counts <- lapply(counts, function(x){
    dplyr::filter(x, Freq > 0)
  })

  col2 <- lapply(counts, varinlength) %>%
          lapply(as.data.frame)

  col2 <- bind_rows(col2)
  col2$Locus <- rep(locNames(x), each = length(popNames(x)))
  col2 <- arrange(col2, Var2, Locus)
  col2$varinlength <- round(col2$varinlength, 1)

  # Number of haplotypes - Col 3
  haplotypes <- genind2df(x)
  nhaplotypes <- function(x){
    haplotypes %>%
      group_by(pop) %>%
      summarise(across(where(is.character), ~ length(unique(na.omit(.x)))))
  }

  col3 <- nhaplotypes(haplotypes) %>%
    pivot_longer(cols = -pop,
                 names_to = "Locus",
                 values_to = "nhapl")

  treem_table <- list(col1, col2, col3)

# Put them all together!
  treem_table <- bind_cols(treem_table[[1]][ ,-3],
                           treem_table[[2]]$varinlength,
                           treem_table[[3]]$nhapl)
  colnames(treem_table) <- c("pop", "mean", "var", "nhapl")

  treem_table <- treem_table %>%
    group_split(pop)

  treem_table <- lapply(treem_table, function(x){
                        x[ ,-1]
  })

  treem_table <- lapply(treem_table, function(x){
    unite(x, 'concat', colnames(x), sep = ",")
  })

  treem_table <- bind_cols(treem_table)

  treem_table <- treem_table %>%
    unite('concat', colnames(treem_table), sep = " ")

  popnames <- paste(popNames(x), sep = " ",
                    collapse = " ")

  treem_table <- rbind(popnames, treem_table)

  write.table(x = treem_table,
              file = export_name, sep = ",",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
