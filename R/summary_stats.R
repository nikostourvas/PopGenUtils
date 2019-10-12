## functions

###### basic statistics and their SEs
table_out <- function(obj, variable, name){

        means <- colMeans(variable, na.rm=T)
        out <- c(means, mean(variable, na.rm = T))

        out <- as.data.frame(out)
        Pops <- c(popNames(obj), "Total")
        rownames(out) <- Pops
        colnames(out) <- name

        sem_out <- apply(variable, 2, function(x) sd(x) / sqrt(length(x)))
        ## 2 means work along columns
        sem_out_mean <- sd(variable) / sqrt(length(variable))

        sem_out <- as.data.frame(c(sem_out, sem_out_mean))
        rownames(sem_out) <- Pops
        colnames(sem_out) <- paste("SE", name, sep = "_")

        table_out <- cbind(out, sem_out)

        return(table_out)
}



###### basic statistics reported from poppr and their SEs
poppr2hierfstat_out <- function(obj, variable){

        obj_list <- seppop(obj)

        stats_poppr <- list()
        for(i in 1: length(obj_list)){
                stats_poppr[[i]] <- locus_table(obj_list[[i]], information = F)
        }

        table_out <- list()
        for(i in 1:length(obj_list))
                table_out[[i]] <- stats_poppr[[i]][-nrow(stats_poppr[[1]]), variable]

        table_out <- as.matrix(as.data.frame(table_out))
        colnames(table_out) <- popNames(obj)

        return(table_out)
}




genalex.table <- function(obj){
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



summary_df <- cbind(N, na[,1], ne_Hs[,1], Ho[,1], uHe[,1], Fis[,1])
rownames(summary_df) <- c(popNames(obj), "Total")
colnames(summary_df) <- c("N", "na", "ne", "Ho", "uHe", "Fis")
summary_df <- round(as.data.frame(summary_df), digits = 3)

return(summary_df)

}

# knitr::kable(summary_df, "html", digits = 3) %>%
#         kable_styling(bootstrap_options = "striped", full_width = F)

