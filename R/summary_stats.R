## functions

###### basic statistics and their SEs
table_out <- function(obj, variable, name){

        means <- colMeans(variable, na.rm = TRUE)
        out <- c(means, mean(variable, na.rm = TRUE))

        out <- as.data.frame(out)
        Pops <- c(popNames(obj), "Mean")
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
                stats_poppr[[i]] <- locus_table(obj_list[[i]], information = FALSE)
        }

        table_out <- list()
        for(i in 1:length(obj_list))
                table_out[[i]] <- stats_poppr[[i]][-nrow(stats_poppr[[1]]), variable]

        table_out <- as.matrix(as.data.frame(table_out))
        colnames(table_out) <- popNames(obj)

        return(table_out)
}



