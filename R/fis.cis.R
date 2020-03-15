fis.cis <- function(obj, nboot){
        res <- boot.ppfis(obj, nboot = nboot)
        res <- data.frame(sapply(res[[2]], c))
        return(res)
}