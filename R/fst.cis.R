fst.cis <- function(obj, nboot){
        res <- boot.ppfst(obj, nboot = nboot)

# data.frame from lower tri matrix
ind <- which(upper.tri(res[["ll"]], diag = TRUE), arr.ind = TRUE)
nn <- dimnames(res[["ll"]])
fst_ll <- data.frame(pop1 = nn[[1]][ind[, 1]],
           pop2 = nn[[2]][ind[, 2]],
           ll = res[["ll"]][ind])

# data.frame from upper tri matrix
ind <- which(upper.tri(res[["ul"]], diag = TRUE), arr.ind = TRUE)
nn <- dimnames(res[["ul"]])
fst_ul <- data.frame(pop1 = nn[[1]][ind[, 1]],
           pop2 = nn[[2]][ind[, 2]],
           ul = res[["ul"]][ind])

fst_cis <- na.omit(
        cbind(fst_ll, ul=fst_ul[,3])
)
return(fst_cis)
}