#' Experimental function: Calculate the standard deviation of one or more
#' parameters by bootstraping
#'
#' @param x a vector or matrix (with each column as one parameter)
#' @param nboot number of bootstraps; default: 1000
#'
#' @return the sd of a vector, or a data.frame with the sd of each parameter
#' in each row
#' @export
#'
#' @author Nikolaos Tourvas
#' @import boot
boot.param.sd <- function(x, nboot=1000){

  if(is.vector(x) == TRUE){
        boot.mean <- function(vec, i){
          return(mean(vec[i]))
        }
    res <- boot(data = x,
                statistic = boot.mean,
                R = nboot)
    return(sd(res$t))

  } else {
        boot.mean.df <- function(df, i){
          df2 <- df[i,]
          return(colMeans(df2))
        }
    res <- boot(data = x,
                statistic = boot.mean.df,
                R=nboot)
    res <- apply(res$t, 2, sd)
    res <- data.frame(sd = res)
    rownames(res) <- colnames(x)
    return(res)
  }
}
