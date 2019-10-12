### bootstraping a df of a parameter to calc sd (bootstrap over loci)
boot.sd.param <- function(x, nboot){
n <- nrow(x)
bootstrapsample <- list()
for(i in 1:ncol(x)){
        tmpdata = sample(x,n*nboot, replace=TRUE)
        bootstrapsample[[i]] = matrix(tmpdata, nrow=n, ncol=nboot)
}

bsmeans <- list()
for(i in 1:ncol(x)){
        bsmeans[[i]] <- colMeans(bootstrapsample[[i]])
}

# compute sd
sdev <- list()
for(i in 1:ncol(x)){
        sdev[[i]] <- sd(bsmeans[[i]])
}

sdev <- do.call(rbind.data.frame, sdev)
colnames(sdev) <- "sd"
rownames(sdev) <- colnames(x)
return(sdev)

}

# compute CIs from bootstraping
# compute d*
# deltastar <- list()
# d <- list()
# for(i in 1:ncol(x)){
#         deltastar[[i]] <- bsmeans[[i]] - mean(x[[i]])
#         d[[i]] <- quantile(deltastar[[i]], c(0.025, 0.975))
# }

