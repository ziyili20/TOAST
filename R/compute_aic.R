## compute aic to decide cell type number
compute_aic <- function(estProp,Y.raw){

     K = ncol(estProp)
     Nsample = dim(Y.raw)[2]
     idx = apply(estProp,2,function(x) sum(x) == 0)
     estProp[,idx] = matrix(runif(Nsample*sum(idx),0.0001,0.0002),Nsample,sum(idx))
     estProf <- t(mycsfit(estProp, t(Y.raw))$ghat)
     tmpmat <- estProf %*% t(estProp)
     rss = norm(Y.raw-tmpmat,type = "F")^2
     nSample = ncol(Y.raw) * nrow(Y.raw)
     nParam = K*(nrow(Y.raw)+ncol(Y.raw))
     aic = nSample*log(rss/nSample)+ 2*nParam + (2*nParam*(nParam+1))/(nSample-nParam-1)

     return(aic)
}
