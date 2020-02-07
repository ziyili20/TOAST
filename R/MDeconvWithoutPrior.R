MDeconvWithoutPrior <- function(Ymat, SelMarker, 
                    epsilon = 1e-3,
                    maxIter = 1000,
                    verbose = TRUE) {
    
    message("Deconvolution without prior information.")
    
    K <- length(SelMarker)
    N <- ncol(Ymat)
    M <- length(unlist(SelMarker))
    
    numMarker <- list()
    for(k in 1:K) {
        numMarker[[k]] <- match(SelMarker[[k]], rownames(Ymat))
    }
    
    Wmat <- matrix(0, M, K)
    Hmat <- matrix(0, K, N)
    
    niter = 1
    diff = 1
    
    ## Initialize Wmat
    for (k in 1:K) {
        Wmat[numMarker[[k]], k] <- rowMeans(Ymat[numMarker[[k]], ])*K
    }
    
    while(diff > epsilon) {
        if(verbose) {
            message("Iteration = ", niter)
        }
        ## update Hmat
        Horig <- Hmat
        Hmat <- LMsolve(Ymat, Wmat, Prop_constraint = TRUE)
        if(is.infinite(max(Hmat)) | is.infinite(min(Hmat)) | is.na(sum(Hmat))) {
            message("H has infinite values!")
            break()
        }
        ## update Wmat
        for(k in 1:K) {
            idx = (Hmat[k,]!=0)
            Wmat[numMarker[[k]],k] <- rowMeans(Ymat[numMarker[[k]],idx])/mean(Hmat[k,idx])
        }
        Wmat[Wmat<0] <- 0
        if(is.infinite(max(Wmat)) | is.infinite(min(Wmat)) | is.na(sum(Wmat))) {
            message("W has infinite values!")
            break()
        }
        
        diff <- max(abs(Horig - Hmat))
        niter <- niter + 1
        if(niter == maxIter) {
            message("Reach maximum number of iterations.")
            break
        }
    }
    
    rownames(Hmat) <- names(SelMarker)
    colnames(Hmat) <- colnames(Ymat)
    
    return(list(H = Hmat,
                W = Wmat))
}