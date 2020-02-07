MDeconv <- function(Ymat, SelMarker,
                    alpha = NULL,
                    sigma = NULL,
                    epsilon = 1e-3,
                    maxIter = 1000,
                    verbose = TRUE) {

    K = length(SelMarker)

    if(is.null(names(SelMarker))) {
        message("Input marker list has empty cell type names.")
        message(paste0("Assign cell type names to CellType1 to Celltype", K))
        names(SelMarker) <- paste0("CellType", 1:K)
    }

    mres1 <- match(unlist(SelMarker), rownames(Ymat))
    if(all(is.na(mres1))) {
        stop("None of the marker genes match row names of input Ymat!")
    } else if(any(is.na(mres1))){
        nna <- sum(is.na(mres1))
        message(paste0("Remove ", nna, " markers that do not match Ymat row names."))
        reduceYmat <- Ymat[na.omit(mres1),]

        for(k in 1:K) {
            mres2 <- match(SelMarker[[k]], rownames(Ymat))
            SelMarker[[k]] <- SelMarker[[k]][!is.na(mres2)]
        }
    } else {
        reduceYmat <- Ymat[mres1,]
    }

    if(is.null(alpha)) {
        res1 <- MDeconvWithoutPrior(Ymat = reduceYmat,
                                    SelMarker = SelMarker,
                                    epsilon = epsilon,
                                    verbose = verbose)
    } else {

        priorres <- GetPrior(alpha = alpha, sigma = sigma)
        alpha_prior <- priorres$alpha_prior
        sigma_prior <- priorres$sigma_prior

        if(length(alpha_prior) != K | length(sigma_prior) != K) {
            stop(paste0("Alpha and sigma need to be of length ", K, "!"))
        } else {
            res1 <- MDeconvWithPrior(Ymat = reduceYmat,
                                        SelMarker = SelMarker,
                                        alpha = alpha_prior,
                                        sigma = sigma_prior,
                                        epsilon = epsilon,
                                        verbose = verbose)

        }
    }

    return(res1)
}
