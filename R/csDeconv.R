csDeconv <- function(Y_raw,
                    K,
                    FUN = RefFreeCellMix_wrapper,
                    nMarker = 1000,
                    InitMarker = NULL,
                    TotalIter = 30) {
    # Y_raw is the high-throughput measurement from complex tissues
    #      (rows for features and columns for samples);
    #      or a SummarizedExperiment object.
    # K is the pre-specified number of pure cell types
    # FUN is the reference-free deconvolution function,
    #    this function should take Y_raw and K,
    #    and the return values should be a N by K proportion matrix.
    #    N is the number of samples and K is the number of cell types.
    # nMarker is the number of marker used in the deconvolution
    # InitMarker is the initial marker used in the deconvolution,
    #          if not specified, the top variable features will be used
    # TotalIter is the total number of iterations specified

    if (is(Y_raw, "SummarizedExperiment")) {
         se <- Y_raw
         Y_raw <- assays(se)$counts
    } else if (!is(Y_raw, "matrix")) {
         stop("Y_raw should be a matrix or a SummarizedExperiment object!")
    }

    if (is.null(rownames(Y_raw))) {
        row.names(Y_raw) <- seq(nrow(Y_raw))
    }
    if (is.null(InitMarker)) {
        InitMarker <- findRefinx(Y_raw, nmarker = nMarker)
    } else {
        if (sum(!(InitMarker %in% rownames(Y_raw))) > 0) {
                stop("Discrepancy between
                    InitMarker and the row names of Y_raw!")
        }
    }

    allProp <- list()
    allRMSE <- rep(0, TotalIter + 1)

    Y <- Y_raw[InitMarker, ]
    Prop0 <- FUN(Y, K)
    allProp[[1]] <- Prop0

    out_all <- csSAM::csfit(Prop0, t(Y_raw))
    prof <- t(out_all$ghat)
    tmpmat <- prof %*% t(Prop0)
    allRMSE[1] <- sqrt(mean((t(Y_raw) - t(tmpmat)) ^ 2))

    message("+========================================+")
    message("+======= Total iterations = ",
                TotalIter, " ==========+")

    for (i in seq_len(TotalIter)) {
        message("Current iter = ", i)

        updatedInx <- DEVarSelect(Y_raw, Prop0, nMarker)
        Y <- Y_raw[updatedInx, ]
        Prop0 <- FUN(Y, K)
        allProp[[i + 1]] <- Prop0

        out_all <- csSAM::csfit(Prop0, t(Y_raw))
        prof <- t(out_all$ghat)
        tmpmat <- prof %*% t(Prop0)
        allRMSE[i + 1] <- sqrt(mean((t(Y_raw) - t(tmpmat)) ^ 2))
    }

    min_idx <- which.min(allRMSE)
    Prop0 <- allProp[[min_idx]]

    return(list(allRMSE = allRMSE,
        allProp = allProp,
        estProp = Prop0
    ))
}
