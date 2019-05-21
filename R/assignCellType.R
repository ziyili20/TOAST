assignCellType <- function(input, reference) {
    # input: estimated proportions matrix
    # reference: reference proportions matrix,
    #    either from RB deconvolution or experiment

    if (ncol(input) != ncol(reference)) {
        stop("Input matrix should have the
            same dimensions as reference matrix!")
    } else {
        K <- ncol(input)
    }

    colnames(input) <- seq(K)
    colnames(reference) <- seq(K)
    corMat <- cor(input, reference,
                use = "pairwise.complete.obs")
    prop_cor <- rep(0, K)
    tmpmat <- corMat
    for (i in seq(K)) {
        maxind <- which(tmpmat == max(tmpmat),
                        arr.ind = TRUE)
        prop_cor[maxind[1]] <- colnames(corMat)[maxind[2]]
        tmpmat[maxind[1], ] <- rep(-1, K)
        tmpmat[, maxind[2]] <- rep(-1, K)
    }
    colnames(input) <- prop_cor
    trans_input <- input[, colnames(reference)]
    return(trans_input)
}
