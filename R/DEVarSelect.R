DEVarSelect <- function(Y_raw, Prop0, nMarker = 1000, bound_negative = FALSE){

    if (is(Y_raw, "SummarizedExperiment")) {
         se <- Y_raw
         Y_raw <- assays(se)$counts
    } else if (!is(Y_raw, "matrix")) {
         stop("Y_raw should be a matrix or a SummarizedExperiment object!")
    }

    if (nrow(Prop0) < ncol(Prop0)) {
         stop("Prop0 should have dimension N (samples) by K (cell types)!")
    }
    if (!ncol(Y_raw) == nrow(Prop0)) {
         stop("Y_raw should have dimension P (features) by N (samples)!")
    }

    K <- ncol(Prop0)
    N_sample <- nrow(Prop0)

    ## find tissue specific genes
    idx <- NULL
    for(k in seq_len(K)) {
        cvec <- rep(-1/(K-1),K)
        cvec[k] <- 1
        design <- rep(0,N_sample)
        tmp <- DEKTissue(K, Y=Y_raw,
                      Prop=Prop0,
                      design=design,
                      contrast_vec=cvec,
                      bound_negative=bound_negative)
        idx[[k]] <- sort(abs(tmp$t_stat),
                      decreasing=TRUE,
                      index=TRUE)$ix
    }
    nmarker <- nMarker
    ## number of markers per tissue. Consider overlaps
    nmarker_tissue <- nmarker/K * 1.2
    idxMarker <- NULL
    for(k in seq_len(K)) {
        idxMarker <- c(idxMarker,
                    idx[[k]][seq_len(nmarker_tissue)])
    }
    idxMarker <- unique(idxMarker)

    return(idxMarker)
}
