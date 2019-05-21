DEVarSelect <- function(Y_raw, Prop0, nMarker = 1000){

    K <- dim(Prop0)[2]
    N_sample <- dim(Prop0)[1]

    ## find tissue specific genes
    idx <- NULL
    for(k in seq_len(K)) {
        cvec <- rep(-1/(K-1),K)
        cvec[k] <- 1
        design <- rep(0,N_sample)
        tmp <- DEKTissue(K, Y=Y_raw,
                      Prop=Prop0,
                      design=design,
                      contrast_vec=cvec)
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
