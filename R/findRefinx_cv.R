findRefinx_cv <- function(rawdata, nmarker=1000) {

    mm <- Matrix::rowMeans(rawdata)
    vv <- matrixStats::rowVars(rawdata)
    cv <- sqrt(vv) / mm
    cv[is.na(cv)] <- 0
    ix <- sort(cv, decreasing=TRUE, index=TRUE)$ix

    return(ix[seq(nmarker)])
}
