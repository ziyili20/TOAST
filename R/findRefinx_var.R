findRefinx_var <- function(rawdata, nmarker=1000) {

    vv <- matrixStats::rowVars(log(rawdata+1))
    vv[is.na(vv)] <- 0
    ix <- sort(vv, decreasing=TRUE, index=TRUE)$ix

    return(ix[seq(nmarker)])
}
