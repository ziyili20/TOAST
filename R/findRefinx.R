findRefinx <- function(rawdata,
                       nmarker = 1000,
                       sortBy = "var") {

     if (is(rawdata, "SummarizedExperiment")) {
          se <- rawdata
          rawdata <- assays(se)$counts
     } else if (!is(rawdata, "matrix")) {
          stop("rawdata should be a matrix
               or a SummarizedExperiment object!")
     }

     if (nrow(rawdata) < ncol(rawdata)) {
          stop("rawdata matrix should have
               dimension P (features) by N (samples)!")
     }
     if (nmarker > nrow(rawdata)) {
          stop("You have specified nmarker larger
               than the number of original markers!")
     }

     if (sortBy == "cv") {
          mm <- Matrix::rowMeans(rawdata)
          vv <- matrixStats::rowVars(rawdata)
          cv <- sqrt(vv) / mm
          cv[is.na(cv)] <- 0
          final_v <- cv
     } else if(sortBy == "var") {
          vv <- matrixStats::rowVars(log(rawdata+1))
          vv[is.na(vv)] <- 0
          final_v <- vv
     } else {
          stop("sortBy should be either 'cv' or 'var'!")
     }

     ix <- sort(final_v, decreasing=TRUE, index=TRUE)$ix
     return(ix[seq(nmarker)])
}
