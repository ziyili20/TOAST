RefFreeCellMix_wrapper <- function(Y, K) {

     if (is(Y, "SummarizedExperiment")) {
          se <- Y
          Y <- assays(se)$counts
     } else if (!is(Y, "matrix")) {
          stop("Y should be a matrix
               or a SummarizedExperiment object!")
     }

    if (K<0 | K>ncol(Y)) {
         stop("K should be between 0 and N (samples)!")
    }
    outY <- myRefFreeCellMix(Y,
           mu0=myRefFreeCellMixInitialize(Y, K = K))
    Prop0 <- outY$Omega
    return(Prop0)
}
