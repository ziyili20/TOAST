RefFreeCellMix_wrapper <- function(Y, K) {
    outY = RefFreeEWAS::RefFreeCellMix(Y,
           mu0=RefFreeEWAS::RefFreeCellMixInitialize(Y,
                                                     K = K))
    Prop0 = outY$Omega
    return(Prop0)
}
