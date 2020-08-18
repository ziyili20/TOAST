CornerToEstProp <- function(corner){

     N_sample = dim(corner)[1]
     tmp <- nnls(corner, rep(1,N_sample))
     estProp <- diag(as.numeric(tmp$x),length(as.numeric(tmp$x))) %*% t(corner)
     estProp[estProp < 0] = 0
     estProp[estProp > 1] = 1

     return(t(estProp/colSums(estProp)))
}
