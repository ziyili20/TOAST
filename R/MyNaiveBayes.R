### my naive bayes classifier
MyNaiveBayes <- function(selProf, knowRef) {

     nbres <- matrix(0, ncol(knowRef), ncol(selProf))
     for(ii in 1:ncol(selProf)) {
          initrec <- rep(0, ncol(knowRef))
          for(jj in 1:ncol(knowRef)) {
               initrec[jj] <- -sum((knowRef[,jj] - selProf[,ii])^2)
          }
          nbres[,ii] <- exp(initrec)/exp(Rsumlog(initrec))
     }
     return(nbres)
}
