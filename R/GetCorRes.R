GetCorRes <- function(selProf, knowRefAll) {

     ### implement correlation-based assignment
     isc = intersect(rownames(selProf),rownames(knowRefAll))
     selProf = selProf[isc,]
     knowRef <- knowRefAll[isc, ]
     cormat <- cor(knowRef, selProf)

     ### implement naive bayes classifier
     nbmat <- MyNaiveBayes(selProf, knowRef)

     summat <- (cormat + nbmat)/2
     summat_org <- summat
     assignLabel <- rep("Unassigned", ncol(selProf))
     for(i in 1:ncol(knowRefAll)) {
          thisidx <- which(summat == max(summat), arr.ind = TRUE)
          assignLabel[thisidx[2]] <- rownames(summat)[thisidx[1]]
          summat[thisidx[1],] <- -1
          summat[,thisidx[2]] <- -1 ## I add this
     }
     return(list(assignLabel = assignLabel,
                 probMat = summat_org))
}
