Tsisal <- function(Y_raw, K = NULL, knowRef = NULL,
                   possibleCellNumber = 3:15){

     # Y_raw is the DNA methylation 450K array data from complex tissues,
     # rows for CpG sites and columns for samples.
     # K is the number of pure cell types, we allow users to
     # pre-specify or use our method to estimate.
     # knowRef is the external reference panel for cell type label assignment

     ## cell type number estimation
     if(is.null(K)){
          Kres = getCellNumber(Y_raw, possibleCellNumber)
          K = Kres$bestK
     }

     ## feature selection
     outRF <- csDeconv(Y_raw, K = K, TotalIter = 5, bound_negative = TRUE)
     reduceYmat = Y_raw[outRF$updatedInx,]

     ## simplex corner identification and marker selection
     tsisal_res = mysisal(reduceYmat, K = K, topN = 50)
     estProp = tsisal_res$estProp
     selMarker = tsisal_res$selMarker

     ## cell type label assignment
     if(!is.null(knowRef)){
          out_all <- csSAM::csfit(estProp, t(Y_raw))
          prof <- t(out_all$ghat)
          rownames(prof) <- rownames(Y_raw)
          selProf <- prof[unlist(selMarker),]
          labres <- GetCorRes(selProf, knowRefAll = knowRef)
          colnames(estProp) <- labres$assignLabel
     }

     return(list(estProp = estProp,
                 selMarker = selMarker,
                 K = K))
}
