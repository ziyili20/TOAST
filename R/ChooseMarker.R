ChooseMarker <- function(pure_all, CellType,
                         nMarkCT = 10,
                         chooseSig = FALSE,
                         verbose = TRUE) {

    K <- length(CellType)
    SelMarker <- list()

    for(k in 1:K) {
        desn <- rep(0, ncol(pure_all))
        desn[CellType[[k]]] <- 1
        fit <- lmFit(pure_all, design = desn)
        fit <- eBayes(fit)
        res <- topTable(fit, number = nMarkCT*5)
        bestRes2 <- res[res$logFC>0,]
        if(chooseSig) {
            tmpMar <- row.names(bestRes2[which(bestRes2$P.Value<0.05),])
            tmpMar2 <- tmpMar[is.na(match(tmpMar, unlist(SelMarker)))]
            tt <- length(tmpMar2)
            if(tt == 0) {
                if(verbose) {
                    message(paste0("Cell type ", k, " has no significant markers."))
                    message(paste0("Switch to selecting top", nMarkCT, "markers."))
                }
                tmpMar <- row.names(bestRes2)[1:(5*nMarkCT)]
                tmpMar2 <- tmpMar[is.na(match(tmpMar, SelMarker))]
                SelMarker[[k]] <- tmpMar2[1:nMarkCT]
            } else {
                if(tt < nMarkCT) {
                    if(verbose) {
                        message(paste0("Cell type ", k, " has ", tt, " significant markers."))
                        message(paste0("Select all of them for cell type ", k, "."))
                    }
                    SelMarker[[k]] <- tmpMar2
                } else {
                    if(verbose) {
                        message(paste0("Cell type ", k, " has ", tt, " significant markers."))
                        message(paste0("Select the top ", nMarkCT,
                                       " markers for cell type ", k, "."))
                    }
                    SelMarker[[k]] <- tmpMar2[1:nMarkCT]
                }
            }
        } else if(!chooseSig) {
            tmpMar <- row.names(bestRes2)[1:(5*nMarkCT)]
            tmpMar2 <- tmpMar[is.na(match(tmpMar, unlist(SelMarker)))]
            SelMarker[[k]] <- tmpMar2[1:nMarkCT]
        }
    }
    names(SelMarker) <- names(CellType)

    return(SelMarker)
}


