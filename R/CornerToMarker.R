## marker selection by using sisal
CornerToMarker <- function(distances,topN){

     markerList <- apply(distances, 2, function(xx) {
          pure <- rownames(distances)[order(xx)[1:topN]]
          return(pure)
     })
     markerList <-  split(markerList, rep(1:ncol(markerList), each = nrow(markerList)))

     return(markerList)
}
