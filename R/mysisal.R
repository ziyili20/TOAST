mysisal <- function(Y, K, topN){

     Y.norm = simplenormalize(Y)
     sisalres = sisal(t(Y.norm), p = K, iters = 100)
     corner = sisalres$endpoints
     distances = sisalres$distances
     estProp = CornerToEstProp(corner)
     selMarker = CornerToMarker(distances, topN)

     return(list(estProp = estProp,selMarker = selMarker))
}
