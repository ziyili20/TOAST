simplenormalize <- function(Y){
     Y.norm = Y/rowSums(Y)
     return(Y.norm)
}
