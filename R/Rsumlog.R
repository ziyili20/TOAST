Rsumlog <- function(a) {
     s <- a[1]
     for(i in 2:length(a))
          s <- Raddlog(s, a[i])
     s
}
