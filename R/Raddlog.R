Raddlog <- function(a, b) {
     result <- rep(0, length(a))
     idx1 <- a>b+200
     result[idx1] <- a[idx1]

     idx2 <- b>a+200
     result[idx2] <- b[idx2]

     idx0 <- !(idx1|idx2)
     result[idx0] <- a[idx0] + log1p(exp(b[idx0]-a[idx0]))
     result
}
