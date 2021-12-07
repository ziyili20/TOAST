myRefFreeCellMixInitialize <- function (Y, K = 2, Y.Distance = NULL, Y.Cluster = NULL, largeOK = FALSE,
          dist.method = "euclidean", ...)
{
     if (!is.matrix(Y) | !is.numeric(Y)) {
          stop("Y is not a numeric matrix\n")
     }
     n <- dim(Y)[2]
     if (is.null(Y.Cluster)) {
          if (is.null(Y.Distance)) {
               if (n > 2500 & !largeOK) {
                    stop("Y has a large number of subjects!  If this is what you really want, change 'largeOK' to TRUE\n")
               }
               Y.Distance <- dist(t(Y), method = dist.method)
          }
          Y.Cluster <- hclust(Y.Distance, ...)
     }
     classes <- cutree(Y.Cluster, K)
     s <- split(1:n, classes)
     sapply(s, function(u) apply(Y[, u, drop = FALSE], 1, mean,
                                 na.rm = TRUE))
}
