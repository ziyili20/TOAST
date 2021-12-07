myRefFreeCellMix <- function (Y, mu0 = NULL, K = NULL, iters = 10, Yfinal = NULL,
          verbose = TRUE)
{
     if (is.null(mu0)) {
          if (K == 1) {
               if (!is.null(Yfinal))
                    Y <- Yfinal
               n <- dim(Y)[2]
               mu <- matrix(apply(Y, 1, mean, na.rm = TRUE), ncol = 1)
               omega <- matrix(1, n, 1)
               o <- list(Mu = mu, Omega = omega)
               class(o) <- "RefFreeCellMix"
               return(o)
          }
          else mu0 <- myRefFreeCellMixInitialize(Y, K = K, method = "ward")
     }
     incrementalChangeSummary <- list()
     for (i in 1:iters) {
          flag <- !apply(is.na(mu0), 1, any)
          omega <- myprojectMix(Y[flag, ], mu0[flag, ])
          mu <- myprojectMix(t(Y), omega, sumLessThanOne = FALSE)
          incrementalChangeSummary[[i]] <- summary(abs(as.vector(mu -
                                                                      mu0)))
          if (verbose)
               print(incrementalChangeSummary[[i]])
          mu0 <- mu
     }
     if (!is.null(Yfinal)) {
          mu <- myprojectMix(t(Yfinal), omega, sumLessThanOne = FALSE)
     }
     o <- list(Mu = mu, Omega = omega, incrementalChangeSummary = incrementalChangeSummary)
     class(o) <- "RefFreeCellMix"
     o
}
