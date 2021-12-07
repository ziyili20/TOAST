myprojectMix <- function (Y, Xmat, nonnegative = TRUE, sumLessThanOne = TRUE,
          lessThanOne = !sumLessThanOne)
{
     nCol = dim(Xmat)[2]
     nSubj = dim(Y)[2]
     mixCoef = matrix(0, nSubj, nCol)
     rownames(mixCoef) = colnames(Y)
     colnames(mixCoef) = colnames(Xmat)
     if (nonnegative) {
          if (sumLessThanOne) {
               Amat = cbind(rep(-1, nCol), diag(nCol))
               b0vec = c(-1, rep(0, nCol))
          }
          else if (lessThanOne) {
               Amat = cbind(-diag(nCol), diag(nCol))
               b0vec = c(rep(-1, nCol), rep(0, nCol))
          }
          else {
               Amat = diag(nCol)
               b0vec = rep(0, nCol)
          }
          for (i in 1:nSubj) {
               obs = which(!is.na(Y[, i]))
               Dmat = t(Xmat[obs, ]) %*% Xmat[obs, ]
               mixCoef[i, ] = quadprog::solve.QP(Dmat, t(Xmat[obs, ]) %*%
                                            Y[obs, i], Amat, b0vec)$sol
          }
     }
     else {
          for (i in 1:nSubj) {
               obs = which(!is.na(Y[, i]))
               Dmat = t(Xmat[obs, ]) %*% Xmat[obs, ]
               mixCoef[i, ] = solve(Dmat, t(Xmat[obs, ]) %*% Y[obs,
                                                               i])
          }
     }
     return(mixCoef)
}
