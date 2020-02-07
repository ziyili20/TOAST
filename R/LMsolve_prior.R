LMsolve_prior <- function(Ymat, Param_mat, alpha = NULL,
                          sigma = NULL, Ysigma = NULL, 
                          Prop_constraint = TRUE) {
    K <- ncol(Param_mat)
    N <- ncol(Ymat)
    
    Solved_mat <- matrix(0, K, N)
    Dsigma <- diag(c(1/sigma^2))
    for(n in 1:N) {
        tmp1 <- solve(t(Param_mat)%*%Param_mat + Ysigma^2 * Dsigma)
        tmp2 <- t(Param_mat) %*%Ymat[,n] + Ysigma^2 * Dsigma %*% alpha
        Solved_mat[,n] <- tmp1 %*% tmp2
    }
    Solved_mat[Solved_mat < 0] <- 0
    if (Prop_constraint) {
        Solved_mat<- sweep(Solved_mat, 2, colSums(Solved_mat), "/")
    }
    
    return(Solved_mat)
}