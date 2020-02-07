LMsolve <- function(Ymat, Param_mat, Prop_constraint = TRUE) {

    K <- ncol(Param_mat)
    N <- ncol(Ymat)

    Solved_mat <- matrix(0, K, N)
    for(n in 1:N) {
        tmp <- nnls(Param_mat, Ymat[,n])
        Solved_mat[,n] <- tmp$x
    }
    if (Prop_constraint) {
        Solved_mat<- sweep(Solved_mat, 2, colSums(Solved_mat), "/")
    }

    return(Solved_mat)
}
