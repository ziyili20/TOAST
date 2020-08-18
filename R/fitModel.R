fitModel <- function(Design_out, Y) {

    # Design_out is the output from function makeDesign()
    # Y is a G*N matrix,
    #     G is the number of features, N is the number of subjects

    if (is(Y, "SummarizedExperiment")) {
         se <- Y
         Y <- assays(se)$counts
    } else if (!is(Y, "matrix")) {
         stop("Y should be a matrix or a SummarizedExperiment object!")
    }

    N <- ncol(Y)
    Y <- t(Y)
    W <- Design_out$design_matrix

    # model fit
    inv_WW <- solve(t(W)%*%W)
    H <- W%*%inv_WW%*%t(W)
    Ypred <- H%*%Y
    coefs <- inv_WW%*%t(W)%*%Y
    resi <- Y - Ypred

    Ymean <- colMeans(Y,na.rm = TRUE)
    SSR <- colSums( (Ypred-Ymean)^2, na.rm = TRUE)
    SSE <- colSums( resi^2, na.rm = TRUE)
    MSE <- SSE / (N - ncol(W))
    coefs_var <- matrix( diag(solve(t(W)%*%W)), ncol(W), 1 ) %*% MSE

    rownames(coefs) <- colnames(Design_out$design_matrix)
    rownames(coefs_var) <- colnames(Design_out$design_matrix)

    fitted_model <- list(Design_out = Design_out,
                    N = N,
                    coefs = coefs,
                    coefs_var = coefs_var,
                    Y = t(Y),
                    Ypred = Ypred,
                    resi = resi,
                    all_coefs = Design_out$all_coefs,
                    all_cell_types = Design_out$all_cell_types,
                    MSE = MSE,
                    model_names = colnames(Design_out$design_matrix))

    return(fitted_model)
}
