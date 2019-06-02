GetCMatrix_contrast <- function(contrast_matrix, W, verbose) {

    if (is.vector(contrast_matrix)) {
        if (verbose) {
            message("contrast_matrix is a vector.")
        }
        if (length(contrast_matrix) != ncol(W)) {
            stop("But it should have length ", ncol(W), "!")
        } else {
            cmatrix <- matrix(contrast_matrix, nrow = 1)
        }
    } else if (is.matrix(contrast_matrix)) {
        if (verbose) {
            message("contrast_matrix is matrix.")
        }
        if (ncol(contrast_matrix) != ncol(W)) {
            stop("But it should have ", ncol(W), " columns!")
        } else {
            cmatrix = contrast_matrix
        }
    }

    return(cmatrix)
}
