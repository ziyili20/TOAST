GetCMatrix_CTJoint_CoefOne <- function(coef, model_names, W, verbose) {

    if (length(grep(coef, model_names)) > 0) {
        if (verbose) {
            message("Test the joint effect of ",
                coef,
                " in all cell types. \n",
                sep = "")
        }
        param_vec = grep(coef, model_names)
        cmatrix = matrix(rep(0,
                         ncol(W) * length(param_vec)),
                         nrow = length(param_vec))
        for (i in seq_len(length(param_vec))) {
            cmatrix[i, param_vec[i]] = 1
        }
    } else {
        stop("Coef should be a valid phenotype!")
    }

    return(cmatrix)
}
