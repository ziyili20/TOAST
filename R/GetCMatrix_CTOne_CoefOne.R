GetCMatrix_CTOne_CoefOne <- function(coef, cell_type, model_names, W, verbose) {

    newK = length(cell_type)
    celltype_names = grep(cell_type,
                      model_names,
                      value = TRUE)

    if (length(grep(coef, celltype_names)) > 0){
        if (verbose) {
            message("Test the effect of ", coef, " in ",
                cell_type, ". \n", sep = "")
        }
        param_vec = match(grep(coef, celltype_names, value = TRUE),
                       model_names)
        cmatrix = matrix(rep(0, ncol(W) * length(param_vec)),
                      nrow = length(param_vec))
        for (i in seq_len(length(param_vec))) {
            cmatrix[i, param_vec[i]] = 1
        }
    } else {
        stop("Coef should be a valid phenotype!")
    }

    return(cmatrix)
}
