GetCMatrix_argument <- function(coef, cell_type, W, K, model_names,
                                fitted_model, verbose) {

    if(length(cell_type) == 1){
        if (tolower(cell_type) == "joint") {
            if (length(coef) == 1) {
                cmatrix <- GetCMatrix_CTJoint_CoefOne(coef,
                                               model_names,
                                               W,
                                               verbose)
            } else if (length(coef) == 3) {
                cmatrix <- GetCMatrix_CTJoint_CoefThree(coef,
                                                model_names,
                                                W,
                                                K,
                                                verbose)
            } else {
                stop("Coef should be a valid phenotype!")
            }
        } else if (cell_type %in% colnames(fitted_model$Design_out$Prop)) {
            if (length(coef) == 1) {
                cmatrix <- GetCMatrix_CTOne_CoefOne(coef,
                                             cell_type,
                                             model_names,
                                             W,
                                             verbose)
            } else if (length(coef == 3)){
                cmatrix <- GetCMatrix_CTOne_CoefThree(coef,
                                               cell_type,
                                               model_names,
                                               W,
                                               verbose)
            } else {
                stop("Coef should be a valid phenotype!")
            }
        } else {
            stop("Cell_type is invalid!")
        }

    } else if (length(cell_type) == 2) {
        ct1_indx <- grep(cell_type[1], model_names)
        ct2_indx <- grep(cell_type[2], model_names)
        if (length(ct1_indx)>0 & length(ct2_indx)>0) {
            if(length(coef) == 1){
                cmatrix <- GetCMatrix_CTTwo_CoefOne(coef,
                                             cell_type,
                                             model_names,
                                             W,
                                             ct1_indx,
                                             ct2_indx,
                                             verbose)
            } else if (length(coef) == 2 && length(grep(coef[1],
                            model_names))>0) {
                cmatrix <- GetCMatrix_CTOne_CoefTwo(coef,
                                             cell_type,
                                             model_names,
                                             K,
                                             W,
                                             verbose)
            } else {
                stop("Trp to do cross-cell type testing:
                    Invalid coef value!")
            }
        } else {
            stop("Trp to do cross-cell type testing:
                invalid cell_type values!!")
        }

    } else {
        stop("Please specify a valid cell_type value!")
    }

    return(cmatrix)
}
