GetCMatrix_CTJoint_CoefThree <- function(coef, model_names, W, K, verbose) {

    if (length(grep(coef[1], model_names)) >0) {
        if (verbose) {
            message("Test the joint effect of ", coef[1],
                " level ", coef[2], " vs. level ", coef[3],
                " in all cell types. \n", sep = "")
        }
        param_names = grep(coef[1], model_names, value = TRUE)
        if (length(grep(paste0(coef[1],
                           coef[3]),
                           param_names)) >0) {
            control_indx = grep(paste0(coef[1],coef[3]),
                            model_names)
            if (length(grep(paste0(coef[1],
                               coef[2]),
                               param_names)) >0) {
                case_indx = grep(paste0(coef[1],coef[2]),
                              model_names)
                cmatrix = matrix(rep(0, ncol(W) * K),
                              nrow = K)
                for (k in seq_len(K)) {
                    cmatrix[k, control_indx[k]] = -1
                    cmatrix[k, case_indx[k]] = 1
                }
            } else if (length(grep(paste0(coef[1],
                                    coef[2]),
                                    param_names)) == 0) {
                cmatrix = matrix(rep(0, ncol(W) * K), nrow = K)
                for(k in seq_len(K)) {
                    cmatrix[k, control_indx[k]] = -1
                }
            }
        } else if (length(grep(paste0(coef[1],
                                coef[3]),
                                param_names)) == 0) {
            if (length(grep(paste0(coef[1],
                               coef[2]),
                               param_names)) == 0) {
                stop("Contrast levels are not valid!")
            }
            case_indx = grep(paste0(coef[1], coef[2]),
                          model_names)
            cmatrix = matrix(rep(0, ncol(W) * K), nrow = K)
            for(k in seq_len(K)){
                cmatrix[k, case_indx[k]] = 1
            }
        } else {
            stop("Coef should be a valid phenotype!")
        }
    } else {
        stop("Coef should be a valid phenotype!")
    }

    return(cmatrix)
}
