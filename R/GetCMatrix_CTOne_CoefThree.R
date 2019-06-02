GetCMatrix_CTOne_CoefThree <- function(coef,
                                       cell_type,
                                       model_names,
                                       W,
                                       verbose) {

    newK = length(cell_type)
    celltype_names = grep(cell_type,
                          model_names,
                          value = TRUE)

    if (length(grep(coef[1], model_names)) >0) {
        if (verbose) {
            message("Test the effect of ", coef[1],
                " level ", coef[2], " vs. level ", coef[3],
                " in ", cell_type,". \n", sep = "")
        }
        param_names = grep(coef[1], celltype_names, value = TRUE)
        if (length(grep(paste0(coef[1],coef[3]),
                        param_names)) >0) {
            control_indx = match(grep(paste0(coef[1],coef[3]),
                                 param_names, value = TRUE),
                             model_names)
            if (length(grep(paste0(coef[1],coef[2]),
                            param_names)) >0) {
                case_indx = match(grep(paste0(coef[1],coef[2]),
                                   param_names, value = TRUE),
                               model_names)
                cmatrix = matrix(rep(0, ncol(W) * length(case_indx)),
                              nrow = length(case_indx))
                for(k in seq_len(length(case_indx))) {
                    cmatrix[k, control_indx[k]] = -1
                    cmatrix[k, case_indx[k]] = 1
                }
            } else if (length(grep(paste0(coef[1],coef[2]),
                                   param_names)) == 0) {
                cmatrix = matrix(rep(0, ncol(W) * length(control_indx)),
                              nrow = length(control_indx))
                for(k in seq_len(length(control_indx))) {
                    cmatrix[k, control_indx[k]] = -1
                }
            }
        } else if (length(grep(paste0(coef[1],coef[3]),
                               param_names)) == 0) {
            if (length(grep(paste0(coef[1],coef[2]),
                            param_names)) == 0) {
                stop("Contrast levels are not valid!")
            }
            case_indx = match(grep(paste0(coef[1],coef[2]),
                               param_names, value = TRUE),
                           model_names)
            cmatrix = matrix(rep(0, ncol(W) * length(case_indx)),
                          nrow = length(case_indx))
            for (k in seq_len(length(case_indx))) {
                cmatrix[k, case_indx] = 1
            }
        }
    } else {
        stop("Coef should be a valid phenotype!")
    }

    return(cmatrix)
}
