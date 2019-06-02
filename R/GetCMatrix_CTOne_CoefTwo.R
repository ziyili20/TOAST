GetCMatrix_CTOne_CoefTwo <- function(coef,
                              cell_type,
                              model_names,
                              K,
                              W,
                              verbose) {

    if (verbose) {
        message("Test the differences of ", cell_type[1], " vs. ",
            cell_type[2], " in ", coef[1], ":", coef[2], ". \n", sep = "")
    }
    to_test1 <- paste0(cell_type[1], ":", coef[1], coef[2])
    to_test2 <- paste0(cell_type[2], ":", coef[1], coef[2])
    to_test1_indx <- grep(to_test1, model_names)
    to_test2_indx <- grep(to_test2, model_names)
    ct1_indx <- grep(cell_type[1], model_names[seq_len(K)])
    ct2_indx <- grep(cell_type[2], model_names[seq_len(K)])
    if (length(to_test1_indx)>0) {
        cmatrix <- matrix(0, nrow = 1, ncol = ncol(W))
        cmatrix[1, c(ct1_indx, to_test1_indx)] <- 1
        cmatrix[1, c(ct2_indx, to_test2_indx)] <- -1
    } else {
        cmatrix <- matrix(0, nrow = 1, ncol = ncol(W))
        cmatrix[1, ct1_indx] <- 1
        cmatrix[1, ct2_indx] <- -1
    }

    return(cmatrix)
}
