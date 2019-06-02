GetCMatrix_CTTwo_CoefOne <- function(coef, cell_type, model_names,
                              W, ct1_indx, ct2_indx, verbose) {
    if (tolower(coef) == "joint") {
        if (verbose) {
            message("Test the joint effect of ", cell_type[1],
                " vs. ", cell_type[2], ". \n", sep = "")
        }
        cmatrix <- matrix(0, nrow = 1, ncol = ncol(W))
        cmatrix[1, ct1_indx] <- 1
        cmatrix[1, ct2_indx] <- -1
    } else if (length(grep(coef, model_names))>0) {
        if (verbose) {
            message("Test the difference of ", cell_type[1], " vs. ",
                cell_type[2], " in different values of ",
                coef, ". \n", sep = "")
        }
        cmatrix <- matrix(0, nrow = 1, ncol = ncol(W))
        ct1_names <- grep(cell_type[1],
                       model_names,
                       value = TRUE)
        ct2_names <- grep(cell_type[2],
                       model_names,
                       value = TRUE)
        ct1_indx <- match(grep(coef, ct1_names, value = TRUE),
                       model_names)
        ct2_indx <- match(grep(coef, ct2_names, value = TRUE),
                       model_names)
        cmatrix[1, ct1_indx] <- 1
        cmatrix[1, ct2_indx] <- -1
    }

    return(cmatrix)
}
