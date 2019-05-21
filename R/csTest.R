csTest <- function(fitted_model,
             coef = NULL,
             cell_type = NULL,
             contrast_matrix = NULL,
             var_shrinkage = TRUE,
             verbose = TRUE,
             sort = TRUE) {

    # fitted_model is the output from fitModel().
    # coef is a phenotype name, e.g. "disease",
    #    or a vector of contrast terms, e.g. c("disease", "case", "control").
    # cell_type is a cell type name, e.g. "celltype1", or "neuron".
    #    If cell_type is NULL
    #    or specified as "ALL",
    #    compound effect of coef in all cell types will be tested.
    # contrast_matrix is a matrix (or a vector) to specify contrast, e.g.,
    #    cmat <- matrix(0, 2, 6); cmat[1,3] <- 1: cmat[2,4] <- 1
    #    is to test whether the 3rd parameter
    #    and 4th parameter are zero simultaneously
    #    i.e. \beta_{3} = \beta_{4} = 0.
    #    If contrast_matrix is specified, coef and cell_type will be ignored!
    # var_shrinkage is to apply shrinkage on estimated MSE.  Applying shrinkage
    #          helps remove extremely small variance estimation
    #          and stablize statistics.
    # verbose is a boolean parameter. Testing information will be printed,
    #          if verbose = TRUE.
    # sort is a boolean parameter.
    #       The output results will be sorted by p value
    #       if sort = TRUE.


    if (!is.null(contrast_matrix)) {
      out <- csTest_separate(fitted_model = fitted_model,
                    coef = coef,
                    cell_type = cell_type,
                    contrast_matrix = contrast_matrix,
                    var_shrinkage = var_shrinkage,
                    verbose = verbose,
                    sort = sort)
      res_table <- out$res_table
    } else {
      if (!is.null(coef)) {
        if(!is.null(cell_type)) {
           out <- csTest_separate(fitted_model = fitted_model,
                          coef = coef,
                          cell_type = cell_type,
                          contrast_matrix = contrast_matrix,
                          var_shrinkage = var_shrinkage,
                          verbose = verbose,
                          sort = sort)
           res_table <- out$res_table
        } else {
           res_table <- list()
           cell_type_all <- c(fitted_model$Design_out$all_cell_types,
                       "joint")
           for (c in seq_len(length(cell_type_all))) {
             out <- csTest_separate(fitted_model = fitted_model,
                            coef = coef,
                            cell_type = cell_type_all[c],
                            contrast_matrix = contrast_matrix,
                            var_shrinkage = var_shrinkage,
                            verbose = verbose,
                                   sort = sort)
                res_table[[c]] <- out$res_table
             }
             names(res_table) <- cell_type_all
          }
       } else {
           res_table <- list()
           coef_all <- c(fitted_model$Design_out$all_coefs, "joint")
           for (c in seq_len(length(coef_all))) {
               out <- csTest_separate(fitted_model = fitted_model,
                                 coef = coef_all[c],
                                 cell_type = cell_type,
                                 contrast_matrix = contrast_matrix,
                                 var_shrinkage = var_shrinkage,
                                 verbose = verbose,
                                 sort = sort)
               res_table[[c]] <- out$res_table
           }
           names(res_table) <- coef_all
       }
    }

    return(res_table)
}
