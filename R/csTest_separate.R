csTest_separate <- function(fitted_model,
                        coef = NULL,
                        cell_type = NULL,
                        contrast_matrix = NULL,
                        var_shrinkage = TRUE,
                        verbose = TRUE,
                        sort = TRUE) {

    N <- fitted_model$N
    K <- ncol(fitted_model$Design_out$Prop)
    W <- fitted_model$Design_out$design_matrix
    beta <- fitted_model$coefs
    beta_var <- fitted_model$coefs_var
    MSE <- fitted_model$MSE
    model_names <- fitted_model$model_names
    G <- length(MSE)

    if (is.null(contrast_matrix)) {
        cmatrix <- GetCMatrix_argument(coef, cell_type,
                                W, K, model_names,
                                fitted_model, verbose)
    } else {
        print("contrast_matrix is specified.")
        print("Coef and cell_type will be ignored.")
        cmatrix <- GetCMatrix_contrast(contrast_matrix, W, verbose)
    }

    ## apply shrinkage on estimated MSE:
    if (var_shrinkage) {
        MSE_threshold <- quantile(c(MSE), 0.1, na.rm = TRUE)
        MSE[MSE < MSE_threshold] <- MSE_threshold
    }

    ## calculate F statistics
    L = cmatrix
    c = rep(0, nrow(L))
    tmp1 <- solve(L %*% solve(t(W) %*% W) %*% t(L))
    inv_sigma <- 1 / MSE
    Lb <-  L %*% beta
    tmp2 <- Lb - matrix(rep(c, ncol(beta)), ncol = ncol(beta))
    a1 <- t(tmp2) %*% tmp1

    res_table = data.frame(
        f_statistics = rep(0, G),
        p_value = rep(0, G),
        fdr = rep(0, G)
    )
    res_table$f_statistics <- rowSums(a1 * t(tmp2)) * inv_sigma
    res_table$p_value <- pf(abs(res_table$f_stat),
                        df1 = 1,
                        df2 = N - ncol(W),
                        lower.tail = FALSE)
    res_table$fdr <- p.adjust(res_table$p_value, method = 'fdr')

    ## If it is to test one parameter, add beta, beta_var, mu and effect size
    if (nrow(cmatrix) == 1 & sum(cmatrix != 0) == 1 & length(coef) == 1) {
        param_num = which(cmatrix != 0)
        res_table$beta <- beta[param_num, ]
        res_table$beta_var <- beta_var[param_num, ]
        if (param_num > K) {
            res_table$mu <- beta[param_num - K, ]
            res_table$effect_size = res_table$beta /
                                (res_table$mu + res_table$beta / 2)
            res_table <- res_table[, c(4, 5, 6, 7, seq_len(3))]
        }
    }

    rownames(res_table) = rownames(fitted_model$Y)

    if (sort) {
        res_table = res_table[order(res_table$p_value), ]
    }

    return(list(
        res_table = res_table,
        design = fitted_model$Design_out$design
    ))
}



