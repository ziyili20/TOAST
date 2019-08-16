DEKTissue <- function(K, Y, Prop, design,
                      WhichPar=NULL,
                      contrast_vec=NULL,
                      sort=FALSE,
                      var_threshold=0.1,
                      logged = FALSE, 
                      bound_negative = FALSE) {
    
    if (is(Y, "SummarizedExperiment")) {
        se <- Y
        Y <- assays(se)$counts
    } else if (!is(Y, "matrix")) {
        stop("Y should be a matrix or a SummarizedExperiment object!")
    }
    
    N <- ncol(Y)
    if(nrow(Prop)!=N | ncol(Prop)!=K){
        stop("Dimension of proportion input is not correct!")
    }
    
    Y <- t(na.omit(Y))
    G <- dim(Y)[2]
    
    if(!all(design==0)){
        W <- cbind(Prop,Prop*design)
    }else{
        W <- Prop
    }
    H <- solve(t(W)%*%W)%*%t(W)
    coefs <- H%*%Y
    Ypred <- W%*%coefs
    resi <- Y-Ypred
    
    s2.case <- colSums(resi^2) / (N - ncol(W))
    varBeta <- matrix(diag(solve(t(W)%*%W)),ncol(W),1)%*%s2.case
    
    if(!logged){
        ## bound varBeta a bit
        varBeta[varBeta<var_threshold] <- var_threshold
    }
    
    res_table <- data.frame(t_stat=rep(0,G), t_pval=rep(0,G), t_fdr=rep(0,G))
    rownames(res_table) <- colnames(Y)
    
    if(!all(design==0)){
        if(is.null(contrast_vec)){
            res_table <- data.frame(beta=rep(0,G),
                                    mu=rep(0,G),
                                    effect_size=rep(0,G),
                                    t_stat=rep(0,G),
                                    t_pval=rep(0,G),
                                    t_fdr=rep(0,G))
            rownames(res_table) <- colnames(Y)
            res_table$beta <- coefs[WhichPar,]
            res_table$mu <- coefs[WhichPar-K,]
            res_table$effect_size <- res_table$beta/
                (res_table$mu + res_table$beta/2)
            res_table$t_stat <- coefs[WhichPar,]/sqrt(varBeta[WhichPar,])
            res_table$t_pval <- 2*pt(-abs(res_table$t_stat),df=N-ncol(W))
            res_table$t_fdr <- p.adjust(res_table$t_pval, method = 'fdr')
        }else{
            res_table <- data.frame(t_stat=rep(0,G),
                                    t_pval=rep(0,G),
                                    t_fdr=rep(0,G))
            rownames(res_table) <- colnames(Y)
            cc = matrix(contrast_vec, nrow = 1)
            tmp1 <- cc %*% solve(t(W) %*% W) %*% t(cc)
            demtmp <- as.numeric(sqrt(tmp1 * s2.case))
            res_table$t_stat <- as.numeric(contrast_vec%*%coefs)/demtmp
            res_table$t_pval <- 2*pt(-abs(res_table$t_stat),df=N-ncol(W))
            res_table$t_fdr <- p.adjust(res_table$t_pval, method = 'fdr')
        }
    }else{
        
        if(is.null(contrast_vec)){
            res_table <- data.frame(t_stat=rep(0,G),
                                    t_pval=rep(0,G),
                                    t_fdr=rep(0,G))
            rownames(res_table) <- colnames(Y)
            res_table$t_stat <- coefs[WhichPar,]/sqrt(varBeta[WhichPar,])
            res_table$t_pval <- 2*pt(-abs(res_table$t_stat),df=N-ncol(W))
            res_table$t_fdr <- p.adjust(res_table$t_pval, method = 'fdr')
        }else{
            res_table <- data.frame(muA=rep(0,G),
                                    muB=rep(0,G),
                                    t_stat=rep(0,G),
                                    t_pval=rep(0,G),
                                    t_fdr=rep(0,G))
            rownames(res_table) <- colnames(Y)
            i <- which(contrast_vec!=0)[1]
            j <- which(contrast_vec!=0)[2]
            if (bound_negative) {
                coefs[coefs<0] <- 0
            }
            res_table$muA <- coefs[i,]
            res_table$muB <- coefs[j,]
            cc = matrix(contrast_vec, nrow = 1)
            tmp1 <- cc %*% solve(t(W) %*% W) %*% t(cc)
            demtmp <- as.numeric(sqrt(as.numeric(tmp1) * s2.case))
            res_table$t_stat <- as.numeric(contrast_vec%*%coefs)/demtmp
            res_table$t_pval <- 2*pt(-abs(res_table$t_stat),df=N-ncol(W))
            res_table$t_fdr <- p.adjust(res_table$t_pval, method = 'fdr')
        }
    }
    
    if(sort){
        return(res_table[order(res_table$t_pval),])
    }else{
        return(res_table)
    }
}
