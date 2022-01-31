cedar <- function(Y_raw,  # bulk observed data row:feature, col:sample
                  prop,   # cell type proportion: row:sample, col:cell type
                  design.1, # covariates with cell type specific (cs-) effects
                            # row: sample, col: covariate
                  design.2=NULL, # covariates without cs-effects
                                 # row: sample, col: covariate
                  factor.to.test=NULL, # covariate to be  tested
                                       # form1: only name of covariate
                                       # form2: covariate name + contrast levels
                                       #        (ref level at last)
                  pval = NULL, # independent csTest inference result "p-value"
                               # could come from TOAST, TCA or any other methods
                               # row: feature, col: cell type
                  p.adj = NULL, # independent csTest inference result "fdr"
                                # row: feature, col: cell type
                  tree = NULL, # tree used to account cell type correlation
                               # col: cell type, row: tree layer
                               # in same row, different numbers represent different nodes
                               # cell types with same number in same row
                               # means they have same internal node at this level.
                               # example:
                               #   1,1,1,1
                               #   1,1,2,2
                               #   1,2,3,4
                  de.state = NULL, # de.state of each feature in each cell type
                                   # 0: non-DE; 1: DE
                  cutoff.tree = c('fdr', 0.01), # cut off used to define DE
                                                # state to estimate tree
                                                # could be 'fdr' or 'pval'
                                                # default it 'fdr'=0.01
                  cutoff.prior.prob = c('pval', 0.01), # cut off used to
                                                       # define DE state to
                                                       # estimate prior prob
                                                       # could be 'fdr' or 'pval'
                                                       # default is 'pval'=0.01
                  parallel.core = NULL, # integer to specificy number of
                                        # cores to use
                  corr.fig = FALSE,     # whether to plot pval correlation
                  tree.type = c('single','full') # two tree structures as input
){

  ### Basic info
  gene.num <- dim(Y_raw)[1]
  cell.num <- dim(prop)[2]

  ### Step 1: first round inference with any possible packages.
  ###         If choose to use TOAST (pval is null), then could use a function
  ###         called 'toast.first.round'
  ###         This will give run toast for each cell type and put all results in
  ###         a list

  ### no matter whether we use toast for first round inference, the design matrix
  ### is needed for posterior probability calculation in later steps
  Design_out <- makeDesign_tree(design.1 = design.1, design.2=design.2, Prop=prop,
                                factor.to.test=factor.to.test)

  Design_matrix <- Design_out$design_matrix

  toast_res <- NULL
  if( is.null(pval) & is.null(p.adj) ){
    message('No prior inference information, run TOAST for first round inference \n')

    fitted_model <- fitModel( Design_out, Y_raw )
    toast_res <- toast.first.round(fitted_model = fitted_model,
                                   celltypes = fitted_model$all_cell_types,
                                   coef = factor.to.test)

    ### extract pval and fdr information from first round toast analysis

    cell.types <-  fitted_model$all_cell_types

    pval = p.adj <- matrix(NA, nrow = gene.num, ncol = cell.num )
    colnames(pval) = colnames(p.adj) <- cell.types
    rownames(pval) = rownames(p.adj) <- rownames(Y_raw)

    for( i in 1:cell.num ){
      pval[,i] <- toast_res[[i]]$p_value
      p.adj[,i] <- toast_res[[i]]$fdr
    }

  }else if( !is.null(pval) & is.null(p.adj) ){
    cell.types <- colnames(pval)
    p.adj <- matrix(NA, nrow = gene.num, ncol = cell.num )
    colnames(p.adj) <- cell.types
    rownames(p.adj) <- rownames(pval)
    ### calculate fdr for each cell type p-value.
    for( i in 1:cell.num ){
      p.adj[,i] <- p.adjust(pval[,i],method='fdr')
    }

  }else if( is.null(pval) & !is.null(p.adj) ){
    stop( 'p-value is needed')
  }
  
  ### determine DE state for tree estimation if not provided
  de.res <- matrix(NA,ncol=cell.num, nrow=gene.num)
  if(is.null(de.state)){
    # if de.state not provided then
    # use fdr as cutoff
    if(cutoff.tree[1]=='fdr'){
      if(length(cutoff.tree) == (cell.num + 1) ){
        # in this case, each cell type has different cutoffs
        for(cell.ix in 1:cell.num){
          de.res[,cell.ix] <- (p.adj[,cell.ix] <
                                 as.numeric(cutoff.tree[,(cell.ix+1)]))*1
        }
      }else if(length(cutoff.tree) == 2){
        # in this case, all cell types have same cutoffs
        de.res <- (p.adj < as.numeric(cutoff.tree[2])*1)
      }else{
        stop("length of cutoff.tree is incorrect.")
      }

    }else if(cutoff.tree[1]=='pval'){
      # use pvalue as cutoff
      if(length(cutoff.tree) == (cell.num + 1) ){
        for(cell.ix in 1:cell.num){
          de.res[,cell.ix] <- (pval[,cell.ix] <
                                 as.numeric(cutoff.tree[,(cell.ix+1)]))*1
        }
      }else if(length(cutoff.tree) == 2){
        de.res <- (pval < as.numeric(cutoff.tree[2]))*1
      }else{
        stop("length of cutoff.tree is incorrect.")
      }
    }else{
      stop('Invalid input of cutoff.tree variable')
    }
  }else{
    de.res <- de.state
  }

  ### transformation of pvalue = 0 to pvalue = min.pval* 0.001
  ### in case log10 transformation can not work correctly
  for(cell.ix in 1:cell.num){
    min.tmp <- min(pval[pval[,cell.ix]>0 ,cell.ix])
    if(-log10(min.tmp) > 300 ){
      pval[ pval[,cell.ix]==0, cell.ix] <- min.tmp
    }else{
      pval[ pval[,cell.ix]==0, cell.ix] <- min.tmp*0.001
    }
  }


  ### Step 1.5: plotCorr show correlation between cell types based on first round
  ###           independent test result
  fig.res <- NULL
  if(corr.fig == TRUE){
    message('Generating scatter plot of -log10(pval) among cell types \n')
    fig.res <- plotCorr(pval = data.frame(pval),
                        de.state = data.frame(de.res))
  }

  ### Step 2: Estimate a tree structure based on -log10(pval)
  ###         Features used should be filtered by users
  ###         This includes two modes: 1. user gives features to build up trees
  ###                                  2. Estimate a tree by FDR/pval cutoff

  if( is.null(tree) ){
    tree.input <- tree.est(pval, de.res)
  }else{
    tree.type <- 'custom'
    tree.input <- list()
    tree.input[['custom']] <- tree
  }

  ### Step 3: Estimate prior probabilities based on given tree structure
  ### Step 4: Estimate posterior probability
  ### Input includes:
  ### 1. tree structure
  ### 2. design matrix
  ### 3. observed bulk data
  ### 4. first round pvalue result (could come from other packages)

  if(cutoff.prior.prob[1]=='fdr'){
    cutoff.fdr <- as.numeric(cutoff.prior.prob[2])
    cutoff.pval <- NULL
    inference_input <- p.adj
  }else if(cutoff.prior.prob[1] == 'pval'){
    cutoff.fdr <- NULL
    cutoff.pval <- as.numeric(cutoff.prior.prob[2])
    inference_input <- pval
  }

  tree_res <- list()

  for(tree.ix in tree.type){
    message(paste0('inference with tree: ',tree.ix, '\n'))
    tree_res[[tree.ix]]<- post_calc_tree(Design_matrix = Design_matrix,
                                         Y_raw = Y_raw,
                                         tree.input = tree.input[[tree.ix]],
                                         toast_res = inference_input,
                                         design.1 = design.1,
                                         design.2 = design.2,
                                         factor.to.test = factor.to.test,
                                         fdr=cutoff.fdr,
                                         p.value = cutoff.pval,
                                         core.num = parallel.core
    )
  }


  all.res <- list('toast_res'=toast_res, 'tree_res'=tree_res, 'fig'=fig.res)
  return(all.res)
}


### function name: makeDesign_tree
### function usage: create design matrix
### parameter explanation:
### design.1: covariates have cell type specific effect (N * p1)
### design.2: covariates have same effect for all cell types (N * p2)
### Prop: proportion for each cell type (N * K)
### factor.to.test: the covariate going to be tested (covariate name or
###                                                   covariate name with
###                                                   two levels to compare )

makeDesign_tree <- function(design.1, design.2, Prop, factor.to.test){
  if (!is.vector(Prop)) {
    # cell type number K
    K <- ncol(Prop)
    if (is.null(colnames(Prop))) {
      # create cell type name if it does not exist
      colnames(Prop) <- paste0("celltype", seq_len(K))
    }
  }
  else {
    K <- 1
  }

  covariate.name.1 <- colnames(design.1)
  covariate.name.2 <- colnames(design.2)
  celltype.name <- colnames(Prop)

  if(ncol(design.1) > 1){
    if(length(factor.to.test) == 1){
      # In this situation, user wants to test the covariate globally if
      # factor.to.test contains multiple levels.

      # put the covariate to last column for convenience
      covariate.name.1 <- c(covariate.name.1[covariate.name.1 != factor.to.test],
                            factor.to.test)
      design.1 <- data.frame(design.1[,covariate.name.1])
    }else if(length(factor.to.test) == 3){
      # In this situation, user wants to compare two levels
      # make the second level as reference for convenience

      covariate.name.1 <- c(covariate.name.1[covariate.name.1 != factor.to.test[1]],
                            factor.to.test[1])
      design.1 <- data.frame(design.1[,covariate.name.1])
      # set the third element in 'factor.to.test' as reference level
      design.1[,factor.to.test[1]] <- relevel(design.1[,factor.to.test[1]], factor.to.test[3])
    }

  }

  # in our model, we also include cell type composition as main terms
  # in addition, design.2 provide covariates have NO cell type specific effects
  # so they are also modeled as main terms.
  # Put them in front of cell type specific covariates for convenience
  if(is.null(design.2)){
    dd <- cbind(Prop, design.1)
  }else{
    dd <- cbind(design.2, Prop, design.1)
  }

  # create formula for main terms
  formul <- paste("~", paste(c(covariate.name.2, celltype.name), collapse = "+"))

  # create interaction term
  for (i in 1:ncol(design.1)) {
    tmp <- paste(celltype.name, covariate.name.1[i],
                 sep = ":")
    intterms <- paste(tmp, collapse = "+")
    formul <- paste(formul, intterms, sep = "+")
  }

  design_matrix <- model.matrix(as.formula(formul), dd)[,-1] # remove intercept
  # remove column with all zeros
  design_matrix <- design_matrix[, !colSums(design_matrix == 0) == nrow(Prop)]
  # remove replicated column
  design_matrix <- unique(design_matrix, MARGIN = 2)

  formul <- paste0("~ ", paste(colnames(design_matrix), sep = "+",
                               collapse = "+"))
  return(list(design_matrix = design_matrix, Prop = Prop, design.1 = design.1,
              design.2 = design.2, all_coefs = c(covariate.name.1, covariate.name.2),
              all_cell_types = colnames(Prop), formula = formul))
}


### function name: toast.first.round
### function usage: if the user does not provide any independent inference result
###                 then do csTest with TOAST method
toast.first.round <- function(fitted_model, celltypes, coef, var_shrinkage =T){
  ### used to store result
  res <- list()
  ### test for each cell type
  for(cell in 1:length(celltypes)){
    ### store the TOAST result in res[[cell]] list
    res[[celltypes[cell] ]] <- csTest(fitted_model, coef = coef,
                                      cell_type = celltypes[cell], sort = F,
                                      var_shrinkage = var_shrinkage, verbose = F)
  }

  return(res) ### return list res
}


### function name: tree.est
### function usage: estimate tree structure based on p-values
tree.est <- function(pval, de.res){

    tree.input <- list()
    cell.num <- ncol(pval)

    log.pval <- -log10(pval[rowSums(de.res) > 0, ])
    dist.corr <- as.dist( (1 - cor(log.pval, method = 'pearson' ))/2 )
    hc.tree <- hclust( dist.corr )

    cut.height <- sort(c(0,hc.tree$height),decreasing = T )
    cut.thres <- cut.height

    tree.input[['full']] <- t(cutree(hc.tree, h=cut.thres))
    tree.input[['single']] <- rbind(rep(1,cell.num), seq(1,cell.num,1))

    colnames(tree.input[['full']]) = colnames(tree.input[['single']]) <- colnames(pval)

    tree.input

}



post_calc_tree <- function( Design_matrix,  Y_raw, tree.input, toast_res,
                            factor.to.test, fdr, design.1, design.2,
                            p.value =NULL, no_prior_info =FALSE,
                            p.matrix.input = NULL,core.num=NULL){

  gene.num <- nrow(Y_raw)           ### gene number
  sample.size <- ncol(Y_raw)
  cell.num <- ncol(tree.input)      ### cell number

  # create all DZ combination based on estimated tree structure
  dz.combine.tmp <- DZ_combination_gen(tree.input)
  var.index  <- dz.combine.tmp[[1]]
  dz.combine <- dz.combine.tmp[[2]]

  z.states <- dz.combine[,(max(var.index)-cell.num+1):(max(var.index))]
  all.z.states <- expand.grid(rep(list(0:1), cell.num ))
  colnames(all.z.states) <- colnames(tree.input)

  match.index <- index.match( dz.combine, all.z.states )

  M.reduce <- M.cal(Design_matrix, all.z.states,
                    design.1=design.1, design.2 = design.2,
                    factor.to.test = factor.to.test)

  W.all    <- W.cal(var.index, dz.combine, toast_res, fdr = fdr,
                    p.value = p.value, p.matrix = p.matrix.input)
  weight.all <- W.all[['weight']]
  est.prob <- W.all[['est_prob']]

  ### different combination numbers between Z's and D
  dz.combine.num <- nrow(dz.combine)
  all.z.num <- nrow(all.z.states)

  ### pp is used to store posterior probability:
  pp <- matrix(NA, ncol=(cell.num),nrow=gene.num)

  p.ycz.sum <- matrix(NA, nrow = all.z.num, ncol = gene.num)

  ### log P(Z=1, Y) a vector, in which each element corresponding to a cell type
  p.z1y <- matrix(NA, nrow = gene.num, ncol = cell.num)

  if(is.null(core.num)){
    cores_num <- 1
  }else{
    cores_num <- min(core.num, parallel::detectCores()) #detectCores() - 2
  }

  p.ycz.sum.temp <- Marginal.L.parallel(X=M.reduce, Y=Y_raw, numCores = cores_num)
  p.ycz.sum <- t(p.ycz.sum.temp)

  if(no_prior_info == TRUE){
    p.yzd <- t(p.ycz.sum)
    p.y <- Rsumlog.matrix(p.yzd)

    for( i in 1:cell.num){
      p.z1y[,i] <- Rsumlog.matrix(p.yzd[,all.z.states[,i]==1 ])
    }

  }else{

    p.ycz.sum.dz.states <- p.ycz.sum[match.index, ]

    p.yzd <- t(p.ycz.sum.dz.states + log(weight.all))

    p.y   <- Rsumlog.matrix(p.yzd)

    for( i in 1:cell.num){
      p.z1y[,i] <- Rsumlog.matrix(p.yzd[,z.states[,i]==1 ])
    }

  }

  pp <- exp(p.z1y - p.y )
  colnames(pp) <- colnames(tree.input)
  rownames(pp) <- rownames(Y_raw)

  res.all <- list('est_prob'= est.prob, 'weight'=weight.all, 'dz_combine'=dz.combine,
                  'tree_structure'=tree.input, 'node_index'=var.index ,'pp'=pp)
  return(res.all)

}


DZ_combination_gen <- function(tree.input ){

  tree.level <- dim(tree.input)[1]
  cell.num   <- dim(tree.input)[2]

  var.index <- matrix(NA, ncol = cell.num, nrow = tree.level)
  var.num <- rep(NA, tree.level)
  for( i in 1:tree.level){
    var.num[i] <- length( unique(tree.input[i,]) )
    var.index[i,] <- tree.input[i,]  + sum(var.num[1:(i-1)])*( as.integer(i>1))
  }

  combine.current.tmp <- NULL

  for(tl.ix in 1:tree.level){
    tl.element.num <- length(unique(var.index[tl.ix,]))
    combine.single.tmp <-  expand.grid( rep(list(0:1), tl.element.num))
    combine.current.tmp <-  as.matrix( tidyr::expand_grid( x1=combine.current.tmp, x2 = combine.single.tmp ) )
    ## filter criteria:
    ## 1. splitting node: 0 -> 0, 1 -> 1 or 0;
    ## 2. single node keep previous result
    ## start from second layer of tree to the end of the tree
    if(tl.ix > 1){
      for( cell in 1:cell.num){
        if(tl.ix == 2){
          ix.rm <-  which(  apply(as.matrix(combine.current.tmp[,var.index[1:tl.ix,cell]]),1,diff) == 1 )
        }else{
          ix.rm <-  which( apply( apply(as.matrix(combine.current.tmp[,var.index[1:tl.ix,cell]]),1,diff), 2, max) == 1 )
        }
        if(length(ix.rm) > 0 ){
          combine.current.tmp <- combine.current.tmp[-ix.rm,]
        }

      }

      # criteria 2
      ## 2. single node keep previous result
      ## strategy: a. check whether up node split; b. if not split, we only keep combination that keep same between the two nodes in two layers
      tl.up.ix <- tl.ix - 1
      nodes.up <- var.index[tl.up.ix, ]
      nodes.current <- var.index[tl.ix,]
      nodes.up.unique <- unique(nodes.up)
      for(nodes.up.ix in nodes.up.unique){
        if( length(unique(nodes.current[nodes.up == nodes.up.ix])) == 1 ){
          ix.rm <- which( apply(combine.current.tmp[ ,  c(nodes.up.ix, unique(nodes.current[nodes.up == nodes.up.ix]) ) ], 1, diff) !=0 )
          if(length(ix.rm) > 0 ){
            combine.current.tmp <- combine.current.tmp[-ix.rm,]
          }
        }


      }
    }


  }

  return(list(node.index = var.index, dz.combine = combine.current.tmp))
}


index.match <- function(dz.combine, all.z.states){
  match.ix <- c()
  cell.num <- ncol(all.z.states)
  col.num <- dim(dz.combine)[2]
  z.cols <- (col.num - cell.num + 1):col.num

  for( i in 1:nrow(dz.combine)){

    for( n in 1:nrow(all.z.states)){
      if(all(dz.combine[i, z.cols] == all.z.states[n,])  ){
        match.ix <- c(match.ix, n)
        break
      }
    }
  }
  return(match.ix)
}


M.cal <- function(Design_matrix, all.z.states, design.1, design.2, factor.to.test){
  M.reduce <- list()
  cell.num <- ncol(all.z.states)
  cell.types <- colnames(all.z.states)
  dz.combine.num <- nrow(all.z.states)

  for(dz.ix in 1:dz.combine.num){
    ### remove interaction columns corresponding to Z's = 0
    if(all(all.z.states[dz.ix,]==1)){
      design.reduce <- Design_matrix
    }else{
      celltype.to.remove <- cell.types[all.z.states[dz.ix,] == 0]

      if(length(factor.to.test)==1){

        term.tmp <- paste0(celltype.to.remove,':' ,factor.to.test)
        name.index <- c()
        for(term.ix in 1:length(term.tmp)){
          name.index <- c(name.index, grep(term.tmp[term.ix], colnames(Design_matrix)))
        }
        name.index <- unique(sort(name.index))

      }else if(length(factor.to.test)==3){

        term.tmp <- paste0(celltype.to.remove, ':', factor.to.test[1],
                           factor.to.test[2])
        name.index <- c()
        for(term.ix in 1:length(term.tmp)){
          name.index <- c(name.index, grep(term.tmp[term.ix], colnames(Design_matrix)))
        }
        name.index <- unique(sort(name.index))

      }

      design.reduce <- Design_matrix[,  - name.index ]
    }

    ### projection matrix for column spaced defined on design.reduce
    M.reduce[[dz.ix]] <- design.reduce %*% solve( t(design.reduce) %*% design.reduce ) %*% t(design.reduce)
  }
  return(M.reduce)
}

W.cal <- function(var.index, dz.combine, toast_res, fdr = 0.05, p.value = NULL, p.matrix = NULL){
  dz.combine <- as.matrix(dz.combine)

  if(is.list(toast_res)){
    cell.num <- length(toast_res)
    gene.num <- nrow(toast_res[[1]])
  }else if(is.matrix(toast_res)){
    cell.num <- ncol(toast_res)
    gene.num <- nrow(toast_res)
  }

  layer.num <- nrow(var.index)
  dz.combine.num <- nrow(dz.combine)
  fdr.all.cell <- matrix(NA,ncol=cell.num,nrow = gene.num )

  if(is.list(toast_res)){
    for(i in 1:cell.num){
      if( is.null(p.value)){
        fdr.all.cell[,i] <- toast_res[[i]]$fdr
      }else if( is.null(fdr)){
        fdr.all.cell[,i] <- toast_res[[i]]$p_value
      }
    }
  }else if(is.matrix(toast_res)){
    fdr.all.cell <- toast_res
  }

  if( is.null(p.value)){
    de.all.cell <- (fdr.all.cell < fdr)*1
  }else if( is.null(fdr)){
    de.all.cell <- (fdr.all.cell < p.value)*1
  }

  if( is.null(p.matrix) ){
    p.cond.info <- matrix(NA, ncol = cell.num , nrow = layer.num)
    ### Calculate conditional probability
    for(l in (layer.num -1):1){
      nodes <- unique(var.index[l,])
      node.num <- length(nodes)

      for( node in nodes){

        cell.ix <- which(var.index[l,] == node)
        child.var <- unique(var.index[l+1,cell.ix])

        if(length(child.var) ==1){
          p.cond.info[(l+1),cell.ix] <- 1
        }else{
          de.all.ix <- which(rowSums(de.all.cell[,cell.ix]) > 0)

          for(child in child.var){
            ix <- which(var.index[l+1, ] == child)
            if(length(ix) == 1){
              de.num.tmp <- length(which(de.all.cell[,ix] >0 ) )
              if(de.num.tmp ==0){
                de.num.tmp <- 1
              }
              p.cond.info[l+1, ix ] <- de.num.tmp/ length(de.all.ix)
            }else{
              de.num.tmp <- length(which(rowSums(de.all.cell[,ix]) >0 ) )
              if(de.num.tmp ==0){
                de.num.tmp <- 1
              }
              p.cond.info[l+1, ix ] <- de.num.tmp/ length(de.all.ix)
            }
          }
        }

        if(l == 1){
          if(length(cell.ix) ==1){
            p.cond.info[l, cell.ix] <- mean(de.all.cell[,cell.ix]>0)
          }else{
            p.cond.info[l, cell.ix] <- mean(rowSums(de.all.cell[,cell.ix])>0 )
          }
        }

      }
    }
  }else{
    p.cond.info <- p.matrix
  }

  ###### calculate prior probs (weights) for each DZ combination
  prod.keep <- matrix(0, ncol= cell.num, nrow= layer.num )
  for( i in 1:layer.num){
    nodes <- unique(var.index[i,])
    for(node in nodes){
      prod.keep[ i, which(var.index[i, ] == node )[1] ] <- 1
    }

  }

  weights.all <- rep(NA, dz.combine.num )
  for(dz.combine.ix in 1:dz.combine.num){
    dz.state <- dz.combine[dz.combine.ix, ]
    #cat(dz.state,'\n')
    dz.state.matrix <- matrix(dz.state[c(var.index)], nrow = layer.num, ncol = cell.num, byrow=F)

    dz.state.cond.matrix.1 <- rbind(dz.state.matrix[1,], dz.state.matrix[-layer.num,] )
    dz.state.cond.matrix.2 <- rbind(1 - dz.state.matrix[1,], dz.state.matrix[-layer.num,] )

    dz.state.power.1 <- dz.state.matrix * dz.state.cond.matrix.1
    dz.state.power.2 <- (1 - dz.state.matrix) * dz.state.cond.matrix.2

    p.matrix <- (p.cond.info ^ dz.state.power.1 ) * ( (1 - p.cond.info) ^ dz.state.power.2)

    weights.all[dz.combine.ix] <- prod( p.matrix ^ prod.keep)
  }


  res <- list('weight' = weights.all, 'est_prob'= p.cond.info)
  return(res)

}

Marginal.L.parallel <- function(X, Y, numCores){
  N = length(X)
  sample.size <- ncol(Y)
  foo <- function(i, Y, X, sample.size) {
    mu_est <- Y %*% t(X[[i]])
    resi.temp <- Y - mu_est
    sd_est <- sqrt( rowMeans( (resi.temp - rowMeans(resi.temp))^2 ) * (1 - 1/(sample.size)) )
    p.ycz <- dnorm(Y, mu_est, sd_est, log = T)
    rowSums( p.ycz )
  }
  doParallel::registerDoParallel(cores=numCores)
  cl <- parallel::makeCluster(numCores)
  p.ycz.sum <- parallel::parSapply(cl, 1:N, foo, Y, X, sample.size)
  parallel::stopCluster(cl)
  return(p.ycz.sum)
}


Rsumlog.matrix <- function(a){
  s <- a[,1]

  for( i in 2:dim(a)[2] ){
    s <- Raddlog.matrix(s, a[,i])

  }

  return(s)
}

Raddlog.matrix <- function(a,b){
  result <- rep(0, length(a))
  idx1 <- (a>b+200) | is.infinite(b)
  result[idx1] <- a[idx1]

  idx2 <- (b>a+200) | is.infinite(a)
  result[idx2] <- b[idx2]

  idx0 <- !(idx1|idx2)
  result[idx0] <- a[idx0] + log1p(exp(b[idx0]-a[idx0]))
  return(result)
}


