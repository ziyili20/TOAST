### Required input is a matrix of p-values and one of following:
### 1. a matrix of DE/DM state (column is cell type, row is feature);
### 2. a pval threshold for each cell type
### 3. a fdr threshold for each cell type
### If the three inputs are not specified, then we would use fdr < 0.05 as
### threshold for each cell type.

### For odds ratio calculation, the count of 2 by 2 table could be zero
### To avoid this, we add 1 to each cell of the 2 by 2 table

### For log10 pval transformation, pval = 0 may exist, which would lead to -Inf
### Then, for each cell type, we select the minimum non-zero pval, and time
### it with 0.01, then replace it for features with pval = 0

### Requirement1: the two matrix should match dimensions dataframe
### Requirement2: column name: cell type; row name: feature name
### NA suggested to be removed, if not would be removed by function


plotCorr <- function(pval,
                     de.state=NULL,
                     pval.thres=NULL,
                     fdr.thres=NULL,
                     p.size = 0.2,
                     p.color = grDevices::adjustcolor( "black", alpha.f = 0.2),
                     fig.margin = c(1,1,1,1),
                     fig.margin.unit = 'in',
                     line.type = 'dashed',
                     line.color = 'blue'
){
  cell.names <- colnames(pval)
  if(!is.data.frame(pval)){
    if(is.matrix(pval)){
      pval <- data.frame(pval)
    }else{
      stop('pval should be a data frame or matrix')
    }
  }

  dim.pval <- dim(pval)

  if(is.null(cell.names)){
    message('No cell names detected, new names created for each cell type. \n')
    cell.names <- paste0('cell.',seq(1,dim(pval)[2],1))
  }

  if(is.null(de.state) & is.null(pval.thres) & is.null(fdr.thres)){
    fdr.thres = 0.05
    message('No threshold as input to define DE/DMC. Use fdr = 0.05 as threshold.\n')
  }

  if(!is.null(de.state)){
    message('Detect input of de.state; Use de.state to calculate odds ratio.')
    if(!is.data.frame(de.state)){
      if(is.matrix(de.state)){
        de.state <- data.frame(de.state)
      }else{
        stop('de.state should be a data frame or matrix \n')
      }
    }

    dim.de.state <- dim(de.state)
    if( ! all(dim.pval == dim.de.state) ){
      stop('dimensions of two inputs do not match\n')
    }

  }else if(is.null(de.state) & !is.null(pval.thres) ){
    de.state <- matrix(NA, ncol= ncol(pval), nrow = nrow(pval))
    if(length(pval.thres)==1){
      de.state <- (pval < pval.thres)*1
    }else if(length(pval.thres)==ncol(pval)){
      for( i in 1:ncol(pval)){
        de.state[,i] <- (pval[,i] < pval.thres[i])*1
      }
    }else{
      stop('pval.thres length is not correct: assign one value for all cell types
           or assign different values to each cell type.')
    }

    message('Detect input of pval.thres; Use pval.thres to calculate odds ratio.\n')
  }else if(is.null(de.state) & !is.null(fdr.thres)){
    de.state = fdr.trans <- matrix(NA, ncol= ncol(pval), nrow = nrow(pval))

    for( i in 1:ncol(pval)){
      fdr.trans[,i] <- p.adjust(pval[,i], 'fdr')
    }

    if(length(fdr.thres)==1){
      de.state <- (fdr.trans < fdr.thres)*1

    }else if(length(fdr.thres)==ncol(pval)){
      for( i in 1:ncol(fdr.trans)){
        de.state[,i] <- (fdr.trans[,i] < fdr.thres[i])*1
      }
    }else{
      stop('fdr.thres length is not correct: assign one value for all cell types
           or assign different values to each cell type.')
    }

    message('Detect input of fdr.thres; Use fdr.thres to calculate odds ratio.')
  }


  for( i in 1:ncol(pval)){
    if(min(pval[,i])==0){
      pval[pval[,i]==0,i] <- min(pval[pval[,i]!=0,i])*0.01
    }
  }

  colnames(pval) <- cell.names
  colnames(de.state) <- cell.names
  de.state <- data.frame(de.state)

  pval.log10.trans <- data.frame(-log10(pval))
  colnames(pval.log10.trans)  <- cell.names

  cell.num <- dim.pval[2]

  cell.thres <- rep(NA, cell.num)
  for( cell.ix in 1:cell.num ){
    cell.thres[cell.ix] <- min(pval.log10.trans[de.state[,cell.ix]==1,cell.ix])
  }
  message(paste0(c("-log10(pval) threshold for each cell type:\n", round(cell.thres,digits=3) ), collapse = " " ))

  fig.tmp <- GGally::ggpairs(data = pval.log10.trans,
                             upper = list(continuous = GGally::wrap(my_custom_cor,
                                                            de.info = de.state)),
                             lower = list(continuous = GGally::wrap("points", size = p.size,
                                                            colour = p.color)),
                             axisLabels = 'internal')

  for( i in 2:cell.num){
    for(j in 1:(i-1)){
      fig.tmp[i,j] <- fig.tmp[i,j] +
        ggplot2::geom_vline(xintercept = cell.thres[j], linetype = line.type,
                   colour= line.color ) +
        ggplot2::geom_hline(yintercept = cell.thres[i], linetype = line.type,
                   colour= line.color)
    }
  }

  fig.tmp + ggplot2::theme(plot.margin = ggplot2::unit(fig.margin, fig.margin.unit))

}






my_custom_cor <- function(data, mapping, de.info, color = 'black',
                          sizeRange = c(1,5), ...){

  # get the x and y data to use the other code
  x <-GGally::eval_data_col(data, mapping$x)
  y <-GGally::eval_data_col(data, mapping$y)

  ct <- cor.test(x,y) ### test -log10(pval) correlation

  ct.symbol <- c("***", "**", "*", ".", " ")[(ct$p.value  < c(0.001, 0.01, 0.05, 0.1, 1) )][1]
  r <- unname(ct$estimate)
  rt <- format(r, digits=2)[1]

  ## odds ratio
  de.x <- GGally::eval_data_col(de.info, mapping$x)
  de.y <- GGally::eval_data_col(de.info, mapping$y)

  m.tmp <- cbind(de.x, de.y)

  a <- as.numeric( sum( rowSums(m.tmp) == 2 ) ) + 1
  b <- as.numeric( sum( (de.x - de.y) == 1 )) + 1
  c <- as.numeric( sum( (de.x - de.y) == -1 ) ) + 1
  d <- as.numeric( sum( rowSums(m.tmp) == 0 ) ) + 1

  or.est <- format( d/b/c*a, digits=2)
  or.pval <- fisher.test(x=matrix(c(a,c,b,d),2,2))$p.value

  or.symbol <- c("***", "**", "*", ".", " ")[( or.pval <=  c(0.001, 0.01, 0.05, 0.1, 1) )][1]

  cex <- max(sizeRange)

  GGally::ggally_text(
    label = paste0( 'Corr: ',as.character(rt), ct.symbol,'\n',
                    'OR: ', as.character(or.est), or.symbol),
    mapping = ggplot2::aes(),
    xP = 0.5, yP = 0.5,
    size = 5,
    color = color,
    ...
  ) +
    ggplot2::geom_text(
      ggplot2::aes_string(
        x = 0.8,
        y = 0.8
      ),
      label = " ",
      size = 1,
      color = color,
      ...
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(
        color = color,
        linetype = "longdash"
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank()
    )



}
