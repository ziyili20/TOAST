\name{cedar}
\alias{cedar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing cell type specific differential signals for specified phenotype by
considering DE/DM state corrleation between cell types.
%%  ~~function to do ... ~~
}
\description{
This function provides posterior probability of whether a feature is DE/DM in
certain cell type given observed bulk data.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cedar(Y_raw, prop, design.1, design.2=NULL, factor.to.test=NULL,
            pval = NULL, p.adj = NULL, tree = NULL, p.matrix.input = NULL,
            de.state = NULL, cutoff.tree = c('fdr', 0.01), 
            cutoff.prior.prob = c('pval', 0.01),
            similarity.function = NULL, parallel.core = NULL, corr.fig = FALSE, 
            run.time = TRUE, tree.type = c('single','full'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y_raw}{matrix of observed bulk data, with rows representing features and
  columns representing samples}
  \item{prop}{matrix of cell type composition of samples, with rows representing
  samples and columns representing cell types}
  \item{design.1}{covariates with cell type specific effect, with rows representing
  samples and columns representing covariates}
  \item{design.2}{covariates without cell type sepcific effect, with rows representing
  samples and columns representing covariates}
  \item{factor.to.test}{  A phenotype name, e.g. "disease", or a vector of contrast
  terms, e.g. c("disease", "case", "control").}
  \item{pval}{matrix of p-values, with rows representing features and columns
  representing cell types. colnames must be same as input of prop}
  \item{p.adj}{matrix of adjusted p-values, with rows representing features and columns
  representing cell types. colnames must be same as input of prop}
  \item{tree}{tree structure between cell types, a matrix with row representing 
  layers andcolumn representing cell types (column name is required)}
  \item{p.matrix.input}{prior probability on each node of the tree structure.
  only work when tree structure has been specified. the dimension must be same as tree input.}
  \item{de.state}{DE/DM state of each feature in each cell type, with row representing
  features and column representing cell types (1:DE/DM, 0:non-DE/DM)}
  \item{cutoff.tree}{cut off used to define DE state to estimate tree could be
  'fdr' or 'pval' default it 'fdr'=0.01. suggest to start with restrictive cut off and 
  change to relative loose value when the restrictive cut off is failed}
  \item{cutoff.prior.prob}{cut off used to define DE state to estimate prior probs
  of nodes on tree could be 'fdr' or 'pval' default it 'fdr'=0.01. suggest to start 
  with restrictive cut off and change to relative loose value when the restrictive cut off is failed}
  \item{similarity.function}{custom function used to calculate similarity between cell types 
  that used for tree structure estimation. the input of the custom is assumed to be a matrix
  of log transformed p-value. dimension is: selected gene number * cell number}
  \item{parallel.core}{number of cores for parallel running, default is NULL}
  \item{corr.fig}{a boolean value, whether to plot corrrelation between cell types
  use function plotCorr()}
  \item{run.time}{a boolean value, whether to report running time in seconds}
  \item{tree.type}{tree type for inference, default is c('single','full')}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list
\item{toast_res}{If pval is NULL, then TOAST result by function csTest() is returned}
\item{tree_res}{matrix of posterior probability of each feature for each cell type}
\item{fig}{If corr.fig = TRUE, then figure show DE/DM state correlation between
cell types will be returned}
\item{time_used}{If run.time = TRUE, then running time (seconds) of CeDAR will be 
returned}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Luxiao Chen <luxiao.chen@emory.edu>%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
N <- 300 # simulation a dataset with 300 samples
K <- 3 # 3 cell types
P <- 500 # 500 features

### simulate proportion matrix
Prop <- matrix(runif(N*K, 10,60), ncol=K)
Prop <- sweep(Prop, 1, rowSums(Prop), FUN="/")
colnames(Prop) <- c("Neuron", "Astrocyte", "Microglia")

### simulate phenotype names
design <- data.frame(disease=factor(sample(0:1,size = N,replace=TRUE)),
                     age=round(runif(N, 30,50)),
                     race=factor(sample(1:3, size = N,replace=TRUE)))
Y <- matrix(rnorm(N*P, N, P), ncol = N)
rownames(Y) <- paste0('gene',1:nrow(Y))
d1 <-  data.frame('disease' = factor(sample(0:1,size = N,replace=TRUE)))

### Only provide bulk data, proportion
res <- cedar(Y_raw = Y, prop = Prop,
             design.1 = design[,1:2],
             design.2 = design[,3],
             factor.to.test = 'disease',
             cutoff.tree = c('pval',0.1),
             corr.fig = TRUE,
             cutoff.prior.prob = c('pval',0.1) )

### result of toast (independent test)
str(res$toast_res)
### posterior probability of DE/DM of cedar with single layer tree structure
head(res$tree_res$single$pp)
### posterior probability of DE/DM of cedar with muliple layer tree structure
head(res$tree_res$full$pp)
### estimated tree structure of three cell types
head(res$tree_res$full$tree_structure)
### scatter plot of -log10(pval) showing DE/DM state correlation between cell types
res$fig


### Using custom similarity function to estimate tree structure
### In CeDAR, the input is assumed to be a matrix of log transformed p-values 
### with row representing genes and columns represening cell types

sim.fun <- function(log.pval){
  similarity.res <- sqrt((1 - cor(log.pval, method = 'spearman'))/2)
  return(similarity.res)
}

res <- cedar(Y_raw = Y, prop = Prop,
             design.1 = design[,1:2],
             design.2 = design[,3],
             factor.to.test = 'disease',
             cutoff.tree = c('pval',0.1),
             similarity.function = sim.fun,
             corr.fig = FALSE,
             cutoff.prior.prob = c('pval',0.1) )

### posterior probability of DE/DM of cedar with muliple layer tree structure
head(res$tree_res$full$pp) 


### Using custom tree structure as input
### cell type 1 and cell type 3 are more similar
tree.input <- rbind(c(1,1,1),c(1,2,1),c(1,2,3))
### If column name is provided for the matrix; make sure it is same as variable Prop
colnames(tree.input) <- c("Neuron", "Astrocyte", "Microglia")

res <- cedar(Y_raw = Y, prop = Prop,
             design.1 = design[,1:2],
             design.2 = design[,3],
             factor.to.test = 'disease',
             cutoff.tree = c('pval',0.1),
             tree = tree.input,
             corr.fig = FALSE,
             cutoff.prior.prob = c('pval',0.1) )

### posterior probability of DE/DM of cedar with muliple layer tree structure
head(res$tree_res$custom$pp) 



### Using custom tree structure and prior probability of each node as input
### cell type 1 and cell type 3 are more similar
tree.input <- rbind(c(1,1,1),c(1,2,1),c(1,2,3))
colnames(tree.input) <- c("Neuron", "Astrocyte", "Microglia")

p.matrix.input <- rbind(c(0.2,0.2,0.2), c(0.5,0.25,0.5), c(0.5,1,0.5))
# marginally, each cell type has 0.05 (cell 1: 0.2 * 0.5 * 0.5, cell 2: 0.2 * 0.25 * 1)
# probability to be DE for a randomly picked gene
# there will be about 50\% DE genes in cell type 1 overlaped with cell type 3; 
# while there will be about 25\% DE genes in cell type 1 overlaped with cell type 2

res <- cedar(Y_raw = Y, prop = Prop,
             design.1 = design[,1:2],
             design.2 = design[,3],
             factor.to.test = 'disease',
             cutoff.tree = c('pval',0.1),
             tree = tree.input,
             p.matrix.input = p.matrix.input,
             corr.fig = FALSE,
             cutoff.prior.prob = c('pval',0.1) )

### posterior probability of DE/DM of cedar with muliple layer tree structure
head(res$tree_res$custom$pp) 



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
