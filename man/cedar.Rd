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
            pval = NULL, p.adj = NULL, tree = NULL, de.state = NULL,
            cutoff.tree = c('fdr', 0.01), cutoff.prior.prob = c('pval', 0.01),
            parallel.core = NULL, corr.fig = FALSE, tree.type = c('single','full'))
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
  representing cell types.}
  \item{p.adj}{matrix of adjusted p-values, with rows representing features and columns
  representing cell types.}
  \item{tree}{tree structure between cell types}
  \item{de.state}{DE/DM state of each feature in each cell type, with row representing
  features and column representing cell types (1:DE/DM, 0:non-DE/DM)}
  \item{cutoff.tree}{cut off used to define DE state to estimate tree could be
  'fdr' or 'pval' default it 'fdr'=0.01}
  \item{cutoff.prior.prob}{cut off used to define DE state to estimate prior probs
  of nodes on tree could be 'fdr' or 'pval' default it 'fdr'=0.01}
  \item{parallel.core}{number of cores for parallel running, default is NULL}
  \item{corr.fig}{a boolean value, whether to plot corrrelation between cell types
  use function plotCorr()}
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
