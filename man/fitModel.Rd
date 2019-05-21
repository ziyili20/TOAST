\name{fitModel}
\alias{fitModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit model with proportions and phenotypes.
}
\description{
This function receives design matrix from makeDesign()
and fits the model including all cell types and phenotypes.
}
\usage{
fitModel(Design_out, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Design_out}{
The output from function makeDesign().
}
  \item{Y}{
A G*N matrix, G is the number of features,
N is the number of subjects.
}
}
\value{
  \item{Design_out }{The input Design_out object.}
  \item{N }{Number of samples from matrix Y.}
  \item{coefs }{Estimated coefficients (beta) in the model.}
  \item{coefs_var }{Estimated variance of the
  coefficients (beta variance) in the model.}
  \item{Y }{Observation Y matrix.}
  \item{Ypred }{Predicted Y from the fitted model.}
  \item{all_coefs }{The names of all phenotypes.}
  \item{all_cell_types }{The names of all cell types.}
  \item{MSE }{Estimated mean squared error.}
  \item{model_names }{The names of all terms in the fitted model.}
}
\references{
Ziyi Li, Zhijin Wu, Peng Jin, Hao Wu.
"Dissecting differential signals in high-throughput data from
complex tissues."
}
\author{
Ziyi Li <ziyi.li@emory.edu>
}
\examples{
N = 300 # simulation a dataset with 300 samples
K = 3 # 3 cell types
P <- 500 # 500 features

### simulate proportion matrix
Prop = matrix(runif(N*K, 10,60), ncol=K)
Prop = sweep(Prop, 1, rowSums(Prop), FUN="/")
colnames(Prop) = c("Neuron", "Astrocyte", "Microglia")
Y <- matrix(rnorm(N*P, N, P), ncol = N)

### simulate phenotype names
design <- data.frame(disease=factor(sample(0:1,
                     size = N,replace=TRUE)),
                     age=round(runif(N, 30,50)),
                     race=factor(sample(1:3, size = N,replace=TRUE)))
Design_out <- makeDesign(design, Prop)

### fit model
fitted_model <- fitModel(Design_out, Y)

}
