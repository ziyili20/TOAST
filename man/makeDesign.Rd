\name{makeDesign}
\alias{makeDesign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate design matrix from input phenotypes and proportions.
}
\description{
This function generate design matrix and
make preparations for following fitModel and csTest.
}
\usage{
makeDesign(design, Prop)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{
  A N by P phenotype matrix, with rows as samples and
  columns as phenotypes (e.g. age, gender, disease, etc.).
}
  \item{Prop}{
  A N by K proportion matrix, with rows as samples and
  columns as cell types
}
}
\value{
  \item{design_matrix }{A comprehensive design matrix
  incorporated phenotype and proportion information.}
  \item{Prop }{The input proportion matrix.}
  \item{design }{The input design/phenotype matrix.}
  \item{all_coefs }{The names of all phenotypes.}
  \item{all_cell_types }{The names of all cell types.}
  \item{formula }{The formula of the tested model,
  including all phenotypes, cell types and interaction terms.}
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

### simulate proportion matrix
Prop = matrix(runif(N*K, 10,60), ncol=K)
Prop = sweep(Prop, 1, rowSums(Prop), FUN="/")
colnames(Prop) = c("Neuron", "Astrocyte", "Microglia")

### simulate phenotype names
design <- data.frame(disease=factor(sample(0:1, size = N,replace=TRUE)),
                     age=round(runif(N, 30,50)),
                     race=factor(sample(1:3, size = N,replace=TRUE)))
Design_out <- makeDesign(design, Prop)
}
