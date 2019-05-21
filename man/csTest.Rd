\name{csTest}
\alias{csTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing differential signals for specified phenotype and cell type(s).
}
\description{
This function conducts statistical tests
for specified phenotype and cell type(s).
}
\usage{
csTest(fitted_model, coef = NULL, cell_type = NULL,
       contrast_matrix = NULL, var_shrinkage = TRUE,
       verbose = TRUE, sort = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitted_model}{
  The output from fitModel() function.
  }
  \item{coef}{
  A phenotype name, e.g. "disease",
  or a vector of contrast terms, e.g. c("disease", "case", "control").
  }
  \item{cell_type}{
  A cell type name, e.g. "celltype1", or "neuron".
  If cell_type is NULL or specified as "ALL",
  compound effect of coef in all cell types will be tested.
  }
  \item{contrast_matrix}{
  If contrast_matrix is specified,
  coef and cell_type will be ignored!
  A matrix (or a vector) to specify contrast, e.g.,
  cmat <- matrix(0, 2, 6); cmat[1,3] <- 1: cmat[2,4] <- 1
  is to test whether the 3rd parameter and 4th parameter
  are zero simultaneously i.e. beta3 = beta4 = 0. }
  \item{var_shrinkage}{
  Whether to apply shrinkage on estimated MSE or not.
  Applying shrinkage helps remove extremely small variance
  estimation and stablize statistics.}
  \item{verbose}{
  A boolean parameter. Testing information will be printed
  if verbose = TRUE.
  }
  \item{sort}{
  A boolean parameter.
  The output results will be sorted by p value if sort = TRUE.
  }
}
\value{
A matrix including the results from
testing the phenotype in specified cell type(s).
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
N <- 300 # simulation a dataset with 300 samples
K <- 3 # 3 cell types
P <- 500 # 500 features

### simulate proportion matrix
Prop <- matrix(runif(N*K, 10,60), ncol=K)
Prop <- sweep(Prop, 1, rowSums(Prop), FUN="/")
colnames(Prop) <- c("Neuron", "Astrocyte", "Microglia")

### simulate phenotype names
design <- data.frame(disease=factor(sample(0:1,
                             size = N,replace=TRUE)),
                     age=round(runif(N, 30,50)),
                     race=factor(sample(1:3, size = N,replace=TRUE)))
Y <- matrix(rnorm(N*P, N, P), ncol = N)

### generate design matrix and fit model
Design_out <- makeDesign(design, Prop)
fitted_model <- fitModel(Design_out, Y)

### check the names of cell types and phenotypes
fitted_model$all_cell_types
fitted_model$all_coefs

### detect age effect in neuron
test <- csTest(fitted_model, coef = "age",
cell_type = "Neuron", contrast_matrix = NULL)

## coef can be specified in different ways:
#### jointly test a phenotype:
test <- csTest(fitted_model, coef = "age",
cell_type = "joint", contrast_matrix = NULL)

#### if I do not specify cell_type
test <- csTest(fitted_model, coef = "age",
cell_type = NULL, contrast_matrix = NULL)
## this is exactly the same as
test <- csTest(fitted_model, coef = "age",
contrast_matrix = NULL)

#### other examples
test <- csTest(fitted_model, coef = "race",
cell_type = "Astrocyte", contrast_matrix = NULL)
test <- csTest(fitted_model, coef = "age",
cell_type = "Microglia", contrast_matrix = NULL)

#### specify contrast levels
test <- csTest(fitted_model, coef = c("race", 3, 2),
cell_type = "Neuron", contrast_matrix = NULL)
#### specify contrast levels in all cell types
test <- csTest(fitted_model, coef = c("race", 3, 2),
cell_type = "joint", contrast_matrix = NULL)

#### csTest can tolerate different ways of specifying contrast level
#### note race=1 is used as reference when fitting model
#### we can here specify race=2 as reference
test <- csTest(fitted_model, coef = c("race", 1, 2),
cell_type = "Neuron", contrast_matrix = NULL)
## get exactly the same results as
test <- csTest(fitted_model, coef = c("race", 2, 1),
cell_type = "Neuron", contrast_matrix = NULL)

#### specify a contrast matrix:
cmatrix = rep(0,15)
cmatrix[c(4,5)] = c(1,-1)
test <- csTest(fitted_model, coef = NULL,
cell_type = NULL, contrast_matrix = cmatrix)
#### specific a contrast matrix with two rows:
cmatrix = matrix(rep(0,30),2,15)
cmatrix[1,4] = 1
cmatrix[2,5] = 1
test <- csTest(fitted_model, coef = NULL,
contrast_matrix = cmatrix)
}
