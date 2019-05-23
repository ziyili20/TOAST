# Analyses of high-throughput data from heterogeneous samples with TOAST

`TOAST` is an R package designed for the analyses of high-throughput data from complex, heterogeneous tissues. It is designed for the analyses of high-throughput data from 
  heterogeneous tissues,
  which is a mixture of different cell types.  
  
  TOAST offers functions for detecting cell-type 
  specific differential expression (csDE) or 
  differential methylation (csDM) for microarray data,
  and improving reference-free deconvolution 
  based on cross-cell type differential analysis. 
  TOAST implements a rigorous staitstical framework, 
  based on linear model, which provides great 
  flexibility for csDE/csDM detection and 
  superior computationl performance. 
  
  In this readme file, we briefly present how to install TOAST package through GitHub. For detailed usage of TOAST, please refer to the vignette file.

## Installation and quick start

### Install TOAST


In R, install the `TOAST` package by

```{r install, message=FALSE, warning=FALSE}
library(devtools)
install_github("ziyili20/TOAST", build_vignettes=TRUE)
 
```

To view the package vignette in pdf format, run the following lines in R

```{r vig, message=FALSE, warning=FALSE}
library(TOAST)
vignette("TOAST")
```
The content in this README file is essentially the same as the package vignette.

### How to get help for TOAST

Any TOAST questions should be posted
to the GitHub Issue section of TOAST 
homepage at https://github.com/ziyili20/TOAST/issues.

### Quick start on detecting cell type-specific differential signals

Here we show the key steps for a cell 
type-specific different analysis. This 
code chunk assumes you have an expression
or DNA methylation matrix called `Y_raw`,
a data frame of sample information called
`design`, and a table of cellular composition
information (i.e. mixing proportions) 
called `prop`. If the cellular composition
is not available, our vignette file 
provides discussions about how to obtain mixing 
proportions using reference-free deconvolution 
or reference-based deconvolution.

```{r quick_start, eval = FALSE}
Design_out <- makeDesign(design, Prop)
fitted_model <- fitModel(Design_out, Y_raw)
fitted_model$all_coefs # list all phenotype names
fitted_model$all_cell_types # list all cell type names
# coef should be one of above listed phenotypes
# cell_type should be one of above listed cell types
res_table <- csTest(fitted_model, coef = "age", 
                    cell_type = "Neuron", contrast_matrix = NULL)
head(res_table)
```
**For detailed usage of TOAST, please refer to the vignette file through**

```{r vignette}
vignette("TOAST")
# or
browseVignettes("TOAST")
```
