\name{RA_100samples}
\alias{RA_100samples}
\docType{data}
\title{
An example dataset for cellular proportion
estimation and multiple factor design
}
\description{
The dataset contains normalized beta valutes
for 3000 CpGs from 100 samples (50 Rheumatoid
arthritis patients and 50 controls) and their
phenotypes (disease status, age, and gender).
The dataset also contains a sub-setted blood
reference matrix for the matched 3000 CpGs.
This data was obtained and processed based
on GSE42861.
}
\usage{data("RA_100samples")}
\references{
Liu Y, Aryee MJ, Padyukov L, Fallin MD et al.
Epigenome-wide association data implicate DNA
methylation as an intermediary of genetic risk
in rheumatoid arthritis. Nat Biotechnol
2013 Feb;31(2):142-7. PMID: 23334450
}
\examples{
data(RA_100samples)
RA_100samples$Y_raw[1:5,1:5]
head(RA_100samples$Pheno)
head(RA_100samples$Blood_ref)
}
\keyword{datasets}
