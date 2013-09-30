\name{BJWTct}
\alias{BJWTct}
\docType{data}
\title{BJWTct}
\description{Dataset from the Institute of plant nutrition (LUH). In the
  experiments an insertion has introduced into the genome of different
  rice lines. This insertion changes gene expression from the gene of
  interest. By using the quantitative RT-PCR the changes in gene
  expression of a target gene will analyse between different samples and
  the wild type. 
  Content: ID of biological replicates
  Target: gene of interest and EF (elongation factor) gene of control
  Samples: rice lines with an insertion (different letters) and wild type (W)}
\usage{data(BJWTct)}
\format{
  A data frame with 84 observations on the following 5 variables.
  \describe{
    \item{\code{Well}}{a factor with well IDs}
    \item{\code{Content}}{factor levels denoting biological replicates}
    \item{\code{Target}}{3 different genes with Ef as control}
    \item{\code{Sample}}{3 treatments B, J, WT}
    \item{\code{Ct}}{a numeric response variable with estimated threshold cycles}
  }
}
\author{Eline Biedermann <eline.biedermann@gmx.de>}

\examples{
## only designs can be analyzed with qpcr()
## when all combinations
##   of factor levels and gene levels are available
##
## Here:
## - gene 06g16350 is not observed for Sample B
## - gene 01g67540 is not observed for Sample J
##
## Therefore the dataset is split into 2 parts:

data(BJWTct)

## B vs WT
BWT <- subset(BJWTct, Target != "06g16350" & Sample != "J")
BWT <- droplevels(BWT)
qpcr(data=BWT, response="Ct", gene="Target", control_gene="EF",
     fixed1="Sample", rep_id="Content", contrasts=TRUE)

## J vs WT
JWT <- subset(BJWTct, Target != "01g67540" & Sample != "B")
JWT <- droplevels(JWT)
qpcr(data=JWT, response="Ct", gene="Target", control_gene="EF",
     fixed1="Sample", rep_id="Content", contrasts=TRUE)
}

\keyword{datasets}
