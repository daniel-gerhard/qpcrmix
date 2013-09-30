\name{PRWTct}
\alias{PRWTct}
\docType{data}
\title{PRWTct}
\description{Dataset from the Institute of plant nutrition (LUH). In the
  experiments an insertion has introduced into the genome of different
  rice lines. This insertion changes gene expression from the gene of
  interest. By using the quantitative RT-PCR the changes in gene
  expression of a target gene will analyse between different samples and
  the wild type. 
  Content: ID of biological replicates
  Target: gene of interest and EF (elongation factor) gene of control
  Samples: rice lines with an insertion (different letters) and wild type (W)}
\usage{data(PRWTct)}
\format{
  A data frame with 84 observations on the following 5 variables.
  \describe{
    \item{\code{Well}}{a factor with well IDs}
    \item{\code{Content}}{factor levels denoting biological replicates}
    \item{\code{Target}}{3 different genes with EF as control}
    \item{\code{Sample}}{3 treatments P, R, WT}
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
## - gene 10g30610 is not observed for Sample P
## - gene 08g02110 is not observed for Sample R
##
## Therefore the dataset is split into 2 parts:

data(PRWTct)

## P vs WT
PWT <- subset(PRWTct, Target != "10g30610" & Sample != "R")
PWT <- droplevels(PWT)
qpcr(data=PWT, response="Ct", gene="Target", control_gene="EF",
     fixed1="Sample", rep_id="Content", contrasts=TRUE)

## R vs WT
RWT <- subset(PRWTct, Target != "08g02110" & Sample != "P")
RWT <- droplevels(RWT)
qpcr(data=RWT, response="Ct", gene="Target", control_gene="EF",
     fixed1="Sample", rep_id="Content", contrasts=TRUE)
}

\keyword{datasets}
