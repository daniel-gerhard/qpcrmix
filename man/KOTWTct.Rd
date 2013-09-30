\name{KOTWTct}
\alias{KOTWTct}
\docType{data}
\title{KOTWTct}
\description{Dataset from the Institute of plant nutrition (LUH). In the
  experiments an insertion has introduced into the genome of different
  rice lines. This insertion changes gene expression from the gene of
  interest. By using the quantitative RT-PCR the changes in gene
  expression of a target gene will analyse between different samples and
  the wild type. 
  Content: ID of biological replicates
  Target: gene of interest and EF (elongation factor) gene of control
  Samples: rice lines with an insertion (different letters) and wild type (W)}
\usage{data(KOTWTct)}
\format{
  A data frame with 84 observations on the following 5 variables.
  \describe{
    \item{\code{Well}}{a factor with well IDs}
    \item{\code{Content}}{factor levels denoting biological replicates}
    \item{\code{Target}}{2 different genes with Ef as control}
    \item{\code{Sample}}{3 treatments KO, T, WT}
    \item{\code{Ct}}{a numeric response variable with estimated threshold cycles}
  }
}
\author{Eline Biedermann <eline.biedermann@gmx.de>}

\examples{
data(KOTWTct)
qpcr(data=KOTWTct, response="Ct", gene="Target", control_gene="EF",
     fixed1="Sample", rep_id="Content", contrasts=TRUE)
}

\keyword{datasets}
