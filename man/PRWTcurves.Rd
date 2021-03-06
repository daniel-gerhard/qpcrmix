\name{PRWTcurves}
\alias{PRWTcurves}
\docType{data}
\title{PRWTcurves}
\description{Dataset from the Institute of plant nutrition (LUH). In the
  experiments an insertion has introduced into the genome of different
  rice lines. This insertion changes gene expression from the gene of
  interest. By using the quantitative RT-PCR the changes in gene
  expression of a target gene will analyse between different samples and
  the wild type. 
  Content: ID of biological replicates
  Target: gene of interest and EF (elongation factor) gene of control
  Samples: rice lines with an insertion (different letters) and wild type (W)}
\usage{data(PRWTcurves)}
\format{
  A data frame with 3780 observations on the following 6 variables.
  \describe{
    \item{\code{Well}}{a factor with well IDs}
    \item{\code{Content}}{factor levels denoting biological replicates}
    \item{\code{Target}}{3 different genes with EF as control}
    \item{\code{Sample}}{treatments P, R, WT}   
    \item{\code{Cycle}}{cycle number} 
    \item{\code{Fluorescence}}{observed fluorescence (proportional to PCR material)}
  }
}
\author{Eline Biedermann <eline.biedermann@gmx.de>}

\examples{
data(PRWTcurves)
ggplot(PRWTcurves, aes(x=Cycle, y=Fluorescence, colour=Sample,
group=Well)) + geom_line() + facet_wrap(~Target)
}

\keyword{datasets}
