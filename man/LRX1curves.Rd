\name{LRX1curves}
\alias{LRX1curves}
\docType{data}
\title{Fluorescence curves for LRX1 in transgenic roots under phosphorus deficiency.}
\description{Roots of a wildtype as well as transgenic roots are
  observed under standard nutritional conditions and under phosphorus
  deficiency. The relative expression of LRX1 in with K7GWIWG2(II)::LRX1
generated roots are compared to the wildtype and to a standard gene AtUBC-2.}
\usage{data(LRX1curves)}
\format{
  A data frame with 3240 observations on the following 7 variables.
  \describe{
    \item{\code{Well}}{a factor with well IDs}
    \item{\code{Content}}{factor levels denoting biological replicates}
    \item{\code{Target}}{2 different genes LRX1 100bp and At-UBC2 as control}
    \item{\code{P}}{+/- Phosphorus}
    \item{\code{WT}}{Wildtype yea/no} 
    \item{\code{Cycle}}{cycle number} 
    \item{\code{Fluorescence}}{observed fluorescence (proportional to PCR material)}
  }
}
\author{Martin Hinrichs <martin_hinrichs@t-online.de>}

\examples{
data(LRX1curves)
ggplot(LRX1curves, aes(x=Cycle, y=Fluorescence, colour=P:WT,
group=Well)) + geom_line() + facet_wrap(~Target)
}

\keyword{datasets}
