\name{HRGPcurves}
\alias{HRGPcurves}
\docType{data}
\title{Fluorescence curves for HRGP in transgenic roots under phosphorus deficiency.}
\description{Roots of a wildtype as well as transgenic roots are
  observed under standard nutritional conditions and under phosphorus
  deficiency. The relative expression of HRGP in with pK7GWIWG2(II)::HRGP
generated roots are compared to the wildtype and to a standard gene AtUBC-2.}
\usage{data(HRGPcurves)}
\format{
  A data frame with 3240 observations on the following 7 variables.
  \describe{
    \item{\code{Well}}{a factor with well IDs}
    \item{\code{Content}}{factor levels denoting biological replicates}
    \item{\code{Target}}{2 different genes HRGP and At-UBC2 as control}
    \item{\code{P}}{+/- Phosphorus}
    \item{\code{WT}}{Wildtype yea/no} 
    \item{\code{Cycle}}{cycle number} 
    \item{\code{Fluorescence}}{observed fluorescence (proportional to PCR material)}
  }
}
\author{Martin Hinrichs <martin_hinrichs@t-online.de>}

\examples{
data(HRGPcurves)
ggplot(HRGPcurves, aes(x=Cycle, y=Fluorescence, colour=P:WT,
group=Well)) + geom_line() + facet_wrap(~Target)
}

\keyword{datasets}
