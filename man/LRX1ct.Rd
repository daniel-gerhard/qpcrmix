\name{LRX1ct}
\alias{LRX1ct}
\docType{data}
\title{Relative expression of LRX1 in transgenic roots under phosphorus deficiency.}
\description{Roots of a wildtype as well as transgenic roots are
  observed under standard nutritional conditions and under phosphorus
  deficiency. The relative expression of LRX1 in with K7GWIWG2(II)::LRX1
generated roots are compared to the wildtype and to a standard gene AtUBC-2.}
\usage{data(LRX1ct)}
\format{
  A data frame with 72 observations on the following 6 variables.
  \describe{
    \item{\code{Well}}{a factor with well IDs}
    \item{\code{Content}}{factor levels denoting biological replicates}
    \item{\code{Target}}{2 different genes LRX1 100bp and At-UBC2 as control}
    \item{\code{P}}{+/- Phosphorus}
    \item{\code{WT}}{Wildtype yes/no}
    \item{\code{ct}}{a numeric response variable with estimated threshold cycles}
  }
}
\author{Martin Hinrichs <martin_hinrichs@t-online.de>}


\examples{
data(LRX1ct)
qpcr(data=LRX1ct, response="ct",
     gene="Target", control_gene="At-UBC2",
     fixed1="P", fixed2 ="WT",
     rep_id="Content", contrasts=TRUE)
}

\keyword{datasets}
