\name{IPE}
\alias{IPE}
\docType{data}
\title{IPE dataset}
\description{Dataset from the Institute of plant nutrition (LUH). Gene
  XTH should be upregulated at phosphate deficiency.}
\usage{data(IPE)}
\format{
  A data frame with 72 observations on the following 5 variables.
  \describe{
    \item{\code{Well}}{a factor with well IDs}
    \item{\code{Content}}{ID of biological replicates}
    \item{\code{Target}}{2 different genes: UBQ (control) and XTH (gene
      of interest)}
    \item{\code{Sample}}{treatment factor}
    \item{\code{Ct}}{a numeric response variable with estimated threshold cycles}
  }
}
\author{Melanie Bremer <bremer@pflern.uni-hannover.de>}

\examples{
data(IPE)
qipe <- qpcr(data=IPE, response="Ct",
             gene="Target", control_gene="UBQ",
             fixed1="Sample",
             rep_id="Content",
             block=FALSE, splitplot=FALSE,
             contrasts=TRUE, interaction=FALSE, adjusted=TRUE)

# model
qipe$model
# testing deltaCt contrasts
qipe$pvalues 
}

\keyword{datasets}
