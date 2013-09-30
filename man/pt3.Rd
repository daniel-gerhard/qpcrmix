\name{pt3}
\alias{pt3}
\docType{data}
\title{PT3 expression at phosphate deficiency}
\description{RT-PCR analysis of the phosphate transporter PT3 in rice at
  several time periods of phosphate deficiency.}
\usage{data(pt3)}
\format{
  A data frame with 90 observations on the following 5 variables.
  \describe{
    \item{\code{Well}}{a factor with well IDs}
    \item{\code{Content}}{factor levels denoting biological replicates}
    \item{\code{Target}}{2 different genes PT3 and eEf (control)}
    \item{\code{Time}}{0,2,4,6, and 8 hours in a nutrient solution at
      low phosphate content}
    \item{\code{Ct}}{a numeric response variable with estimated threshold cycles}
  }
}
\author{Melanie Bremer <bremer@pflern.uni-hannover.de>}

\examples{
data(pt3)
ptct <- qpcr(data=pt3, response="Ct",
             gene="Target", control_gene="eEf",
             fixed1="Time",
             rep_id="Content",
             block=FALSE, splitplot=FALSE,
             contrasts=TRUE, interaction=FALSE, adjusted=TRUE)

# model
ptct$model
# testing deltaCt contrasts
ptct$pvalues 
}

\keyword{datasets}
