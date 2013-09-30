\name{ptb}
\alias{ptb}
\docType{data}
\title{Paratuberculosis of Cows}
\description{Two-way factorial in a randomized Split-Plot design.}
\usage{data(ptb)}
\format{
  A data frame with 96 observations on the following 7 variables.
  \describe{
    \item{\code{Well}}{a factor}
    \item{\code{Type}}{a factor denoting the paratuberculosis status
      (Clinical, Subclinical, or Control)}
    \item{\code{Name}}{a factor indicating the experimental unit (an
      individual cow)}
    \item{\code{group}}{a factor referring to PTB challange (or not) of
      blood cells from each cow}
    \item{\code{Detector}}{a factor}
    \item{\code{Ct}}{a numeric response variable with cycle numbers as
      observations}
    \item{\code{StdDev}}{Standard Deviation}
  }
}
\references{
Coussens P. M., Colvin C. J., Rosa G. J., Perez Laspiur J., and Elftman
M. D., Evidence for a novel gene expression program in peripheral blood
mononuclear cells from Mycobacterium avium subsp. paratuberculosis-infected cattle, Infect Immun 71 (2003) 6487-98.
}

\examples{
data(ptb)
qptb <- qpcr(data=ptb, response="Ct",
             gene="Detector", control_gene="BetaActin",
             fixed1="Type", fixed2="group",
             rep_id="Name",
             block=FALSE, splitplot=TRUE,
             contrasts=TRUE, interaction=FALSE, adjusted=FALSE)

# model
qptb$model
# testing deltaCt contrasts
qptb$pvalues 
}

\keyword{datasets}
