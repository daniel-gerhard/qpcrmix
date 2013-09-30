\name{qpcr}
\alias{qpcr}
\title{Mixed model analysis of qpcr data}
\description{Fits a linear mixed model to qpcr data, allowing for
  several genes, two-factorial block and split-plot designs. Delta-Ct
  and DeltaDelta-Ct values are calculated by multiple contrasts.}
\usage{
qpcr(data, response, gene, control_gene, fixed1, fixed2 = NULL, rep_id,
block = FALSE, splitplot = FALSE, contrasts=TRUE, interaction=FALSE,
adjusted = FALSE, df=NULL)
}
\arguments{
  \item{data}{A data.frame object.}
  \item{response}{Name of a column in data, denoting a numeric response
    vector, usually estimated Ct values.}
  \item{gene}{Name of a column in data, denoting a factor coding for
    different genes.}
  \item{control_gene}{Name of a single control gene in the vector of
    gene names.}
  \item{fixed1}{Name of a column in data, denoting a factor coding for
    several treatment groups}
  \item{fixed2}{Name of a column in data, denoting a second factor coding for
    several treatment groups, default is NULL}
  \item{rep_id}{Name of a column in data, denoting a factor which
    identifies a biological replication or block.}
  \item{block}{Logical; are the rep_ids denoting a block variable?
    default is FALSE.}
  \item{splitplot}{Logical; Is the experiment designed as a split-plot
    with fixed1 as mainplots and fixed2 as subplot? default is FALSE.}
  \item{contrasts}{Logical; If TRUE: Tests for general linear hypotheses
    are performed.}
  \item{interaction}{Logical; If TRUE: Tests for interaction between fixed1 and fixed2 are performed if contrasts==TRUE.}
  \item{adjusted}{Logical; Should multiplicity adjusted p-values be
    reported? default is FALSE.}
  \item{df}{Providing degrees of freedom as a numeric value, an input of
  NULL corresponds to the assumption of a Normal- instead of a t-distribution.}
}
\value{An object of class qpcrmacro.}
\references{
Steibel, JP, Poletto, R, Coussens, PM, Rosa, GJM (2009): A powerful and
flexible linear mixed model framework for the analysis of relative
quantification RT-PCR data. Genomics 94:146--152.
}
\author{Daniel Gerhard <gerhard@biostat.uni-hannover.de>}


\seealso{\code{\link{lmer}}, \code{\link{glht}}}

\keyword{ htest }
