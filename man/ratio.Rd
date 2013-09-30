\name{ratio}
\alias{ratio}
\docType{data}
\title{Expression of the gene diazepam binding inhibitor (DBI)}
\description{Two-way ANOVA in Randomized complete block
  design. Quantitative RT-PCR was used to study expression of the gene
  diazepam binding inhibitor (DBI) in the brain of piglets subject to
  weaning and social isolation treatments. The experimental layout
followed a randomized complete block design (n = 3 litters) and the
treatments consisted of a 2 x 2 factorial combination of weaning
(early-weaned or non-weaned) and social isolation (isolated or
control). Preliminary assays indicated that Sus scrofa 18S ribosomal RNA
(18S) was suitable for use as an endogenous control gene and that the
amplification efficiency for primers of the two genes (18S and DBI)
was close to 2. All reactions were performed in triplicate but
some observations were excluded from the analysis because of
evidence of non-specific amplifications (as revealed by dissociation
curve analyses).}
\usage{data(ratio)}
\format{
  A data frame with 134 observations on the following 8 variables.
  \describe{
    \item{\code{Well}}{a factor}
    \item{\code{Sample}}{a factor}
    \item{\code{gene}}{factor with gene IDs: DBI or 18S as a control.}
    \item{\code{y}}{efficiency adjusted Ct levels}
    \item{\code{TRT}}{treatment factor (combination of weaning
(early-weaned or non-weaned) and social isolation (isolated or
control))}
    \item{\code{litter}}{the experimental unit (blocks)}
    \item{\code{wean}}{weaning treatment}
    \item{\code{iso}}{isolation treatment}
  }
}
\references{
Steibel, JP, Poletto, R, Coussens, PM, Rosa, GJM (2009): A powerful and
flexible linear mixed model framework for the analysis of relative
quantification RT-PCR data. Genomics 94:146--152.
}

\examples{
data(ratio)
qratio <- qpcr(data=ratio, response="y",
               gene="gene", control_gene="18S",
               fixed1="TRT", fixed2=NULL,
               rep_id="litter", block=TRUE,
               splitplot=FALSE,
               contrasts=TRUE, interaction=FALSE, adjusted=FALSE)

# model
qratio$model
# testing deltaCt contrasts
qratio$pvalues 
}


\keyword{datasets}
