\name{apStats}
\alias{apStats}
\title{Compute statistics used to construct antiProfile}
\usage{
  apStats(e, cl, tiss = NULL, minL = 10, cutoff = 5,
    OnCutoff = 2.54)
}
\arguments{
  \item{e}{matrix of gene expression, with one column per
  sample}

  \item{cl}{vector of normal/cancer indicators as 0/1}

  \item{tiss}{vector of tissue types for each sample}

  \item{minL}{minimum number of samples of a given
  tissue/class to compute stats}

  \item{cutoff}{median absolute deviation multiplier used
  to determine proportion of samples within normal range of
  expression}

  \item{OnCutoff}{gene expression barcode z-score to
  determine if a gene is expressed}
}
\value{
  An object of class SuccsStats
}
\description{
  This function calculates normal ranges of expressions and
  variance ratios for all probesets. To create an
  anti-profile, call buildAntiProfile on the output of this
  function
}
\examples{
if (require(antiProfilesData)) {
   data(apColonData)
   colonStats = apStats(exprs(apColonData), pData(apColonData)$Status)
 }
}
\author{
  Hector Corrada Bravo \email{hcorrada@gmail.com}
}
\seealso{
  \code{\linkS4class{AntiProfileStats}} for the type of
  object returned. \code{\link{buildAntiProfile}} to
  construct anti-profiles with objects returned by this
  function.
}

