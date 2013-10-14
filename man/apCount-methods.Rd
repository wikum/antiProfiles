\docType{methods}
\name{apCount}
\alias{apCount}
\alias{apCount,AntiProfile,matrix-method}
\alias{apCount,TissueSpecAntiProfile,matrix-method}
\title{Obtain the anti-profile score for a set of samples}
\arguments{
  \item{fit}{an object of class AntiProfile as produced by
  the buildAntiProfile method}

  \item{expr}{a matrix of gene expression, rownames are
  used as identifiers}
}
\value{
  a numeric vector of anti-profile scores
}
\description{
  This function applies the given anti-profile to a new set
  of samples. Rownames in the expression matrix are used to
  match probenames in the AntiProfile object.
}
\examples{
if (require(antiProfilesData)) {
    data(apColonData)

    # compute statistics
    colonStats = apStats(exprs(apColonData), pData(apColonData)$Status)

    # create an anti-profile, ignoring tissue-specificity of probesets, with 10 probesets
    ap = buildAntiProfile(colonStats, tissueSpec=FALSE, sigsize=10)

    # get counts for the original dataset
    counts =apCount(ap, exprs(apColonData))
 }
}
\author{
  Hector Corrada Bravo \email{hcorrada@gmail.com}
}

