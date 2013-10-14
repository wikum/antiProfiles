\name{TissueSpecAntiProfile-class}
\alias{TissueSpecAntiProfile-class}
\alias{getNormalRegions,TissueSpecAntiProfile-method}
\alias{getNormalTissueRegions}
\alias{getNormalTissueRegions,TissueSpecAntiProfile-method}
\alias{getProbesetIds,TissueSpecAntiProfile-method}
\title{A gene expression anti-profile using tissue-specific regions}
\usage{
  ## Accessors

  \S4method{getProbesetIds}{TissueSpecAntiProfile}(object) 
  \S4method{getNormalRegions}{TissueSpecAntiProfile}(object)
  getNormalTissueRegions(object)
}

\arguments{
  \item{object}{Object of class
  \code{\linkS4class{TissueSpecAntiProfile}}}
}
\description{
  This class encapsulates anti-profiles with
  tissue-specific normal expression regions. Objects of
  this class should be built from
  \code{\linkS4class{AntiProfileStats}} objects using the
  \code{\link{buildAntiProfile}} method. Anti-profile
  scores can be computed for new samples using the
  \code{\link{apCount}} method.
}

\section{Accessors}{
  In the following code \code{object} is a
  \code{TissueSpecAntiProfile} object. \describe{
  \item{\code{getProbesetIds}:}{vector of Affymetrix ids
  for probesets included in the anti-profile (character)}
  \item{\code{getNormalRegions}:}{median and upper boundary
  of normal expression regions (numeric matrix of dimension
  s-by-2, where s is the size of the anti-profile)}
  \item{\code{getNormalTissueRegions}:}{median and upper
  boundary of normal expression regions (numeric array of
  dimension s-by-2-by-t, where s is the size of the
  anti-profile, and t the number of normal tissues used in
  the anti-profile)} }
}
\examples{
if (require(antiProfilesData)) {
    data(apColonData)
    # fake tissues
    tissue=rep(c("colon","lung"), len=length(sampleNames(apColonData)))
    tissStats=apStats(exprs(apColonData), pData(apColonData)$Status, tiss=tissue, minL=3)
    tissAP=buildAntiProfile(tissStats, sigsize=10)
    show(tissAP)

    head(getProbesetIds(tissAP))
    head(getNormalRegions(tissAP))
    head(getNormalTissueRegions(tissAP))
  }
}
\author{
  Hector Corrada Bravo \email{hcorrada@gmail.com}
}
\seealso{
  \code{\linkS4class{AntiProfileStats}} for the class of
  objects from which anti-profiles are built.
  \cite{\link{buildAntiProfile}} for the method used to
  construct objects of this class. \code{\link{apCount}}
  for the function used to calculate anti-profile scores
  from objects of this class.  #'
}

