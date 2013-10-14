\name{AntiProfile-class}
\alias{AntiProfile-class}
\alias{getNormalRegions}
\alias{getNormalRegions,AntiProfile-method}
\alias{getProbesetIds}
\alias{getProbesetIds,AntiProfile-method}
\title{A gene expression anti-profile}
\usage{
  ## Accessors

  getProbesetIds(object) 
  getNormalRegions(object)
}

\arguments{
  \item{object}{Object of class (or inheriting from)
  \code{\linkS4class{AntiProfile}}}
}
\description{
  This class encapsulates anti-profiles. Objects of this
  class should be built from
  \code{\linkS4class{AntiProfileStats}} objects using the
  \code{\link{buildAntiProfile}} method. Anti-profile
  scores can be computed for new samples using the
  \code{\link{apCount}} method.
}
\section{Accessors}{
  In the following code \code{object} is a
  \code{AntiProfile} object. \describe{
  \item{\code{getProbesetIds}:}{vector of Affymetrix ids
  for probesets included in the anti-profile (character)}
  \item{\code{getNormalRegions}:}{median and upper boundary
  of normal expression regions (numeric matrix of dimension
  s-by-2, where s is the size of the anti-profile)} }
}
\examples{
if (require(antiProfilesData)) {
    data(apColonData)
    colonStats=apStats(exprs(apColonData), pData(apColonData)$Status)
    colonAP = buildAntiProfile(colonStats, tissueSpec=FALSE, sigsize=10)
    show(colonAP)

    head(getProbesetIds(colonAP))
    head(getNormalRegions(colonAP))
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
  from objects of this class.
}

