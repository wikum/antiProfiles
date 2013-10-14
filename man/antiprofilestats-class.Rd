\name{AntiProfileStats-class}
\alias{AntiProfileStats-class}
\alias{getProbeStats}
\alias{getProbeStats,AntiProfileStats-method}
\title{Statistics used to build anti-profiles}
\usage{
  ## Accessors 
  getProbeStats(object)
}

\arguments{
  \item{object}{Object of class
  \code{\linkS4class{AntiProfileStats}}}
}
\description{
  This class stores statistics required to build
  anti-profiles. Objects of this class should be built
  using the \code{\link{apStats}} function. To build
  anti-profiles from objects of this class, use the
  \code{\link{buildAntiProfile}} function.
}

\section{Accessors}{
  In the following code \code{object} is a
  \code{AntiProfileStats} object. \describe{
  \item{\code{getProbeStats}:}{Gets probeset statistics. A
  \code{data.frame} with columns \describe{
  \item{\code{affyid}:}{Affymetrix probeset id (character)}
  \item{\code{SD0}:}{Normal expression standard deviation
  aggregated over tissue types (numeric)}
  \item{\code{SD1}:}{Cancer expression standard deviation
  aggregated over tumor types (numeric)}
  \item{\code{stat}:}{The log2-variance ratio statistic
  (numeric)} \item{\code{meds0}:}{Median normal expression
  aggregated over tissue types (numeric)}
  \item{\code{mads0}:}{Median absolute deviation of normal
  expression aggregate over tissue types (numeric)} } } }
}
\examples{
if (require(antiProfilesData)) {
   data(apColonData)
   colonStats = apStats(exprs(apColonData), pData(apColonData)$Status)
   show(getProbeStats(colonStats))
 }
}
\author{
  Hector Corrada Bravo \email{hcorrada@gmail.com}
}
\seealso{
  \code{\link{apStats}} to construct objects of this class,
  \code{\link{buildAntiProfile}} to build anti-profiles
  from objects of this class.
}

