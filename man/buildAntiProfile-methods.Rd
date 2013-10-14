\docType{methods}
\name{buildAntiProfile}
\alias{buildAntiProfile}
\alias{buildAntiProfile,AntiProfileStats-method}
\title{Create an anti-profile from a AntiProfileStats object}
\arguments{
  \item{stats}{an object of class AntiProfileStats as
  produced by the apStats function}

  \item{tissueSpec}{use tissue-specific regions of normal
  expression}

  \item{tissueFilter}{use only tissue-specific genes in the
  anti-profile}

  \item{sigsize}{desired size of signature, if NULL,
  computed from statCutoff}

  \item{cutoff}{median absolute deviation multiplier used
  to define normal regions of expression}

  \item{statCutoff}{cutoff used to include probesets in
  anti-profile}
}
\value{
  an object of class AntiProfilesSig or
  AntiProfilesTissueSig depending on the tissueSpec
  argument
}
\description{
  This function creates anti-profile using statistics
  stored in a AntiProfileStats object
}
\examples{
if (require(antiProfilesData)) {
  # create an anti-profile, ignoring tissue-specificity of probesets, with 10 probesets
  data(apColonData)
  colonStats = apStats(exprs(apColonData), pData(apColonData)$Status)
  ap = buildAntiProfile(colonStats, tissueSpec=FALSE, sigsize=10)
}
}
\author{
  Hector Corrada Bravo \email{hcorrada@gmail.com}
}

