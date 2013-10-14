\docType{methods}
\name{apReorder}
\alias{apReorder}
\alias{apReorder,AntiProfileStats,numeric-method}
\title{Reorder an AntiProfileStats object}
\arguments{
  \item{stats}{An object of class AntiProfileStats}

  \item{o}{A numeric vector giving new probe ordering}
}
\value{
  A reordered AntiProfileStats object
}
\description{
  Reorders given AntiProfileStats object using provided
  ordering o
}
\examples{
if (require(antiProfilesData)) {
  data(apColonData)
  colonStats = apStats(exprs(apColonData), pData(apColonData)$Status)
  o = sample(seq(len=nrow(slot(colonStats,"probes"))))
  newStats = apReorder(colonStats, o)
}
}
\author{
  Hector Corrada Bravo \email{hcorrada@gmail.com}
}

