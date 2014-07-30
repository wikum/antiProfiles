#' Compute expression variability measure
#' 
#' This function computes expression variability in a way that
#' removes dependence on mean expression.
#' It uses a local polynomial likelihood method to estimate variance
#' as gamma distributed around given mean expression for each probeset.
#' This function makes this calculation using all samples in argument. To
#' calculate expression variability for samples in different groups, call this
#' function for each subset of columns separately.
#'
#' @param x matrix of gene expression, with one column per sample
#' @param cutoff minimum expression value to be included in computation (for \code{frma} normalized data, we find 2.54 to be a good value for determining if a probeset is expressed in a given sample (default NULL)
#' @param plot make a plot of local likelihood model using \code{smoothScatter} (default=FALSE)
#' @param ... arguments passed to \code{smoothScatter}
#'
#' @return numeric vector of length equal to number of rows of \code{x}
#'
#' @examples
#' if (require(antiProfilesData)) {
#'   data(apColonData)
#'   e <- exprs(apColonData)[,pData(apColonData)$Status==1]
#'   ev <- ev(e, cutoff=2.54)
#' }
#'
#' @author Hector Corrada Bravo \email{hcorrada@@gmail.com}
#' @seealso \code{frma} for normalization
#'
#' @references E. Alemu, H. Corrada Bravo, S. Hannenhalli (2014). Determinants of Expression Variability. Nucleic Acids Research, 42 (6), 3503-14.
#' @export
ev <- function(x, cutoff=NULL, plot=FALSE, ...) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("argument 'x' must be a numeric matrix") 
  }

  if (!is.null(cutoff)) {
    if (!is.numeric(cutoff)) {
      stop("argument 'cutoff' must be numeric")
    }

    x[x<cutoff] <- NA
  }
  
  mns <- rowMeans(x, na.rm=TRUE)
  sds <- matrixStats::rowSds(x, na.rm=TRUE)

  drop <- is.na(sds) & is.na(mns)
  mns <- mns[!drop]
  sds <- sds[!drop]
  
  fit <- locfit(sds^2 ~ lp(mns), family="gamma")
  expSd <- sqrt(predict(fit, mns))

  ev <- rep(NA, nrow(x))
  ev[!drop] <- log2(sds) - log2(expSd)

  if (plot) {
    smoothScatter(mns, sds, xlab="mean expression", ylab="std. dev. expression", ...)
    f1 <- function(x) sqrt(predict(fit,x))
    curve(f1, from=min(mns), to=max(mns), col="red", add=TRUE)
  }
  ev
}
