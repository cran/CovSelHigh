\name{cov.sel.high.lasso}
\alias{cov.sel.high.lasso}
\title{cov.sel.high.lasso}
\description{
Function called by cov.sel.high if type="lasso". Not meant to be used on its own.
}
\usage{cov.sel.high.lasso(Y, X, minscreen = 2, ...)
}
\arguments{

  \item{Y  }{Outcome variable or treatment variable.}
  \item{X  }{A matrix or data frame containing columns of covariates and all functions of covariates e.g. interactions that should be included in the lasso model.}
\item{minscreen}{The minimum number of columns in X that should be selected.}
 \item{\dots  }{Additional arguments passed on to \code{cv.glmnet}.}
}
\details{See \code{cv.glmnet}.
}

\value{
 \code{cov.sel.high.lasso} returns a logical vector of the same length as the number of columns in X. The positions of values in the vector refers to the (functions of) covariates in the corresponding X columns. Value TRUE implies that (the function of) the covariate has a corresponding coefficient not equal to zero.}


\author{Jenny Häggström, <jenny.haggstrom@umu.se> }

\seealso{\code{\link[glmnet]{cv.glmnet}}}

