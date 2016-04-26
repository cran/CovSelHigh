\name{cov.sel.high.rf}
\alias{cov.sel.high.rf}
\title{cov.sel.high.rf}
\description{
Function called by cov.sel.high if type="rf". Not meant to be used on its own.
}
\usage{cov.sel.high.rf(Y, X, threshold = 0.25, ntree = 1000, ...)
}
\arguments{

  \item{Y  }{Outcome variable or treatment variable.}
  \item{X  }{A matrix or data frame containing columns of covariates.}
\item{threshold}{Variable importance threshold, see Value.}
\item{ntree}{Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.}
 \item{\dots  }{Additional arguments passed on to \code{randomForest}.}
}
\details{See \code{randomForest}.
}

\value{
 \code{cov.sel.high.rf} returns a logical vector of the same length as the number of columns in X. The positions of values in the vector refers to the covariates in the corresponding X columns. Value TRUE implies that the covariate has a variable importance value of more than \code{threshold}*the largest observed variable importance value.
}


\author{Jenny Häggström, <jenny.haggstrom@umu.se> }

\seealso{\code{\link[randomForest]{randomForest}}}

