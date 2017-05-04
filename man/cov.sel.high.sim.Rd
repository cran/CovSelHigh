\name{cov.sel.high.sim}
\alias{cov.sel.high.sim}
\title{Simulate Example Data for CovSelHigh }
\description{
Function used internally by cov.sel.high to simulates example data used in Häggström (2016).
}
\usage{cov.sel.high.sim(N, Setting, rep, Models)
	
}
\arguments{
\item{N}{The number of observations to be simulated.}
\item{Setting}{The simulation setting to be used. Unconfoundedness holds given X if Setting=1. M-bias given X if Setting=2.}
\item{rep}{The number of replications to be simulated.}
\item{Models}{The type of outcome models to be used, options are \code{"Linear"}, \code{"Nonlinear"} and \code{"Binary"}.}
}

\value{
 \code{cov.sel.high} returns a list with the following content:
  \item{dat }{A data frame with simulated data.}
}
\references{ 

Häggström, J. (2017). Data-Driven Confounder Selection via Markov and Bayesian Networks. \emph{ArXiv e-prints}.

}

\author{Jenny Häggström, <jenny.haggstrom@umu.se> }
\note{ \code{cov.sel.high.sim} calls the functions \code{\link[bindata]{rmvbin}} and \code{\link[MASS]{mvrnorm}}.
}
