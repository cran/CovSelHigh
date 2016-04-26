\name{cov.sel.high.sim.res}
\alias{cov.sel.high.sim.res}
\title{Summarize Simulation Results for CovSelHigh }
\description{
Function used to summarize results from cov.sel.high when simulated data is used.
}
\usage{cov.sel.high.sim.res(object)
	
}
\arguments{
\item{object}{A list returned from cov.sel.high.}

}

\value{
 \code{cov.sel.high.sim.res} returns a list with the following content:
  \item{resmat }{A matrix with the following columns:
    \code{XTuc}
    \code{Quc}  
    \code{XYuc}  
    \code{Zuc}  
    \code{XTYuc}  
    \code{SinXT}  
    \code{SinQ}  
    \code{SinXY}
    \code{SinZ}
    \code{SinXTY}
    \code{XTeqS}
    \code{QeqS}
    \code{XYeqS}
    \code{ZeqS}
    \code{XTYeqS}
    \code{cards}
    \code{betahatest}
    \code{betahatse}
    \code{betahat_cicov}
    \code{ciL}
    \code{ciU}
  }
    \item{summary_resmat }{A list with the following columns:
    \code{Subset_selection}
    \code{Median_cardinality}  
    \code{Betahat_bias}  
    \code{Betahat_sd}  
    \code{Betahat_MSE}  
    \code{Betahat_CI_coverag}  
    \code{Betahat_mean_lower_C}  
    \code{Betahat_mean_upper_C}  
  }
  \item{xtable}{A LaTeX table}
}
\references{ 

Häggström, J. (2016). Data-Driven Confounder Selection via Markov and Bayesian Networks. \emph{ArXiv e-prints}.

}

\author{Jenny Häggström, <jenny.haggstrom@umu.se> }

