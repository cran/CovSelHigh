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
  \item{resmat }{A matrix with the following columns (columns related to average causal effect estimation are only present when \code{betahat=TRUE}):
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
    \code{betahatest_psm}
    \code{betahatse_psm}
    \code{betahat_cicov_psm}
    \code{ciL_psm}
    \code{ciU_psm}
    \code{betahatest_tmle}
    \code{betahatse_tmle}
    \code{betahat_cicov_tmle}
    \code{ciL_tmle}
    \code{ciU_tmle}
    \code{ciwidth_psm}
    \code{ciwidth_tmle}

  }
  
  
    \item{summary_resmat }{A list with the following columns (columns related to average causal effect estimation are only present when \code{betahat=TRUE}):
    \code{Subset_selection}
    \code{Median_cardinality}  
    \code{Betahat_bias_psm}  
    \code{Betahat_sd_psm}  
    \code{Betahat_mse_psm}  
    \code{Betahat_CI_coverage_psm}  
    \code{Betahat_CI_width_psm}  
    \code{Betahat_mean_lower_CI_psm}  
    \code{Betahat_mean_upper_CI_psm}
    \code{Betahat_bias_psm}  
    \code{Betahat_sd_psm}  
    \code{Betahat_mse_psm}  
    \code{Betahat_CI_coverage_psm}  
    \code{Betahat_CI_width_psm}  
    \code{Betahat_mean_lower_CI_psm}  
    \code{Betahat_mean_upper_CI_psm}  
  }
  
   \item{xtable1}{LaTeX table summarizing the results.}
  \item{xtable2}{LaTeX table summarizing the results.}
}


\references{ 

Häggström, J. (2017). Data-Driven Confounder Selection via Markov and Bayesian Networks. \emph{ArXiv e-prints}.

}

\author{Jenny Häggström, <jenny.haggstrom@umu.se> }

