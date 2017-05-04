\name{cov.sel.high}
\alias{cov.sel.high}
\title{Model-Free Covariate Selection in High Dimensions}
\description{
Model-free selection of covariates in high dimensions under unconfoundedness for situations where the parameter of interest is an average causal effect. This package is based on  model-free backward elimination algorithms proposed in de Luna, Waernbaum and Richardson (2011) and VanderWeele and Shpitser (2011). Confounder selection can be performed via either Markov/Bayesian networks, random forests or LASSO.
}
\usage{cov.sel.high(T=NULL, Y=NULL, X=NULL,type=c("mmpc","mmhc","rf","lasso"), 
                    betahat=TRUE, parallel=FALSE, Simulate=TRUE,N=NULL, Setting=1,
                    rep=1, Models=c("Linear", "Nonlinear", "Binary"), 
                    alpha=0.05, mmhc_score=c("aic","bic"))
	
}
\arguments{
  \item{T  }{A vector, containing \code{0} and \code{1}, indicating a binary treatment variable.}
  \item{Y  }{A vector of observed outcomes.}
  \item{X  }{A matrix or data frame containing columns of covariates. The covariates may be a mix of continuous, unordered discrete
(to be specified in the data frame using \code{factor}), and ordered discrete (to be specified in the data frame using \code{ordered}).}
  \item{type }{The type of method used for selection. The networks algorithms are \code{"mmpc"} for min-max parents and children (Markov network) and \code{"mmhc"} for max-min hill climbing (Bayesian network). Other available methods are random forests, \code{"rf"}, and LASSO, \code{"lasso"}.} 
 \item{betahat}{If \code{betahat=TRUE} the average treatment effect for each selected subset and the full covariate vector is estimated using propensity score matching (PSM) via the function \code{\link[Matching]{Match}} and using targeted maximum likelihood estimation (TMLE) via the function \code{\link[tmle]{tmle}}. }
    \item{parallel}{If \code{parallel=TRUE} and there is a registered parallel backend then the computation will be parallelized. Default is \code{parallel=FALSE}. }
\item{Simulate}{If data is to be simulated according to one of the designs in Häggström (2017) then \code{Simulate} should be set to \code{TRUE}.}
\item{N}{If Simulate=TRUE, \code{N} is the number of observations to be simulated.}
\item{Setting}{If Simulate=TRUE, \code{Setting} is the simulation setting to be used. Unconfoundedness holds given X if Setting=1. M-bias given X if Setting=2. }
\item{rep}{If Simulate=TRUE, \code{rep} is the number of replications to be simulated.}
\item{Models}{If Simulate=TRUE, \code{Models} is the type of outcome models to be used, options are \code{"Linear"}, \code{"Nonlinear"} and \code{"Binary"}.}
\item{alpha}{A numeric value, the target nominal type I error rate (tuning parameter) for \code{"mmpc"} and \code{"mmhc"}.}
\item{mmhc_score}{The score to use for \code{"mmhc"}.}

}
\details{See Häggström (2017).
}

\value{
 \code{cov.sel.high} returns a list with the following content:
  \item{X.T  }{The set of covariates targeting the subset containing all causes of \code{T}.}
  \item{Q.0  }{The set of covariates targeting the subset of \code{X.T} which is also associated with \code{Y} given \code{T}=0, the response in the control group.}
  \item{Q.1  }{The set of covariates targeting the subset of \code{X.T} which is also associated with \code{Y} given \code{T}=1, the response in the treatment group.}
  \item{Q }{Union of Q.0 and Q.1.}
  \item{X.0  }{The set of covariates targeting the subset containing all causes of \code{Y} given \code{T}=0.}
  \item{X.1  }{The set of covariates targeting the subset containing all causes of \code{Y} given \code{T}=1.}
  \item{X.Y }{Union of X.0 and X.1.}
  \item{Z.0  }{The set of covariates targeting the subset of \code{X.0} which is also associated with \code{T}.}
  \item{Z.1  }{The set of covariates targeting the subset of \code{X.1} which is also associated with \code{T}.}
  \item{Z }{Union of Z.0 and Z.1.}
  \item{X.TY  }{Union of X.T and X.Y, the set of covariates targeting the subset containing all causes of \code{T} and \code{Y}.} 
  \item{cardinalities}{The cardinalities of each selected subset.}
  \item{est_psm}{The PSM estimate of the average causal effect, for the full covariate vector and each selected subset.}
  \item{se_psm}{The Abadie-Imbens standard error for the PSM estimate of the average causal effect, for the full covariate vector and each selected subset.}
  \item{est_tmle}{The TMLE estimate of the average causal effect, for the full covariate vector and each selected subset.}
  \item{se_psm}{The influence-curve based standard error for the TMLE estimate of the average causal effect, for the full covariate vector and each selected subset.}
  \item{N}{The number of observations.}
  \item{Setting}{The Setting used.}
  \item{rep}{The number of replications.}
  \item{Models}{Models used.}
  \item{type}{type used.}
  \item{alpha}{alpha used.}
  \item{mmhc_score}{score used.}
 \item{varnames}{Variable names of the used data.}
}
\references{ 

de Luna, X., I. Waernbaum, and T. S. Richardson (2011). Covariate selection for the nonparametric estimation of an average treatment effect. \emph{Biometrika 98}. 861-875

Häggström, J. (2017). Data-Driven Confounder Selection via Markov and Bayesian Networks. \emph{ArXiv e-prints}.

Nagarajan, R., M. Scutari and S. Lebre. (2013) \emph{Bayesian Networks in R with Applications in Systems
  Biology}. Springer, New York. ISBN 978-1461464457.

Scutari, M. (2010). Learning Bayesian Networks with the bnlearn R Package. \emph{Journal of Statistical Software,
  35}, 1-22. URL http://www.jstatsoft.org/v35/i03/.

Sekhon, J.S. (2011). Multivariate and Propensity Score Matching Software with Automated Balance
  Optimization: The Matching Package for R. \emph{Journal of Statistical Software, 42}, 1-52. URL
  http://www.jstatsoft.org/v42/i07/.
}

\author{Jenny Häggström, <jenny.haggstrom@umu.se> }
\note{ Depending on the method type specified \code{cov.sel.high} calls one of the functions \code{mmpc}, \code{mmhc}, \code{randomForest}, \code{cv.glmnet} and, if betahat=TRUE, \code{Match} and \code{tmle}, therefore the packages \code{bnlearn}, \code{randomForest}, \code{glmnet}, \code{Matching} and \code{tmle} are required.
}
\seealso{\code{\link[bnlearn]{bnlearn-package}}, \code{\link[randomForest]{randomForest}}, \code{\link[glmnet]{cv.glmnet}}, \code{\link[Matching]{Match}} and \code{\link[tmle]{tmle}}}
\examples{
##Use simulated data, select subsets using mmpc 
ans<-cov.sel.high(type="mmpc",N=1000, rep=2, Models="Linear", betahat=FALSE, mmhc_score="aic")


##Use simulated data, select subsets using mmpc and estimate ACEs, parallell version
#library(doParallel)
#library(doRNG)
#cl <- makeCluster(4)
#registerDoParallel(cl)
#ans<-cov.sel.high(type="mmpc", parallel=TRUE,  N=500, rep=10, Models="Linear", mmhc_score="aic")
#stopCluster(cl)

}