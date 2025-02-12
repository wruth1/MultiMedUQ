% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Reg_Par_Covs.R
\name{all_pars_cov_mat}
\alias{all_pars_cov_mat}
\title{Joint covariance matrix of fixed and random effects' parameters from both models}
\usage{
all_pars_cov_mat(fit_Y, fit_M)
}
\arguments{
\item{fit_Y}{Model for \code{Y}}

\item{fit_M}{Model for \code{M}}
}
\value{
A matrix of the joint covariance matrix of fixed and random effects' parameters from both models
}
\description{
Joint covariance matrix of fixed and random effects' parameters from both models
}
\details{
Note: We assume that parameter estimators from the two models are independent. This may not be true in practice. See, e.g., Bauer, Preacher, & Gil (2006) for a method to incorporate inter-model dependence.

Note: I thought that the actual covariance matrix of the parameter estimators requires further division by K. However, numerical experimentation showed that this is not the case. The output of this function scales (approximately) linearly with the number of groups.

Note: It is unnecessary to specify which random effects are included because we extract all present REs from the fitted \code{lme4} model.
}
