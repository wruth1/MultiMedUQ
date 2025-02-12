% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Generic_Helpers.R
\name{check_theta}
\alias{check_theta}
\title{Check for invalid RE parameters}
\usage{
check_theta(theta)
}
\arguments{
\item{theta}{A vector of SDs and correlations arranged according to our parameterization}
}
\value{
1 if all parameters are in the interior of the parameter space. 0 if any are on the boundary (and gives a warning). Gives an error if any parameter values are invalid (e.g. SD < 0, abs(corr) > 1).
}
\description{
Check for invalid RE parameters
}
