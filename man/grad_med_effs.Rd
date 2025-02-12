% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Med_Effs.R
\name{grad_med_effs}
\alias{grad_med_effs}
\alias{grad_TE_diff}
\alias{grad_TE_rat}
\alias{grad_TE_or}
\alias{grad_DE_diff}
\alias{grad_DE_rat}
\alias{grad_DE_or}
\alias{grad_IE_diff}
\alias{grad_IE_rat}
\alias{grad_IE_or}
\title{Title}
\usage{
grad_TE_diff(ENC_11, ENC_10, ENC_01, ENC_00)

grad_TE_rat(ENC_11, ENC_10, ENC_01, ENC_00)

grad_TE_or(ENC_11, ENC_10, ENC_01, ENC_00)

grad_DE_diff(ENC_11, ENC_10, ENC_01, ENC_00)

grad_DE_rat(ENC_11, ENC_10, ENC_01, ENC_00)

grad_DE_or(ENC_11, ENC_10, ENC_01, ENC_00)

grad_IE_diff(ENC_11, ENC_10, ENC_01, ENC_00)

grad_IE_rat(ENC_11, ENC_10, ENC_01, ENC_00)

grad_IE_or(ENC_11, ENC_10, ENC_01, ENC_00)
}
\arguments{
\item{ENC_11, ENC_10, ENC_01, ENC_00}{Expected nested counterfactuals at specified levels of x and x_m respectively.}
}
\value{
Gradient of the mediation effect of specified flavour and scale.
}
\description{
Title
}
