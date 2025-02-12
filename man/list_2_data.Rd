% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Generate_Sample_Data.R
\name{list_2_data}
\alias{list_2_data}
\title{Construct a single data frame by `rbind`-ing elements of a list. Add/Create group labels.}
\usage{
list_2_data(X_list, group_labels = NULL)
}
\arguments{
\item{X_list}{A list of data frames. Must all have the same number of columns.}

\item{group_labels}{Optionally, a vector of labels for each element of `X_list`. Must be either the same length as `X_list` or equal to the total number of rows among elements in `X_list`. If `NULL`, labels are `G1`, `G2`,...}
}
\value{
A data frame containing all elements of `X_list` with an extra column for labels.
}
\description{
Construct a single data frame by `rbind`-ing elements of a list. Add/Create group labels.
}
\examples{
data = as.data.frame(matrix(c(1,0,0,1), nrow = 2))
more_data = as.data.frame(matrix(c(0,1,1,0), nrow = 2))
data_list = list(data, more_data)

list_2_data(data_list)

data_names = c("Alice", "Bob")
list_2_data(data_list, data_names)
}
