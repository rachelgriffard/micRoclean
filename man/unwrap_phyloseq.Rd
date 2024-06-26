% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/management_functions.R
\name{unwrap_phyloseq}
\alias{unwrap_phyloseq}
\title{unwrap_phyloseq}
\usage{
unwrap_phyloseq(phyloseq)
}
\arguments{
\item{phyloseq}{Phyloseq object to unwrap into required count and meta matrices for pipeline functions}
}
\value{
List containing counts and metadata data frames for input into pipeline functions
}
\description{
Take input phyloseq object and "unwrap" data into matrices expected as input for functions in
micRoclean package.
}
\seealso{
Other management: 
\code{\link{wrap_phyloseq}()}
}
\concept{management}
