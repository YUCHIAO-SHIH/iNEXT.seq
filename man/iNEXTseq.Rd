% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MainFun.R
\name{iNEXTseq}
\alias{iNEXTseq}
\title{Function to calculate phylogenetic gamma, alpha, beta diversity and dissimilarity measures}
\usage{
iNEXTseq(
  data,
  q = c(0, 1, 2),
  base = "coverage",
  datatype = "abundance",
  level = NULL,
  nboot = 10,
  conf = 0.95,
  PDtree = NULL,
  PDreftime = NULL,
  PDtype = "meanPD"
)
}
\arguments{
\item{data}{OTU data can be input as a \code{matrix/data.frame} (species by assemblages), or a \code{list} of \code{matrices/data.frames}, each matrix represents a species-by-assemblages abundance or incidence matrix.\cr}

\item{q}{a numerical vector specifying the diversity orders. Default is \code{c(0, 1, 2)}.}

\item{base}{sample-size-based rarefaction and extrapolation for gamma and alpha diversity (\code{base = "size"}) or coverage-based rarefaction and extrapolation for gamma, alpha, and beta diversity (\code{base = "coverage"}). Default is \code{"coverage"}.}

\item{datatype}{type of input data: \code{"abundance"} (default) or \code{"incidence_raw"}. If \code{"incidence_raw"}, the first column of each matrix must be the total number of sampling units in each assemblage.}

\item{level}{A numerical vector specifying the particular value of sample coverage (between 0 and 1 when \code{base = "coverage"}) or sample size (\code{base = "size"}). See Details.}

\item{nboot}{a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Use \code{0} to skip the bootstrap procedure. Default is \code{10}.}

\item{conf}{a number < 1 specifying the level of confidence interval. Default is \code{0.95}.}

\item{PDtree}{a phylogenetic tree in Newick format for all observed species in the pooled assemblage.}

\item{PDreftime}{a numerical value specifying the reference time for PD. Default is \code{NULL} (i.e., the age of the root of PDtree).}

\item{PDtype}{Type of PD to report: \code{"meanPD"} (default)  or \code{"PD"}}
}
\value{
If \code{base = "coverage"}, returns a list of seven data frames with three coverage-based diversity (gamma, alpha, and beta) and four dissimilarity measures. If \code{base = "size"}, returns a list of two data frames with gamma and alpha diversity.
}
\description{
\code{iNEXTseq}: This function calculates interpolated and extrapolated phylogenetic gamma, alpha, beta diversity, and dissimilarity measures based on abundance or incidence data, using the framework developed by Chiu et al.
}
\examples{
# Abundance data example
data("fungi")
data("fungi_tree")
output <- iNEXTseq(fungi[1], q = c(0,1,2), level = seq(0.5, 1, 0.05), nboot = 10,
                   conf = 0.95, PDtree = fungi_tree, PDreftime = NULL, PDtype = 'meanPD')
output

# Incidence data example (assuming incidence_data is formatted correctly)
# output <- iNEXTseq(incidence_data, q = c(0,1,2), datatype = "incidence_raw",
#                    level = seq(0.5, 1, 0.05), PDtree = incidence_tree, PDtype = 'meanPD')

}
