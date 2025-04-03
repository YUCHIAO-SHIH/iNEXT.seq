\name{fungi_tree}
\alias{fungi_tree}
\docType{data}
\title{Phylogenetic tree for \code{fungi}}

\description{
  This dataset includes a phylogenetic tree spanning 351 OTUs listed in the dataset \code{fungi}. \cr
  The tree was transformed into an ultrametric tree with a fixed height of 1 using the \code{chronos} function from the \code{ape} package.
}

\usage{data(fungi_tree)}

\format{
An object of class \code{phylo}, structured as a list with the following components: \cr
\itemize{
  \item \code{edge}: Integer matrix [653, 2] indicating the tree edges \cr
  \item \code{edge.length}: Numeric vector of edge lengths \cr
  \item \code{Nnode}: Integer; number of internal nodes \cr
  \item \code{tip.label}: Character vector of OTU labels (length 351) \cr
  \item \code{node.label}: Character vector of internal node labels (length 303) \cr
}
}

\keyword{datasets}

