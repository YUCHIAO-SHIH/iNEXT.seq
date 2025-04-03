\name{global_tree}
\alias{global_tree}
\docType{data}
\title{
Phylogenetic Tree for Global OTU Data
}
\description{
  A rooted phylogenetic tree corresponding to the OTUs included in the \code{global} dataset. This tree provides phylogenetic relationships among microbial taxa sampled from sediment, ocean, and freshwater environments.
}
\usage{data(global_tree)}
\format{
  An object of class \code{phylo} (from the \code{ape} package), containing a rooted phylogenetic tree with the following components: \cr
  \code{$edge}: an integer matrix (1052 × 2) defining the tree structure by parent and child node relationships. \cr
  \code{$Nnode}: number of internal nodes (526). \cr
  \code{$tip.label}: a character vector of 527 OTU names corresponding to the tips of the tree (e.g., \code{"OTU_252853"}, \code{"OTU_512943"}, ...). \cr
  \code{$edge.length}: a numeric vector of 1052 branch lengths representing evolutionary distances between nodes.
}
\details{
  This tree is derived from the full dataset originally published by Caporaso et al. (2011), and pruned to match the 5\% OTU subset used in the \code{global} dataset. It is intended for use in phylogenetic diversity and beta diversity analyses.
}
\references{
  Caporaso, J. G., Lauber, C. L., Walters, W. A., Berg-Lyons, D., Lozupone, C. A., Turnbaugh, P. J., ... & Knight, R. (2011). Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample. \emph{Proceedings of the National Academy of Sciences}, \strong{108}(Supplement 1), 4516–4522.
}
\keyword{datasets}
