\name{global}
\alias{global}
\docType{data}
\title{
Global OTU Count Data (5\% Subset)
}
\description{
  This dataset is a 5\% subset of the full global OTU count data derived from the package \code{microbiome}, originally collected and analyzed by Caporaso et al. (2011) to study microbial diversity across a wide range of environments. This demo version focuses on samples from sediment, ocean, and freshwater habitats, selected due to their ecological importance and representation in global surveys.
}
\usage{data(global)}
\format{
  A species-by-assemblages \code{data.frame}, where rows represent OTUs (Operational Taxonomic Units), and columns represent sampling sites from sediment, ocean, and freshwater environments. \cr
  Example columns include: \cr
  \code{$ TRRsed1}: num 0 2 0 0 9 ... \cr
  \code{$ TRRsed2}: num 0 2 2 6 82 ... \cr
  \code{$ TRRsed3}: num 0 5 0 4 27 ... \cr
  \code{$ NP2}: num 0 4 0 0 0 ... \cr
  \code{$ NP3}: num 0 0 0 0 0 ... \cr
  \code{$ SLEpi20M}: num 0 1 0 0 0 ... \cr
  \code{$ LMEpi24M}: num 1 0 0 0 0 ... \cr
}
\details{
  The original dataset contains a comprehensive OTU table with high sequencing depth across multiple global environments. To improve computational efficiency and enable quick demonstrations, this version includes only \strong{5\% of the original data}, filtered to include only samples from \strong{sediment, ocean, and freshwater ecosystems}.
}
\references{
  Caporaso, J. G., Lauber, C. L., Walters, W. A., Berg-Lyons, D., Lozupone, C. A., Turnbaugh, P. J., ... & Knight, R. (2011). Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample. \emph{Proceedings of the National Academy of Sciences}, \strong{108}(Supplement 1), 4516–4522.
}
\keyword{datasets}
