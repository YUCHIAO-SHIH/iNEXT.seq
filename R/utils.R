#' Internal fetch for non-exported functions
#'
#' Dynamically fetch a non-exported function from another namespace.
#' Used to access internal iNEXT.3D and iNEXT.beta3D functions safely.
#' 
#' @param fname Function name (string)
#' @param pkg Package name (string)
#' @return The function object
#' @export
.internal_fetch <- function(fname, pkg) {
  getFromNamespace(fname, pkg)
}