#' Norm Estimation
#' 
#' Estimates the 1-norm of a power of a matrix.
#' 
#' @param x
#' A numeric matrix.
#' @param pow
#' An integer power.
#' 
#' @return
#' exp(x)
#' 
#' @export
normest <- function(x, pow)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call(R_libexpm_normest, x, as.integer(pow))
}
