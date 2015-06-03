#' Matrix Exponentiation
#' 
#' Matrix exponentiation via Pade' approximations.
#' 
#' @param x
#' A numeric matrix.
#' @param p
#' Order of the Pade' approximation.  Should be between 1 and 13.
#' 6 is typical.
#' 
#' @return
#' exp(x)
#' 
#' @export
expm <- function(x, p=6)
{
  if (nrow(x) != ncol(x))
    stop("Matrix exponentiation is only defined for square matrices.")
  
  if (p < 1L || p > 13L)
    stop("argument 'p' must be between 1 and 13.")
   
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  exp <- .Call(R_libexpm_expm_5_1, x, as.integer(p))
  
  return( exp )
}


