# p is the order of the Pade' approximation. expokit uses 6. For us, 0<p<14
expm <- function(x, p=6)
{
  if (nrow(x) != ncol(x))
    stop("Matrix exponentiation is only defined for square matrices.")
  
  if (p < 1L || p > 13L)
    stop("argument 'p' must be between 1 and 13.")
   
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  exp <- .Call("R_expm", x, as.integer(p))
  
  return( exp )
}


