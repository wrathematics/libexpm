normest <- function(x, pow)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call("R_normest", x, as.integer(pow))
}
