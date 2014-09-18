dyn.load("expm_wrap.so")
source("normest.r")

n <- 1000
x <- matrix(rnorm(n*n), n, n)


matpow <- function(x, pow)
{
  ret <- x
  if (pow == 1) return(ret)
  
  for (i in 2:pow)
  {
    ret <- ret %*% x
  }
  
  ret
}


pow <- 2


t1 <- system.time(test1 <- normest(x, pow))[3]
t2 <- system.time(test2 <- norm(matpow(x, pow), type="O"))[3]

t1
t2

test1
test2
