dyn.load("expm_wrap.so")
source("expm.r")

n <- 500
x <- matrix(rnorm(n*n), n, n)

t1 <- system.time(exp1 <- expm(x))[3]

library(Matrix)
t2 <- system.time(exp2 <- Matrix::expm(x))[3]

t1
t2

all.equal(exp1, as.matrix(exp2), check.attributes=FALSE)
