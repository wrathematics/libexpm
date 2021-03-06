~~~
 _     ___ ____  _______  ______  __  __ 
| |   |_ _| __ )| ____\ \/ /  _ \|  \/  |
| |    | ||  _ \|  _|  \  /| |_) | |\/| |
| |___ | || |_) | |___ /  \|  __/| |  | |
|_____|___|____/|_____/_/\_\_|   |_|  |_|

__     __            _             0.2.0
\ \   / ___ _ __ ___(_) ___  _ __  
 \ \ / / _ | '__/ __| |/ _ \| '_ \ 
  \ V |  __| |  \__ | | (_) | | | |
   \_/ \___|_|  |___|_|\___/|_| |_|
~~~

Copyright 2013-2014, Drew Schmidt.  All rights reserved.



# Introduction 

libexpm is a portable, high-performance C library for computing
the matrix exponential of a dense matrix.

Formally, a matrix exponential is the power series:

  expm(X) = I + X/1! + X^2/2! + X^3/3! + ...



# Building

To build the library, you will need:

  * Cmake 2.6
  * BLAS and LAPACK libraries.
  * R >= 2.15.0 if building the R/ subtree.

Build with make.  Both shared and static libraries will be produced
in the build/ tree.



# License 

The library is released under the 2-clause BSD license.  See the 
file LICENSE for more details.



# Routines and Documentation



