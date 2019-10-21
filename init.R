# run this script first 
setwd("/axiom/home/wconsagra/Workspace/git/FEFA")
# make all functions in module available to one another
# probably want to be more precise with who can use what, but this rough approach works for testing 
source("basis.R")
source("boundary.R")
source("create.basis.R")
source("eval.fd.R")
source("eval.inprod.R")
source("error_polynomial.R")
source("get.basis.matrix.R")
source("innerprod.R")
source("repr.R")
source("representation.basis.R")
source("seeder.R")
source("shape.functions.R")
source("smooth.basis.R")
source("utility.R")
# external libraries used in code
library(geometry)
library(Matrix)
library(gtools) # for permutation in seeder.R, consider implementing in house solution