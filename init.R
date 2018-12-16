# run this script first 
setwd("/axiom/home/wconsagra/Workspace/git/FEFA")
# make all functions in module available to one another
# probably want to be more precise with who can use what, but this rough approach works for testing 
source("basis.R")
source("create.basis.R")
source("eval.fd.R")
source("get.basis.matrix.R")
source("utility.R")
# external libraries used in code
library(geometry)
library(Matrix)