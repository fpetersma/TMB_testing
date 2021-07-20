# ============================================================================ #
# Trial of funciton from chapter 7
# ============================================================================ #

library(TMB) # load TMB

compile("linreg.cpp") # compile the C++ file

dyn.load(dynlib("linreg")) # turn the executable into a dynamic library, and load

set.seed(123) # set seed

data <- list(Y = rnorm(100) + 1:100, x=1:100) # create data

parameters <- list(a=0, b=0, logSigma=0) # define starting values for parameters

obj <- MakeADFun(data, parameters, DLL="linreg")
obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
obj$he()    ## <-- Analytical hessian
sdreport(obj)