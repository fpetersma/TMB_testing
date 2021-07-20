## =============================================================================
## First attempt at a simple ASCR model in TMB                    
## =============================================================================

## Mark Bravington post on why gdb no longer works in Rtools4.0, and how to fix it
## https://r.789695.n4.nabble.com/problem-adding-gdb-to-RTOOLS40-on-Windows-td4768577.html#a4768604

## Libraries
library(TMB)        # load TMB
# library(BuildSys)   # For better debugging
# library(TMBdebug)   # avoids R crashing when TMB crashes
# library(TMBhelper)  # Includes: AICTMB -- calculate AIC based on model output; and Check_Identifiable -- automatically check for non-identiable fixed effects

## Compile the code and load the dll
compile("ascr.cpp", "-O1 -g", DLLFLAGS="") # To make gdbsource() work (https://github.com/kaskr/adcomp/issues/67)
                                           
# also, see documentation
dyn.load(dynlib("ascr")) # turn the executable into a dynamic library, and load

## Set a seed for reproducibility
set.seed(123) # set seed

## Load RData file with data from main.R, and extract data and parameters for TMB
load("test_data.RData")

data <- list(Y_rec = dat$bearings_rad,
             Y_grid = dat$grid_bearings,
             X = dat$distances,
             W = dat$det_hist,
             R = dat$received_levels,
             A = dat$A_x$area,
             trunc_level = dat$trunc_level) # create data


parameters <- list(log_D = 0.01 , 
                   logit_g0 = 0.6 , 
                   # log_kappa = 3, 
                   log_beta_r = 2.9 , 
                   log_sd_r = 0.2 , 
                   log_mu_s = 5.1 ,
                   log_kappa_low = 0.0 ,
                   log_kappa_high = 3.0 ,
                   logit_mix_bear = -0.4 ) # define starting values for parameters

# parameters <- list(log_D = -5.5,
#                    logit_g0 = -5.5 ,
#                    # log_kappa = 3,
#                    log_beta_r = 5.5 ,
#                    log_sd_r = -5.5,
#                    log_mu_s = -5.5,
#                    log_kappa_low = 5.5,
#                    log_kappa_high = 5.5,
#                    logit_mix_bear = -5.5) # define starting values for parameters

rm("dat", "par")

## Do TMB things
obj <- MakeADFun(data, parameters, DLL="ascr", silent = FALSE)
# obj <- normalize(obj, flag) # not sure what flat = "flag" does <- read online, its for 
opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1))#, lower = -5.5, upper = 5.5)
sd_report <- sdreport(obj, getJointPrecision = TRUE)

## Start the bootstrap =========================================================
boot_fits <- list()

B <- 1000 # number of bootstraps
for(i in 1:B) {
  seed <- 18012021
  set.seed(18012021 + i)
  boot_rows <- sample(1:nrow(data$W), nrow(data$W), replace = TRUE)
  
  boot_data <- list(Y_rec = data$Y_rec[boot_rows, ],
                    Y_grid = data$Y_grid,
                    X = data$X,
                    W = data$W[boot_rows, ],
                    R = data$R[boot_rows, ],
                    A = data$A,
                    trunc_level = data$trunc_level) # create data
  
  ## Do bootstrap TMB things
  obj_boot <- MakeADFun(boot_data, parameters, DLL="ascr", silent = FALSE)
  # obj <- normalize(obj, flag) # not sure what flat = "flag" does <- read online, its for 
  opt_boot <- stats::nlminb(obj_boot$par, obj_boot$fn, obj_boot$gr, control = list(trace = 1))
  
  boot_fits[[i]] <- opt_boot
}
opt_boot$par

pars <- sapply(boot_fits, function(fit) fit$par)

Fview

## Do bootstrap TMB things
obj_boot <- MakeADFun(boot_data, parameters, DLL="ascr", silent = FALSE)
# obj <- normalize(obj, flag) # not sure what flat = "flag" does <- read online, its for 
opt_boot <- stats::nlminb(obj_boot$par, obj_boot$fn, obj_boot$gr, control = list(trace = 1))#, lower = -5.5, upper = 5.5)
sdreport(obj, getJointPrecision = TRUE)









## Do old TMB things
obj <- MakeADFun(data, parameters, DLL="ascr", silent = FALSE)
obj$hessian <- TRUE
obj$method <- "BFGS"
obj$control <- list(trace = 6)
# obj$lower <- c(10, ) #rep(-10, length(parameters))
# obj$upper <- 5.5 #rep( 10, length(parameters))

opt <- do.call("optim", obj)

# opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1))

opt
opt$hessian ## <-- FD hessian from optim
obj$he()    ## <-- Analytical hessian
sdreport(obj)
## 