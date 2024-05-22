#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Name: eft_examples.R
# Author: John Yannotty (yannotty.1@buckeyemail.osu.edu)
# Desc: R code to perform model mixing using BART on three EFT examples.
#       The tuning parameters that are currently set were used to produce the plots,
#       which appear in the file "EFT-Model-Mixing-Using-BART.html" 
# Version: 1.0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("/home/johnyannotty/Documents/openbt/src")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/openbt/R/honda_eft.R")
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Optimization/optimization_helpers.R")

#----------------------------------------------------------
#Get training data
ex1_data = get_data(20, 300, 0.005, 4, 4, 0.03, 0.5, 321)

f_train = ex1_data$f_train
f_test = ex1_data$f_test
y_train = ex1_data$y_train
x_train = ex1_data$x_train
x_test = ex1_data$x_test
f_train = ex1_data$f_train
f_test = ex1_data$f_test
f0_test = ex1_data$fg_test

pwwts_train = (1/ex1_data$f_train_dsd^2)/rowSums(1/ex1_data$f_train_dsd^2)
pwwts_test = (1/ex1_data$f_test_dsd^2)/rowSums(1/ex1_data$f_test_dsd^2)

#----------------------------------------------------------
nu = 30
rho = 1
sig2_hat = max(apply(apply(f_train, 2, function(x) (x-y_train)^2),2,min))
lam = rho*sig2_hat*(nu+2)/nu
q0 = 4
fit=train.openbtmixing(x_train,y_train,f_train,pbd=c(1.0,0),ntree = 10,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="physics_model",
                       ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
                       power = 2.0, base = 0.95, minnumbot = 3, overallsd = sqrt(sig2_hat), k = 3, overallnu = nu,
                       summarystats = FALSE, rpath = TRUE, q = q0, rshp1 = 2, rshp2 = 10,
                       stepwpert = 0.1, probchv = 0.1, batchsize = 1000)


#Get mixed mean function
fitp=predict.openbtmixing(fit,x.test = x_test, f.test = f_test,tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "all")

plot(x_test, f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(1,3.6))
points(x_train, y_train, pch = 3)
lines(x_test, f_test[,1], col = 'red', lty = 2)
lines(x_test, f_test[,2], col = 'blue', lty = 2)
lines(x_test, fitp$mmean, col = 'purple', lwd = 2)
lines(x_test, fitp$m.lower, col = 'orange', lwd = 2, cex = 0.5)
lines(x_test, fitp$m.upper, col = 'orange', lwd = 2, cex = 0.5)

# Get wsum
wsum = 0*fitp$wdraws[[1]]
for(i in 1:ncol(f_test)){
  wsum = wsum + fitp$wdraws[[i]]
}
wsum_mean = apply(wsum,2,mean)
wsum_lb = apply(wsum,2,quantile, 0.025)
wsum_ub = apply(wsum,2,quantile, 0.975)


#Plot model weights
plot(x_test, fitp$wmean[,1], pch = 16, col = 'red', type = 'l', ylim = c(-1,2.0), lwd = 2, 
     panel.first = {grid(col = 'lightgrey')})
lines(x_test, fitp$wmean[,2], col = 'blue', lwd = 2)
lines(x_test, fitp$w.upper[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, fitp$w.lower[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, fitp$w.upper[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(x_test, fitp$w.lower[,2], col = 'blue', lty = 'dashed', lwd = 1)

abline(h = 1, col = 'grey', lty = 'dashed')
abline(h = 0, col = 'grey', lty = 'dashed')

# Plot sigma
hist(fitp$sdraws[,1])

plot(x_test, f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(1,3.6))
points(x_train, y_train, pch = 3)
lines(x_test, f_test[,1], col = 'red', lty = 2)
lines(x_test, f_test[,2], col = 'blue', lty = 2)
lines(x_test, fitp$pmmean, col = 'purple', lwd = 2)
lines(x_test, fitp$pm.lower, col = 'orange', lwd = 2, cex = 0.5)
lines(x_test, fitp$pm.upper, col = 'orange', lwd = 2, cex = 0.5)



# Save results
fit_data = list(
  pred_mean = fitp$mmean,
  pred_ub = fitp$m.upper,
  pred_lb = fitp$m.lower,
  wts_mean = fitp$wmean,
  wts_ub = fitp$w.upper,
  wts_lb = fitp$w.lower,
  wsum_mean = wsum_mean,
  wsum_lb = wsum_lb,
  wsum_ub = wsum_ub,
  proj_mean = fitp$pmmean,
  proj_ub = fitp$pm.upper,
  proj_lb = fitp$pm.lower,
  pwts_mean = fitp$pwmean,
  pwts_ub = fitp$pw.upper,
  pwts_lb = fitp$pw.lower,
  delta_mean = fitp$dmean,
  delta_lb = fitp$d.lower,
  delta_ub = fitp$d.upper,
  wsum_mean = wsum_mean,
  x_train = x_train,
  y_train = y_train,
  m = fit$m,
  k = fit$k,
  shp1 = fit$rshp1,
  shp2 = fit$rshp2,
  minnodesz = fit$minnumbot,
  q = q0,
  base = fit$base,
  power = fit$power,
  nu = fit$overallnu,
  lam = fit$overalllambda
)
filedir = "/home/johnyannotty/Documents/RandomPathBART/EFT_Results/"
saveRDS(fit_data, paste0(filedir,"honda_sg4_lg4_rpath_12_21_23.rds"))





#---------------------------------------------------
# Simplex HARD CODE
#---------------------------------------------------

# Simplex Wts
swts = t(apply(fitp$wmean,1,simplex_l2))

plot(x_test, swts[,1], pch = 16, col = 'red', type = 'l', ylim = c(-1,2.0), lwd = 2, 
     panel.first = {grid(col = 'lightgrey')})
lines(x_test, swts[,2], col = 'blue', lwd = 2)

# Delta 
delta = fitp$mmean - rowSums(swts*f_test)
plot(x_test,delta, type = "l")

plot(x_test, f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(1,2.8))
points(x_train, y_train, pch = 3)
lines(x_test, f_test[,1], col = 'red', lty = 2)
lines(x_test, f_test[,2], col = 'blue', lty = 2)
#lines(x_test, fitp$mmean, col = 'purple', lwd = 2)
lines(x_test, rowSums(swts*f_test), col = "green", lwd = 2)


# Project results
N = nrow(fitp$mdraws)
K = ncol(f_test)
n_test = nrow(f_test)
swmean = swub = swlb = matrix(0, nrow = n_test, ncol = K)
spmean = spub = splb = 0
for(i in 1:n_test){
  wts_matrix = matrix(0, nrow = N, ncol = K)
  for(j in 1:K){
    wts_matrix[,j] = fitp$wdraws[[j]][,i]  
  }
  # Weights
  s_wts = t(apply(wts_matrix,1,simplex_l2))
  swmean[i,] = apply(s_wts, 2, mean)
  swlb[i,] = apply(s_wts, 2, quantile, 0.025)
  swub[i,] = apply(s_wts, 2, quantile, 0.975)
  
  # Predictions
  spdraws = rowSums(t(t(s_wts)*f_test[i,]))
  spmean[i] = mean(spdraws)
  splb[i] = quantile(spdraws,0.025)
  spub[i] = quantile(spdraws,0.975)
  cat("Progress: ", round(i/n_test,4),"\r")
}

plot(x_test, spmean, type = "l")
lines(x_test, spub, col = "green2")
lines(x_test, splb, col = "green2")


plot(x_test, swmean[,1], type = "l", ylim = c(-1,2), col = "red")
lines(x_test, swub[,1], col = "red")
lines(x_test, swlb[,1], col = "red")
lines(x_test, swmean[,2], col = "blue")
lines(x_test, swub[,2], col = "blue")
lines(x_test, swlb[,2], col = "blue")
