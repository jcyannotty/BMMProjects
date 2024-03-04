#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("/home/johnyannotty/Documents/openbt/src")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")

source("/home/johnyannotty/Documents/SciData/NuclearPhysics/R/heat_2d_ising.R")
#source("/home/johnyannotty/Documents/BayesToolBox/bayestb/LinearRegression/bayesian_regression.R")
#source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/kernels.R")
#source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/gp_utils.R")

#----------------------------------------------------------
# Load Data 
filedir = "/home/johnyannotty/Documents/Dissertation/results/EFT/spheat_model_fit/"
ms = readRDS(paste0(filedir,"ms_n20_01_23_24.rds"))

f_train = cbind(fsg = ms$fsg_train, flg = ms$flg_train)
f_test = cbind(fsg = ms$fsg_test, flg = ms$flg_test)

nu = 30
rho = 1
sig2_hat = max(apply(apply(f_train, 2, function(x) (x-ms$y_train)^2),2,min))
q0 = 4

fit=train.openbtmixing(ms$g_train,ms$y_train,f_train,pbd=c(1.0,0),ntree = 10,ntreeh=1,numcut=100,tc=2,model="mixbart",modelname="physics_model",
                       ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
                       power = 2.0, base = 0.95, minnumbot = 2, overallsd = sqrt(sig2_hat), k = 1.0, overallnu = nu,
                       summarystats = FALSE, rpath = FALSE, q = q0, rshp1 = 4, rshp2 = 20,
                       stepwpert = 0.1, probchv = 0.1, batchsize = 1000)


#Get mixed mean function
fitp=predict.openbtmixing(fit,x.test = ms$g_test, f.test = f_test,
                          tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "all", proj_type = "euclidean")

plot(ms$g_test, ms$f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(-2,7))
points(ms$g_train, ms$y_train, pch = 3)
lines(ms$g_test, f_test[,1], col = 'red', lty = 2)
lines(ms$g_test, f_test[,2], col = 'blue', lty = 2)
lines(ms$g_test, fitp$mmean, col = 'purple', lwd = 2)
lines(ms$g_test, fitp$m.lower, col = 'purple', lwd = 2, cex = 0.5, lty = 'dashed')
lines(ms$g_test, fitp$m.upper, col = 'purple', lwd = 2, cex = 0.5, lty = 'dashed')

#Plot model weights
plot(ms$g_test, fitp$wmean[,1], pch = 16, col = 'red', type = 'l', ylim = c(-0.2,1.2), lwd = 2, 
     panel.first = {grid(col = 'lightgrey')})
lines(ms$g_test, fitp$wmean[,2], col = 'blue', lwd = 2)
lines(ms$g_test, fitp$w.upper[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(ms$g_test, fitp$w.lower[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(ms$g_test, fitp$w.upper[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(ms$g_test, fitp$w.lower[,2], col = 'blue', lty = 'dashed', lwd = 1)

# Get wsum
wsum = 0*fitp$wdraws[[1]]
for(i in 1:length(fitp$wdraws)){
  wsum = wsum + fitp$wdraws[[i]]
}
wsum_mean = apply(wsum,2,mean)
wsum_lb = apply(wsum,2,quantile, 0.025)
wsum_ub = apply(wsum,2,quantile, 0.975)

#Get projected mixed mean function
plot(ms$g_test, ms$f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(-15,15))
points(ms$g_train, ms$y_train, pch = 3)
lines(ms$g_test, f_test[,1], col = 'red', lty = 2)
lines(ms$g_test, f_test[,2], col = 'blue', lty = 2)
lines(ms$g_test, fitp$pmmean, col = 'purple', lwd = 2)
lines(ms$g_test, fitp$pm.lower, col = 'orange', lwd = 2, cex = 0.5)
lines(ms$g_test, fitp$pm.upper, col = 'orange', lwd = 2, cex = 0.5)

#Plot Projected model weights
plot(ms$g_test, fitp$pwmean[,1], pch = 16, col = 'red', type = 'l', ylim = c(-1,2.0), lwd = 2, 
     panel.first = {grid(col = 'lightgrey')})
lines(ms$g_test, fitp$pwmean[,2], col = 'blue', lwd = 2)
lines(ms$g_test, fitp$pw.upper[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(ms$g_test, fitp$pw.lower[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(ms$g_test, fitp$pw.upper[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(ms$g_test, fitp$pw.lower[,2], col = 'blue', lty = 'dashed', lwd = 1)

#Get estimated delta
plot(ms$g_test, ms$f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(-2,2))
lines(ms$g_test, f_test[,1], col = 'red', lty = 2)
lines(ms$g_test, f_test[,2], col = 'blue', lty = 2)
lines(ms$g_test, fitp$dmean, col = 'purple', lwd = 2)
lines(ms$g_test, fitp$d.upper, col = 'orange', lwd = 2, cex = 0.5)
lines(ms$g_test, fitp$d.lower, col = 'orange', lwd = 2, cex = 0.5)


#------------------------------------------------
# Manual Softmax projection
#------------------------------------------------
library(rBayesianOptimization)
softmax_l2 = function(tmp){
  fita = rowSums(f_test*exp(fitp$wmean/tmp)/rowSums(exp(fitp$wmean/tmp)))
  score = -sum((fitp$mmean - fita)^2) # max the negative loss
  return(list(Score = score, Pred = 0))
} 

tmp_bounds = list(tmp = c(0.01,1))
init_grid_dt = data.frame(tmp = seq(0.015,0.95, length = 10))
bayes_temp = BayesianOptimization(FUN = softmax_l2, acq = "ei",
                                  bounds = tmp_bounds, init_grid_dt = init_grid_dt,
                                  init_points = 0.8, n_iter = 10)


# Plot the loss
tmp_grid = seq(0.01,0.99,length = 100)
plot(tmp_grid, sapply(tmp_grid, function(x) softmax_l2(x)$Score),ylim = c(-10,0))
points(bayes_temp$Best_Par,bayes_temp$Best_Value,pch = 3, col = "red")

fitp_sm=predict.openbtmixing(fit,x.test = ms$g_test, f.test = f_test,
                          tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "all", proj_type = "softmax", temperature = bayes_temp$Best_Par)


#Get projected mixed mean function
plot(ms$g_test, ms$f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(-15,15))
points(ms$g_train, ms$y_train, pch = 3)
lines(ms$g_test, f_test[,1], col = 'red', lty = 2)
lines(ms$g_test, f_test[,2], col = 'blue', lty = 2)
lines(ms$g_test, fitp_sm$pmmean, col = 'purple', lwd = 2)
lines(ms$g_test, fitp_sm$pm.lower, col = 'orange', lwd = 2, cex = 0.5)
lines(ms$g_test, fitp_sm$pm.upper, col = 'orange', lwd = 2, cex = 0.5)

#Plot Projected model weights
plot(ms$g_test, fitp_sm$pwmean[,1], pch = 16, col = 'red', type = 'l', ylim = c(-1,2.0), lwd = 2, 
     panel.first = {grid(col = 'lightgrey')})
lines(ms$g_test, fitp_sm$pwmean[,2], col = 'blue', lwd = 2)
lines(ms$g_test, fitp_sm$pw.upper[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(ms$g_test, fitp_sm$pw.lower[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(ms$g_test, fitp_sm$pw.upper[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(ms$g_test, fitp_sm$pw.lower[,2], col = 'blue', lty = 'dashed', lwd = 1)

#Get estimated delta
plot(ms$g_test, ms$f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(-2,2))
lines(ms$g_test, f_test[,1], col = 'red', lty = 2)
lines(ms$g_test, f_test[,2], col = 'blue', lty = 2)
lines(ms$g_test, fitp_sm$dmean, col = 'purple', lwd = 2)
lines(ms$g_test, fitp_sm$d.upper, col = 'orange', lwd = 2, cex = 0.5)
lines(ms$g_test, fitp_sm$d.lower, col = 'orange', lwd = 2, cex = 0.5)



#------------------------------------------------
# Save Results
#------------------------------------------------
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
  x_train = ms$g_train,
  y_train = ms$y_train,
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

# Euclidean Projection Results
wts_l2_proj = list(
  pred_mean = fitp$pmmean,
  pred_ub = fitp$pm.upper,
  pred_lb = fitp$pm.lower,
  wts_mean = fitp$pwmean,
  wts_ub = fitp$pw.upper,
  wts_lb = fitp$pw.lower,
  delta_mean = fitp$dmean,
  delta_lb = fitp$d.lower,
  delta_ub = fitp$d.upper
)

# Softmax Projection Results
softmax_res = list(
  pred_mean = fitp_sm$pmmean,
  pred_ub = fitp_sm$pm.upper,
  pred_lb = fitp_sm$pm.lower,
  wts_mean = fitp_sm$pwmean,
  wts_ub = fitp_sm$pw.upper,
  wts_lb = fitp_sm$pw.lower,
  delta_mean = fitp_sm$dmean,
  delta_lb = fitp_sm$d.lower,
  delta_ub = fitp_sm$d.upper,
  tmp = bayes_temp$Best_Par 
)

filedir = "/home/johnyannotty/Documents/Dissertation/results/EFT/spheat_model_fit/"
saveRDS(fit_data, paste0(filedir,"bmm_dpath_res_02_05_24.rds"))
saveRDS(wts_l2_proj, paste0(filedir,"bmm_dpath_wl2_proj_02_05_24.rds"))
saveRDS(softmax_res, paste0(filedir,"bmm_dpath_softmax_proj_02_05_24.rds"))
