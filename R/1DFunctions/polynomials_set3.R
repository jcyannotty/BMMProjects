#-------------------------------------------------
# Mixing polynomials - 3 Models in 1D
#-------------------------------------------------
setwd("/home/johnyannotty/Documents/openbt/src")

# Load the R wrapper functions to the OpenBT library.
source("/home/johnyannotty/Documents/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/openbt/R/polynomials.R")
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Samplers/sampler_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Optimization/optimization_helpers.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Plotting/computer_expt_plots.R")

library(plotly)
library(viridis)
library(latex2exp)
library(geoR)

# Set outfile directory
filedir = '/home/johnyannotty/Documents/RandomPathBART/VariogramResults/1d_functions/'

#-------------------------------------------------
set.seed(2)
n_train = 35 #15
n_test = 200
xoff = 0.10
xmin = -4*pi
xmax = 4*pi
s = 0.1

# generate x's for field data
x_train = round(seq(xmin+xoff,xmax-xoff, length = n_train) + runif(n_train, min = -xoff, max = xoff),4)
x_test = round(seq(xmin-2*pi,xmax+2*pi, length = n_test),4)

# Generate model set
f1_train = 2*sin(x_train-pi)-0.5*(x_train+2*pi) + 4
f2_train = 3.5*cos(1.1*(x_train)) - 2.5
f3_train = 2*cos(0.5*(x_train-pi))-6-3*pi/2 + sapply(x_train,function(x) max(0,pi-x)) #0.2*(x_test-3*pi)^2-12

f1_test = 2*sin(x_test-pi)-0.5*(x_test+2*pi) + 4
f2_test = 3.5*cos(1.1*(x_test)) - 2.5
f3_test = 2*cos(0.5*(x_test-pi))-6-3*pi/2 + sapply(x_test,function(x) max(0,pi-x)) #0.2*(x_test-3*pi)^2-12

# True function
fdagger = function(x){
  f = ifelse(x< (-2)*pi, 2.1*sin(x-pi)-0.55*(x+2*pi)+4, 
    ifelse(x<pi, 4*cos(x) -0.5*(x+2*pi),
           2*cos(0.5*(x-pi))-6-3*pi/2
  ))
  return(f)
}

# Get f_train, f_test, and fdagger
f_train = cbind(f1=f1_train, f2=f2_train, f3=f3_train)
f_test = cbind(f1=f1_test, f2=f2_test, f3=f3_test)
f0_train = fdagger(x_train)
f0_test = fdagger(x_test)

# Get y's
set.seed(1)
y_train = f0_train + rnorm(n_train, 0, s)

# Plot the field obs data
plot(x_train, y_train, pch = 16, cex = 1.1, main = 'Example Model Set', panel.first = {grid(col = 'lightgrey')},
     ylim = c(-13,10), ylab = "F(x)", xlab = "X")
lines(x_test, f0_test)
lines(x_test, f_test[,1], col = 'red',lwd = 2)
lines(x_test, f_test[,2], col = 'blue',lwd = 2)
lines(x_test, f_test[,3], col = 'green',lwd = 2)

#-------------------------------------------------
# Fit Model
nu = 30
rho = 1
sig2_hat = max(apply(apply(f_train, 2, function(x) (x-y_train)^2),2,min))
lam = rho*sig2_hat*(nu+2)/nu
#sig2_hat = nu/(nu+2)*0.01
q0 = 4
fit=openbt(x_train,y_train,f_train,pbd=c(1.0,0),ntree = 10,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="physics_model",
           ndpost = 10000, nskip = 2000, nadapt = 4000, adaptevery = 500, printevery = 500,
           power = 2.0, base = 0.95, minnumbot = 2, overallsd = sqrt(sig2_hat), k = 1, overallnu = nu,
           summarystats = FALSE, selectp = FALSE, rpath = TRUE, q = q0, rshp1 = 2, rshp2 = 5,
           stepwpert = 0.1, probchv = 0.1)

#Get mixed mean function
fitp=predict.openbt(fit,x.test = x_test, f.test = f_test,tc=4, q.lower = 0.025, q.upper = 0.975)

plot(x_test, f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(-15,15))
points(x_train, y_train, pch = 3)
lines(x_test, f_test[,1], col = 'red', lty = 2)
lines(x_test, f_test[,2], col = 'blue', lty = 2)
lines(x_test, f_test[,3], col = 'green4', lty = 2)
lines(x_test, fitp$mmean, col = 'purple', lwd = 2)
lines(x_test, fitp$m.lower, col = 'orange', lwd = 2, cex = 0.5)
lines(x_test, fitp$m.upper, col = 'orange', lwd = 2, cex = 0.5)

#Get the model weights
K = ncol(f_train)
fitw = mixingwts.openbt(fit, x.test = x_test, numwts = K, tc = 4)

#Plot model weights
plot(x_test, fitw$wmean[,1], pch = 16, col = 'red', type = 'l', ylim = c(-1,2.0), lwd = 2, 
     panel.first = {grid(col = 'lightgrey')})
lines(x_test, fitw$wmean[,2], col = 'blue', lwd = 2)
lines(x_test, fitw$wmean[,3], col = 'green3', lwd = 2)
lines(x_test, fitw$w.upper[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, fitw$w.lower[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, fitw$w.upper[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(x_test, fitw$w.lower[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(x_test, fitw$w.upper[,3], col = 'green3', lty = 'dashed', lwd = 1)
lines(x_test, fitw$w.lower[,3], col = 'green3', lty = 'dashed', lwd = 1)
abline(h = 1, col = 'grey', lty = 'dashed')
abline(h = 0, col = 'grey', lty = 'dashed')


# Sum of the weights
wsum = 0*fitw$wdraws[[1]]
for(j in 1:K){
  wsum = wsum + fitw$wdraws[[j]]
}
wsum_mean = apply(wsum,2,mean)
wsum_lb = apply(wsum,2,quantile, 0.025)
wsum_ub = apply(wsum,2,quantile, 0.975)


# Difference of the weights
wdiff = vector("list",length = 3)
for(j in 1:K){
  wd_mean = matrix(0,nrow = length(x_test),ncol = K-1)
  wd_lb = matrix(0,nrow = length(x_test),ncol = K-1)
  wd_ub = matrix(0,nrow = length(x_test),ncol = K-1)
  iter = 1
  for(l in 1:K){
    if(l!=j){
      wdtemp = fitw$wdraws[[j]] - fitw$wdraws[[l]]
      wd_mean[,iter] = apply(wdtemp,2,mean) 
      wd_ub[,iter] = apply(wdtemp,2,quantile, 0.975)
      wd_lb[,iter] = apply(wdtemp,2,quantile, 0.025)
      iter = iter + 1
    }
    wdiff[[j]]$mean = wd_mean
    wdiff[[j]]$lb = wd_lb
    wdiff[[j]]$ub = wd_ub
  }
}

# Plot sum of the weights
plot(x_test, wsum_mean, type = 'l', ylim = c(0,2))
lines(x_test, wsum_lb)
lines(x_test, wsum_ub)

# Plot difference in the weights
plot(x_test, wdiff[[1]]$mean[,2], type = 'l', ylim = c(-2,2))
lines(x_test, wdiff[[1]]$lb[,2])
lines(x_test, wdiff[[1]]$ub[,2])

# Save results
fit_data = list(
  pred_mean = fitp$mmean,
  pred_ub = fitp$m.upper,
  pred_lb = fitp$m.lower,
  wts_mean = fitw$wmean,
  wts_ub = fitw$w.upper,
  wts_lb = fitw$w.lower,
  wsum_mean = wsum_mean,
  wsum_lb = wsum_lb,
  wsum_ub = wsum_ub,
  #wdiff = wdiff,
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

saveRDS(fit_data, paste0(filedir,"pwtrig3_results4_11_09_23.rds"))

#-------------------------------------------------
# Theoretical Variogram
#-------------------------------------------------
# Theoretical 
xbnds = c(xmin,xmax)
h_grid = seq(0.1,24,by = 0.1)
vg = variogram.openbtmixing(xbnds,h_grid,10000,1,
                            k=0.75,
                            0.95,
                            power = 1,
                            a1 = 2,
                            a2 = 20,
                            4)

plot(h_grid,vg$vmean, type = "l", ylim = c(0,1))
points(geor_vghat$u,vghat_mean)

# Theoretical grid search
grid = expand.grid(a1=c(5,10),a2=c(10,15,20),beta=c(1,2),k=c(0.75,0.85,0.95,1))
vg_mat = matrix(0,nrow = nrow(grid), ncol = length(h_grid))
xbnds = c(xmin,xmax)
for(i in 1:nrow(grid)){
  vg = variogram.openbtmixing(xbnds, h_grid,10000,1,
                              k=grid[i,"k"],
                              0.95,
                              power = grid[i,"beta"],
                              a1 = grid[i,"a1"],
                              a2 = grid[i,"a2"],
                              4)
  vg_mat[i,] = vg$vmean
  cat("Progress: i = ", round(i/nrow(grid),4),"\r")
}

# Plot variogram
plot(h_grid,vg_mat[1,], type = "l", ylim = c(0,1))
lines(h_grid,vg_mat[2,], col = 'red')
lines(h_grid,vg_mat[3,], col = 'blue')
points(geor_vghat$u,vghat_mean)

plot(h_grid,vg_mat[40,], type = "l", ylim = c(0,1))
lines(h_grid,vg_mat[41,], col = 'red')
lines(h_grid,vg_mat[42,], col = 'blue')
points(geor_vghat2$u,vghat_mean)


# Save results
out = list(vg = vg_mat,x = x_train,y = y_train,h = h_grid, pgrid = grid)
saveRDS(out, paste0(filedir,"pwtrig3_thr_variogram_wt_11_15_23.rds"))


#-------------------------------------------------
# Precision Weights for Variogram
#-------------------------------------------------
# Precision Weight calculation
pw_se = apply(f_train,2,function(x) (y_train-x)^2)
pw_denom = rowSums(1/pw_se)
pw = (1/pw_se)/pw_denom
rowSums(pw)

# Get nearest neightbors
knn = 2
d = abs(outer(x_train,x_train, "-"))
#d = outer(x_train[,1],x_train[,1], "-")^2 + outer(x_train[,2],x_train[,2], "-")^2
nn = t(matrix(apply(d,1,function(x) order(x)[2:(knn+1)]),nrow = knn))
pw_nn = 0*pw

for(i in 1:n_train){
  ind = c(i,nn[i,])
  pw_se_nn = apply(pw_se[ind,],2,sum)
  pw_nn[i,] = (1/pw_se_nn)/sum(1/pw_se_nn)
}

rowSums(pw_nn)
plot(pw[,2])
plot(pw_nn[,1])

xdata = cbind(x1 = x_train,x2 = rep(0,n_train))
vghat = c() 
vghat_nn = c()
#h_grid = seq(0.9,max(x_train)-min(x_train),by = 1)
ugrid = seq(0.5,max(x_train)-min(x_train),by = 1)
K = ncol(f_train)
for(j in 1:K){
  geor_vghat = variog(coords = xdata, data = pw[,j],uvec = ugrid)
  vgh = 2*geor_vghat$v
  vghat = cbind(vghat,vgh)
  
  geor_vghat = variog(coords = xdata, data = pw_nn[,j],uvec = ugrid)
  vgh = 2*geor_vghat$v
  vghat_nn = cbind(vghat_nn,vgh)
}

colnames(vghat) = paste0("vw",1:K)
colnames(vghat_nn) = paste0("vw",1:K)


# Plot the empirical variogram for the mean weight
vghat_mean = rowMeans(vghat)
plot(geor_vghat$u,vghat_mean)

plot(geor_vghat$u,vghat[,2])
abline(h = 2*var(pw[,2]))

#geor_vghat = variog(coords = xdata, data = f_train[,2],uvec = h_grid)
#geor_vghat = variog(coords = xdata, data = f_train[,1], uvec=h_grid[seq(1,55,by = 2)])
#geor_vghat = variog(coords = xdata, data = wresi, uvec=h_grid)
#plot(geor_vghat$u,2*geor_vghat$v)

# Save results
out = list(vghat = vghat,vghat_nn = vghat_nn,hgrid = geor_vghat$u, 
           what = pw, what_nn = pw_nn, x = x_train,y = y_train)
saveRDS(out, paste0(filedir,"pwtrig3_emp_variogram_precwts_11_15_23.rds"))


#-------------------------------------------------
# Batch Fit
#-------------------------------------------------
nu = 30
rho = 1
sig2_hat = max(apply(apply(f_train, 2, function(x) (x-y_train)^2),2,min))
lam = rho*sig2_hat*(nu+2)/nu
q0 = 4

#f_train = matrix(rep(1,n_train), ncol = 1)
#f_test = matrix(rep(1,n_test), ncol = 1)
fit=train.openbtmixing(x_train,y_train,f_train,pbd=c(1.0,0),ntree = 10,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="physics_model",
                       ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
                       power = 1.0, base = 0.95, minnumbot = 2, overallsd = sqrt(sig2_hat), k = 1.5, overallnu = nu,
                       summarystats = FALSE, rpath = TRUE, q = q0, rshp1 = 2, rshp2 = 10,
                       stepwpert = 0.1, probchv = 0.1, batchsize = 100)

#openbt.save(fit,fname = "/home/johnyannotty/Documents/temp_res.obt")


#Get mixed mean function
fitp=predict.openbtmixing(fit,x.test = x_test, f.test = f_test,tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "mean_and_proj", proj_type = "euclidean", temperature = 0.8)

max(abs(fitp$mmean - rowSums(fitp$wmean*f_test)))
ind = 8452
max(abs(fitp$mdraws[ind,] - (fitp$wdraws[[1]][ind,]*f_test[,1]+fitp$wdraws[[2]][ind,]*f_test[,2]+
                             fitp$wdraws[[3]][ind,]*f_test[,3])))


max(abs(fitp$pmdraws[ind,] - (fitp$pwdraws[[1]][ind,]*f_test[,1]+fitp$pwdraws[[2]][ind,]*f_test[,2]+
                               fitp$pwdraws[[3]][ind,]*f_test[,3])))


# Sum of the weights
K= ncol(f_train)
wsum = 0*fitp$wdraws[[1]]
for(j in 1:K){
  wsum = wsum + fitp$wdraws[[j]]
}
wsum_mean = apply(wsum,2,mean)
wsum_lb = apply(wsum,2,quantile, 0.025)
wsum_ub = apply(wsum,2,quantile, 0.975)

# Predictions
plot(x_test, f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(-15,15))
points(x_train, y_train, pch = 3)
lines(x_test, f_test[,1], col = 'red', lty = 2)
lines(x_test, f_test[,2], col = 'blue', lty = 2)
lines(x_test, f_test[,3], col = 'green3', lty = 2)
lines(x_test, fitp$mmean, col = 'purple', lwd = 2)
lines(x_test, fitp$m.lower, col = 'orange', lwd = 2, cex = 0.5)
lines(x_test, fitp$m.upper, col = 'orange', lwd = 2, cex = 0.5)

#Plot model weights
plot(x_test, fitp$wmean[,1], pch = 16, col = 'red', type = 'l', ylim = c(-1,2.0), lwd = 2, 
     panel.first = {grid(col = 'lightgrey')})
lines(x_test, fitp$wmean[,2], col = 'blue', lwd = 2)
lines(x_test, fitp$wmean[,3], col = 'green3', lwd = 2)
lines(x_test, fitp$w.upper[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, fitp$w.lower[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, fitp$w.upper[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(x_test, fitp$w.lower[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(x_test, fitp$w.upper[,3], col = 'green3', lty = 'dashed', lwd = 1)
lines(x_test, fitp$w.lower[,3], col = 'green3', lty = 'dashed', lwd = 1)

abline(h = 1, col = 'grey', lty = 'dashed')
abline(h = 0, col = 'grey', lty = 'dashed')

# Plot sigma
hist(fitp$sdraws[,1])

#Get projected mixed mean function
plot(x_test, f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(-15,15))
points(x_train, y_train, pch = 3)
lines(x_test, f_test[,1], col = 'red', lty = 2)
lines(x_test, f_test[,2], col = 'blue', lty = 2)
lines(x_test, f_test[,3], col = 'green3', lty = 2)
lines(x_test, fitp$pmmean, col = 'purple', lwd = 2)
lines(x_test, fitp$pm.lower, col = 'orange', lwd = 2, cex = 0.5)
lines(x_test, fitp$pm.upper, col = 'orange', lwd = 2, cex = 0.5)

#Plot model weights
plot(x_test, fitp$pwmean[,1], pch = 16, col = 'red', type = 'l', ylim = c(-1,2.0), lwd = 2, 
     panel.first = {grid(col = 'lightgrey')})
lines(x_test, fitp$pw.5[,2], col = 'blue', lwd = 2)
lines(x_test, fitp$pw.5[,3], col = 'green3', lwd = 2)
lines(x_test, fitp$pw.upper[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, fitp$pw.lower[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, fitp$pw.upper[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(x_test, fitp$pw.lower[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(x_test, fitp$pw.upper[,3], col = 'green3', lty = 'dashed', lwd = 1)
lines(x_test, fitp$pw.lower[,3], col = 'green3', lty = 'dashed', lwd = 1)

rowSums(fitp$pwmean)


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

# Projection Results
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

filedir = "/home/johnyannotty/Documents/Dissertation/results/1d_functions/pwtrig3/"
saveRDS(fit_data, paste0(filedir,"pwtrig3_results4_12_12_23.rds"))
saveRDS(wts_l2_proj, paste0(filedir,"bmm_wl2_proj_01_17_24.rds"))


#------------------------------------------------
# Manual Softmax projection
#------------------------------------------------
library(rBayesianOptimization)
softmax_l2 = function(tmp){
  fita = rowSums(f_test*exp(fitp$wmean/tmp)/rowSums(exp(fitp$wmean/tmp)))
  score = -sum((fitp$mmean - fita)^2) # max the negative loss
  return(list(Score = score, Pred = 0))
} 

tmp_bounds = list(tmp = c(0,1))
init_grid_dt = data.frame(tmp = seq(0.01,0.95, length = 10))
bayes_temp = BayesianOptimization(FUN = softmax_l2, acq = "ei",
                                  bounds = tmp_bounds, init_grid_dt = init_grid_dt,
                                  init_points = 0.8, n_iter = 10)


# Plot the loss
tmp_grid = seq(0.01,0.99,length = 100)
plot(tmp_grid, sapply(tmp_grid, function(x) softmax_l2(x)$Score),ylim = c(-45,-30))
points(bayes_temp$Best_Par,bayes_temp$Best_Value,pch = 3, col = "red")

# Set the temperature
tmp = 0.20 #bayes_temp$Best_Par # from bayes opt
K = ncol(f_test)
pwts = list()
wts_array = array(0,dim = c(dim(fitp$wdraws[[1]]),K))

# Create new array
for(j in 1:K){
  wts_array[,,j] = fitp$wdraws[[j]]/tmp
}

# Get the max
wtsmax_tmp = apply(wts_array,c(1,2),max)

# Shift wts_array by the max's
for(j in 1:K){
  wts_array[,,j] = wts_array[,,j] - wtsmax_tmp
}

# Project Project Project !!
pwts_draws = wts_array*0
for(j in 1:K){
  pwts_draws[,,j] = exp(wts_array[,,j])
}

# Normalize
denom = apply(pwts_draws, c(1,2), sum)
for(j in 1:K){
  pwts_draws[,,j] = pwts_draws[,,j]/denom
}

# Get the predictions
pf_draws = fitp$mdraws*0
for(j in 1:K){
  pf_draws = pf_draws + t(t(pwts_draws[,,j])*f_test[,j])
}

# Get the discrepancy
delta_draws = fitp$mdraws - pf_draws

# Summary Stats
pf_mean = apply(pf_draws,2,mean)
pf_ub = apply(pf_draws,2,quantile,0.025)
pf_lb = apply(pf_draws,2,quantile,0.975)

delta_mean = apply(delta_draws,2,mean)
delta_lb = apply(delta_draws,2,quantile,0.025)
delta_ub = apply(delta_draws,2,quantile,0.975)

pwts_mean = matrix(0, nrow = nrow(f_test), ncol = K)
pwts_lb = matrix(0, nrow = nrow(f_test), ncol = K)
pwts_ub = matrix(0, nrow = nrow(f_test), ncol = K)
for(j in 1:K){
  pwts_mean[,j] = apply(pwts_draws[,,j],2,mean)
  pwts_lb[,j] = apply(pwts_draws[,,j],2,quantile, 0.025)
  pwts_ub[,j] = apply(pwts_draws[,,j],2,quantile, 0.975)
}


#Get projected mixed mean function
plot(x_test, f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(-15,15))
points(x_train, y_train, pch = 3)
lines(x_test, f_test[,1], col = 'red', lty = 2)
lines(x_test, f_test[,2], col = 'blue', lty = 2)
lines(x_test, f_test[,3], col = 'green3', lty = 2)
lines(x_test, pf_mean, col = 'purple', lwd = 2)
lines(x_test, pf_lb, col = 'orange', lwd = 2, cex = 0.5)
lines(x_test, pf_ub, col = 'orange', lwd = 2, cex = 0.5)

#Plot model weights
plot(x_test, pwts_mean[,1], pch = 16, col = 'red', type = 'l', ylim = c(-1,2.0), lwd = 2, 
     panel.first = {grid(col = 'lightgrey')})
lines(x_test, pwts_mean[,2], col = 'blue', lwd = 2)
lines(x_test, pwts_mean[,3], col = 'green3', lwd = 2)
lines(x_test, pwts_ub[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, pwts_lb[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, pwts_ub[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(x_test, pwts_lb[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(x_test, pwts_ub[,3], col = 'green3', lty = 'dashed', lwd = 1)
lines(x_test, pwts_lb[,3], col = 'green3', lty = 'dashed', lwd = 1)


# Plot the discrepancy
plot(x_test, delta_mean, pch = 16, cex = 0.8, main = 'Delta', type = 'l',ylim = c(-15,15))
lines(x_test, delta_lb, col = 'orange', lwd = 2, cex = 0.5)
lines(x_test, delta_ub, col = 'orange', lwd = 2, cex = 0.5)

softmax_res = list(
  pred_mean = pf_mean,
  pred_lb = pf_lb,
  pred_ub = pf_ub,
  wts_mean = pwts_mean,
  wts_lb = pwts_lb,
  wts_ub = pwts_ub,
  delta_mean = delta_mean,
  delta_lb = delta_lb,
  delta_ub = delta_ub,
  tmp = tmp 
)

filedir = "/home/johnyannotty/Documents/Dissertation/results/1d_functions/pwtrig3/"
saveRDS(softmax_res, paste0(filedir,"bmm_softmax_proj6_01_17_24.rds"))

#------------------------------------------------
# Non-negative Projection
#------------------------------------------------
# Draws
# Nonnegative Porjection.....!
pwts_draws = array(0, dim = c(dim(fitp$wdraws[[1]]), ncol(f_test)))
for(j in 1:K){
  pwts_draws[,,j] = (abs(fitp$wdraws[[j]]) + fitp$wdraws[[j]])/2
}

# Get the predictions
pf_draws = fitp$mdraws*0
for(j in 1:K){
  pf_draws = pf_draws + t(t(pwts_draws[,,j])*f_test[,j])
}

# Get the discrepancy
delta_draws = fitp$mdraws - pf_draws

# Summary Stats
pf_mean = apply(pf_draws,2,mean)
pf_ub = apply(pf_draws,2,quantile,0.025)
pf_lb = apply(pf_draws,2,quantile,0.975)

delta_mean = apply(delta_draws,2,mean)
delta_lb = apply(delta_draws,2,quantile,0.025)
delta_ub = apply(delta_draws,2,quantile,0.975)

pwts_mean = matrix(0, nrow = nrow(f_test), ncol = K)
pwts_lb = matrix(0, nrow = nrow(f_test), ncol = K)
pwts_ub = matrix(0, nrow = nrow(f_test), ncol = K)
for(j in 1:K){
  pwts_mean[,j] = apply(pwts_draws[,,j],2,mean)
  pwts_lb[,j] = apply(pwts_draws[,,j],2,quantile, 0.025)
  pwts_ub[,j] = apply(pwts_draws[,,j],2,quantile, 0.975)
}


# Save the non-negative projection results
nonneg_res = list(
  pred_mean = pf_mean,
  pred_lb = pf_lb,
  pred_ub = pf_ub,
  wts_mean = pwts_mean,
  wts_lb = pwts_lb,
  wts_ub = pwts_ub,
  delta_mean = delta_mean,
  delta_lb = delta_lb,
  delta_ub = delta_ub
)

filedir = "/home/johnyannotty/Documents/Dissertation/results/1d_functions/pwtrig3/"
saveRDS(softmax_res, paste0(filedir,"bmm_nonneg_proj_01_17_24.rds"))


#------------------------------------------------
# Bounds on the sum of wts
#------------------------------------------------
e_train = matrix(y_train, ncol = 3, nrow = n_train, byrow = FALSE) - f_train
xind = c(5:6)
mbest = which.min(colSums(e_train[xind,]^2))

Delta_bnd = sqrt(sum((3*e_train[xind,mbest])^2))
err_bnd = sqrt(sum(e_train[xind,mbest]^2))
ynorm = sqrt(sum(y_train[xind]^2))
B = (Delta_bnd + err_bnd)/ynorm
c(1-B,1+B)


#------------------------------------------------
# Discrepancy and sum-of-wts
#------------------------------------------------
delta_test = matrix(f0_test, ncol = 3, nrow = n_test, byrow = FALSE) - f_test
wgrid = expand.grid(seq(-0.4,1.4, length = 60),
                    seq(-0.4,1.4, length = 60))

# Fix a wt 
wgrid = cbind(rep(0.4,nrow(wgrid)),wgrid)
wgrid = wgrid[which(rowSums(wgrid)<1.4),]
colnames(wgrid) = c("w1","w2","w3")

bias = 0
for(i in 1:nrow(wgrid)){
  if(length(xind)>1){
    bias[i] = mean(t(t(delta_test[xind,])*wgrid[i,])) + mean(f0_test[xind]*(1-sum(wgrid[i,])))  
  }else{
    bias[i] = sum(delta_test[xind,]*wgrid[i,]) + f0_test[xind]*(1-sum(wgrid[i,]))
  }
}
  
pbias = plot_mean2d_viridis(wgrid[,2:3],bias^2, 
                                  viridis_opt = "viridis",
                                  scale_limit = c(0,25), title = "Bias^2") 

hp = which(wgrid[,2]>0 & wgrid[,2]<0.6 & wgrid[,3]>0 & wgrid[,3]<0.6)
pbias_sub = plot_mean2d_viridis(wgrid[hp,2:3],bias[hp]^2, 
                            viridis_opt = "viridis",
                            scale_limit = c(0,15), title = "Bias^2") 
hist(bias^2)
max(bias[hp]^2)

# Fix a wt 
agrid = expand.grid(seq(0,1, length = 60),
                    seq(0,1, length = 60))

agrid = cbind(rep(0.4,nrow(agrid)),agrid)
agrid = agrid[which(rowSums(agrid)<=1 & rowSums(agrid)>0),]
colnames(agrid) = c("w1","w2","w3")

abias = 0
for(i in 1:nrow(agrid)){
  if(length(xind)>1){
    abias[i] = mean(t(t(delta_test)*agrid[i,])) + mean(f0_test*(1-sum(agrid[i,])))    
  }else{
    abias[i] = sum(delta_test[xind,]*agrid[i,]) + f0_test[xind]*(1-sum(agrid[i,]))
  }
}

pabias = plot_mean2d_viridis(agrid[,2:3],abias^2, 
                            viridis_opt = "viridis",
                            scale_limit = c(0,15), title = "Bias^2") 

hist(abias^2)
max(abias^2)
