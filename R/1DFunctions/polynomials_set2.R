#-------------------------------------------------
# Mixing polynomials - 2 Models in 1D
#-------------------------------------------------
setwd("/home/johnyannotty/Documents/openbt/src")

# Load the R wrapper functions to the OpenBT library.
source("/home/johnyannotty/Documents/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/openbt/R/polynomials.R")
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Samplers/sampler_functions.R")

library(plotly)
library(viridis)
library(latex2exp)
library(geoR)

# Set outfile directory
filedir = '/home/johnyannotty/Documents/RandomPathBART/VariogramResults/1d_functions/'

#-------------------------------------------------
# Control setup
set.seed(2)
n_train = 35 #15
n_test = 100
xoff = 0.10
xmin = 0
xmax = 8
s = 0.1

# generate x's for field data
x_train = round(seq(xmin+xoff,xmax-xoff, length = n_train) + runif(n_train, min = -xoff, max = xoff),4)
x_test = round(seq(xmin,xmax, length = n_test),4)

# get simulator parameters
h_list = list(model1 = list(a = 2, b = 0, c = 1.5, p = 2), model2 = list(a = 6, b = 0, c = 1.5, p = 2))
xc1_train = x_train
xc2_train = x_train

# get weight functions via softmax and set the bias
wts_basis = matrix(1, nrow = n_train, ncol = 2)
wts_basis[,2] = fp(x_train, a = 3,b = 0.75,p = 1, c = 3) 
wts_matrix = softmax(wts_basis)
bias_vec = rep(0,n_train)
#bias_vec = fp(x_train, a = 0, b = 0.25, p = 1, c = 0.5)

# get true system
f_train = wts_basis # easy initialization
f_train[,1] = fp(x_train, a=h_list$model1$a, b=h_list$model1$b, c=h_list$model1$c, p=h_list$model1$p)
f_train[,2] = fp(x_train, a=h_list$model2$a, b=h_list$model2$b, c=h_list$model2$c ,p=h_list$model2$p)
ftrue = get_ftrue(f_train,wts_matrix,bias_vec)
set.seed(1)
y_train = ftrue + rnorm(n_train, 0, s)

# Get test data for plotting
wts_basis_test = matrix(1, nrow = n_test, ncol = 2)
wts_basis_test[,2] = fp(x_test, a = 3,b = 0.75,p = 1, c = 3) 
wts_matrix_test = softmax(wts_basis_test)
bias_vec_test = rep(0,n_test)
#bias_vec_test = fp(x_test, a = 0, b = 0.25, p = 1, c = 0.5)
f_test = matrix(0,nrow = n_test, ncol = 2)
f_test[,1] = fp(x_test, a=h_list$model1$a, b=h_list$model1$b, c=h_list$model1$c, p=h_list$model1$p)
f_test[,2] = fp(x_test, a=h_list$model2$a, b=h_list$model2$b, c=h_list$model2$c ,p=h_list$model2$p)
ftrue_test = get_ftrue(f_test,wts_matrix_test,bias_vec_test)

# Plot the field obs data
plot(x_train, y_train, pch = 16, cex = 0.8, main = 'Field Obs', panel.first = {grid(col = 'lightgrey')},ylim = c(0,10))
lines(x_train, f_train[,1], col = 'red')
lines(x_train, f_train[,2], col = 'blue')
lines(x_train, ftrue, col = 'black')

# Setup objects
x_train = matrix(x_train, ncol = 1)
x_test = matrix(x_test, ncol = 1)

#-----------------------------------------------------
# Fit the model
#-----------------------------------------------------
nu = 30
rho = 1
sig2_hat = max(apply(apply(f_train, 2, function(x) (x-y_train)^2),2,min))
lam = rho*sig2_hat*(nu+2)/nu
q0 = 4
fit=openbt(x_train,y_train,f_train,pbd=c(1.0,0),ntree = 10,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="physics_model",
           ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
           power = 2.0, base = 0.95, minnumbot = 2, overallsd = sqrt(sig2_hat), k = 1, overallnu = nu,
           summarystats = FALSE, selectp = FALSE, rpath = TRUE, q = q0, rshp1 = 2, rshp2 = 10,
           stepwpert = 0.1, probchv = 0.1)

#Get mixed mean function
fitp=predict.openbt(fit,x.test = x_test, f.test = f_test,tc=4, q.lower = 0.025, q.upper = 0.975)

plot(x_test, ftrue_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(-1,10))
points(x_train, y_train, pch = 3)
lines(x_test, f_test[,1], col = 'red', lty = 2)
lines(x_test, f_test[,2], col = 'blue', lty = 2)
lines(x_test, fitp$mmean, col = 'green4', lwd = 2)
lines(x_test, fitp$m.lower, col = 'orange', lwd = 2, cex = 0.5)
lines(x_test, fitp$m.upper, col = 'orange', lwd = 2, cex = 0.5)

#Get the model weights
K = 2
fitw = mixingwts.openbt(fit, x.test = x_test, numwts = K, tc = 4)

#Plot model weights
plot(x_test, fitw$wmean[,1], pch = 16, col = 'red', type = 'l', ylim = c(-1,2.0), lwd = 2, 
     panel.first = {grid(col = 'lightgrey')})
lines(x_test, fitw$wmean[,2], col = 'blue', lwd = 2)
lines(x_test, fitw$w.upper[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, fitw$w.lower[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, fitw$w.upper[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(x_test, fitw$w.lower[,2], col = 'blue', lty = 'dashed', lwd = 1)
abline(h = 1, col = 'grey', lty = 'dashed')
abline(h = 0, col = 'grey', lty = 'dashed')

# Get posterior of wsum
wsum = 0*fitw$wdraws[[1]]
for(j in 1:K){
  wsum = wsum + fitw$wdraws[[j]]
}
wsum_mean = apply(wsum,2,mean)
wsum_lb = apply(wsum,2,quantile, 0.025)
wsum_ub = apply(wsum,2,quantile, 0.975)

# Plot sum of the weights
plot(x_test, wsum_mean, type = 'l', ylim = c(0,2))
lines(x_test, wsum_lb)
lines(x_test, wsum_ub)

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

saveRDS(fit_data, paste0(filedir,"quad2_results_11_15_23.rds"))


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
ugrid = seq(0.2,max(x_train)-min(x_train),by = 0.3)
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

plot(geor_vghat$u,vghat_nn[,2])
abline(h = 2*var(pw[,1]))

# Save results
out = list(vghat = vghat,vghat_nn = vghat_nn,hgrid = geor_vghat$u, 
           what = pw, what_nn = pw_nn, x = x_train,y = y_train)
saveRDS(out, paste0(filedir,"quad2_emp_variogram_precwts_11_15_23.rds"))


# Data variogram
y_vghat = variog(coords = xdata, data = f_train[,2],uvec = ugrid)
plot(geor_vghat$u,2*y_vghat$v)

#-------------------------------------------------
# Theoretical Variogram
#-------------------------------------------------
# Theoretical 
xbnds = c(xmin,xmax)
h_grid = seq(0.1,8,by = 0.1)
vg = variogram.openbtmixing(xbnds,h_grid,5000,1,
                            k=0.75,
                            0.95,
                            power = 1,
                            a1 = 2,
                            a2 = 20,
                            4)

plot(h_grid,vg$vmean, type = "l", ylim = c(0,1))
points(geor_vghat$u,vghat_mean)

# Theoretical grid search
grid = expand.grid(a1=c(2,10),a2=c(10,15,20),beta=c(2),k=seq(0.5,1,by=0.1))
h_grid = seq(0.1,8,by = 0.1)
vg_mat = matrix(0,nrow = nrow(grid), ncol = length(h_grid))
xbnds = c(xmin,xmax)
for(i in 1:nrow(grid)){
  vg = variogram.openbtmixing(xbnds, h_grid,5000,1,
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

plot(h_grid,vg_mat[29,], type = "l", ylim = c(0,1))
lines(h_grid,vg_mat[35,], col = 'red')
lines(h_grid,vg_mat[23,], col = 'blue')
lines(h_grid,vg_mat[17,], col = 'green3')
lines(h_grid,vg_mat[11,], col = 'yellow')
lines(h_grid,vg_mat[5,], col = 'purple')
points(geor_vghat$u,vghat_mean)

# Save results
out = list(vg = vg_mat,x = x_train,y = y_train,h = h_grid, pgrid = grid)
saveRDS(out, paste0(filedir,"quad2_thr_variogram_wt_11_15_23.rds"))


#-------------------------------------------------
# Batch Fit
#-------------------------------------------------
fit=train.openbtmixing(x_train,y_train,f_train,pbd=c(1.0,0),ntree = 10,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="physics_model",
           ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
           power = 2.0, base = 0.95, minnumbot = 2, overallsd = sqrt(sig2_hat), k = 1, overallnu = nu,
           summarystats = FALSE, rpath = TRUE, q = q0, rshp1 = 2, rshp2 = 10,
           stepwpert = 0.1, probchv = 0.1, batchsize = 1000)

#Get mixed mean function
fitp=predict.openbtmixing(fit,x.test = x_test, f.test = f_test,tc=4, q.lower = 0.025, q.upper = 0.975)

plot(x_test, ftrue_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(-1,10))
points(x_train, y_train, pch = 3)
lines(x_test, f_test[,1], col = 'red', lty = 2)
lines(x_test, f_test[,2], col = 'blue', lty = 2)
lines(x_test, fitp$mmean, col = 'green4', lwd = 2)
lines(x_test, fitp$m.lower, col = 'orange', lwd = 2, cex = 0.5)
lines(x_test, fitp$m.upper, col = 'orange', lwd = 2, cex = 0.5)

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
