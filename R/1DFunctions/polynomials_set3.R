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
           summarystats = FALSE, selectp = FALSE, rpath = TRUE, q = q0, rshp1 = 2, rshp2 = 10,
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



# De-trend it
outlm = lm(log(pw[,1]/(1-pw[,1]))~x_train)
wfit = exp(outlm$fitted.values)/(1 +  exp(outlm$fitted.values))
wresi = pw[,1] - wfit
plot(x_train, pw[,1], col = "black", pch = 16)
lines(x_train, wfit, col = "red", lwd = 2)

geor_vghat = variog(coords = xdata, data = wresi, uvec=h_grid)
plot(geor_vghat$u,2*geor_vghat$v)

# For climate models
plot(x_train[,1],x_train[,2], 
     col = ifelse(pw[,1]>0.75,"black",ifelse(pw[,1]>0.5,"darkred",
                                               ifelse(pw[,1]>0.25,"orange2","yellow"))), pch = 16)

plot(x_train[,1],x_train[,2], 
     col = ifelse(pw_nn[,1]>0.75,"black",ifelse(pw_nn[,1]>0.5,"darkred",
                                             ifelse(pw_nn[,1]>0.25,"orange2","yellow"))), pch = 16)


#-------------------------------------------------
# Old - Pointwise Weights via regularization for Variogram
#-------------------------------------------------
tau2 = (1/(2*1))^2
h_grid = seq(0.1,16,by = 0.1)

wt_ptw = matrix(0,nrow = n_train,ncol = 3)
for(i in 1:n_train){
  wt_ptw[i,] = t(solve((f_train[i,]%*%t(f_train[i,]))/sig2_hat + 
                         1/tau2*diag(3))%*%(f_train[i,]*y_train[i]/sig2_hat + rep(0.5,3)/tau2))   
}

# Plot the pointwise weights
plot(x_train,wt_ptw[,1], col = 'red', pch = 16)
points(x_train,wt_ptw[,2], col = 'blue', pch = 16)
points(x_train,wt_ptw[,3], col = 'green3', pch = 16)


# Plot the empirical variogram for w1
xdata = cbind(x1 = x_train,x2 = rep(0,n_train))
vghat = c() 
for(j in 1:K){
  geor_vghat = variog(coords = xdata, data = wt_ptw[,j],uvec = h_grid)
  vgh = 2*geor_vghat$v
  vghat = cbind(vghat,vgh)
}

colnames(vghat) = paste0("vw",1:K)

# Plot the empirical variogram for the mean weight
vghat_mean = rowMeans(vghat)
plot(geor_vghat$u,vghat_mean)

# Save results
out = list(vghat = vghat, hgrid = geor_vghat$u, what = wt_ptw, x = x_train,y = y_train)
saveRDS(out, paste0(filedir,"pwtrig3_emp_variogram_wt_10_30_23.rds"))


#-------------------------------------------------
# Batch Fit
#-------------------------------------------------
nu = 30
rho = 1
sig2_hat = max(apply(apply(f_train, 2, function(x) (x-y_train)^2),2,min))
lam = rho*sig2_hat*(nu+2)/nu
q0 = 4

fit=train.openbtmixing(x_train,y_train,f_train,pbd=c(1.0,0),ntree = 10,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="physics_model",
                       ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
                       power = 1.0, base = 0.95, minnumbot = 5, overallsd = sqrt(sig2_hat), k = 2, overallnu = nu,
                       summarystats = FALSE, rpath = TRUE, q = q0, rshp1 = 2, rshp2 = 10,
                       stepwpert = 0.1, probchv = 0.1, batchsize = 1000)

#Get mixed mean function
fitp=predict.openbtmixing(fit,x.test = x_test, f.test = f_test,tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "mean_and_sigma")

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
plot(x_test, fitp$pw.5[,1], pch = 16, col = 'red', type = 'l', ylim = c(-1,2.0), lwd = 2, 
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

xx = cbind(fitp$wdraws[[1]][1:500,1],fitp$wdraws[[2]][1:500,1],fitp$wdraws[[3]][1:500,1])
cbind(fitp$pwdraws[[1]][1:5,1],fitp$pwdraws[[2]][1:5,1],fitp$pwdraws[[3]][1:5,1])
rowSums(cbind(fitp$pwdraws[[1]][1:5,1],fitp$pwdraws[[2]][1:5,1],fitp$pwdraws[[3]][1:5,1]))

