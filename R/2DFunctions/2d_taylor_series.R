#-----------------------------------------------------
# 2D-Taylor Series Expansions
#-----------------------------------------------------
setwd("/home/johnyannotty/Documents/openbt/src")

# Load the R wrapper functions to the OpenBT library.
source("/home/johnyannotty/Documents/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/openbt/R/polynomials.R")
source("/home/johnyannotty/Documents/openbt/R/mixing_priors.R")
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Samplers/sampler_functions.R")

library(plotly)
library(viridis)
library(latex2exp)
library(geoR)

# Out file directory
filedir = '/home/johnyannotty/Documents/RandomPathBART/VariogramResults/2d_functions/'

#-----------------------------------------------------
#-----------------------------------------------------
set.seed(2)
n_train = 80 
n_test = 50
xoff = 0.10
xmin = -pi
xmax = pi
s = 0.1

# Generate data
x_train = grid_2d_design(n1=10,n2=8, xmin = c(xmin,xmin), xmax = c(xmax,xmax))
#x_train = grid_2d_design(n1=25,n2=20, xmin = c(xmin,xmin), xmax = c(xmax,xmax))
xt_min = apply(x_train,2,min);  xt_max = apply(x_train,2,max); 
x1_test = round(seq(xt_min[1]-0.1,xt_max[1]+0.1,length = n_test),4)
x2_test = round(seq(xt_min[1]-0.1,xt_max[1]+0.1,length = n_test),4)
x_test = expand.grid(x1_test,x2_test)

#x_test = grid_2d_design(n1=40,n2=40,n=n_test, xmin = xt_min-0.2, xmax = xt_max + 0.2)
#x_test = grid_2d_testset(x_train,d = 0.8,nbd = 20)

plot(x_train[,1],x_train[,2])
plot(x_test[,1],x_test[,2], pch = 16, col = 'grey70')
points(x_train[,1],x_train[,2], pch = 16, col = 'blue')

# Generate True function
f0_train = sin(x_train[,1]) + cos(x_train[,2])
f0_test = sin(x_test[,1]) + cos(x_test[,2])

set.seed(99)
y_train = f0_train + rnorm(n_train,0,s)

ks1 = 7
kc1 = 10 #6
ks2 = 13 #9
kc2 = 6

h1_sinexp = sapply(x_train[,1],function(x) sin_exp(k=ks1,x0=pi,x))
h1_cosexp = sapply(x_train[,2],function(x) cos_exp(k=kc1,x0=pi,x))
h1_train = h1_sinexp + h1_cosexp
h1_test = sapply(x_test[,1],function(x) sin_exp(k=ks1,x0=pi,x)) + sapply(x_test[,2],function(x) cos_exp(k=kc1,x0=pi,x))

h2_sinexp = sapply(x_train[,1],function(x) sin_exp(k=ks2,x0=-pi,x))
h2_cosexp = sapply(x_train[,2],function(x) cos_exp(k=kc2,x0=-pi,x))
h2_train = h2_sinexp + h2_cosexp
h2_test = sapply(x_test[,1],function(x) sin_exp(k=ks2,x0=-pi,x)) + sapply(x_test[,2],function(x) cos_exp(k=kc2,x0=-pi,x))

f_train = cbind(h1_train, h2_train)
f_test = cbind(h1_test, h2_test)


#-----------------------------------------------------
# Model Training 
#-----------------------------------------------------
nu = 20
rho = 1
sig2_hat = max(apply(apply(f_train, 2, function(x) (x-y_train)^2),2,min))
lam = rho*sig2_hat*(nu+2)/nu
q0 = 4
fit=openbt(x_train,y_train,f_train,pbd=c(1.0,0),ntree = 15,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="physics_model",
           ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
           power = 2.0, base = 0.95, minnumbot = 3, overallsd = sqrt(sig2_hat), k = 3, overallnu = nu,
           summarystats = FALSE, selectp = FALSE, rpath = TRUE, q = 4.0, rshp1 = 15, rshp2 = 15,
           stepwpert = 0.1, probchv = 0.1)


#Get mixed mean function
fitp=predict.openbt(fit,x.test = x_test, f.test = f_test,tc=4, q.lower = 0.025, q.upper = 0.975)
p1 = plot_pred_2d_gg2(x_test,matrix(fitp$mmean),scale_vals = c(-2.7,-0.25,2.05), title = "Mean Prediction from BMM")

#Get the model weights
K = 2
fitw = mixingwts.openbt(fit, x.test = x_test, numwts = K, tc = 4)
w1 = plot_wts_2d_gg2(x_test,fitw$wmean,wnum = 1,xcols = c(1,2), title="Mean Weights", scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.1,0.5,1.1))
w2 = plot_wts_2d_gg2(x_test,fitw$wmean,wnum = 2,xcols = c(1,2), title="Mean Weights", scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.1,0.5,1.1))

# Get wsum
wsum = 0*fitw$wdraws[[1]]
for(i in 1:K){
  wsum = wsum + fitw$wdraws[[i]]
}
wsum_mean = apply(wsum,2,mean)
wsum_lb = apply(wsum,2,quantile, 0.025)
wsum_ub = apply(wsum,2,quantile, 0.975)

# Difference of the weights
wdiff = vector("list",length = K)
for(j in 1:K){
  wd_mean = matrix(0,nrow = nrow(x_test),ncol = K-1)
  wd_lb = matrix(0,nrow = nrow(x_test),ncol = K-1)
  wd_ub = matrix(0,nrow = nrow(x_test),ncol = K-1)
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
  wdiff = wdiff,
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

saveRDS(fit_data, paste0(filedir,"taylor_sincos2_results3_10_30_23.rds"))


#-----------------------------------------------------
# Precision Weights Empirical Variogram
#-----------------------------------------------------
# Precision Weight calculation
pw_se = apply(f_train,2,function(x) (y_train-x)^2)
pw_denom = rowSums(1/pw_se)
pw = (1/pw_se)/pw_denom
rowSums(pw)

# Get nearest neightbors
knn = 3
d = outer(x_train[,1],x_train[,1], "-")^2 + outer(x_train[,2],x_train[,2], "-")^2
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

#xdata = cbind(x1 = x_train[,1],x2 = rep(0,n_train))
xdata = cbind(x1 = x_train[,1],x2 = x_train[,2])
vghat = c() 
vghat_nn = c()
ugrid = seq(0.3,7.5,by = 0.3)
K = ncol(f_train)
for(j in 1:K){
  geor_vghat = variog(coords = xdata, data = pw[,j])#,uvec = ugrid)
  vgh = 2*geor_vghat$v
  vghat = cbind(vghat,vgh)
  
  geor_vghat = variog(coords = xdata, data = pw_nn[,j])#,uvec = ugrid)
  vgh = 2*geor_vghat$v
  vghat_nn = cbind(vghat_nn,vgh)
}

colnames(vghat) = paste0("vw",1:K)
colnames(vghat_nn) = paste0("vw",1:K)


# Plot the empirical variogram for the mean weight
vghat_mean = rowMeans(vghat)
plot(geor_vghat$u,vghat_mean)

plot(geor_vghat$u,vghat[,1])
abline(h = 2*var(pw_nn[,1]))

# Save results
out = list(vghat = vghat,vghat_nn = vghat_nn,hgrid = geor_vghat$u, 
           what = pw, what_nn = pw_nn, x = x_train,y = y_train)
saveRDS(out, paste0(filedir,"taylor_sincos2_emp_variogram_precwts_11_15_23.rds"))


#-----------------------------------------------------
# Theoretical Variogram
#-----------------------------------------------------
xbnds = matrix(c(xmin,xmax),nrow = 2,ncol = 2, byrow = TRUE)
h_grid = seq(0.1,6.5,by = 0.1)
vg = variogram.openbtmixing(xbnds,h_grid,10000,1,
                            k=0.8,
                            0.95,
                            power = 2,
                            a1 = 10,
                            a2 = 10,
                            4)

plot(h_grid,vg$vmean, type = "l", ylim = c(0,0.8))
points(geor_vghat$u,vghat_mean)

# Theoretical grid search
grid = expand.grid(a1=c(5,10),a2=c(10,15,20),beta=c(2),k=c(0.6,0.7,0.8,1))
h_grid = seq(0.1,7.5,by = 0.1)
xbnds = matrix(c(xmin,xmax),nrow = 2,ncol = 2, byrow = TRUE)
vg_mat = matrix(0,nrow = nrow(grid), ncol = length(h_grid))
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
plot(h_grid,vg_mat[19,], type = "l", ylim = c(0,0.3))
lines(h_grid,vg_mat[20,], col = 'red')
lines(h_grid,vg_mat[21,], col = 'blue')
points(geor_vghat$u,vghat_mean)

plot(h_grid,vg_mat[1,], type = "l", ylim = c(0,0.3))
lines(h_grid,vg_mat[2,], col = 'red')
lines(h_grid,vg_mat[3,], col = 'blue')
points(geor_vghat$u,vghat_mean)

# Get the best fitting variogram
which.min(apply(vg_mat,1,function(x) sqrt(mean((x-vghat_mean)^2))))
order(apply(vg_mat,1,function(x) sqrt(mean((x-vghat_mean)^2))))

# Save results
out = list(vg = vg_mat,x = x_train,y = y_train,h = h_grid, pgrid = grid)
saveRDS(out, paste0(filedir,"taylor_sincos2_thr_variogram_wt_11_15_23.rds"))


#-----------------------------------------------------
# Original Pointwise variogram idea
#-----------------------------------------------------
tau2 = (1/(2*1))^2
h_grid = seq(0.1,6,by = 0.1)

wt_ptw = matrix(0,nrow = n_train,ncol = 2)
for(i in 1:n_train){
  wt_ptw[i,] = t(solve((f_train[i,]%*%t(f_train[i,]))/sig2_hat + 
                         1/tau2*diag(2))%*%(f_train[i,]*y_train[i]/sig2_hat + rep(0.5,2)/tau2))   
}

# Plot the pointwise weights
plot(x_train[,1],wt_ptw[,1], col = 'red', pch = 16)
points(x_train[,1],wt_ptw[,2], col = 'blue', pch = 16)

# Plot the empirical variogram for w2
xdata = x_train
vghat = c()
K = ncol(f_train)
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
filedir = '/home/johnyannotty/Documents/Dissertation/results/2d_functions/'
out = list(vghat = vghat, hgrid = geor_vghat$u, what = wt_ptw, x = x_train,y = y_train, h = h_grid)
saveRDS(out, paste0(filedir,"taylor_sincos2_emp_variogram_wt_10_30_23.rds"))


#-----------------------------------------------------
# Deterministic Fit
#-----------------------------------------------------
fit=openbt(x_train,y_train,f_train,pbd=c(0.7,0),ntree = 20,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="physics_model",
           ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
           power = 2.0, base = 0.95, minnumbot = 3, overallsd = sqrt(sig2_hat), k = 2.5, overallnu = nu,
           summarystats = FALSE, selectp = FALSE, rpath = FALSE, q = 4.0, rshp1 = 2, rshp2 = 15,
           stepwpert = 0.1, probchv = 0.1)


#Get mixed mean function
fitp=predict.openbt(fit,x.test = x_test, f.test = f_test,tc=4, q.lower = 0.025, q.upper = 0.975)

#Get the model weights
K = 2
fitw = mixingwts.openbt(fit, x.test = x_test, numwts = K, tc = 4)

# Get wsum
wsum = 0*fitw$wdraws[[1]]
for(i in 1:K){
  wsum = wsum + fitw$wdraws[[i]]
}
wsum_mean = apply(wsum,2,mean)
wsum_lb = apply(wsum,2,quantile, 0.025)
wsum_ub = apply(wsum,2,quantile, 0.975)
