#------------------------------------------------
# Specifc Heat Model Fitting
#------------------------------------------------
source("/home/johnyannotty/Documents/SciData/NuclearPhysics/R/heat_2d_ising.R")
#source("/home/johnyannotty/Documents/BayesToolBox/bayestb/LinearRegression/bayesian_regression.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/kernels.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/gp_utils.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/gp_regression.R")

#------------------------------------------------
# Generate Data
#------------------------------------------------
L = 8
n_train = 15
n_test = 30
s = 0.02
g_train = c(0.02,0.1,0.3,seq(0.02,5,length = n_train-3))
g_test = seq(0.01,5, length = n_test)
gvec = seq(0.01,5, length = 200)
x_train = 0.5*log(g_train+1)
x_test = 0.5*log(g_test+1)

fsg_train = heat_sg_exp(g_train,10,L)
flg_train = heat_lg_exp(g_train,10,L)

fsg_test = heat_sg_exp(g_test,10,L)
flg_test = heat_lg_exp(g_test,10,L)

fsg_grid = heat_sg_exp(gvec,10,L)
flg_grid = heat_lg_exp(gvec,10,L)


dz2_train = sapply(x_train,function(x) d2_logz(x,L,0.1))

y_train = dz2_train/L^2 + rnorm(n_train,0,s)
f0_test = sapply(x_test,function(x) d2_logz(x,L,0.1)/L^2)


#------------------------------------------------
# Define covariance
#------------------------------------------------
eft_sg_sqr_exp_wrapper = function(x1,x2,params,sig2 = TRUE){
  # EFT squared exponential kernel
  a = 0.5
  fxgrid = cbind(x1,x2)
  fxgrid[,1] = ifelse(fxgrid[,1] < a,fxgrid[,1],a)
  fxgrid[,2] = ifelse(fxgrid[,2] < a,fxgrid[,2],a)
  
  lscale = 1
  sscale = 0.05
  Rfx = apply(fxgrid, 1, function(x) sqr_exp_kernel(x[1],x[2],sscale,lscale))
  
  # EFT add in scalar functions
  yref = 1
  qx = 0.5
  N = 10
  Rfx = ((yref^2)*(1 - qx^2)^(N+1)/(1-qx^2))*Rfx
  
  # Intermediate discrepancy model
  gxgrid = cbind(x1,x2)
  Rgx = apply(gxgrid, 1, function(x) sqr_exp_kernel(x[1],x[2],params[1],params[2]))
  Rgx = sapply(1:nrow(gxgrid), function(i) if(gxgrid[i,1] < a || gxgrid[i,2] < a){return(0)}else{return(Rgx[i])}   )
  
  if(sig2){
    Rsx = apply(gxgrid, 1, function(x) ifelse(x[1] == x[2],params[3],0))
    Rx = Rsx + Rfx + Rgx
  }else{
    Rx = Rfx + Rgx
  }
  return(Rx)
}



eft_lg_sqr_exp_wrapper = function(x1,x2,params,sig2 = TRUE){
  # EFT squared exponential kernel
  b = 1.9
  fxgrid = cbind(x1,x2)
  fxgrid[,1] = ifelse(fxgrid[,1] > b,fxgrid[,1],b)
  fxgrid[,2] = ifelse(fxgrid[,2] > b,fxgrid[,2],b)
  
  lscale = 2
  sscale = 0.05
  Rfx = apply(fxgrid, 1, function(x) sqr_exp_kernel(x[1],x[2],sscale,lscale))
  
  # EFT add in scalar functions
  yref = 1
  qx = 0.5
  N = 10
  Rfx = ((yref^2)*(1 - qx^2)^(N+1)/(1-qx^2))*Rfx
  
  # Intermediate discrepancy model
  gxgrid = cbind(x1,x2)
  Rgx = apply(gxgrid, 1, function(x) sqr_exp_kernel(x[1],x[2],params[1],params[2]))
  Rgx = sapply(1:nrow(gxgrid), function(i) if(gxgrid[i,1] < a || gxgrid[i,2] < a){return(0)}else{return(Rgx[i])}   )
  
  if(sig2){
    Rsx = apply(gxgrid, 1, function(x) ifelse(x[1] == x[2],params[3],0))
    Rx = Rsx + Rfx + Rgx
  }else{
    Rx = Rfx + Rgx
  }
  return(Rx)
}


#------------------------------------------------
# Model Fitting
#------------------------------------------------
hplist = list(c(3,0.5),c(5,5),c(31,1))
mhp = c(1,1,0.2)

xind = expand.grid(1:n_train,1:n_train)
rx0 = eft_sg_sqr_exp_wrapper(x_train[xind[,1]],x_train[xind[,2]],c(0.5,1,0.2))

fsg_train = ifelse(fsg_train>3.2,3.2,fsg_train)
fsg_test = ifelse(fsg_test>3.2,3.2,fsg_test)

fit = gp_train(y_train - fsg_train, g_train, cov_function = eft_sg_sqr_exp_wrapper, 
               priornames = c("invgamma","invgamma","invgamma"), hp = hplist, 
               mh_proposal = mhp, nd = 5000, nadapt = 2000, nburn = 1000, 
               adaptevery = 500, nug = 1e-6)


plot(fit$post[,1],type = 'l')
plot(fit$post[,2],type = 'l')
plot(fit$post[,3],type = 'l')

m1 = fsg_train
m2 = fsg_test
pred = gp_predict(fit,y_train,g_train,g_test,m1,m2,cov_function = eft_sg_sqr_exp_wrapper)

plot(g_test,f0_test, type = 'l')
points(g_train, y_train, pch = 3)
lines(g_test,pred$pred_mean, col = "purple")
lines(g_test,pred$pred_mean + 2*pred$pred_sd , col = "purple", lty = 'dashed')
lines(g_test,pred$pred_mean - 2*pred$pred_sd , col = "purple", lty = 'dashed')
lines(gvec,fsg_grid, col = "red")


hplist = list(c(3,0.2),c(5,4),c(31,1))
mhp = c(1,1,0.2)

sh_train = which(g_train<gvec[76])
sh_test = which(g_test<gvec[76])
flg_train[sh_train] = flg_grid[76]
flg_test[sh_test] = flg_grid[76]

fit = gp_train(y_train - flg_train, g_train, cov_function = eft_lg_sqr_exp_wrapper, 
               priornames = c("invgamma","invgamma","invgamma"), hp = hplist, 
               mh_proposal = mhp, nd = 5000, nadapt = 2000, nburn = 1000, 
               adaptevery = 500, nug = 1e-6)


plot(fit$post[,1],type = 'l')
plot(fit$post[,2],type = 'l')
plot(fit$post[,3],type = 'l')

m1 = flg_train
m2 = flg_test
pred = gp_predict(fit,y_train,g_train,g_test,m1,m2,cov_function = eft_lg_sqr_exp_wrapper)

plot(g_test,f0_test, type = 'l')
points(g_train, y_train, pch = 3)
lines(g_test,pred$pred_mean, col = "purple")
lines(g_test,pred$pred_mean + 2*pred$pred_sd , col = "purple", lty = 'dashed')
lines(g_test,pred$pred_mean - 2*pred$pred_sd , col = "purple", lty = 'dashed')
lines(gvec,flg_grid, col = "blue")

#------------------------------------------------
# Code Example
#------------------------------------------------
sqr_exp_wrapper = function(x1,x2,params,sig2 = FALSE){
  Rx = apply(cbind(x1,x2), 1, function(x) sqr_exp_kernel(x[1],x[2],params[1],params[2]))
  return(Rx)
}

hplist = list(c(3,0.5),c(5,5))
mhp = c(1,1)

x = seq(1,5,length = 10)
xg = seq(1,5,length = 25)
y = x^2
y = y - mean(y)

fit = gp_train(y, x, cov_function = sqr_exp_wrapper, 
               priornames = c("invgamma","invgamma"), hp = hplist, 
               mh_proposal = mhp, nd = 2000, nadapt = 2000, nburn = 1000, 
               adaptevery = 500, nug = 1e-6)


plot(fit$post[,1])
plot(fit$post[,2])

m1 = rep(0,10)
m2 = rep(0,25)
pred = gp_predict(fit,y,x,xg,m1,m2,cov_function = sqr_exp_wrapper)

plot(xg,xg^2)
lines(xg,pred$pred_mean + mean(x^2))
