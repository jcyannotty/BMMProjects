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
n_train = 20
n_test = 50
s = 0.02
g_train = c(0.02,0.1,0.3,seq(0.4,3.0,length = n_train-7), seq(3.25,5,length=4))
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
# Assume x's are ordered
eft_sg_sqr_exp_wrapper = function(x1,x2,params,nrow,ncol,sig2 = TRUE){
  # EFT squared exponential kernel
  a = gvec[14]
  fxgrid = cbind(x1,x2)
  fxgrid[,1] = ifelse(fxgrid[,1] < a,fxgrid[,1],a)
  fxgrid[,2] = ifelse(fxgrid[,2] < a,fxgrid[,2],a)
  
  lscale = 1
  sscale = 0.05
  Rfx = apply(fxgrid, 1, function(x) sqr_exp_kernel(x[1],x[2],sscale,lscale))
  Rfx = matrix(Rfx, nrow = nrow, ncol = ncol, byrow = TRUE)
  
  # EFT add in scalar functions
  yref = 1
  qx = 0.5
  N = 10
  Rfx = ((yref^2)*(1 - qx^(N+1))/(1-qx^2))*Rfx
  
  # Intermediate discrepancy model
  gxgrid = cbind(x1,x2)
  Rgx = apply(gxgrid, 1, function(x) sqr_exp_kernel(x[1],x[2],params[1],params[2]))
  Rgx = sapply(1:nrow(gxgrid), function(i) if(gxgrid[i,1] < a || gxgrid[i,2] < a){return(0)}else{return(Rgx[i])})
  
  # Shape Rgx into matrix then get the non-zero block
  Rgx = matrix(Rgx, nrow = nrow, ncol = ncol, byrow = FALSE)
  Rgx_rows = which(rowSums(Rgx)>0)
  Rgx_cols = which(colSums(Rgx)>0)
  
  # Get the cov between g(a) and g(x)
  xg1 = unique(x1)
  xg1 = xg1[which(xg1>=a)]
  Rgxa1 = sapply(xg1, function(x) sqr_exp_kernel(x,a,params[1],params[2]))
  
  xg2 = unique(x2)
  xg2 = xg2[which(xg2>=a)]
  Rgxa2 = sapply(xg2, function(x) sqr_exp_kernel(x,a,params[1],params[2]))
  
  Raa = sqr_exp_kernel(a,a,params[1],params[2])
  
  # Get the conditional distribution covariance
  Rgx[Rgx_rows,Rgx_cols] = Rgx[Rgx_rows,Rgx_cols] - Rgxa1%*%t(Rgxa2)/Raa
  
  if(sig2){
    Rsx = apply(gxgrid, 1, function(x) ifelse(x[1] == x[2],params[3],0))
    Rx = Rsx + Rfx + Rgx
  }else{
    Rx = Rfx + Rgx
  }
  return(Rx)
}



eft_lg_sqr_exp_wrapper = function(x1,x2,params,nrow,ncol,sig2 = TRUE){
  # EFT squared exponential kernel
  b = gvec[84]
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
  Rfx = ((yref^2)*(1 - qx^(N+1))/(1-qx^2))*Rfx
  
  # Intermediate discrepancy model
  gxgrid = cbind(x1,x2)
  Rgx = apply(gxgrid, 1, function(x) sqr_exp_kernel(x[1],x[2],params[1],params[2]))
  Rgx = sapply(1:nrow(gxgrid), function(i) if(gxgrid[i,1] > b || gxgrid[i,2] > b){return(0)}else{return(Rgx[i])}   )
  
  # Shape Rgx into matrix then get the non-zero block
  Rgx = matrix(Rgx, nrow = nrow, ncol = ncol, byrow = FALSE)
  Rgx_rows = which(rowSums(Rgx)>0)
  Rgx_cols = which(colSums(Rgx)>0)
  
  # Get the cov between g(b) and g(x)
  xg1 = unique(x1)
  xg1 = xg1[which(xg1<=b)]
  Rgxb1 = sapply(xg1, function(x) sqr_exp_kernel(x,b,params[1],params[2]))
  
  xg2 = unique(x2)
  xg2 = xg2[which(xg2<=b)]
  Rgxb2 = sapply(xg2, function(x) sqr_exp_kernel(x,b,params[1],params[2]))
  
  Rbb = sqr_exp_kernel(b,b,params[1],params[2])
  
  Rgx[Rgx_rows,Rgx_cols] = Rgx[Rgx_rows,Rgx_cols] - Rgxb1%*%t(Rgxb2)/Rbb
  
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
hplist = list(c(3,15),c(5,5),c(31,1))
mhp = c(1,1,0.2)

xind = expand.grid(1:n_train,1:n_train)
rx0 = eft_sg_sqr_exp_wrapper(x_train[xind[,1]],x_train[xind[,2]],c(0.5,1,0.2),20,20)

sh_train = which(g_train>gvec[14])
sh_test = which(g_test>gvec[14])
fsg_train[sh_train] = fsg_grid[14]
fsg_test[sh_test] = fsg_grid[14]

ysc_train = y_train - fsg_train
#odds = 2*(1:ceiling(n_train/2))-1
set.seed(12)
odds = sample(1:5,4,replace = TRUE) + seq(0,15, by = 5)

fit_sg = gp_train(ysc_train[odds], g_train[odds], cov_function = eft_sg_sqr_exp_wrapper, 
               priornames = c("invgamma","invgamma","invgamma"), hp = hplist, 
               mh_proposal = mhp, nd = 5000, nadapt = 2000, nburn = 1000, 
               adaptevery = 500, nug = 1e-6)


plot(fit_sg$post[,1],type = 'l')
plot(fit_sg$post[,2],type = 'l')
plot(fit_sg$post[,3],type = 'l')

m1 = fsg_train
m2 = fsg_test

pred_sg = gp_predict(fit_sg,y_train[odds],g_train[odds],g_test,m1[odds],m2,cov_function = eft_sg_sqr_exp_wrapper)

plot(g_test,f0_test, type = 'l')
points(g_train, y_train, pch = 3)
lines(g_test,pred_sg$pred_mean, col = "purple")
lines(g_test,pred_sg$pred_mean + 2*pred_sg$pred_sd , col = "purple", lty = 'dashed')
lines(g_test,pred_sg$pred_mean - 2*pred_sg$pred_sd , col = "purple", lty = 'dashed')
lines(gvec,fsg_grid, col = "red")
lines(gvec,c(fsg_grid[1:14],rep(fsg_grid[14],186)) , col = "red", lty = "dashed")


hplist = list(c(3,15),c(5,4),c(31,1))
mhp = c(1,1,0.2)

lh_train = which(g_train<gvec[84])
lh_test = which(g_test<gvec[84])
flg_train[lh_train] = flg_grid[84]
flg_test[lh_test] = flg_grid[84]
#odds = c(2,8,16,19)

ylc_train = y_train - flg_train
fit_lg = gp_train(ylc_train[odds], g_train[odds], cov_function = eft_lg_sqr_exp_wrapper, 
               priornames = c("invgamma","invgamma","invgamma"), hp = hplist, 
               mh_proposal = mhp, nd = 5000, nadapt = 2000, nburn = 1000, 
               adaptevery = 500, nug = 1e-6)


plot(fit_lg$post[,1],type = 'l')
plot(fit_lg$post[,2],type = 'l')
plot(fit_lg$post[,3],type = 'l')

m1 = flg_train
m2 = flg_test
pred_lg = gp_predict(fit_lg,y_train[odds],g_train[odds],g_test,m1[odds],m2,cov_function = eft_lg_sqr_exp_wrapper)

plot(g_test,f0_test, type = 'l')
points(g_train, y_train, pch = 3)
lines(g_test,pred_lg$pred_mean, col = "purple")
lines(g_test,pred_lg$pred_mean + 2*pred_lg$pred_sd , col = "purple", lty = 'dashed')
lines(g_test,pred_lg$pred_mean - 2*pred_lg$pred_sd , col = "purple", lty = 'dashed')
lines(gvec,flg_grid, col = "blue")
lines(g_test,flg_test, col = "blue", lty = "dashed")
points(g_train[odds], flg_train[odds])


#------------------------------------------------
# Build Model Set
#------------------------------------------------
pred_sg_train = gp_predict(fit_sg,y_train[odds],
                           g_train[odds],g_train[-odds],fsg_train[odds],fsg_train[-odds],
                           cov_function = eft_sg_sqr_exp_wrapper)

pred_lg_train = gp_predict(fit_lg,y_train[odds],
                           g_train[odds],g_train[-odds],flg_train[odds],flg_train[-odds],
                           cov_function = eft_lg_sqr_exp_wrapper)

plot(g_test,f0_test, type = 'l')
points(g_train, y_train, pch = 3)
lines(g_test,pred_sg$pred_mean, col = "purple")
lines(g_test,pred_sg$pred_mean + 2*pred_sg$pred_sd , col = "purple", lty = 'dashed')
lines(g_test,pred_sg$pred_mean - 2*pred_sg$pred_sd , col = "purple", lty = 'dashed')
lines(gvec,fsg_grid, col = "red")
lines(gvec,c(fsg_grid[1:14],rep(fsg_grid[14],186)) , col = "red", lty = "dashed")
points(g_train[-odds], pred_sg_train$pred_mean, col = "red", pch = 12)

plot(g_test,f0_test, type = 'l')
points(g_train, y_train, pch = 3)
lines(g_test,pred_lg$pred_mean, col = "purple")
lines(g_test,pred_lg$pred_mean + 2*pred_lg$pred_sd , col = "purple", lty = 'dashed')
lines(g_test,pred_lg$pred_mean - 2*pred_lg$pred_sd , col = "purple", lty = 'dashed')
lines(gvec,flg_grid, col = "blue")
lines(g_test,flg_test, col = "blue", lty = "dashed")
points(g_train[-odds], pred_lg_train$pred_mean, col = "blue", pch = 12)

ms = list(
  fsg_train = pred_sg_train$pred_mean,
  fsg_test = pred_sg$pred_mean,
  fsg_test_lb = pred_sg$pred_mean - 2*pred_sg$pred_sd,
  fsg_test_ub = pred_sg$pred_mean + 2*pred_sg$pred_sd,
  fsg_grid = fsg_grid,
  a = gvec[14],
  aind = 14,
  
  flg_train = pred_lg_train$pred_mean,
  flg_test = pred_lg$pred_mean,
  flg_test_lb = pred_lg$pred_mean - 2*pred_lg$pred_sd,
  flg_test_ub = pred_lg$pred_mean + 2*pred_lg$pred_sd,
  flg_grid = flg_grid,
  b = gvec[84],
  bind = 84,
  
  g_interp = g_train[odds],
  g_train = g_train[-odds],
  g_test = g_test,
  gvec = gvec,

  y_train = y_train[-odds],
  y_interp = y_train[odds],
  f0_test = f0_test
)

filedir = "/home/johnyannotty/Documents/Dissertation/results/EFT/spheat_model_fit/"
saveRDS(ms, paste0(filedir,"ms_n20_01_23_24.rds"))
saveRDS(fit_sg, paste0(filedir,"fsg_fit_n4_01_23_24.rds"))
saveRDS(fit_lg, paste0(filedir,"flg_fit_n4_01_23_24.rds"))

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
