#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("/home/johnyannotty/Documents/openbt/src")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")

source("/home/johnyannotty/Documents/SciData/NuclearPhysics/R/heat_2d_ising.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/LinearRegression/bayesian_regression.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/kernels.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/gp_utils.R")

#----------------------------------------------------------
#Get training data for interpolations
gvec = seq(0.01,5,length=200)
gset1 = c(0.02,0.1,0.15,0.25)
gset2 = c(3,3.5,4,4.5,5)
gset = c(gset1,gset2)

fsg = heat_c2x_sg(gset1)
flg = heat_c2x_lg(gset2)
fg = c(fsg,flg)

L = 2
xvec = 0.5*log(gvec+1)
dz2 = sapply(xvec,function(x) d2_logz(x,L,0.1))


#----------------------------------------------------------
# Train interpolators
#----------------------------------------------------------
# Polynomial - order 8
# Set train and test basis
pbasis = cbind(gset^0,gset,gset^2,gset^3,gset^4,gset^5,gset^6,gset^7)
colnames(pbasis) = paste0("g",0:(ncol(pbasis)-1))

ptest = cbind(gvec^0,gvec,gvec^2,gvec^3,gvec^4,gvec^5,gvec^6,gvec^7)
colnames(ptest) = paste0("g",0:(ncol(ptest)-1))

# Set the priors and fit
mu_vec = c(mean(fg),rep(0,ncol(pbasis)-1))
V = diag(5,ncol(pbasis))
fpoly_post = bayes_reg(pbasis,fg,mu_vec,V,30,0.05,5000,1000)
fpoly_mean = predict_bayes(fpoly_post,ptest)

# Get predictions
pm = apply(fpoly_mean$post_dist,2,mean)
pm_lb = apply(fpoly_mean$post_dist,2,quantile, 0.025)
pm_ub = apply(fpoly_mean$post_dist,2,quantile, 0.975)

# GP with squared exp
gp_modelset = list() 
lam = c(0.4,1)
sc = c(0.2,0.4)
#gvec_gp = c(0,0.005,seq(0.021,5.2,length=70))
gvec_gp = seq(0,5.2,length = 100)
gset_gp = gset[-c(2,3)]
fg_gp = fg[-c(2,3)]

for(j in 1:length(sc)){
  g1grid = expand.grid(gset_gp,gset_gp)
  R11 = apply(g1grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc[j],lam[j]))
  R11 = matrix(R11, nrow = length(gset_gp), ncol = length(gset_gp))
  dim(R11)
  #Rx = Rx + diag(sig2,n)
  
  g2grid = expand.grid(gvec_gp,gvec_gp)
  R22 = apply(g2grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc[j],lam[j]))
  R22 = matrix(R22, nrow = length(gvec_gp), ncol = length(gvec_gp))
  dim(R22)
  
  g12grid = expand.grid(gset_gp,gvec_gp)
  R12 = apply(g12grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc[j],lam[j]))
  R12 = matrix(R12, nrow = length(gset_gp), ncol = length(gvec_gp), byrow = FALSE)
  dim(R12)
  
  m1 = rep(mean(fg_gp),length(fg_gp))
  m2 = rep(mean(fg_gp),length(gvec_gp))
  pdgp = predict_dist_gp(fg_gp,m1,m2,R11,R22,R12)
  pdgp_lb = pdgp$mp - 2*diag(pdgp$Rp)^0.5
  pdgp_ub = pdgp$mp + 2*diag(pdgp$Rp)^0.5

  gp_modelset[[j]] = data.frame(m = pdgp$mp, lb = pdgp_lb, ub = pdgp_ub)
}

# Plots
plot(gvec,dz2/L^2, type = 'l')
lines(gvec,heat_c2x_sg(gvec),col = "red",lty = "dashed")
lines(gvec,heat_c2x_lg(gvec),col = "blue",lty = "dashed")
lines(gvec,pm,col = "green",lty = "dashed")
lines(gvec,pm_lb,col = "green",lty = "dotted")
lines(gvec,pm_ub,col = "green",lty = "dotted")
lines(gvec_gp,gp_modelset[[1]]$m,col = "purple",lty = "dashed")
lines(gvec_gp,gp_modelset[[1]]$ub,col = "purple",lty = "dotted")
lines(gvec_gp,gp_modelset[[1]]$lb,col = "purple",lty = "dotted")
lines(gvec_gp,gp_modelset[[2]]$m,col = "orange",lty = "dashed")
lines(gvec_gp,gp_modelset[[2]]$ub,col = "orange",lty = "dotted")
lines(gvec_gp,gp_modelset[[2]]$lb,col = "orange",lty = "dotted")


#----------------------------------------------------------
# Generate data for mixing
#----------------------------------------------------------
L = 2
n_train = 15
n_test = 100
s = 0.02
x_train = c(0.01,0.10,0.2,seq(0.5,5,length = n_train-3))
x_test = seq(0.01,5, length = n_test)

y_train = sapply(0.5*log(x_train+1),function(x) d2_logz(x,L,0.1))/L^2 + rnorm(n_train,0,s)
f0_test = sapply(0.5*log(x_test+1),function(x) d2_logz(x,L,0.1)/L^2)

# Build f sets
f_train = matrix(0, nrow = n_train, ncol = 4)
f_test = matrix(0, nrow = n_test, ncol = 4)

f_train[,1] = heat_c2x_sg(x_train)
f_train[,2] = heat_c2x_lg(x_train)
f_test[,1] = heat_c2x_sg(x_test)
f_test[,2] = heat_c2x_lg(x_test)

pb_train = cbind(x_train^0,x_train,x_train^2,x_train^3,x_train^4,x_train^5,x_train^6,x_train^7)
pb_test = cbind(x_test^0,x_test,x_test^2,x_test^3,x_test^4,x_test^5,x_test^6,x_test^7)

xtr_grid = expand.grid(x_train,x_train)
Rtr = apply(xtr_grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc,lam))
Rtr = matrix(Rtr, nrow = n_train, ncol = n_train)

x1trgrid = expand.grid(gset_gp,x_train)
R1tr = apply(x1trgrid,1,function(x) sqr_exp_kernel(x[1],x[2],sc,lam))
R1tr = matrix(R1tr, nrow = length(gset_gp), ncol = n_train, byrow = FALSE)

mtr = rep(mean(fg_gp),n_train)

xts_grid = expand.grid(x_test,x_test)
Rts = apply(xts_grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc,lam))
Rts = matrix(Rts, nrow = n_test, ncol = n_test)

x1tsgrid = expand.grid(gset_gp,x_test)
R1ts = apply(x1tsgrid,1,function(x) sqr_exp_kernel(x[1],x[2],sc,lam))
R1ts = matrix(R1ts, nrow = length(gset_gp), ncol = n_test, byrow = FALSE)

mts = rep(mean(fg_gp),n_test)

# Get the polynomial and GP predictions
f_train[,3] = apply(predict_bayes(fpoly_post,pb_train)$post_dist,2,mean)
f_train[,4] = predict_dist_gp(fg_gp,m1,mtr,R11,Rtr,R1tr)$mp

f_test[,3] = apply(predict_bayes(fpoly_post,pb_test)$post_dist,2,mean)
f_test[,4] = predict_dist_gp(fg_gp,m1,mts,R11,Rts,R1ts)$mp

#----------------------------------------------------------
# Mixing with Rpath
nu = 30
rho = 1
sig2_hat = max(apply(apply(f_train, 2, function(x) (x-y_train)^2),2,min))
lam = rho*sig2_hat*(nu+2)/nu
q0 = 4
#f_train[,2] = ifelse(f_train[,2]>10^5,10^5,f_train[,2])
#f_test[,2] = ifelse(f_test[,2]>10^5,10^5,f_test[,2])

#f_train[,1] = ifelse(f_train[,1]<3.9,3.9,f_train[,1])
#f_test[,1] = ifelse(f_test[,1]<3.9,3.9,f_test[,1])
#f_train[,2] = ifelse(f_train[,2]>2,2,f_train[,2])
#f_test[,2] = ifelse(f_test[,2]>2,2,f_test[,2])

fit=train.openbtmixing(g_train,y_train,f_train[,c(1,2)],pbd=c(1.0,0),ntree = 10,ntreeh=1,numcut=100,tc=2,model="mixbart",modelname="physics_model",
                       ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
                       power = 2.0, base = 0.95, minnumbot = 2, overallsd = sqrt(sig2_hat), k = 1.0, overallnu = nu,
                       summarystats = FALSE, rpath = TRUE, q = q0, rshp1 = 4, rshp2 = 25,
                       stepwpert = 0.1, probchv = 0.1, batchsize = 1000)


#Get mixed mean function
fitp=predict.openbtmixing(fit,x.test = g_test, f.test = f_test[,c(1,2)],
                          tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "mean_and_sigma")

plot(g_test, f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(-2,7))
points(g_train, y_train, pch = 3)
lines(g_test, f_test[,1], col = 'red', lty = 2)
lines(g_test, f_test[,2], col = 'blue', lty = 2)
lines(g_test, f_test[,3], col = 'green', lty = 2)
lines(g_test, f_test[,4], col = 'gold', lty = 2)
lines(g_test, fitp$mmean, col = 'purple', lwd = 2)
lines(g_test, fitp$m.lower, col = 'purple', lwd = 2, cex = 0.5, lty = 'dashed')
lines(g_test, fitp$m.upper, col = 'purple', lwd = 2, cex = 0.5, lty = 'dashed')

# Get wsum
wsum = 0*fitp$wdraws[[1]]
for(i in 1:length(fitp$wdraws)){
  wsum = wsum + fitp$wdraws[[i]]
}
wsum_mean = apply(wsum,2,mean)
wsum_lb = apply(wsum,2,quantile, 0.025)
wsum_ub = apply(wsum,2,quantile, 0.975)


#Plot model weights
plot(g_test, fitp$wmean[,1], pch = 16, col = 'red', type = 'l', ylim = c(-5,5.0), lwd = 2, 
     panel.first = {grid(col = 'lightgrey')})
lines(g_test, fitp$wmean[,2], col = 'blue', lwd = 2)
lines(g_test, fitp$wmean[,3], col = 'green', lwd = 2)
#lines(g_test, fitp$wmean[,4], col = 'gold', lwd = 2)
lines(g_test, fitp$w.upper[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(g_test, fitp$w.lower[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(g_test, fitp$w.upper[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(g_test, fitp$w.lower[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(g_test, fitp$w.upper[,3], col = 'green', lty = 'dashed', lwd = 1)
lines(g_test, fitp$w.lower[,3], col = 'green', lty = 'dashed', lwd = 1)
#lines(g_test, fitp$w.upper[,4], col = 'gold', lty = 'dashed', lwd = 1)
#lines(g_test, fitp$w.lower[,4], col = 'gold', lty = 'dashed', lwd = 1)

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


# Save model set
ms = list(f_train = f_train,
          f_test = f_test,
          x_train = x_train,
          x_test = x_test,
          y_train = y_train,
          fg_test = f0_test,
          gp_modelset = gp_modelset,
          gp_params = list(kernel = klist, sc = sc, ls = lam, shp = shp)
    )

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
#  proj_mean = fitp$pmmean,
#  proj_ub = fitp$pm.upper,
#  proj_lb = fitp$pm.lower,
#  pwts_mean = fitp$pwmean,
#  pwts_ub = fitp$pw.upper,
#  pwts_lb = fitp$pw.lower,
#  delta_mean = fitp$dmean,
#  delta_lb = fitp$d.lower,
#  delta_ub = fitp$d.upper,
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
saveRDS(fit_data, paste0(filedir,"spheat_gp_set12_rpath_01_18_24.rds"))
saveRDS(ms, paste0(filedir,"spheat_ms2_01_18_24.rds"))



#----------------------------------------------------------
# General - Multiple GPs over multiple training sets 
#----------------------------------------------------------
L = 8
n_train = 20
n_test = 100
s = 0.02
#g_train = c(0.01,0.10,0.2,seq(0.5,5,length = n_train-3))
g_train = seq(0.02,5,length = n_train)
g_test = seq(0.01,5, length = n_test)
x_train = 0.5*log(g_train+1)
x_test = 0.5*log(g_test+1)

fsg_grid = heat_sg_exp(g_test,10,L)
flg_grid = heat_lg_exp(g_test,10,L)

L = 8
xvec = 0.5*log(gvec+1)
dz2_train = sapply(x_train,function(x) d2_logz(x,L,0.1))

y_train = dz2_train/L^2 + rnorm(n_train,0,s)
f0_test = sapply(x_test,function(x) d2_logz(x,L,0.1)/L^2)


# Build f sets
f_train = matrix(0, nrow = n_train, ncol = 4)
f_test = matrix(0, nrow = n_test, ncol = 4)

# GP with squared exp
gp_modelset = list() 
lam = c(0.5,0.5,1,0.5)
sc = c(0.5,0.5,0.5,1)
shp = c(1.5,1.5,1.5,1.5)
klist = c("m","m","m","m")

gslist = list(c(0.02, 0.25),c(0.1, 0.2, 0.4),c(0.05, 0.2),c(0.05,0.15,0.45))
gllist = list(c(3,3.5,4,4.5,5),c(2.5,3,3.5,4,4.5),c(1.6,2.5,3.5,4.5,5),c(1.5,2.5,3.5,4.5,5))
fslist = list()
fllist = list()
for(j in 1:length(sc)){
  fsg0 = heat_sg_exp(gslist[[j]], 10, L)
  flg0 = heat_lg_exp(gllist[[j]], 10, L)
  fslist[[j]] = fsg0
  fllist[[j]] = flg0
}

for(j in 1:length(sc)){
  # Evaluate Taylor series
  gset_gp = c(gslist[[j]],gllist[[j]])
  fg_gp = c(fslist[[j]], fllist[[j]])
  
  # GP Train
  g1grid = expand.grid(gset_gp,gset_gp)
  g2grid = expand.grid(g_train,g_train)
  g12grid = expand.grid(gset_gp,g_train)
  if(klist[j] == "s"){
    R11 = apply(g1grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc[j],lam[j])) 
    R22 = apply(g2grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc[j],lam[j]))
    R12 = apply(g12grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc[j],lam[j]))
  }else{
    R11 = apply(g1grid,1,function(x) matern_kernel(x[1],x[2],sc[j],lam[j],shp[j]))
    R22 = apply(g2grid,1,function(x) matern_kernel(x[1],x[2],sc[j],lam[j],shp[j]))
    R12 = apply(g12grid,1,function(x) matern_kernel(x[1],x[2],sc[j],lam[j],shp[j]))
  }
  
  # Reshape
  R11 = matrix(R11, nrow = length(gset_gp), ncol = length(gset_gp))
  R22 = matrix(R22, nrow = length(g_train), ncol = length(g_train))
  R12 = matrix(R12, nrow = length(gset_gp), ncol = length(g_train), byrow = FALSE)
  
  m1 = rep(mean(fg_gp),length(fg_gp))
  mtr = rep(mean(fg_gp),n_train)
  
  # Predict
  f_train[,j] = predict_dist_gp(fg_gp,m1,mtr,R11,R22,R12)$mp
  
  # Obs Test Set
  g2grid = expand.grid(g_test,g_test)
  g12grid = expand.grid(gset_gp,g_test)
  if(klist[j] == "s"){
    R22 = apply(g2grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc[j],lam[j]))
    R12 = apply(g12grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc[j],lam[j]))
  }else{
    R22 = apply(g2grid,1,function(x) matern_kernel(x[1],x[2],sc[j],lam[j],shp[j]))
    R12 = apply(g12grid,1,function(x) matern_kernel(x[1],x[2],sc[j],lam[j],shp[j]))
  }
  
  # Reshape
  R11 = matrix(R11, nrow = length(gset_gp), ncol = length(gset_gp))
  R22 = matrix(R22, nrow = length(g_test), ncol = length(g_test))
  R12 = matrix(R12, nrow = length(gset_gp), ncol = length(g_test), byrow = FALSE)
  
  mts = rep(mean(fg_gp),n_test)
  
  # Predictions
  pgp = predict_dist_gp(fg_gp,m1,mts,R11,R22,R12)
  pgp_lb = pgp$mp - 2*diag(pgp$Rp)^0.5
  pgp_ub = pgp$mp + 2*diag(pgp$Rp)^0.5
  f_test[,j] = pgp$mp
  gp_modelset[[j]] = data.frame(m = pgp$mp, lb = pgp_lb, ub = pgp_ub)
}


par(mfrow = c(2,2))
for(j in 1:length(klist)){
  plot(g_test,f0_test, type = 'l', col="black", main = paste(klist[j],j), ylim=c(-5,8))
  points(gslist[[j]],fslist[[j]], col = "red", pch = 3)
  points(gllist[[j]],fllist[[j]], col = "blue", pch = 3)
  lines(g_test,fsg_grid,col = "red",lty = "dashed")
  lines(g_test,flg_grid,col = "blue",lty = "dashed")
  lines(g_test,gp_modelset[[j]]$m,col = "orange",lty = "dashed")
  lines(g_test,gp_modelset[[j]]$lb,col = "orange",lty = "dotted")
  lines(g_test,gp_modelset[[j]]$ub,col = "orange",lty = "dotted")
}

