#-------------------------------------------------
# BART Variogram
#-------------------------------------------------
setwd("/home/johnyannotty/Documents/openbt/src")

# Load the R wrapper functions to the OpenBT library.
source("/home/johnyannotty/Documents/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/kernels.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Plotting/computer_expt_plots.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/gp_utils.R")

library(plotly)
library(viridis)
library(latex2exp)
library(geoR)

#library(sp)
#data(s100)
#plot(s100, lowess = TRUE)

#------------------------------------------------
# Get test grid
#------------------------------------------------
set.seed(101)
n_train = 30
n_test = 300
s = 0.1
xmin = -1; xmax = 1
x_train = seq(xmin+0.1,xmax-0.1, length = n_train) + runif(n_train, -0.1, 0.1)
x_test = seq(xmin,xmax,length = n_test)

sc = 1.5; lam = 0.5; pwr = 1

xtr_ind = expand.grid(1:n_train,1:n_train)
R11 = apply(xtr_ind,1,function(x) power_exp_kernel(x_train[x[1]],x_train[x[2]],sc,lam,pwr))
R11 = matrix(R11, nrow = n_train, ncol = n_train, byrow = FALSE)
R11 = R11 + diag(s^2,n_train)
dim(R11)

xts_ind = expand.grid(1:n_test,1:n_test)
R22 = apply(xts_ind,1,function(x) power_exp_kernel(x_test[x[1]],x_test[x[2]],sc,lam,pwr))
R22 = matrix(R22, nrow = n_test, ncol = n_test, byrow = FALSE)

xtrs_ind = expand.grid(1:n_train,1:n_test)
R12 = apply(xtrs_ind,1,function(x) power_exp_kernel(x_train[x[1]],x_test[x[2]],sc,lam,pwr))
R12 = matrix(R12, nrow = n_train, ncol = n_test, byrow = FALSE)
dim(R12)

m1 = rep(0,n_train)
m2 = rep(0,n_test)

# Generate y
set.seed(1334)
y_train = sample_gp(m1,R11)
plot(x_train, y_train)

# True Predictive Distribution 
pdgp = predict_dist_gp(y_train,m1,m2,R11,R22,R12)
f0_test = pdgp$mp
plot(x_test, f0_test, type = "l")


#------------------------------------------------
# Variogram Time
#------------------------------------------------
ugrid = seq(0.05,1.5, length = 20)
vgyhat = variog(coords = cbind(x_train,rep(0,n_train)), data = y_train ,uvec = ugrid)

plot(vgyhat$u,vgyhat$v)

xbnds = matrix(c(xmin,xmax), nrow = 1, byrow = FALSE)

# Set outfile directory
h_grid = ugrid
vg = variogram.openbtmixing(xbnds,h_grid,10000,1,
                            k=1.25,
                            0.95,
                            power = 1.25,
                            a1 = 2,
                            a2 = 30,
                            4, 
                            type = 'b',
                            ymin = min(y_train),
                            ymax = max(y_train),
                            sigma2 = 0.005,
                            ncut = 200)

# Semi-variogram
plot(h_grid,vg$vmean/2, type = "l", ylim = c(0,3.0))
points(vgyhat$u,vgyhat$v)
abline(h = var(y_train), col = "grey")
#points(h_grid,sapply(h_grid,function(h) (1-power_exp_kernel(0,h,1,0.3,1))), pch = 3, col = 'red')


#-------------------------------------------------
# Batch Fit
#-------------------------------------------------
nu = 15
q0 = 4

# Train a bart model
fit=train.openbtmixing(x_train,y_train,as.matrix(rep(1,n_train)),pbd=c(1.0,0),ntree = 30,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="physics_model",
                       ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
                       power = 1.25, base = 0.95, minnumbot = 1, overallsd = sqrt(0.005), k = 1.15, overallnu = nu,
                       summarystats = FALSE, rpath = FALSE, q = q0, rshp1 = 2, rshp2 = 100,
                       stepwpert = 0.1, probchv = 0.1, batchsize = 10000)


#Get mixed mean function
fitp=predict.openbtmixing(fit,x.test = x_test, f.test = as.matrix(rep(1,n_test)),tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "mean_and_sigma")

# Sigma
hist(fitp$sdraws[,1])
plot(fitp$sdraws[,1], type = 'l')

# compare vs true curve
sqrt(mean((fitp$mmean - f0_test)^2))

plot_mean1d(x_test, pmean = fitp$mmean, plb = fitp$m.lower, pub = fitp$m.upper,
            amean = c(f0_test), colors = c("black","purple2"),
            apts_x = x_train, apts_y = y_train)

#plot(sapply(h_grid,function(h) 2*(1-power_exp_kernel(0,h,1,0.3,1))))

#------------------------------------------------
# Prediction with gp
#------------------------------------------------
#geodata = list(data = y_train, coords = cbind(x_train,runif(n_train,-0.001,0.001)))
geodata = list(data = y_train, coords = cbind(x_train,rep(0,n_train)))
geor_fit = likfit(geodata, cov.model = "gaussian", ini.cov.pars = c(0.5,0.1),fix.nugget = FALSE,
                  nugget = 0.01) 

xtr_ind = expand.grid(1:n_train,1:n_train)
xts_ind = expand.grid(1:n_test,1:n_test)
xtrs_ind = expand.grid(1:n_train,1:n_test)

if(geor_fit$cov.model == "gaussian"){
  gps2 = geor_fit$tausq
  gpsc = sqrt(geor_fit$cov.pars[1])
  gpls = geor_fit$cov.pars[2]^2/2

  Rp11 = apply(xtr_ind,1,function(x) sqr_exp_kernel(x_train[x[1]],x_train[x[2]],gpsc,gpls))
  Rp22 = apply(xts_ind,1,function(x) sqr_exp_kernel(x_test[x[1]],x_test[x[2]],gpsc,gpls))
  Rp12 = apply(xtrs_ind,1,function(x) sqr_exp_kernel(x_train[x[1]],x_test[x[2]],gpsc,gpls))
}else{
  gps2 = geor_fit$tausq
  gpsc = sqrt(geor_fit$cov.pars[1])
  gpls = geor_fit$cov.pars[2]
  
  Rp11 = apply(xtr_ind,1,function(x) power_exp_kernel(x_train[x[1]],x_train[x[2]],gpsc,gpls,1))
  Rp22 = apply(xts_ind,1,function(x) power_exp_kernel(x_test[x[1]],x_test[x[2]],gpsc,gpls,1))
  Rp12 = apply(xtrs_ind,1,function(x) power_exp_kernel(x_train[x[1]],x_test[x[2]],gpsc,gpls,1))
}

Rp11 = matrix(Rp11, nrow = n_train, ncol = n_train, byrow = FALSE)
Rp11 = Rp11 + diag(gps2,n_train)

Rp22 = matrix(Rp22, nrow = n_test, ncol = n_test, byrow = FALSE)

Rp12 = matrix(Rp12, nrow = n_train, ncol = n_test, byrow = FALSE)

mp1 = geor_fit$trend.matrix*geor_fit$beta
mp2 = rep(geor_fit$beta,n_test)
pdgp = predict_dist_gp(y_train,mp1,mp2,Rp11,Rp22,Rp12)
fhat = pdgp$mp
fshat = sqrt(diag(pdgp$Rp))
fhat_lb = fhat - 1.96*fshat
fhat_ub = fhat + 1.96*fshat

sqrt(mean((fhat - f0_test)^2))

plot_mean1d(x_test, pmean = fhat, plb = fhat_lb, pub = fhat_ub,
            amean = c(f0_test), colors = c("black","purple2"),
            apts_x = x_train, apts_y = y_train)


#------------------------------------------------
# Compare methods
#------------------------------------------------
plot_mean1d(x_test, pmean = cbind(fitp$mmean,fhat), 
            plb = cbind(fitp$m.lower,fhat_lb), 
            pub = cbind(fitp$m.upper,fhat_ub),
            amean = c(f0_test), colors = c("black","red","blue"))

plot(h_grid,vg$vmean/2, type = "l", ylim = c(0,3.5), col = "blue")
points(vgyhat$u,vgyhat$v)
lines(h_grid,sapply(h_grid,function(h) (gpsc^2-sqr_exp_kernel(0,h,gpsc,gpls))+gps2),
      col = 'red')
abline(h = var(y_train), col = "grey")
lines(geor_fit)

#------------------------------------------------
# Save Results
#------------------------------------------------
fit_data = list(
  pred_mean = fitp$mmean,
  pred_ub = fitp$m.upper,
  pred_lb = fitp$m.lower,
  m = fit$m,
  k = fit$k,
  shp1 = fit$rshp1,
  shp2 = fit$rshp2,
  minnodesz = fit$minnumbot,
  q = q0,
  base = fit$base,
  power = fit$power,
  maxd = fit$maxd,
  nu = fit$overallnu,
  lam = fit$overalllambda,
  sigma = fitp$smean[1]
)

# Model set
ms = list(
  x_train = x_train,
  y_train = y_train,
  x_test = x_test,
  f0_test = f0_test,
  cov = "power_exp",
  sc = sc,
  ls = lam,
  pwr = pwr,
  sigma = s
)

# GP Results
gp_data = list(
  fhat = fhat,
  fshat = fshat,
  fhat_lb = fhat_lb,
  fhat_ub = fhat_ub,
  sc = gpsc,
  ls = gpls,
  beta = geor_fit$beta,
  model = geor_fit$cov.model,
  sigma = sqrt(gps2)
)

# Save results
filedir = "/home/johnyannotty/Documents/RandomPathBART/BART_Results/Simulations/ms1/"
saveRDS(fit_data,paste0(filedir,"dpath_res1.rds"))
saveRDS(gp_data,paste0(filedir,"gp_res2.rds"))
saveRDS(fitp$sdraws[,1], paste0(filedir,"rpath_sdraws1.rds"))
saveRDS(ms,paste0(filedir,"ms.rds"))


