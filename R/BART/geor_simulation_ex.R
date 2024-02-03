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

library(sp)
data(s100)
#plot(s100, lowess = TRUE)

#------------------------------------------------
# Get test grid
#------------------------------------------------
set.seed(100)
n_test = 300
x_train = s100$coords
n_train = nrow(x_train)
y_train = s100$data
x_test = cbind(runif(n_test),runif(n_test))

xtr_ind = expand.grid(1:n_train,1:n_train)
R11 = apply(xtr_ind,1,function(x) power_exp_kernel(x_train[x[1],],x_train[x[2],],1,0.3,1))
R11 = matrix(R11, nrow = n_train, ncol = n_train, byrow = FALSE)
dim(R11)

xts_ind = expand.grid(1:n_test,1:n_test)
R22 = apply(xts_ind,1,function(x) power_exp_kernel(x_test[x[1],],x_test[x[2],],1,0.3,1))
R22 = matrix(R22, nrow = n_test, ncol = n_test, byrow = FALSE)

xtrs_ind = expand.grid(1:n_train,1:n_test)
R12 = apply(xtrs_ind,1,function(x) power_exp_kernel(x_train[x[1],],x_test[x[2],],1,0.3,1))
R12 = matrix(R12, nrow = n_train, ncol = n_test, byrow = FALSE)
dim(R12)

m1 = rep(0,n_train)
m2 = rep(0,n_test)
pdgp = predict_dist_gp(y_train,m1,m2,R11,R22,R12)

f0_test = pdgp$mp


#------------------------------------------------
# Variogram Time
#------------------------------------------------
ugrid = seq(0.05,1, length = 20)
vgyhat = variog(coords = s100$coords, data = s100$data ,uvec = ugrid)

plot(vgyhat$u,vgyhat$v)

xbnds = matrix(c(0,0,1,1), nrow = 2, byrow = FALSE)

# Set outfile directory
h_grid = seq(0.02,1,by = 0.02)
vg = variogram.openbtmixing(xbnds,h_grid,10000,1,
                            k=1.3,
                            0.95,
                            power = 1.5,
                            a1 = 2,
                            a2 = 20,
                            4, 
                            type = 'b',
                            ymin = min(s100$data),
                            ymax = max(s100$data),
                            sigma2 = 0.05)

# Semi-variogram
plot(h_grid,vg$vmean/2, type = "l", ylim = c(0,2.5))
points(vgyhat$u,vgyhat$v)
points(h_grid,sapply(h_grid,function(h) (1-power_exp_kernel(0,h,1,0.3,1))), pch = 3, col = 'red')
lines(geor_fit)

#-------------------------------------------------
# Batch Fit
#-------------------------------------------------
nu = 30
rho = 1
sig2_hat = 0.1  #0.05*var(y_train)
q0 = 4

# Train a bart model
fit=train.openbtmixing(s100$coords,s100$data-mean(s100$data),as.matrix(rep(1,100)),pbd=c(1.0,0),ntree = 30,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="physics_model",
                       ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
                       power = 1.5, base = 0.95, minnumbot = 2, overallsd = 0.05^2, k = 1.3, overallnu = nu,
                       summarystats = FALSE, rpath = TRUE, q = q0, rshp1 = 2, rshp2 = 20,
                       stepwpert = 0.1, probchv = 0.1, batchsize = 1000)


#Get mixed mean function
fitp=predict.openbtmixing(fit,x.test = x_test, f.test = as.matrix(rep(1,n_test)),tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "mean_and_sigma", proj_type = "softmax", temperature = 0.2)

# Sigma
hist(fitp$sdraws[,1])

# compare vs true curve
sqrt(mean((fitp$mmean - f0_test)^2))

#plot(sapply(h_grid,function(h) 2*(1-power_exp_kernel(0,h,1,0.3,1))))

#------------------------------------------------
# Prediction with gp
#------------------------------------------------
geor_fit = likfit(s100, cov.model = "gaussian", ini.cov.pars = c(0.5,0.3)) 

xtr_ind = expand.grid(1:n_train,1:n_train)
Rp11 = apply(xtr_ind,1,function(x) sqr_exp_kernel(x_train[x[1],],x_train[x[2],],
                                                  geor_fit$cov.pars[1],geor_fit$cov.pars[2]^2))
Rp11 = matrix(Rp11, nrow = n_train, ncol = n_train, byrow = FALSE)
dim(Rp11)

xts_ind = expand.grid(1:n_test,1:n_test)
Rp22 = apply(xts_ind,1,function(x) sqr_exp_kernel(x_test[x[1],],x_test[x[2],],
                                                  geor_fit$cov.pars[1],geor_fit$cov.pars[2]^2))
Rp22 = matrix(Rp22, nrow = n_test, ncol = n_test, byrow = FALSE)

xtrs_ind = expand.grid(1:n_train,1:n_test)
Rp12 = apply(xtrs_ind,1,function(x) sqr_exp_kernel(x_train[x[1],],x_test[x[2],],
                                                   geor_fit$cov.pars[1],geor_fit$cov.pars[2]^2))
Rp12 = matrix(Rp12, nrow = n_train, ncol = n_test, byrow = FALSE)
dim(Rp12)

mp1 = geor_fit$trend.matrix
mp2 = rep(1,n_test)
pdgp = predict_dist_gp(y_train,mp1,mp2,Rp11,Rp22,Rp12)
fhat = pdgp$mp
sqrt(mean((fhat - f0_test)^2))

sqrt(mean(geor_fit$model.components$residuals^2))


#------------------------------------------------
# Generate date
#------------------------------------------------
branin <- function(xx, a=1, b=5.1/(4*pi^2), c=5/pi, r=6, s=10, t=1/(8*pi)){
  x1 <- xx[1]
  x2 <- xx[2]
  term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
  term2 <- s*(1-t)*cos(x1)
  y <- term1 + term2 + s
  return(y)
}

gstat_sine = function(x){
  return(sin((x[1]*2*pi)*(x[2])*2*pi))
}

n_train = 50
n_test = 100
s = 0.05
xbnds = matrix(c(-5,-5,5,5), nrow = 2, byrow = FALSE)

set.seed(44)
x1_train = runif(n_train, xbnds[1,1],xbnds[1,2])
x2_train = runif(n_train, xbnds[2,1],xbnds[2,2])
x_train = cbind(x1 = x1_train, x2 = x2_train)

x1_test = seq(xbnds[1,1], xbnds[1,2], length = n_test)
x2_test = seq(xbnds[2,1], xbnds[2,2], length = n_test)
x_test = cbind(x1 = x1_test, x2 = x2_test)

set.seed(44)
y_train = apply(x_train,1,gstat_sine) + rnorm(n_train, 0, s)
f0_test = apply(x_test,1,gstat_sine)
hist(y_train)
