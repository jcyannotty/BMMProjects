#-------------------------------------------------
# Mixing polynomials - 3 Models in 1D
# Hierarchical Stacking
#-------------------------------------------------
source("/home/johnyannotty/Documents/openbt/R/polynomials.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Samplers/sampler_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Samplers/density_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/BMM/fwls.R")

library(plotly)
library(viridis)
library(latex2exp)
library(rstan)
library(loo)

# Set directory
filedir = '/home/johnyannotty/Documents/BMMProjects/R/1DFunctions/'
standir = '/home/johnyannotty/Documents/BMMProjects/R/1DFunctions/'

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
# Fit each individual model (Closed form available)
# Ml: Y ~ N(fl(x),sig^2), sig2 ~ Lam*nu/Chi_nu^2
#-------------------------------------------------
nu = 5; lam = 3
#s2 = rscinvchi2(10000,nu,lam)
#hist(sqrt(s2), main = "Prior draws of sigma")
nupost = lampost = 0
r2matrix = matrix(0,nrow = n_train, ncol = K)
y_loopd = matrix(0,nrow = n_train, ncol = K)
for(l in 1:K){
  # Data info
  r2matrix[,l] = (y_train - f_train[,l])^2
  ssr = sum(r2matrix[,l])
  
  # Posterior info (stored for informative purposes)
  nupost[l] = n_train + nu
  lampost[l] = (nu*lam + ssr)/nupost[l]
  
  # LOO PD
  nuloo = nupost[l] - 1
  lamloo = (nu*lam + ssr-r2matrix[,l])/nuloo
  y_loopd[,l] = dtscaled(y_train, df=nuloo, mean=f_train[,l], scale = sqrt(lamloo)) 
}


#-------------------------------------------------
# Fit the HS Model using stan
#-------------------------------------------------
lpd_point = log(y_loopd)
stan_data = list(X = matrix(x_train), N=n_train, d=1, d_discrete=0,
                     lpd_point = lpd_point, K=ncol(lpd_point),tau_mu = 1,
                     tau_sigma = 1, tau_discrete = 0.5 , tau_con = 1)
fiths = stan(paste0(standir,"polynomials_set3_hs_linear.stan"), data = stan_data)
w_fit = extract(fiths, pars=w)$w

head(extract(fiths)$mu)
head(extract(fiths)$mu_0)
head(extract(fiths)$sigma)
head(extract(fiths)$beta_con)
head(extract(fiths)$sigma_con)


#str(fiths)
names(fiths)
w_fit = extract(fiths)$w

w1 = apply(w_fit[,,1],2,mean)
w2 = apply(w_fit[,,2],2,mean)
w3 = apply(w_fit[,,3],2,mean)

w1_lb = apply(w_fit[,,1],2,quantile, 0.025)
w2_lb = apply(w_fit[,,2],2,quantile, 0.025)
w3_lb = apply(w_fit[,,3],2,quantile, 0.025)

w1_ub = apply(w_fit[,,1],2,quantile, 0.975)
w2_ub = apply(w_fit[,,2],2,quantile, 0.975)
w3_ub = apply(w_fit[,,3],2,quantile, 0.975)

plot(x_train,w1*f_train[,1] + w2*f_train[,2] + w3*f_train[,3])
lines(x_test, f0_test)

plot(x_train,w1, type = "l", col = "red", ylim = c(0,1))
lines(x_train,w2,col='blue')
lines(x_train,w3,col='green3')
lines(x_train,w1_lb,col='red', lty = "dashed")
lines(x_train,w1_ub,col='red', lty = "dashed")
lines(x_train,w2_lb,col='blue', lty = "dashed")
lines(x_train,w2_ub,col='blue', lty = "dashed")
lines(x_train,w3_lb,col='green3', lty = "dashed")
lines(x_train,w3_ub,col='green3', lty = "dashed")

