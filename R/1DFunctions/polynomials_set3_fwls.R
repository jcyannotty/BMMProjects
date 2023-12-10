#-------------------------------------------------
# Mixing polynomials - 3 Models in 1D
# Feature Weighted Linear Stacking
#-------------------------------------------------
source("/home/johnyannotty/Documents/openbt/R/polynomials.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Samplers/sampler_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/BMM/fwls.R")

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
#-------------------------------------------------
# Get wts basis for 1D input
basis = c("linear","quad","linear")
g_train = fwls_construct_basis(x_train, K, basis)
g_test = fwls_construct_basis(x_test, K, basis)

fit_fwls = fwls(y_train, f_train, g_train, lambda = 3)
fitp_fwls = fwls_predict(fit_fwls,f_test,g_test)

plot(x_test, f0_test, pch = 16, cex = 0.8, main = 'Fits', type = 'l',ylim = c(-15,15))
points(x_train, y_train, pch = 3)
lines(x_test, f_test[,1], col = 'red', lty = 2)
lines(x_test, f_test[,2], col = 'blue', lty = 2)
lines(x_test, f_test[,3], col = 'green4', lty = 2)
lines(x_test, fitp_fwls$fx, col = 'purple', lwd = 2)
lines(x_test, fitp_fwls$fx_lb, col = 'orange', lwd = 2, cex = 0.5)
lines(x_test, fitp_fwls$fx_ub, col = 'orange', lwd = 2, cex = 0.5)


plot(x_test, fitp_fwls$wx[,1], pch = 16, col = 'red', type = 'l', ylim = c(-5,5), lwd = 2, 
     panel.first = {grid(col = 'lightgrey')})
lines(x_test, fitp_fwls$wx[,2], col = 'blue', lwd = 2)
lines(x_test, fitp_fwls$wx[,3], col = 'green3', lwd = 2)
lines(x_test, fitp_fwls$wx_ub[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, fitp_fwls$wx_lb[,1], col = 'red', lty = 'dashed', lwd = 1)
lines(x_test, fitp_fwls$wx_ub[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(x_test, fitp_fwls$wx_lb[,2], col = 'blue', lty = 'dashed', lwd = 1)
lines(x_test, fitp_fwls$wx_ub[,3], col = 'green4', lty = 'dashed', lwd = 1)
lines(x_test, fitp_fwls$wx_lb[,3], col = 'green4', lty = 'dashed', lwd = 1)

abline(h = 1, col = 'grey', lty = 'dashed')
abline(h = 0, col = 'grey', lty = 'dashed')

