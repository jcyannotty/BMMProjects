setwd("/home/johnyannotty/Documents/openbt/src")

# Load the R wrapper functions to the OpenBT library.
source("/home/johnyannotty/Documents/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/openbt/R/polynomials.R")
source("/home/johnyannotty/Documents/openbt/R/mixing_priors.R")
#source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")
source("C:/Users/johny/Documents/packages/OpenBT/R/eft_mixing_helper_functions.R")
source("C:/Users/johny/Documents/packages/BayesToolBox/bayestb/Samplers/sampler_functions.R")

library(plotly)
library(viridis)
library(latex2exp)
library(geoR)

# Out file directory
filedir = '/home/johnyannotty/Documents/RandomPathBART/VariogramResults/2d_functions/'

#-----------------------------------------------------
#-----------------------------------------------------
set.seed(23)
n_train = 100 
n_test = 25
xoff = 0.10
xmin = -pi
xmax = pi
s = 0.1

# Generate data
x_train = grid_2d_design(n1=10,n2=10, xmin = c(xmin,xmin), xmax = c(xmax,xmax))
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

#-----------------------------------------------------
# Model Training 
#-----------------------------------------------------
nu = 40
rho = 1
sig2_hat = 0.1
lam = rho*sig2_hat*(nu+2)/nu
q0 = 4


fit=train.openbtmixing(x_train,y_train,as.matrix(rep(1,n_train)),pbd=c(1.0,0),ntree = 20,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="mixmodel",
                       ndpost = 5000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
                       power = 1.0, base = 0.95, minnumbot = 3, overallsd = sqrt(sig2_hat), k = 1.0, overallnu = nu,
                       summarystats = FALSE, rpath = TRUE, q = q0, rshp1 = 2, rshp2 = 20,
                       stepwpert = 0.1, probchv = 0.1, batchsize = 5000)

#openbt.save(fit,fname = "/home/johnyannotty/Documents/temp_res.obt")


#Get mixed mean function
fitp=predict.openbtmixing(fit,x.test = x_test, f.test = as.matrix(rep(1,n_test^2)),tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "mean_and_sigma", proj_type = "euclidean", temperature = 0.8)

#Get mixed mean function
p1 = plot_pred_2d_gg2(x_test,matrix(fitp$mmean),scale_vals = c(-2.0,0,2), 
                      title = "Mean Prediction from BMM", viridis_opt = "turbo")

p1 = plot_pred_2d_gg2(x_test,matrix(f0_test),scale_vals = c(-2.7,-0.25,2), 
                      title = "Mean Prediction from BMM", viridis_opt = "turbo")

p1 = plot_pred_2d_gg2(x_test,matrix(unlist(fitp$mdraws[1,])),scale_vals = c(-2,0,2), 
                      title = "Mean Prediction from BMM", viridis_opt = "turbo")


sqrt(mean((fitp$mmean - f0_test)^2))

hist(unlist(fitp$sdraws[,1]))

# Save results
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
  lam = fit$overalllambda
)

# Save results
filedir = "/home/johnyannotty/Documents/Dissertation/results/rpath_bart/"
saveRDS(fit_data,paste0(filedir,"2d_taylor_bart_res4_04_11_24.rds"))
saveRDS(fitp$sdraws[,1], paste0(filedir,"2d_taylor_bart_sdraws4_04_11_24.rds"))


#-----------------------------------------------------
# Plots
#-----------------------------------------------------
#filedir = "/home/johnyannotty/Documents/Dissertation/results/rpath_bart/2dtaylor/"
filedir = "D:/Dissertation/results/rpath_bart/2dtaylor/"
rp1 = readRDS(paste0(filedir,"2d_taylor_res1_04_11_24.rds"))
rp2 = readRDS(paste0(filedir,"2d_taylor_res2_04_11_24.rds"))
rp4 = readRDS(paste0(filedir,"2d_taylor_res3_04_11_24.rds"))
rp3 = readRDS(paste0(filedir,"2d_taylor_res4_04_11_24.rds"))


c(rp1$shp1,rp1$shp2)
c(rp2$shp1,rp2$shp2)
c(rp3$shp1,rp3$shp2)
c(rp4$shp1,rp4$shp2)


str(rp1)

xt_min = apply(x_train,2,min);  xt_max = apply(x_train,2,max); 
x1_test = round(seq(xt_min[1]-0.1,xt_max[1]+0.1,length = n_test),4)
x2_test = round(seq(xt_min[1]-0.1,xt_max[1]+0.1,length = n_test),4)
x_test = expand.grid(x1_test,x2_test)


pd = plot_pred_2d_gg2(x_test,matrix(f0_test),scale_vals = c(-2,0,2), 
                      title = "True System", viridis_opt = "turbo")

p1 = plot_pred_2d_gg2(x_test,matrix(rp1$pred_mean),scale_vals = c(-2,0,2), 
                      title = TeX("$\\alpha_1 = 2, \\alpha_2 = 20$"), viridis_opt = "turbo")

p2 = plot_pred_2d_gg2(x_test,rp2$pred_mean,scale_vals = c(-2,0,2), 
                      title = TeX("$\\alpha_1 = 2, \\alpha_2 = 10$"), viridis_opt = "turbo")

p3 = plot_pred_2d_gg2(x_test,rp3$pred_mean,scale_vals = c(-2,0,2), 
                      title = TeX("$\\alpha_1 = 2, \\alpha_2 = 5$"), viridis_opt = "turbo")

p4 = plot_pred_2d_gg2(x_test,rp4$pred_mean,scale_vals = c(-2,0,2), 
                      title = TeX("$\\alpha_1 = 10, \\alpha_2 = 10$"), viridis_opt = "turbo")

pleg = g_legend(p1 + theme(legend.position = "bottom", legend.key.width = unit(2.2,"cm"))+
                  labs(fill = ""))

grid.arrange(arrangeGrob(pd + theme(legend.position = "none",axis.text = element_text(size = 14),
                                    axis.title = element_text(size = 14)) + 
                           labs(x = TeX("$x_1"),y = TeX("$x_2")),
                         p1 + theme(legend.position = "none", axis.text = element_text(size = 14),
                                    axis.title = element_text(size = 14)) + 
                           labs(x = TeX("$x_1"),y = TeX("$x_2")),
                         p3 + theme(legend.position = "none", axis.text = element_text(size = 14),
                                    axis.title = element_text(size = 14)) + 
                           labs(x = TeX("$x_1"),y = TeX("$x_2")),
                         p4 + theme(legend.position = "none", axis.text = element_text(size = 14),
                                    axis.title = element_text(size = 14)) + 
                           labs(x = TeX("$x_1"),y = TeX("$x_2")), 
                         ncol  = 2), 
             nrow = 2, pleg, heights = c(10,1))


y_train = as.vector(read.csv("/home/johnyannotty/Downloads/tempy.txt", header = FALSE))$V1
x_train = as.matrix(read.csv("/home/johnyannotty/Downloads/tempx.txt", header = FALSE))
