#-----------------------------------------------------
# 2D-Taylor Series Expansions
#-----------------------------------------------------
setwd("/home/johnyannotty/Documents/openbt/src")

# Load the R wrapper functions to the OpenBT library.
source("/home/johnyannotty/Documents/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/openbt/R/polynomials.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Plotting/computer_expt_plots.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Samplers/sampler_functions.R")

library(plotly)
library(viridis)
library(latex2exp)
library(geoR)

# Out file directory
filedir = '/home/johnyannotty/Documents/Dissertation'

#-----------------------------------------------------
#-----------------------------------------------------
set.seed(2)
n_train = 100 
n_test = 25
xoff = 0.10
xmin = -pi
xmax = pi
s = 0.05

# Generate data
x1_train = seq(xmin + xoff, xmax - xoff, length = sqrt(n_train))
x_train = as.matrix(expand.grid(x1_train,x1_train))

x1_test = round(seq(xmin,xmax,length = n_test),4)
x_test = as.matrix(expand.grid(x1_test,x1_test))

plot(x_train[,1],x_train[,2])
plot(x_test[,1],x_test[,2], pch = 16, col = 'grey70')
points(x_train[,1],x_train[,2], pch = 16, col = 'blue')

# Generate True function
f0_train = sin(x_train[,1]) + cos(x_train[,2])
f0_test = sin(x_test[,1]) + cos(x_test[,2])

set.seed(99)
y_train = f0_train + rnorm(n_train,0,s)


f1 = function(x,xc1,xc2){
  out = 1.05*ifelse(x[1]>xc1,sin(x[1]),sin(xc1)) + 1*ifelse(x[2]>xc2,cos(x[2]),cos(xc2))
  return(out)
}

f2 = function(x,xc1,xc2){
  out = 1*ifelse(x[1]<xc1,sin(x[1]),sin(xc1)) + 1.05*ifelse(x[2]>xc2,cos(x[2]),cos(xc2))
  return(out)
}

f3 = function(x,xc1,xc2){
  out = 0.95*ifelse(x[1]<xc1,sin(x[1]),sin(xc1)) + 1.05*ifelse(x[2]<xc2,cos(x[2]),cos(xc2))
  return(out)
}

f4 = function(x,xc1,xc2){
  out = 1.1*ifelse(x[1]>xc1,sin(x[1]),sin(xc1)) + 1.1*ifelse(x[2]<xc2,cos(x[2]),cos(xc2))
  return(out)
}

f_train = matrix(0, nrow = nrow(x_train), ncol = 4)
f_test = matrix(0, nrow = nrow(x_test), ncol = 4)

xc1 = c(0,0)
xc2 = c(0,0)
xc3 = c(0,0)
xc4 = c(0,0)

f_train[,1] = apply(x_train, 1, function(x) f1(x,xc1[1],xc1[2]))
f_train[,2] = apply(x_train, 1, function(x) f2(x,xc2[1],xc2[2]))
f_train[,3] = apply(x_train, 1, function(x) f3(x,xc3[1],xc3[2]))
f_train[,4] = apply(x_train, 1, function(x) f4(x,xc4[1],xc4[2]))


f_test[,1] = apply(x_test, 1, function(x) f1(x,xc1[1],xc1[2]))
f_test[,2] = apply(x_test, 1, function(x) f2(x,xc2[1],xc2[2]))
f_test[,3] = apply(x_test, 1, function(x) f3(x,xc3[1],xc3[2]))
f_test[,4] = apply(x_test, 1, function(x) f4(x,xc4[1],xc4[2]))


pdagger = plot_mean2d_viridis(x_test,matrix(f0_test), 
                         viridis_opt = "turbo",
                         scale_limit = c(-2.2,2.2), title = TeX("$f_\\dagger(x)$")) 
pdagger = pdagger + labs(fill = TeX("$f_\\dagger(x)$"), x = bquote(x[1]), y = element_blank()) + theme(axis.title = element_text(size=18))
pdagger = pdagger + theme(axis.text=element_text(size=12),axis.title=element_text(size=18), plot.title = element_text(size = 15),legend.title = element_text(size = 16))


# Function path
fp_list = list()
K = ncol(f_test)
wk_lab = bquote(hat(w)[k]*"(x)")
for(i in 1:K){
  fk_lab = bquote(f[.(i)]*"(x)")
  p1 = plot_mean2d_viridis(x_test,matrix(f_test[,i]), 
                                viridis_opt = "turbo",
                                scale_limit = c(-2.2,2.2), title = fk_lab) 
  p1 = p1 + labs(fill = fk_lab, x = bquote(x[1]), y = element_blank()) + theme(axis.title = element_text(size=18))
  p1 = p1 + theme(axis.text=element_text(size=12),axis.title=element_text(size=18), plot.title = element_text(size = 15),legend.title = element_text(size = 16))
  
  if(i == 1){pleg = g_legend(p1 + theme(legend.position = "bottom", 
                                         legend.key.width = unit(2.2,"cm"))+ 
                               labs(fill =  bquote(f[k]*"(x)")))}
  fp_list[[i]] = p1 + theme(legend.position = "none")
}


grid.arrange(arrangeGrob(pdagger + theme(legend.position = "none"),
                         fp_list[[1]],fp_list[[2]],fp_list[[3]],fp_list[[4]],
                         nrow = 1,left = grid.text(bquote(x[2]), rot = 90, 
                                                   gp = gpar(fontsize = 18))), 
             nrow = 2, pleg, heights = c(1,0.1))


grid.arrange(arrangeGrob(pdagger + theme(legend.position = "none"),
                         fp_list[[1]],fp_list[[2]],
                         nrow = 1,left = grid.text(bquote(x[2]), rot = 90, 
                                                   gp = gpar(fontsize = 18)),
             layout_matrix = rbind(c(1, 2, 3))),
             arrangeGrob(fp_list[[3]]+labs(y = bquote(x[2])),fp_list[[4]],
                         nrow = 1,
                         layout_matrix = rbind(c(NA,1,2,NA)), widths = c(0.5,1,1,0.5)),
             nrow = 3, pleg, heights = c(1,1,0.18))


# Residuals
rp_list = list()
for(i in 1:K){
  fk_lab = bquote(f[.(i)]*"(x)")
  p1 = plot_mean2d_gradient(x_test,matrix(f0_test - f_test[,i]), 
                            scale_vals = c(-1,0,1), title = fk_lab) 
  p1 = p1 + labs(fill = wk_lab, x = bquote(x[1]), y = element_blank()) + theme(axis.title = element_text(size=18))
  p1 = p1 + theme(axis.text=element_text(size=12),axis.title=element_text(size=18), plot.title = element_text(size = 15),legend.title = element_text(size = 16))
  
  if(i == 1){pleg = g_legend(p1 + theme(legend.position = "bottom", 
                                        legend.key.width = unit(2.2,"cm"))+ 
                               labs(fill =  bquote(hat(w)[k]*"(x)")))}
  rp_list[[i]] = p1 + theme(legend.position = "none")
}


grid.arrange(arrangeGrob(rp_list[[1]],rp_list[[2]],rp_list[[3]],rp_list[[4]],
                         nrow = 1,left = grid.text(bquote(x[2]), rot = 90, 
                                                   gp = gpar(fontsize = 18))), 
             nrow = 2, pleg, heights = c(1,0.1))

#-----------------------------------------------------
# Model Training 
#-----------------------------------------------------
nu = 20
rho = 1
sig2_hat = max(apply(apply(f_train, 2, function(x) (x-y_train)^2),2,min))
lam = rho*sig2_hat*(nu+2)/nu
q0 = 4
fit=openbt(x_train,y_train,f_train,pbd=c(1.0,0),ntree = 10,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="mixmodel",
           ndpost = 10000, nskip = 2000, nadapt = 4000, adaptevery = 500, printevery = 500,
           power = 1.0, base = 0.95, minnumbot = 3, overallsd = sqrt(sig2_hat), k = 1.0, overallnu = nu,
           summarystats = FALSE, selectp = FALSE, rpath = TRUE, q = 4.0, rshp1 = 2, rshp2 = 10,
           stepwpert = 0.1, probchv = 0.1)


#Get mixed mean function
fitp=predict.openbt(fit,x.test = x_test, f.test = f_test,tc=4, q.lower = 0.025, q.upper = 0.975)
pb1 = plot_mean2d_viridis(x_test,matrix(fitp$mmean), 
                         viridis_opt = "turbo",
                         scale_limit = c(-2.2,2.2), title = "RPBART-BMM\n") 
pb1 = pb1 + labs(fill = wk_lab, x = bquote(x[1]), y = element_blank()) + theme(axis.title = element_text(size=18))
pb1 = pb1 + theme(axis.text=element_text(size=12),axis.title=element_text(size=18), plot.title = element_text(size = 15),legend.title = element_text(size = 16))

rb1 = plot_mean2d_gradient(x_test,matrix(f0_test - fitp$mmean), 
                           scale_vals = c(-0.1,0,0.1), title = TeX("$\\hat{r}(x) = f_\\dagger(x) - \\hat{f}_\\dagger(x)$") ) 
rb1 = rb1 + labs(fill = wk_lab, x = bquote(x[1]), y = element_blank()) + theme(axis.title = element_text(size=18))
rb1 = rb1 + theme(axis.text=element_text(size=12),axis.title=element_text(size=18), plot.title = element_text(size = 15),legend.title = element_text(size = 16))


pbw1 = plot_mean2d_viridis(x_test,matrix(fitp$m.upper - fitp$m.lower), 
                          viridis_opt = "inferno",
                          scale_limit = c(0.05,0.45), title = "Credible Interval Width \n") 
pbw1 = pbw1 + labs(fill = wk_lab, x = bquote(x[1]), y = element_blank()) + theme(axis.title = element_text(size=18))
pbw1 = pbw1 + theme(axis.text=element_text(size=12),axis.title=element_text(size=18), plot.title = element_text(size = 15),legend.title = element_text(size = 16))


grid.arrange(arrangeGrob(pb1 + theme(legend.position = "bottom", legend.key.width = unit(1.5,"cm"),
                                     legend.text = element_text(size = 15)) + 
                           labs(fill = TeX("$\\hat{f}_\\dagger(x)$  ")),
                         rb1 + theme(legend.position = "bottom", legend.key.width = unit(1.5,"cm"),
                                     legend.text = element_text(size = 15)) + 
                           labs(fill = TeX("$\\hat{r}(x)$  ")),
                         pbw1 + theme(legend.position = "bottom", legend.key.width = unit(1.5,"cm"),
                                      legend.text = element_text(size = 15)) + 
                           labs(fill = "Width"),nrow = 1,
                         left = grid.text(bquote(x[2]), rot = 90, 
                                                   gp = gpar(fontsize = 18))), 
             nrow = 2, heights = c(1,0.05))

#Get the model weights
K = 4
fitw = mixingwts.openbt(fit, x.test = x_test, numwts = K, tc = 4)

wtpal = plasma(6)
wtpal[6] = "yellow2"
wtpal = c("black",wtpal)
wsv = seq(-0.4,0.8,length = 7)
wk_lab = bquote(hat(w)[k]*"(x)")

# BART-BMM Wts
# Wt1
wp_list = list()
for(i in 1:K){
  wk_lab = bquote(hat(w)[.(i)]*"(x)")
  w1b = plot_mean2d_map_gradient(x_test,fitw$wmean[,i],xcols = c(1,2), 
                                 title=wk_lab, 
                                 scale_colors = wtpal, scale_vals = wsv)
  w1b = w1b + labs(fill = wk_lab, x = bquote(x[1]), y = element_blank()) + theme(axis.title = element_text(size=18))
  w1b = w1b + theme(axis.text=element_text(size=12),axis.title=element_text(size=18), plot.title = element_text(size = 15),legend.title = element_text(size = 16))
  
  if(i == 1){wleg = g_legend(w1b + theme(legend.position = "bottom", 
                                         legend.key.width = unit(2.2,"cm"))+ 
                               labs(fill =  bquote(hat(w)[k]*"(x)")))}
  wp_list[[i]] = w1b + theme(legend.position = "none")
}

grid.arrange(arrangeGrob(wp_list[[1]],wp_list[[2]],wp_list[[3]],wp_list[[4]],
                         nrow = 1,left = grid.text(bquote(x[2]), rot = 90, 
                                                   gp = gpar(fontsize = 18))), 
             nrow = 2, wleg, heights = c(1,0.1))

hist(unlist(fitp$sdraws[,1]))
plot(unlist(fitp$sdraws[,1]))
mean(unlist(fitp$sdraws[,1]))

# Get wsum
wsum = 0*fitw$wdraws[[1]]
for(i in 1:K){
  wsum = wsum + fitw$wdraws[[i]]
}
wsum_mean = apply(wsum,2,mean)
wsum_lb = apply(wsum,2,quantile, 0.025)
wsum_ub = apply(wsum,2,quantile, 0.975)


# Gamma
gp = gammapost.openbtmixing(fit)
hist(gp[,1])
hist(gp[,2])
summary(gp)


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

filedir = "/home/johnyannotty/Documents/Dissertation/results/2d_functions/"
saveRDS(fit_data, paste0(filedir,"sincos_4k_res_rpath_05_07_24.rds"))

ms = list(
  y_train = y_train,
  x_train = x_train,
  x_test = x_test,
  f_train = f_train,
  f_test = f_test,
  f0_test = f0_test,
  s = s
)

filedir = '/home/johnyannotty/Documents/Dissertation/results/2d_functions/'
saveRDS(ms, paste0(filedir,"ms_sincos_n100.rds"))



write.csv(ms$f_train, paste0(filedir,"sincos_4k_ftrain_rpath_05_07_24.txt"), row.names = FALSE)
write.csv(ms$f_test, paste0(filedir,"sincos_4k_ftest_rpath_05_07_24.txt"), row.names = FALSE)
write.csv(ms$f0_test, paste0(filedir,"sincos_4k_f0test_rpath_05_07_24.txt"), row.names = FALSE)
write.csv(ms$x_train, paste0(filedir,"sincos_4k_xtrain_rpath_05_07_24.txt"), row.names = FALSE)
write.csv(ms$x_test, paste0(filedir,"sincos_4k_xtest_rpath_05_07_24.txt"), row.names = FALSE)
write.csv(ms$y_train, paste0(filedir,"sincos_4k_ytrain_rpath_05_07_24.txt"), row.names = FALSE)
