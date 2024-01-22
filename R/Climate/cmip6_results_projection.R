#------------------------------------------------
# Climate CMIP6 results loader - from Unity
#------------------------------------------------
setwd("/home/johnyannotty/Documents/openbt/src")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")

library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)
library(plotly)
library(viridis)
library(latex2exp)

#------------------------------------------------
#------------------------------------------------
filedir = "/home/johnyannotty/Documents/CMIP6_mixing/"
datadir = "Data/North_Hemisphere/"
resdir = "Results/North_Hemisphere/"

#------------------------------------------------
# Batch Loader
#------------------------------------------------
ffold = "nh_6month_2014_n30k_hp1/"
dataname = "ACC_BCC_CESM2_CNRM_NH_6M14_01_06_24_n30000.rds"
fls = system(paste0("ls ",filedir,resdir,ffold),intern = TRUE)
batch = 0
ms = readRDS(paste0(filedir,datadir,dataname))
for(i in 1:length(fls)){
  if(grepl("batch",fls[i])){
    batch = batch + 1
    temp_fit = readRDS(paste0(filedir,resdir,ffold,fls[i]))
    if(batch == 1){
      fit = temp_fit
    }else{
      fit = list(
        pred_mean = c(fit$pred_mean,temp_fit$pred_mean),
        pred_ub = c(fit$pred_ub,temp_fit$pred_ub),
        pred_lb = c(fit$pred_lb,temp_fit$pred_lb),
        wts_mean = rbind(fit$wts_mean,temp_fit$wts_mean),
        wts_ub = rbind(fit$wts_ub,temp_fit$wts_ub),
        wts_lb = rbind(fit$wts_lb,temp_fit$wts_lb),
        wsum_mean = c(fit$wsum_mean,temp_fit$wsum_mean),
        wsum_lb = c(fit$wsum_lb,temp_fit$wsum_lb),
        wsum_ub = c(fit$wsum_ub,temp_fit$wsum_ub)
        #proj_mean = temp_fit$pmmean,
        #proj_ub = temp_fit$pm.upper,
        #proj_lb = temp_fit$pm.lower,
        #pwts_mean = temp_fit$pwmean,
        #pwts_ub = temp_fit$pw.upper,
        #pwts_lb = temp_fit$pw.lower,
        #delta_mean = temp_fit$dmean,
        #delta_lb = temp_fit$d.lower,
        #delta_ub = temp_fit$d.upper
      )
    }
  }else if(grepl("sdraws",fls[i])){
    sfit = readRDS(paste0(filedir,resdir,ffold,fls[i]))
  }
}


resid = fit$pred_mean - ms$y_test
sqrt(mean(resid^2))

# Resid Plot
rb = plot_residuals_hm_gg2(ms$x_test,resid,xcols = c(1,2), title="BMM Mean Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-10,0,10)) #c(-3.5,0,3.5)
rb = rb + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")


# Projection
library(rBayesianOptimization)
h = sample(1:nrow(ms$f_test),5000)
softmax_l2 = function(tmp){
  fita = rowSums(ms$f_test[h,]*exp(fit$wts_mean[h,]/tmp)/rowSums(exp(fit$wts_mean[h,]/tmp)))
  score = -sum((fit$pred_mean[h] - fita)^2) # max the negative loss
  return(list(Score = score, Pred = 0))
} 

tmp_bounds = list(tmp = c(0.02,1))
tmp_grid = seq(0.01,0.95, length = 200)
init_grid_dt = data.frame(tmp = seq(0.01,0.95, length = 20))
bayes_temp = BayesianOptimization(FUN = softmax_l2, acq = "ei",
                                  bounds = tmp_bounds, init_grid_dt = init_grid_dt,
                                  init_points = 0.8, n_iter = 30)
loss_grid = sapply(tmp_grid, function(x) softmax_l2(x)$Score)
plot(tmp_grid,loss_grid)
points(bayes_temp$Best_Par,bayes_temp$Best_Value, pch = 3, col = 'red')


