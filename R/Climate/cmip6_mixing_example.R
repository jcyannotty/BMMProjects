#------------------------------------------------
# CMIP6 Run for Unity
# sudo openconnect --user=yannotty.1  vpn.asc.ohio-state.edu
# ssh yannotty.1@unity.asc.ohio-state.edu
#------------------------------------------------
filedir = "/home/johnyannotty/Documents/CMIP6_mixing/"
datadir = "Data/North_Hemisphere/"
out_data = readRDS(paste0(filedir,datadir,"/CESM2_CNRM_NWH_June_14_04_18_24_n300.rds"))

setwd("/home/johnyannotty/Documents/openbt/src")
source("/home/johnyannotty/Documents/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Plotting/computer_expt_plots.R")

library(rBayesianOptimization)
library(lattice)
library(RColorBrewer)
library(plotly)
library(viridis)
library(latex2exp)
library(colorRamps)


# Format the data in n_train x 2 matrix
nu = 20
rho = 1
q0 = 4

sig2_hat = 0.05*min(apply(apply(out_data$f_train, 2, function(x) (x-out_data$y_train)^2),2,mean))

# Train the model
fit=train.openbtmixing(as.matrix(out_data$x_train),out_data$y_train,out_data$f_train,
                       pbd=c(1.0,0),ntree = 20,ntreeh=1,numcut=500,tc=4,model="mixbart",modelname="cmip6_mix",
                       ndpost = 5000, nskip = 2000, nadapt = 3000, adaptevery = 300, printevery = 500,
                       power = 0.5, base = 0.95, minnumbot = 2, overallsd = sqrt(sig2_hat), k = 1.0, overallnu = nu,
                       summarystats = FALSE, rpath = TRUE, q = q0, rshp1 = 2, rshp2 = 30,
                       stepwpert = 0.1, probchv = 0.1, batchsize = 500, maxd = 50)


#Get mixed mean function
batchsz = 5000
nb = ceiling(nrow(out_data$f_test)/batchsz)
tmp = 0.1
for(i in 2:nb){
  indst = (i-1)*batchsz + 1
  indend = min(i*batchsz,nrow(out_data$x_test))
  fitp=predict.openbtmixing(fit,x.test = out_data$x_test[indst:indend,], 
                            f.test = out_data$f_test[indst:indend,],
                            tc=4, q.lower = 0.025, q.upper = 0.975,
                            ptype = "mean_and_proj", proj_type = "euclidean", temperature = tmp)
  
  # Get wsum
  wsum = 0*fitp$wdraws[[1]]
  for(j in 1:ncol(out_data$f_test)){
    wsum = wsum + fitp$wdraws[[j]]
  }
  wsum_mean = apply(wsum,2,mean)
  wsum_lb = apply(wsum,2,quantile, 0.025)
  wsum_ub = apply(wsum,2,quantile, 0.975)
  
  
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
    proj_mean = fitp$pmmean,
    proj_ub = fitp$pm.upper,
    proj_lb = fitp$pm.lower,
    pwts_mean = fitp$pwmean,
    pwts_ub = fitp$pw.upper,
    pwts_lb = fitp$pw.lower,
    delta_mean = fitp$dmean,
    delta_lb = fitp$d.lower,
    delta_ub = fitp$d.upper,
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
    sims = out_data$sims,
    lat_bnds = out_data$lat_bnds,
    lon_bnds = out_data$lon_bnds,
    tmp = tmp
  )
  
  # Save results
  dt = gsub("/","_",format(Sys.time(), "%D_%H:%M"))
  sn = paste(unlist(sapply(out_data$sims,
                           function(x) strsplit(strsplit(x,"_Amon_")[[1]][2],"_historical")[[1]][1])
  ), collapse = "_")
  sn = gsub("ACCESS-CM2","ACC",sn)
  sn = gsub("BCC-CSM2-MR","BCC",sn)
  sn = gsub("MIROC-ES2L","MIROC",sn)
  sn = gsub("CMCC-CM2-SR5","CMCC",sn)
  sn = gsub("CNRM-CM6-1-HR","CNRM",sn)
  sn = gsub("KIOST-ESM","KIOST",sn)
  
  desc = "NorthAmerica_June_Dec_2014"
  
  if(i<10){
    bnstr = paste0("0",i)
  }else{
    bnstr = paste0(i)
  }
  saveRDS(fit_data,paste0(filedir,"Results/North_Hemisphere/CESM2_CNRM_North_America_June14_n300_hp3/",
                          sn,"_batch",bnstr,".rds"))
}

fitp=predict.openbtmixing(fit,x.test = out_data$x_test[1:100,], 
                          f.test = out_data$f_test[1:100,],
                          tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "sigma", proj_type = "euclidean", temperature = 0.5)

hist(fitp$sdraws[,1])
saveRDS(fit_data,paste0(filedir,"Results/North_Hemisphere/CESM2_CNRM_North_America_June14_n300_hp3/", sn,"_sdraws",".rds"))
                        

# Helper to order result names
order_fnames = function(fls){
  fls_split = sapply(fls,function(x) strsplit(x,"batch"))
  fls_new = fls
  for(i in 1:length(fls_split)){
    if(length(fls_split[[i]])>1){
      if(nchar(fls_split[[i]][2])< 7){
        fls_new[i] = paste0(fls_split[[i]][1], "batch0",fls_split[[i]][2])
      }  
    }
  }
  fls = fls[order(fls_new)]
  return(fls)
}

resdir = "Results/North_Hemisphere/CESM2_CNRM_North_America_June14_n300_hp1/"
fls = system(paste0("ls ",filedir,resdir),intern = TRUE)
fls = order_fnames(fls)
batch = 0
for(i in 1:length(fls)){
  if(grepl("batch",fls[i])){
    batch = batch + 1
    temp_fit = readRDS(paste0(filedir,resdir,fls[i]))
    if(batch == 1){
      fit1 = temp_fit
    }else{
      fit1 = list(
        pred_mean = c(fit1$pred_mean,temp_fit$pred_mean),
        pred_ub = c(fit1$pred_ub,temp_fit$pred_ub),
        pred_lb = c(fit1$pred_lb,temp_fit$pred_lb),
        wts_mean = rbind(fit1$wts_mean,temp_fit$wts_mean),
        wts_ub = rbind(fit1$wts_ub,temp_fit$wts_ub),
        wts_lb = rbind(fit1$wts_lb,temp_fit$wts_lb),
        wsum_mean = c(fit1$wsum_mean,temp_fit$wsum_mean),
        wsum_lb = c(fit1$wsum_lb,temp_fit$wsum_lb),
        wsum_ub = c(fit1$wsum_ub,temp_fit$wsum_ub),
        proj_mean = c(fit1$proj_mean,temp_fit$proj_mean),
        proj_ub = c(fit1$proj_ub,temp_fit$proj_ub),
        proj_lb = c(fit1$proj_lb,temp_fit$proj_lb),
        pwts_mean = rbind(fit1$pwts_mean,temp_fit$pwts_mean),
        pwts_ub = rbind(fit1$pwts_ub,temp_fit$pwts_ub),
        pwts_lb = rbind(fit1$pwts_lb,temp_fit$pwts_lb),
        delta_mean = c(fit1$delta_mean,temp_fit$delta_mean),
        delta_lb = c(fit1$delta_lb,temp_fit$delta_lb),
        delta_ub = c(fit1$delta_ub,temp_fit$delta_ub)
      )
    }
  }else if(grepl("sdraws",fls[i])){
    sfit1 = readRDS(paste0(filedir,resdir,fls[i]))
  }
}

usa = map_data("world",region = "USA")
world = map_data("world")
states = map_data("state")

# Convert long to (-180,180) scale
xs_test = out_data$x_test
xs_test[,"lon"] = ifelse(xs_test[,"lon"] > 180,xs_test[,"lon"]-360,xs_test[,"lon"])
xs_test[,1] = ifelse(xs_test[,1] < 0,xs_test[,1]-0.25,xs_test[,1])

scmin = -10; scmax = 45

pdagger = plot_mean2d_map_viridis(xs_test[,1:2],out_data$y_test,xcols = c(1,2), 
                                  viridis_opt = "viridis",
                                  scale_limit = c(scmin,scmax), title = "ERA5",
                                  maps_list = list(data.frame(world),data.frame(states)),
                                  maps_cols = c("grey30","grey30"),
                                  lat_bnds = out_data$lat_bnds,
                                  lon_bnds = out_data$lon_bnds
) 
pdagger = pdagger + labs(fill = wk_lab, x = "Longitude", y = "Latitude") + theme(axis.title = element_text(size=18))
pdagger = pdagger + theme(axis.text=element_text(size=12),axis.title=element_text(size=16), plot.title = element_text(size = 15),legend.title = element_text(size = 16))
pdagger = pdagger + theme(axis.title = element_blank(), axis.ticks = element_blank(),
                          axis.text = element_blank())


pb = plot_mean2d_map_viridis(xs_test[,1:2],fit1$pred_mean,xcols = c(1,2), 
                             viridis_opt = "viridis",
                             scale_limit = c(scmin,scmax), title = "BART-BMM",
                             maps_list = list(data.frame(world),data.frame(states)),
                             maps_cols = c("grey30","grey30"),
                             lat_bnds = out_data$lat_bnds,
                             lon_bnds = out_data$lon_bnds
) 
pb = pb + labs(fill = wk_lab, x = "Longitude", y = "Latitude") + theme(axis.title = element_text(size=18))
pb = pb + theme(axis.text=element_text(size=12),axis.title=element_text(size=16), plot.title = element_text(size = 15),legend.title = element_text(size = 16))
pb = pb + theme(axis.title = element_blank(),axis.text = element_blank(),
                axis.ticks = element_blank())
pb + theme(plot.margin=unit(c(-0.5,1,1,1),"cm"))

p1 = plot_mean2d_map_viridis(xs_test[,1:2],out_data$f_test[,1],xcols = c(1,2), 
                                  viridis_opt = "viridis",
                                  scale_limit = c(scmin,scmax), title = "BART-BMM",
                                  maps_list = list(data.frame(world),data.frame(states)),
                                  maps_cols = c("grey30","grey30"),
                                  lat_bnds = out_data$lat_bnds,
                                  lon_bnds = out_data$lon_bnds
) 


p2 = plot_mean2d_map_viridis(xs_test[,1:2],out_data$f_test[,2],xcols = c(1,2), 
                             viridis_opt = "viridis",
                             scale_limit = c(scmin,scmax), title = "BART-BMM",
                             maps_list = list(data.frame(world),data.frame(states)),
                             maps_cols = c("grey30","grey30"),
                             lat_bnds = out_data$lat_bnds,
                             lon_bnds = out_data$lon_bnds
) 


p_leg = g_legend(pb+theme(legend.position = "bottom", legend.key.size = unit(0.8,'cm'),
                           legend.key.width = unit(2.8,'cm'), 
                           legend.title = element_text(size=16))+
                   labs(fill=""))


grid.arrange(arrangeGrob(pdagger+theme(legend.position = "none",
                                       plot.margin=unit(c(0.125,0.5,0.5,0.1),"cm")),
                         pb+theme(legend.position = "none",
                                  plot.margin=unit(c(0.125,0.5,0.5,0.1),"cm")),
                         nrow = 1),
                         nrow=2, heights = c(1,0.1), p_leg)


grid.arrange(arrangeGrob(p1+theme(legend.position = "none"),
                         p2+theme(legend.position = "none"),
                         nrow = 1), nrow=2, heights = c(5,1), p_leg)


# RMSE
sqrt(min((out_data$f_test[,1] - out_data$y_test)^2))
sqrt(min((out_data$f_test[,2] - out_data$y_test)^2))
sqrt(mean((out_data$f_test[,2] - out_data$y_test)^2))
sqrt(mean((fit1$pred_mean - out_data$y_test)^2))


max(fit1$wts_mean)
min(fit1$wts_mean)

wtpal = plasma(6)
wtpal[6] = "yellow2"
wtpal[7] = "ivory2"
wtpal = c("black",wtpal)
wsv = seq(-0.1,1.3,length = 7)
wk_lab = bquote(hat(w)[k]*"(x)")

# BART-BMM Wts
# Wt1
w1b = plot_mean2d_map_gradient(xs_test,fit1$wts_mean[,1],xcols = c(1,2), 
                               title="CESM2", 
                               scale_colors = wtpal, scale_vals = wsv,
                               maps_list = list(data.frame(world),data.frame(states)),
                               maps_cols = c("white","white"),
                               lat_bnds = out_data$lat_bnds,
                               lon_bnds = c(min(xs_test[,'lon']),max(xs_test[,'lon'])))
w1b = w1b + labs(fill = wk_lab, x = "Longitude", y = "Latitude") + theme(axis.title = element_text(size=18))
w1b = w1b + theme(axis.text=element_text(size=12),axis.title=element_text(size=18), plot.title = element_text(size = 15),legend.title = element_text(size = 16))

#wsv = c(-0.1,0,0.2,0.3,0.4,0.5,1.1)

# wt2 
w2b = plot_mean2d_map_gradient(xs_test,fit1$wts_mean[,2],xcols = c(1,2), 
                               title="CNRM", 
                               scale_colors = wtpal, scale_vals = wsv,
                               maps_list = list(data.frame(world),data.frame(states)),
                               maps_cols = c("white","white"),
                               lat_bnds = out_data$lat_bnds,
                               lon_bnds = c(min(xs_test[,'lon']),max(xs_test[,'lon'])))
w2b = w2b + labs(fill = wk_lab, x = "Longitude", y = "Latitude") + theme(axis.title = element_text(size=18))
w2b = w2b + theme(axis.text=element_text(size=12),axis.title=element_text(size=18), plot.title = element_text(size = 15),legend.title = element_text(size = 16))


w_leg = g_legend(w1b+theme(legend.position = "bottom", legend.key.size = unit(0.8,'cm'),
                           legend.key.width = unit(2.8,'cm'), 
                           legend.title = element_text(size=16) )+
                   labs(fill=wk_lab))



grid.arrange(arrangeGrob(w1b+theme(legend.position = "none"),
                         w2b+theme(legend.position = "none"),
                         nrow = 1), nrow=2, heights = c(5,1), w_leg)


h = which(xs_test[,1] == -145)
plot(fit1$wts_mean[h,2], ylim = c(-0.5,1.5))
lines(fit1$wts_lb[h,2])
lines(fit1$wts_ub[h,2])


sum(fit1$wts_mean[5,]*out_data$f_test[5,])
fit1$pred_mean[1:5]


r1 = plot_mean2d_map_gradient(xs_test[,1:2],out_data$y_test - out_data$f_test[,1],xcols = c(1,2), 
                            scale_vals = c(-20,0,6), title = "BART-BMM",
                             maps_list = list(data.frame(world),data.frame(states)),
                             maps_cols = c("grey30","grey30"),
                             lat_bnds = out_data$lat_bnds,
                             lon_bnds = out_data$lon_bnds
) 


r2 = plot_mean2d_map_gradient(xs_test[,1:2],out_data$y_test - out_data$f_test[,2],xcols = c(1,2), 
                             scale_vals = c(-20,0,20), title = "BART-BMM",
                             maps_list = list(data.frame(world),data.frame(states)),
                             maps_cols = c("grey30","grey30"),
                             lat_bnds = out_data$lat_bnds,
                             lon_bnds = out_data$lon_bnds
) 

