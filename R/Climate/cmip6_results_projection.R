#------------------------------------------------
# Climate CMIP6 results loader - from Unity
#------------------------------------------------
setwd("/home/johnyannotty/Documents/openbt/src")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Plotting/computer_expt_plots.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Optimization/optimization_helpers.R")

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
datadir = "Data/World/"
resdir = "Results/World/"

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


# Mean discrepancy
project_wts = function(pmean,wmean,f_test,tmp){
  ahat = apply(wmean,1,function(x) sparse_gen(x,tmp))
  fita = rowSums(f_test*t(ahat))
  delta = pmean - fita
  return(list(ahat = t(ahat), delta = delta))
} 


#------------------------------------------------
# Batch Loader
#------------------------------------------------
ffold = "world_3m_2014_n45k_hp1/"
dataname = "ACC_BCC_MIROC_CMCC_CESM2_CNRM_CanESM5_KIOST_W_3M_2014_01_26_24_n45000.rds"
fls = system(paste0("ls ",filedir,resdir,ffold),intern = TRUE)
fls = order_fnames(fls)
ms = readRDS(paste0(filedir,datadir,dataname))
batch = 0
for(i in 1:length(fls)){
  if(grepl("batch",fls[i])){
    batch = batch + 1
    temp_fit = readRDS(paste0(filedir,resdir,ffold,fls[i]))
    if(batch == 1){
      fit1 = temp_fit
    }else{
      fit1 = list(
        pred_mean = c(fit1$pred_mean,temp_fit$pred_mean),
        pred_ub = c(fit1$pred_upper,temp_fit$pred_ub),
        pred_lb = c(fit1$pred_lower,temp_fit$pred_lb),
        wts_mean = rbind(fit1$wts_mean,temp_fit$wts_mean),
        wts_ub = rbind(fit1$wts_ub,temp_fit$wts_ub),
        wts_lb = rbind(fit1$wts_lbr,temp_fit$wts_lb)
        #wsum_mean = c(fit1$wsum_mean,temp_fit$wsum_mean),
        #wsum_lb = c(fit1$wsum_lb,temp_fit$wsum_lb),
        #wsum_ub = c(fit1$wsum_ub,temp_fit$wsum_ub),
        #proj_mean = c(fit1$proj_mean,temp_fit$proj_mean),
        #proj_ub = c(fit1$proj_ub,temp_fit$proj_ub),
        #proj_lb = c(fit1$proj_lb,temp_fit$proj_lb),
        #pwts_mean = rbind(fit1$pwts_mean,temp_fit$pwts_mean),
        #pwts_ub = rbind(fit1$pwts_ub,temp_fit$pwts_ub),
        #pwts_lb = rbind(fit1$pwts_lb,temp_fit$pwts_lb),
        #delta_mean = c(fit1$delta_mean,temp_fit$delta_mean),
        #delta_lb = c(fit1$delta_lb,temp_fit$delta_lb),
        #delta_ub = c(fit1$delta_ub,temp_fit$delta_ub)
      )
    }
  }else if(grepl("sdraws",fls[i])){
    sfit1 = readRDS(paste0(filedir,resdir,ffold,fls[i]))
  }
}

# Time period and E/W hemisphere
ntp = 3
tp_len = nrow(ms$f_test)/ntp
for(i in 1:ntp){
  assign(paste0("h",i),(tp_len*(i-1)+1):(tp_len*i))
}

# World Map
usa = map_data("world",region = "USA")
world = map_data("world")
states = map_data("state")

# Simulator names
sn = sapply(ms$sims,function(x) strsplit(strsplit(x,"_Amon_")[[1]][2],"_historical")[[1]][1])
sn = gsub("ACCESS-CM2","ACC",sn)
sn = gsub("BCC-CSM2-MR","BCC",sn)
sn = gsub("MIROC-ES2L","MIROC",sn)
sn = gsub("CMCC-CM2-SR5","CMCC",sn)
sn = gsub("CNRM-CM6-1-HR","CNRM",sn)
sn = gsub("KIOST-ESM","KIOST",sn)
names(sn) = NULL

# Convert long to (-180,180) scale
xs_test = ms$x_test
xs_test[,"lon"] = ifelse(xs_test[,"lon"] > 180,xs_test[,"lon"]-360,xs_test[,"lon"])
xs_test[,1] = ifelse(xs_test[,1] < 0,xs_test[,1]-0.25,xs_test[,1])


# Pointwise Projection shrinkage
tmp_grid = seq(0,0.95, by = 0.05)
proj_res = list()
for(i in 1:length(tmp_grid)){
  proj_res[[i]] = project_wts(fit1$pred_mean[h1], fit1$wts_mean[h1,], ms$f_test[h1,], tmp_grid[i])
  cat("progress: ", round(i/length(tmp_grid),4), "...\r")
}

#for(i in 1:length(tmp_grid)){
#  proj_res[[i]]$ahat = t(proj_res[[i]]$ahat)
#}


# View Projections....
h = h1
wtpal = plasma(6)
wtpal[6] = "yellow2"
wtpal = c("black",wtpal)
wsv = seq(0,1,length = 7)

# Default Labels
wk_lab = bquote(hat(w)[k]*"(x)")
K = ncol(ms$f_test)
wplot_list = list()
j = 5
for(i in 1:length(tmp_grid)){
  w0b = plot_mean2d_map_gradient(xs_test[h,],proj_res[[i]]$ahat[,j],xcols = c(1,2), 
                                 title=paste0(sn[j],": T = ", tmp_grid[i]), 
                                 scale_colors = wtpal, scale_vals = wsv,
                                 maps_list = list(data.frame(world),data.frame(states)),
                                 maps_cols = c("white","white"),
                                 lat_bnds = ms$lat_bnds,
                                 lon_bnds = c(min(xs_test[,'lon']),max(xs_test[,'lon'])))
  w0b = w0b + labs(fill = wk_lab, x = "Longitude", y = "Latitude") + 
    theme(axis.title = element_text(size=18))
  w0b = w0b + theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=18), 
                    plot.title = element_text(size = 15),legend.title = element_text(size = 16))
  #assign(paste0("w",i,"sf"), w0b)
  wplot_list[[i]] = w0b
}

w_leg = g_legend(w0b+theme(legend.position = "bottom", legend.key.size = unit(0.8,'cm'),
                            legend.key.width = unit(2.8,'cm'), 
                            legend.title = element_text(size=16) )+
                   labs(fill=wk_lab))

for(i in 1:(length(tmp_grid)/2) ){
  grid.arrange(arrangeGrob(wplot_list[[2*i-1]],wplot_list[[2*i]]))    
}

# Plot the mean delta
mean_delta = 0
max_delta = 0
maxh_delta = 0
sd_delta = 0
for(i in 1:length(proj_res)){
  mean_delta[i] = mean(proj_res[[i]]$delta)
  max_delta[i] = max(proj_res[[i]]$delta)
  sd_delta[i] = sd(proj_res[[i]]$delta)
  maxh_delta[i] = which.max(proj_res[[i]]$delta)
}

plot(mean_delta)
plot(max_delta)
plot(sd_delta)


#------------------------------------------------
# Shrinkage plots
xh = h[which(xs_test[h,1] < -100 & xs_test[h,1] > -100.75 & xs_test[h,2] < 45.75 & xs_test[h,2] > 44.75)]
shp_mat = matrix(0, nrow = length(tmp_grid), ncol = K)

prz = apply(proj_res[[1]]$ahat,1,function(x) sum(x>0))
prd = 0
#xh = which.max(prz)
xh = which.min(prz)
for(i in 1:length(tmp_grid)){
  shp_mat[i,] = t(proj_res[[i]]$ahat)[xh,]
  prd[i] = proj_res[[i]]$delta[xh]
}

plot(tmp_grid, shp_mat[,1], type = 'l', col = 'red', ylim = c(-0.1,1.1))
lines(tmp_grid, shp_mat[,2], col = 'blue')
lines(tmp_grid, shp_mat[,3], col = 'green')
lines(tmp_grid, shp_mat[,4], col = 'pink')
lines(tmp_grid, shp_mat[,5], col = 'purple')
lines(tmp_grid, shp_mat[,6], col = 'orange')
lines(tmp_grid, shp_mat[,7], col = 'cyan')
lines(tmp_grid, shp_mat[,8], col = 'black')

prd

saveRDS(proj_res, paste0(filedir,resdir,"world_3m_2014_n45k_sparsemax_mean/proj_res_april.rds"))

#------------------------------------------------


#------------------------------------------------
resid = fit$pred_mean - ms$y_test
sqrt(mean(resid^2))

# Resid Plot
rb = plot_residuals_hm_gg2(ms$x_test,resid,xcols = c(1,2), title="BMM Mean Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-10,0,10)) #c(-3.5,0,3.5)
rb = rb + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")



