#------------------------------------------------
# Climate CMIP6 results loader - from Unity
#------------------------------------------------
setwd("/home/johnyannotty/Documents/openbt/src")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Plotting/computer_expt_plots.R")


library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)
library(plotly)
library(viridis)
library(latex2exp)

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


#------------------------------------------------
#------------------------------------------------
filedir = "/home/johnyannotty/Documents/CMIP6_mixing/"
datadir = "Data/World/Emulation/"
resdir = "Results/World/"

dataname = "CESM2_CNRM_SWUSA_GRID_EV_Dec_2014_12_03_23_n3000.rds"
resname = "CESM2_CNRM_SWUSA_GRID_EV_Dec_2014_12_12_23_n3000.rds"
sname = "CESM2_CNRM_SWUSA_GRID_EV_Dec_2014_sdraws_12_12_23_n3000.rds"

dataname = "CESM2_CNRM_NorthAmerica_EV_Dec_2014_12_05_23_n4000.rds"
resname = "CESM2_CNRM_NorthAmerica_EV_Dec_2014_12_05_23_n4000.rds"
sname = "CESM2_CNRM_NorthAmerica_EV_Dec_2014_sdraws_12_05_23_n4000.rds"

dataname = "MIROC_CESM2_CNRM_NorthAmerica_EV_Dec_2014_12_06_23_n4000.rds"
resname = "MIROC_CESM2_CNRM_NorthAmerica_EV_Dec_2014_12_06_23_n4000.rds"
sname = "MIROC_CESM2_CNRM_NorthAmerica_EV_Dec_2014_sdraws_12_06_23_n4000.rds"

#------------------------------------------------
# Batch Loader
#------------------------------------------------
#ffold = "na_Dec2014_n4000_proj1/"
#dataname = "CESM2_CNRM_NorthAmerica_EV_Dec_2014_12_05_23_n4000.rds"

dataname = "MIROC_World_ev_Aug2014_03_10_24.rds"
ffold = "world_miroc_emulation_aug2014/"

fls = system(paste0("ls ",filedir,resdir,ffold),intern = TRUE)
fls = order_fnames(fls)
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
        wsum_ub = c(fit$wsum_ub,temp_fit$wsum_ub),
        proj_mean = c(fit$proj_mean,temp_fit$proj_mean),
        proj_ub = c(fit$proj_ub,temp_fit$proj_ub),
        proj_lb = c(fit$proj_lb,temp_fit$proj_lb),
        pwts_mean = rbind(fit$pwts_mean,temp_fit$pwts_mean),
        pwts_ub = rbind(fit$pwts_ub,temp_fit$pwts_ub),
        pwts_lb = rbind(fit$pwts_lb,temp_fit$pwts_lb),
        delta_mean = c(fit$delta_mean,temp_fit$delta_mean),
        delta_lb = c(fit$delta_lb,temp_fit$delta_lb),
        delta_ub = c(fit$delta_lb,temp_fit$delta_ub)
        )
    }
  }else if(grepl("sdraws",fls[i])){
    sfit = readRDS(paste0(filedir,resdir,ffold,fls[i]))
  }
}


resid = fit$pred_mean - ms$f_test[,3]
sqrt(mean(resid^2))

#ms$y_test = c(ms$y_test,ms$y_train)
#ms$x_test = rbind(ms$x_test,ms$x_train)

#resid = ms$f_test[,2] - ms$y_test

#hj = which(ms$x_test[,3] == 2)
rb = plot_residuals_hm_gg2(ms$x_test,resid,xcols = c(1,2), title="BMM Mean Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-10,0,10)) #c(-3.5,0,3.5)
rb = rb + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")

# Simulator 1 residuals
resid1 = ms$f_test[,1] - ms$y_test
#resid1 = abs(ms$y_test/ms$f_test[,1]) 
r1 = plot_residuals_hm_gg2(ms$x_test,resid1,xcols = c(1,2), title="Sim1 Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-25,0,25)) #c(-3.5,0,3.5)
r1 = r1 + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")

# Simulator 2 residuals
resid2 = ms$f_test[,2] - ms$y_test
#resid2 = abs(ms$y_test/ms$f_test[,2])
r2 = plot_residuals_hm_gg2(ms$x_test,resid2,xcols = c(1,2), title="Sim2 Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-25,0,25)) #c(-3.5,0,3.5)
r2 = r2 + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")

# Simulator 3 residuals
resid3 = ms$f_test[,3] - ms$y_test
r3 = plot_residuals_hm_gg2(ms$x_test,resid3,xcols = c(1,2), title="Sim2 Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-25,0,25)) #c(-3.5,0,3.5)
r3 = r3 + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")


max(fit$wts_mean);min(fit$wts_mean)
w1 = plot_wts_2d_gg2(ms$x_test,fit$wts_mean,wnum = 1,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.8,0.5,1.5), title = "W1")


w2 = plot_wts_2d_gg2(ms$x_test,fit$wts_mean,wnum = 2,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.8,0.5,1.5), title = "W2")


w3 = plot_wts_2d_gg2(ms$x_test,fit$wts_mean,wnum = 3,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.8,0.5,1.5),title = "W3")


wsum = plot_wts_2d_gg2(ms$x_test,as.matrix(fit$wsum_mean),wnum = 1,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(0.5,1,1.5), title = "Wsum")

grid.arrange(w1,w2,w3,nrow = 1)
grid.arrange(r1,r2,r3,nrow = 1)
grid.arrange(w1,w2,w3,r1,r2,r3,nrow = 2)
grid.arrange(w1,w2,wsum,r1,r2,rb,nrow = 2)

# Error Standard Deviation
hist(unlist(sfit))
plot(unlist(sfit))

# Prediciton
usa = map_data("world",region = "USA")
world = map_data("world")
states = map_data("state")


# Convert long to (-180,180) scale
ntp = 3
tp_len = nrow(ms$f_test)/ntp
for(i in 1:ntp){
  assign(paste0("h",i),(tp_len*(i-1)+1):(tp_len*i))
}
hlist = list(h1,h2,h3)

xs_test = ms$x_test
xs_test[,"lon"] = ifelse(xs_test[,"lon"] > 180,xs_test[,"lon"]-360,xs_test[,"lon"])
xs_test[,1] = ifelse(xs_test[,1] < 0,xs_test[,1]-0.25,xs_test[,1])

plot_mean2d_map_viridis(xs_test[h1,], fit$pred_mean,xcols = c(1,2), 
                        viridis_opt = "viridis",
                        scale_limit = c(-30,40), title = "ERA5",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(-90,90),
                        lon_bnds = c(-180,180)
)

plot_mean2d_map_viridis(xs_test[h2,], ms$f_test[h2,3],xcols = c(1,2), 
                        viridis_opt = "viridis",
                        scale_limit = c(-30,40), title = "ERA5",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(-90,90),
                        lon_bnds = c(-180,180)
)

plot_mean2d_map_gradient(xs_test, ms$f_test[h2,3]-fit$pred_mean,xcols = c(1,2), 
                        scale_vals = c(-5,0,5), title = "ERA5",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(0,90),
                        lon_bnds = c(0,180)
)


p2 = plot_pred_2d_gg2(ms$x_test, ms$f_test[,2],title = "Simulator 2", 
                      scale_vals = c(-32,0,22) #scale_vals = c(-50,0,40)
) + labs(x = "Longitude", y = "Latitude", 
         fill = "Values")


# Interval width
pred_width = fit$pred_ub - fit$pred_lb
hist(pred_width)
pw1 = plot_residuals_hm_gg2(ms$x_test,pred_width,xcols = c(1,2), title="BMM Mean Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(0,10,20)) #c(-3.5,0,3.5)


#------------------------------------------------
# Weight Projections
#------------------------------------------------
pbmm = plot_pred_2d_gg2(ms$x_test, fit$proj_mean,title = "BMM", 
                        scale_vals = c(-45,0,30)) + 
  labs(x = "Longitude", y = "Latitude",fill = "Values")
         

dbmm = plot_pred_2d_gg2(ms$x_test, fit$delta_mean,title = "BMM", 
                        scale_vals = c(-10,0,10)) + 
  labs(x = "Longitude", y = "Latitude",fill = "Values")


w1 = plot_wts_2d_gg2(ms$x_test,fit$pwts_mean,wnum = 1,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(0,0.5,1), title = "W1")


w2 = plot_wts_2d_gg2(ms$x_test,fit$pwts_mean,wnum = 2,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(0,0.5,1), title = "W2")




#------------------------------------------------
# Multiple time periods
#------------------------------------------------
#h = which(ms$x_test[,3] == 6 & ms$x_test[,1] > (0))
resid = fit$pred_mean[h] - ms$y_test[h]
h = 1:(length(ms$y_test)/10)
r1 = plot_residuals_hm_gg2(ms$x_test[h,],resid,xcols = c(1,2), title="BMM Mean Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-6,0,6)) #c(-3.5,0,3.5)
r1 = r1 + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")


max(fit$wts_mean);min(fit$wts_mean)

w1 = plot_wts_2d_gg2(ms$x_test[h,],fit$wts_mean[h,],wnum = 1,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-1,0.5,2.5))

w2 = plot_wts_2d_gg2(ms$x_test[h,],fit$wts_mean[h,],wnum = 2,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-1,0.5,2.5))

w3 = plot_wts_2d_gg2(ms$x_test[h,],fit$wts_mean[h,],wnum = 3,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(0.5,0.5,1.5))

w4 = plot_wts_2d_gg2(ms$x_test[h,],fit$wts_mean[h,],wnum = 4,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(0.5,0.5,1.5))

grid.arrange(w1+theme(legend.position = "bottom"),
             w2+theme(legend.position = "bottom"),
             w3+theme(legend.position = "bottom"),
             w4+theme(legend.position = "bottom"),nrow = 1)

# Weight Slices
hw = which(ms$x_test[,3] == 6 & ms$x_test[,1] == 100)
plot(ms$x_test[hw,2],fit$wts_mean[hw,1], col = 'red',type = 'l', ylim = c(-0.8,2))
lines(ms$x_test[hw,2],fit$wts_mean[hw,2], col = 'blue')
lines(ms$x_test[hw,2],fit$wts_mean[hw,3], col = 'green3')
lines(ms$x_test[hw,2],fit$wts_mean[hw,4], col = 'orange')
lines(ms$x_test[hw,2],fit$wts_ub[hw,1], col = 'red',lty = 'dashed')
lines(ms$x_test[hw,2],fit$wts_ub[hw,2], col = 'blue',lty = 'dashed')
lines(ms$x_test[hw,2],fit$wts_ub[hw,3], col = 'green3',lty = 'dashed')
lines(ms$x_test[hw,2],fit$wts_ub[hw,4], col = 'orange',lty = 'dashed')
lines(ms$x_test[hw,2],fit$wts_lb[hw,1], col = 'red',lty = 'dashed')
lines(ms$x_test[hw,2],fit$wts_lb[hw,2], col = 'blue',lty = 'dashed')
lines(ms$x_test[hw,2],fit$wts_lb[hw,3], col = 'green3',lty = 'dashed')
lines(ms$x_test[hw,2],fit$wts_lb[hw,4], col = 'orange',lty = 'dashed')


# Predictions
pbmm = plot_pred_2d_gg2(ms$x_test[h,c(1,2)], fit$pred_mean[h],title = "BMM", 
                        scale_vals = c(-45,0,30) #scale_vals = c(-50,0,40)
) + labs(x = "Longitude", y = "Latitude", 
         fill = "Values")

pdag = plot_pred_2d_gg2(ms$x_test[h,c(1,2)], ms$y_test[h],title = "ERA5", 
                        scale_vals = c(-45,0,30) ,#scale_vals = c(-50,0,40)
) + labs(x = "Longitude", y = "Latitude", fill = "Values")


grid.arrange(pbmm,pdag,nrow = 1)

p1 = plot_pred_2d_gg2(ms$x_test[h,c(1,2)], ms$f_test[h,1],title = "Sim1", 
                        scale_vals = c(-45,0,30) ,#scale_vals = c(-50,0,40)
) + labs(x = "Longitude", y = "Latitude", fill = "Values")


# Interval width
pred_width = fit$pred_ub - fit$pred_lb
hist(pred_width)
pw1 = plot_residuals_hm_gg2(ms$x_test[h,],pred_width[h],xcols = c(1,2), 
                            title="BMM Mean Residuals", 
                            scale_colors = c("white","blue","midnightblue"),
                            scale_vals = c(0,3,6)) #c(-3.5,0,3.5)
hist(pred_width[h])

# Accurate regions
arh1 = which(abs(resid1)<abs(resid2) & abs(resid1)<abs(resid3))
arh2 = which(abs(resid2)<abs(resid1) & abs(resid2)<abs(resid3))

ar = rep(0.5,nrow(x_test))
ar[arh1] = 1 
ar[arh2] = 0 

ar = plot_residuals_hm_gg2(x_test, ar,xcols = c(1,2), 
                           title="Accurate Regions (Red = Access, Blue = BCC)", 
                           scale_colors = c("darkblue","green3","darkred"),
                           scale_vals = c(0,0.5,1)) #c(-3.5,0,3.5)
ar = ar + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")
ar = ar + theme(axis.text=element_text(size=12),axis.title=element_text(size=13), 
                plot.title = element_text(size = 14),
                legend.title = element_text(size = 16))
ar = ar + theme(legend.position = "none")

grid.arrange(w1,w2,w3,ar,nrow = 2)

#-----------------------------------------------------
# Precision Weights Empirical Variogram
#-----------------------------------------------------
# Precision Weight calculation
x_train = ms$x_train
f_train = ms$f_train
y_train = ms$y_train
n_train = length(y_train)

pw_se = apply(f_train,2,function(x) (y_train-x)^2)
pw_denom = rowSums(1/pw_se)
pw = (1/pw_se)/pw_denom
rowSums(pw)

# Get nearest neightbors
knn = 5
d = outer(x_train[,1],x_train[,1], "-")^2 + outer(x_train[,2],x_train[,2], "-")^2
nn = t(matrix(apply(d,1,function(x) order(x)[2:(knn+1)]),nrow = knn))
pw_nn = 0*pw

for(i in 1:n_train){
  ind = c(i,nn[i,])
  pw_se_nn = apply(pw_se[ind,],2,sum)
  pw_nn[i,] = (1/pw_se_nn)/sum(1/pw_se_nn)
}

rowSums(pw_nn)

#xdata = cbind(x1 = x_train[,1],x2 = rep(0,n_train))
xdata = cbind(x1 = x_train[,1],x2 = x_train[,2])
vghat = c() 
vghat_nn = c()
ugrid = seq(1,80,by = 2)
K = ncol(f_train)
for(j in 1:K){
  geor_vghat = variog(coords = xdata, data = pw[,j],uvec = ugrid)
  vgh = 2*geor_vghat$v
  vghat = cbind(vghat,vgh)
  
  geor_vghat = variog(coords = xdata, data = pw_nn[,j],uvec = ugrid)
  vgh = 2*geor_vghat$v
  vghat_nn = cbind(vghat_nn,vgh)
}

colnames(vghat) = paste0("vw",1:K)
colnames(vghat_nn) = paste0("vw",1:K)


# Plot the empirical variogram for the mean weight
vghat_mean = rowMeans(vghat)
plot(geor_vghat$u,vghat_mean)

plot(geor_vghat$u,vghat[,1])
abline(h = 2*var(pw_nn[,1]))

# Save results
outdir = '/home/johnyannotty/Documents/RandomPathBART/VariogramResults/Climate/'
out = list(vghat = vghat,vghat_nn = vghat_nn,hgrid = geor_vghat$u, 
           what = pw, what_nn = pw_nn, x = x_train,y = y_train)
saveRDS(out, paste0(outdir,"cmip6_northamerica500_emp_variogram_precwts_11_15_23.rds"))

#-----------------------------------------------------
# Theoretical Variogram
#-----------------------------------------------------
xbnds = matrix(c(ms$lon_bnds[1],ms$lon_bnds[2],
                 ms$lat_bnds[1],ms$lat_bnds[2]),nrow = 2,ncol = 2, byrow = TRUE)
h_grid = seq(0.5,60,by = 0.5)
vg = variogram.openbtmixing(xbnds,h_grid,10000,1,
                            k=1.1,
                            0.95,
                            power = 0.5,
                            a1 = 5,
                            a2 = 30,
                            4)

plot(h_grid,vg$vmean, type = "l", ylim = c(0,0.4))
points(geor_vghat$u,vghat_mean)

# Theoretical grid search
grid = expand.grid(a1=c(5,10),a2=c(10,15,20),beta=c(0.5,1),k=c(0.85,0.95,1,1.05))
vg_mat = matrix(0,nrow = nrow(grid), ncol = length(h_grid))
for(i in 1:nrow(grid)){
  vg = variogram.openbtmixing(xbnds, h_grid,10000,1,
                              k=grid[i,"k"],
                              0.95,
                              power = grid[i,"beta"],
                              a1 = grid[i,"a1"],
                              a2 = grid[i,"a2"],
                              4)
  vg_mat[i,] = vg$vmean
  cat("Progress: i = ", round(i/nrow(grid),4),"\r")
}

# Plot variogram
plot(h_grid,vg_mat[19,], type = "l", ylim = c(0,0.3))
lines(h_grid,vg_mat[20,], col = 'red')
lines(h_grid,vg_mat[21,], col = 'blue')
points(geor_vghat$u,vghat_mean)

plot(h_grid,vg_mat[1,], type = "l", ylim = c(0,0.3))
lines(h_grid,vg_mat[2,], col = 'red')
lines(h_grid,vg_mat[3,], col = 'blue')
points(geor_vghat$u,vghat_mean)

# Get the best fitting variogram
which.min(apply(vg_mat,1,function(x) sqrt(mean((x-vghat_mean)^2))))
order(apply(vg_mat,1,function(x) sqrt(mean((x-vghat_mean)^2))))

# Save results
out = list(vg = vg_mat,x = x_train,y = y_train,h = h_grid, pgrid = grid)
saveRDS(out, paste0(filedir,"cmip6_northamerica500_thr_variogram_wt_11_15_23.rds"))


#----------------------------------------------------------
# Extra Code
#----------------------------------------------------------
h = which(ms$x_test[,3] < 500 & ms$x_test[,3] > 100)
wtm = fit$wts_mean 
wtm[-h,] = -10  
w1 = plot_wts_2d_gg2(ms$x_test,wtm,wnum = 1,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.8,0.5,1.5), title = "W1")


w2 = plot_wts_2d_gg2(ms$x_test,wtm,wnum = 2,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.8,0.5,1.5), title = "W2")


w3 = plot_wts_2d_gg2(ms$x_test,wtm,wnum = 3,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.8,0.5,1.5),title = "W3")

