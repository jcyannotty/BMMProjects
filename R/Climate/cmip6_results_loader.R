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
datadir = "Data/"
resdir = "Results/SWUSA/"
#resname = "ACC_BCC_NorthAmerica_Dec_2014_11_08_23_n1000_m50_d4.rds"
#resname = "ACC_BCC_MIROC_NorthAmerica_Dec_2014_11_08_23_n1000_m50_d4.rds"
#dataname = "ACC_BCC_NorthAmerica_Dec_2014_11_07_23_n1000.rds"
dataname = "CMCC_CESM2_NorthAmerica_Dec_2014_11_30_23_n3000.rds"
resname = "CESM2_CNRM_SWUSA_GRID_EV_Dec_2014_12_03_23.rds"
sname = "CESM2_CNRM_SWUSA_GRID_EV_Dec_2014_12_03_23.rds"

dataname = "CMCC_CESM2_SWUSA_GRID_EV_Dec_2014_12_03_23_n3000.rds"
resname = "CMCC_CESM2_SWUSA_GRID_Dec_2014_12_02_23.rds"
sname = "CMCC_CESM2_SWUSA_GRID_Dec_2014_sdraws_12_02_23.rds"

dataname = "CESM2_CNRM_SWUSA_GRID_EV_Dec_2014_12_03_23_n3000.rds"
resname = "CESM2_CNRM_SWUSA_GRID_EV_Dec_2014_12_03_23.rds"
sname = "CESM2_CNRM_SWUSA_GRID_EV_Dec_2014_sdraws_12_03_23.rds"


ms = readRDS(paste0(filedir,datadir,dataname))
fit = readRDS(paste0(filedir,resdir,resname))
sfit = readRDS(paste0(filedir,resdir,sname))

resid = fit$pred_mean - ms$y_test
sqrt(mean(resid^2))

#resid = ms$f_test[,2] - ms$y_test

#hj = which(ms$x_test[,3] == 2)
rb = plot_residuals_hm_gg2(ms$x_test,resid,xcols = c(1,2), title="BMM Mean Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-10,0,10)) #c(-3.5,0,3.5)
rb = rb + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")

# Simulator 1 residuals
resid1 = ms$f_test[,1] - ms$y_test
r1 = plot_residuals_hm_gg2(ms$x_test,resid1,xcols = c(1,2), title="Sim1 Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-20,0,10)) #c(-3.5,0,3.5)
r1 = r1 + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")

# Simulator 2 residuals
resid2 = ms$f_test[,2] - ms$y_test
r2 = plot_residuals_hm_gg2(ms$x_test,resid2,xcols = c(1,2), title="Sim1 Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-20,0,10)) #c(-3.5,0,3.5)
r2 = r2 + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")


max(fit$wts_mean);min(fit$wts_mean)

w1 = plot_wts_2d_gg2(ms$x_test,fit$wts_mean,wnum = 1,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.1,0,2.5), title = "W1")



w2 = plot_wts_2d_gg2(ms$x_test,fit$wts_mean,wnum = 2,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.5,0,1.5), title = "W2")


w3 = plot_wts_2d_gg2(ms$x_test,fit$wts_mean,wnum = 3,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.5,0.75,2.0))

# Error Standard Deviation
hist(unlist(sfit))
plot(unlist(sfit))

# Prediciton
p1 = plot_pred_2d_gg2(ms$x_test, fit$pred_mean,title = "BMM", 
                      scale_vals = c(-22,0,22) #scale_vals = c(-50,0,40)
                      ) + labs(x = "Longitude", y = "Latitude", 
                                                          fill = "Values")

pdag = plot_pred_2d_gg2(ms$x_test, ms$y_test,title = "ERA5", 
                        scale_vals = c(-22,0,22) ,#scale_vals = c(-50,0,40)
                      ) + labs(x = "Longitude", y = "Latitude", 
                                                       fill = "Values")


# Interval width
pred_width = fit$pred_ub - fit$pred_lb
hist(pred_width)
pw1 = plot_residuals_hm_gg2(ms$x_test,pred_width,xcols = c(1,2), title="BMM Mean Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-0.5,15,30)) #c(-3.5,0,3.5)



# Multiple time periods
h = which(ms$x_test[,3] == 2)
r1 = plot_residuals_hm_gg2(ms$x_test[h,],resid[h],xcols = c(1,2), title="BMM Mean Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-25,0,25)) #c(-3.5,0,3.5)
r1 = r1 + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")


max(fit$wts_mean);min(fit$wts_mean)

w1 = plot_wts_2d_gg2(ms$x_test[h,],fit$wts_mean[h,],wnum = 1,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.1,0.5,1.1))


w2 = plot_wts_2d_gg2(ms$x_test[h,],fit$wts_mean[h,],wnum = 2,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.1,0.5,1.1))


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
