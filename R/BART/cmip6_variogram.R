#------------------------------------------------
# CMIP6 BART Variogram 
#------------------------------------------------
library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)
library(plotly)
library(viridis)
library(latex2exp)
library(colorRamps)
library(geoR)

setwd("/home/johnyannotty/Documents/openbt/src")

source("/home/johnyannotty/Documents/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/kernels.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Plotting/computer_expt_plots.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/gp_utils.R")

library(rBayesianOptimization)

#------------------------------------------------
# Read in Data
#------------------------------------------------
filedir = "/home/johnyannotty/Documents/CMIP6_mixing/"
datadir = "Data/North_Hemisphere/"
dataname = "ACC_BCC_CESM2_CNRM_NH_6M14_01_06_24_n30000.rds"
ms = readRDS(paste0(filedir,datadir,dataname))

xh = which(ms$x_train$lon < 0 & ms$x_train$V3 == 1)
x_train = as.matrix(ms$x_train[xh,])
rownames(x_train) = NULL
y_train = ms$y_train[xh]
n_train = length(xh)

# Empirical Variogram
ugrid = seq(0.5,80, by = 1.5)
vgyhat = variog(coords = x_train, data = y_train ,uvec = ugrid)
plot(vgyhat$u,vgyhat$v)

# x bounds and h grid
xbnds = t(apply(ms$x_test,2,range))
xbnds = xbnds[-3,]
hgrid = seq(0.5,80,by = 1.5)

param_list = list(k = c(1.0,1.25,1.5,1.75,2.0), a1 = c(2,10), a2 = c(5,10,20,40),power = c(0.5,1,2))
param_grid = expand.grid(param_list)
vg_out = matrix(0,nrow = nrow(param_grid), ncol = length(hgrid))
sig2 = 1
for(j in 1:nrow(param_grid)){
  vg = variogram.openbtmixing(xbnds[1:3,],hgrid,10000,1,
                              k=param_grid[j,"k"],
                              0.95,
                              power = param_grid[j,"power"],
                              a1 = param_grid[j,"a1"],
                              a2 = param_grid[j,"a2"],
                              4,
                              ncut = 1000,
                              beta = 0,
                              sigma2 = sig2,
                              maxd = 999,
                              type = "b",
                              ymin = min(y_train),
                              ymax = max(y_train)
                            )
  cat("Progress: ", round(j/nrow(param_grid),4)*100)
  vg_out[j,] = vg$vmean
}

# Semi-variogram
plot(hgrid,vg_out[112,]/2, type = "l", ylim = c(0,700))
points(vgyhat$u,vgyhat$v)
abline(h = var(y_train), col = "grey")
#points(hgrid,sapply(h-grid,function(h) (1-power_exp_kernel(0,h,1,0.3,1))), pch = 3, col = 'red')

out = list(vg = vg_out, xbnds = xbnds, h = hgrid, pgrid = param_grid,
           sigma2 = 1,dataref = paste0(filedir,datadir,dataname),
           vgyhat = 2*vgyhat$v, vgyhat_uvec = vgyhat$u)
saveRDS(out, paste0(filedir,"Variograms/vgb_nwh_n2500.rds"))


#------------------------------------------------
# Roptim for Vg
#------------------------------------------------
# Define the vg loss
bayes_vg_loss = function(k,a1,a2,lam){
  #k = pvec[1]*(ubnds[1] - lbnds[1]) + lbnds[1]
  #a1 = pvec[2]*(ubnds[2] - lbnds[2]) + lbnds[2]
  #a2 = pvec[3]*(ubnds[3] - lbnds[3]) + lbnds[3]
  #pwr = pvec[4]*(ubnds[4] - lbnds[4]) + lbnds[4]
  #lam = pvec[5]*(ubnds[5] - lbnds[5]) + lbnds[5]
  
  k = k*(ubnds[1] - lbnds[1]) + lbnds[1]
  a1 = a1*(ubnds[2] - lbnds[2]) + lbnds[2]
  a2 = a2*(ubnds[3] - lbnds[3]) + lbnds[3]
  pwr = 1
  #pwr = pwr*(ubnds[4] - lbnds[4]) + lbnds[4]
  lam = lam*(ubnds[5] - lbnds[5]) + lbnds[5]
  sig2 = nu*(lam)/(nu + 2)
  
  vg = variogram.openbtmixing(xbnds,hgrid,N,1,
                              k=k,
                              0.95,
                              power = pwr,
                              a1 = a1,
                              a2 = a2,
                              4,
                              ncut = 1000,
                              beta = 0,
                              sigma2 = sig2,
                              maxd = 999,
                              type = "b",
                              ymin = ymin,
                              ymax = ymax
  )
  score = mean((vg$vmean/2-emp_vg)^2) 
  return(list(Score = score, Pred = 0))
}


vg_loss = function(pvec){
  k = pvec[1]*(ubnds[1] - lbnds[1]) + lbnds[1]
  a1 = pvec[2]*(ubnds[2] - lbnds[2]) + lbnds[2]
  a2 = pvec[3]*(ubnds[3] - lbnds[3]) + lbnds[3]
  #pwr = pvec[4]*(ubnds[4] - lbnds[4]) + lbnds[4]
  lam = pvec[4]*(ubnds[4] - lbnds[4]) + lbnds[4]
  pwr = 1
  sig2 = nu*(lam)/(nu + 2)
  
  vg = variogram.openbtmixing(xbnds,hgrid,N,1,
                              k=k,
                              0.95,
                              power = pwr,
                              a1 = a1,
                              a2 = a2,
                              4,
                              ncut = 1000,
                              beta = 0,
                              sigma2 = sig2,
                              maxd = 999,
                              type = "b",
                              ymin = ymin,
                              ymax = ymax
  )
  score = mean((vg$vmean/2-emp_vg)^2) 
  return(score)
}


# Global options
nu = 3
xbnds = t(apply(ms$x_test,2,range))
xbnds = xbnds[-3,]
hgrid = seq(0.5,20,by = 0.75)
N = 10000
ymin = min(y_train)
ymax = max(y_train)
lbnds = c(0.5,2,2,0.1)
ubnds = c(4,10,40,2)

# Empirical
vgyhat = variog(coords = x_train, data = y_train, uvec = hgrid)
emp_vg = vgyhat$v[-1]
plot(vgyhat$u,vgyhat$v)

# Test the loss
#pvec0 = c(1.5,2,10,2,1)
pvec0 = rep(0.5,4)
vg_loss(pvec0)


# Run Optim
vgoptim = optim(pvec0, vg_loss, upper = rep(1,4), method = "L-BFGS-B",
                lower = rep(0,4), control = list(maxit = 50))
vgoptim$par
res = 0 
for(j in 1:4){
  res[j] = vgoptim$par[j]*(ubnds[j]-lbnds[j])+lbnds[j] 
}


# Do bayes opt
par_bounds = list(k = c(0,1),a1 = c(0,1),a2 = c(0,1),
                  #pwr = c(0,1),
                  lam = c(0,1))
init_grid_dt = data.frame(k = seq(0.1,0.9, length = 10),
                          a1 = seq(0.1,0.9, length = 10),
                          a2 = seq(0.1,0.9, length = 10),
                          #pwr = seq(0.1,0.9, length = 10),
                          lam = seq(0.1,0.9, length = 10))
bayes_temp = BayesianOptimization(FUN = vg_loss, acq = "ei",
                                  bounds = par_bounds, init_grid_dt = init_grid_dt,
                                  init_points = pvec0, n_iter = 30)

tmp = bayes_temp$Best_Par # from bayes opt



