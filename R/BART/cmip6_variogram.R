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

source("/home/johnyannotty/Documents/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/kernels.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Plotting/computer_expt_plots.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/gp_utils.R")

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



