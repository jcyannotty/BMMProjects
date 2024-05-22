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
hgrid = seq(0.5,200,by = 1.5)
vgyhat = variog(coords = x_train, data = y_train ,uvec = hgrid)
plot(vgyhat$u,vgyhat$v)

# x bounds and h grid
xbnds = t(apply(ms$x_test,2,range))
xbnds = xbnds[1:2,]

param_list = list(k = c(1.0,1.25,1.5,1.75,2.0), a1 = c(2,10), a2 = c(5,10,20,40),power = c(0.5,1,2))
param_list = list(k = c(2.25), a1 = c(6), a2 = c(42),power = c(0.5))
param_grid = expand.grid(param_list)
vg_out = matrix(0,nrow = nrow(param_grid), ncol = length(hgrid))
sig2 = 1
for(j in 1:nrow(param_grid)){
  vg = variogram.openbtmixing(xbnds,hgrid,10000,1,
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
plot(hgrid,vg_out[1,]/2, type = "l", ylim = c(0,700))
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
  lam = lam*(ubnds[4] - lbnds[4]) + lbnds[4]
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
nu = 80
xbnds = t(apply(ms$x_test,2,range))
xbnds = xbnds[-3,]
#hgrid = seq(0.5,20,by = 0.75)
hgrid = seq(3,25,by = 2.5)
N = 10000
ymin = min(y_train)
ymax = max(y_train)
lbnds = c(0.5,2,2,0.1)
ubnds = c(4,10,80,1)

# Empirical
vgyhat = variog(coords = x_train, data = y_train, uvec = hgrid)
emp_vg = vgyhat$v #vgyhat$v[-1]
plot(vgyhat$u,vgyhat$v)

# Test the loss
#pvec0 = c(1.5,2,10,2,1)
pvec0 = rep(0.5,4)
vg_loss(pvec0)


# Run Optim
vgoptim = optim(pvec0, vg_loss, upper = rep(1,4), method = "L-BFGS-B",
                lower = rep(0,4), control = list(maxit = 100))
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
bayes_temp = BayesianOptimization(FUN = bayes_vg_loss, acq = "ei",
                                  bounds = par_bounds, init_grid_dt = init_grid_dt,
                                  init_points = pvec0, n_iter = 30)

tmp = bayes_temp$Best_Par # from bayes opt

#----------------------------------------------------------------
# Miroc
#----------------------------------------------------------------
cmfiledir = "/home/johnyannotty/NOAA_DATA/CMIP6/historical/"
sim = nc_open(paste0(cmfiledir,"tas_Amon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc"))

lon = ncvar_get(sim,"lon")
lat = ncvar_get(sim,"lat")
lon_merc = ifelse(lon>180,lon-360,lon)
tas = ncvar_get(sim, "tas")
sim_time_days = ncvar_get(sim,"time")

# Convert time units
tunits = ncatt_get(sim,"time","units")

tustr = strsplit(tunits$value, " ")
tdstr = strsplit(unlist(tustr)[3], "-")
tmonth = as.integer(unlist(tdstr)[2])
tday = as.integer(unlist(tdstr)[3])
tyear = as.integer(unlist(tdstr)[1])
sim_time = chron(sim_time_days,origin=c(tmonth, tday, tyear))

nlon = length(lon)
nlat = length(lat)
tas = ncvar_get(sim,"tas",start=c(1,1,1975),count = c(nlon,nlat,1))

tas_vec = as.vector(t(tas)) - 273.15
xgrid = cbind(lon = rep(lon_merc,each = nlat),lat = rep(lat,nlon))
dim(tas)

usa = map_data("world",region = "USA")
world = map_data("world")
states = map_data("state")

# Convert long to (-180,180) scale
xs_test = xgrid
xs_test[,"lon"] = ifelse(xs_test[,"lon"] > 180,xs_test[,"lon"]-360,xs_test[,"lon"])

# Setup the test grid
stepsz = 1
max_lon = 0; min_lon = -180; max_lat = 90; min_lat = 0
x_train_lon = lon_merc[which(lon_merc<=max_lon & lon_merc>=min_lon)]
x_train_lat = lat[which(lat<=max_lat & lat>=min_lat)]
hlon = seq(1,length(x_train_lon),by=stepsz)
hlat = seq(1,length(x_train_lat),by=stepsz)
x_train_grid = expand.grid(lon = x_train_lon[hlon], lat = x_train_lat[hlat])


# Colors
plot_mean2d_map_viridis(xgrid,tas_vec,xcols = c(1,2), 
                        viridis_opt = "viridis",
                        scale_limit = c(0,40), title = "ERA5",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(0,90),
                        lon_bnds = c(-180,0)
) 

# xx = cbind(lon = rep(seq(180,359.75,by = 0.25),each = 721),
#            lat = rep(seq(-90,89.75,by = 0.25),721))
# bytest = bilinear(lon,lat,tas-273.15,xx[,1],xx[,2])
# xx[,1] = ifelse(xx[,1]>180,xx[,1]-360,xx[,1])
# plot_mean2d_map_viridis(xx,bytest$z,xcols = c(1,2),
#                         viridis_opt = "viridis",
#                         scale_limit = c(-20,40), title = "ERA5",
#                         maps_list = list(data.frame(world),data.frame(states)),
#                         maps_cols = c("grey30","grey30"),
#                         lat_bnds = c(0,90),
#                         lon_bnds = c(-180,0)
# )



# Data generation
eo_split = TRUE
xh = which(xgrid[,1] < 0 & xgrid[,2] > 0)
if(eo_split){
  ev = xh[which(xh%%2 == 0)]  
  odd = xh[which(xh%%2 == 1)] 
  y_train = tas_vec[ev]
  x_train = xgrid[ev,]
  n_train = nrow(x_train)
  
  x_test = xgrid[odd,]
  y_test = tas_vec[odd]
  n_test = length(odd)
}else{
  x_test = as.matrix(xgrid[xh,])
  rownames(x_test) = NULL
  y_test = tas_vec[xh]
  n_test = length(xh)
  
  train_data = data.frame(x_test, y = y_test) %>% merge(x_train_grid)
  x_train = train_data[,1:2]
  y_train = train_data[,3]
  n_train = nrow(x_train)
}

# Empirical Variogram
hgrid = seq(0.5,400,by = 5)
vgyhat = variog(coords = x_train, data = y_train ,uvec = hgrid)
plot(vgyhat$u,vgyhat$v, ylim = c(0,400))

#vgm1 = likfit(coords = x_train, data = y_train, cov.model = "gaussian", ini.cov.pars = c(2,4))
#vgm1$parameters.summary

#df = data.frame(coords = x_train, data = y_train)
#df = as.geodata(df)

#vgenv = variog.model.env(geodata = df, obj.variog = vgyhat, nsim= 1000, model.pars = vgm1)
#plot(vgyhat, envelope.obj = vgenv)
#lines(hgrid,sapply(hgrid, function(x) 1000 - power_exp_kernel(x,0,sqrt(1000),550,1)))



# x bounds and h grid
xbnds = t(apply(x_train,2,range))
xbnds = xbnds[1:2,]

param_list = list(k = c(1.0,1.25,1.5,1.75,2.0), a1 = c(2,10), a2 = c(5,10,20,40),power = c(0.5,1,2))
param_list = list(k = c(1.7), a1 = c(2), a2 = c(65),power = c(1.75))
param_grid = expand.grid(param_list)
vg_out = matrix(0,nrow = nrow(param_grid), ncol = length(hgrid))
sig2 = 0.05
for(j in 1:nrow(param_grid)){
  vg = variogram.openbtmixing(xbnds,hgrid,10000,1,
                              k=param_grid[j,"k"],
                              0.95,
                              power = param_grid[j,"power"],
                              a1 = param_grid[j,"a1"],
                              a2 = param_grid[j,"a2"],
                              4,
                              ncut = 500,
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
plot(hgrid,vg_out[1,]/2, type = "l", ylim = c(0,800))
points(vgyhat$u,vgyhat$v)
abline(h = var(y_train), col = "grey")


#-------------------------------------------------
# Batch Fit
#-------------------------------------------------
# Train a bart model
q0 = 4
k = 1.4; a1 = 2; a2 = 100; pwr = 1.0
nu = 50
fit=train.openbtmixing(as.matrix(xgrid[xh,]),tas_vec[xh],as.matrix(rep(1,n_train+n_test)),pbd=c(1.0,0),ntree = 30,ntreeh=1,
                       numcut=300,tc=4,model="mixbart",modelname="miroc",
                       ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
                       power = pwr, base = 0.95, minnumbot = 2, overallsd = sqrt(0.05), k = k, overallnu = nu,
                       summarystats = FALSE, rpath = TRUE, q = q0, rshp1 = a1, rshp2 = a2,
                       stepwpert = 0.1, probchv = 0.1, batchsize = 5000, maxd = 100)

#Get mixed mean function
#x_test = x_train
#n_test = n_train
#y_test = y_train
fitp=predict.openbtmixing(fit,x.test = x_test, f.test = as.matrix(rep(1,n_test)),tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "mean_and_sigma", proj_type = "softmax", temperature = 0.2)

fitp=predict.openbtmixing(fit,x.test = xgrid[xh,], f.test = as.matrix(rep(1,n_test+n_train)),tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "mean_and_sigma", proj_type = "softmax", temperature = 0.2)


plot_mean2d_map_viridis(xgrid[xh,], fitp$mmean,xcols = c(1,2), 
                        viridis_opt = "viridis",
                        scale_limit = c(-20,40), title = "ERA5",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(0,90),
                        lon_bnds = c(-180,0)
)


plot_mean2d_map_viridis(xgrid, tas_vec,xcols = c(1,2), 
                        viridis_opt = "viridis",
                        scale_limit = c(-20,40), title = "ERA5",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(0,90),
                        lon_bnds = c(-180,0)
)

plot_mean2d_map_viridis(xgrid[xh,], fitp$mmean - tas_vec[xh],xcols = c(1,2), 
                        viridis_opt = "inferno",
                        scale_limit = c(-5,5), title = "ERA5",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(0,90),
                        lon_bnds = c(-180,0)
)


sqrt(mean((tas_vec[xh] - fitp$mmean)^2))
hist(fitp$sdraws[,1])
plot(fitp$sdraws[,1])


#------------------------------------------------
# GP Fit
#------------------------------------------------
df = data.frame(coords = x_train, data = y_train)
df = as.geodata(df)
hgrid = seq(0.5,250,by = 5)

vgh = variog(coords = x_train, data = y_train, uvec = hgrid)
plot(vgh)

geor_fit = variofit(vgh, cov.model = "gaussian", ini.cov.pars = c(400,70),fix.nugget = FALSE)#, kappa = 1.5) 
lines(geor_fit)

xtr_ind = expand.grid(1:n_train,1:n_train)
xts_ind = expand.grid(1:(n_train+n_test),1:(n_train+n_test))
xtrs_ind = expand.grid(1:n_train,1:(n_train+n_test))

if(geor_fit$cov.model == "gaussian"){
  gps2 = geor_fit$tausq
  gpsc = sqrt(geor_fit$cov.pars[1])
  gpls = geor_fit$cov.pars[2]^2/2
  
  Rp11 = apply(xtr_ind,1,function(x) sqr_exp_kernel(x_train[x[1],],x_train[x[2],],gpsc,gpls))
  Rp22 = apply(xts_ind,1,function(x) sqr_exp_kernel(xgrid[xh[x[1]],],xgrid[xh[x[2]],],gpsc,gpls))
  Rp12 = apply(xtrs_ind,1,function(x) sqr_exp_kernel(x_train[x[1],],xgrid[xh[x[2]],],gpsc,gpls))
}else{
  gps2 = geor_fit$tausq
  gpsc = sqrt(geor_fit$cov.pars[1])
  gpls = geor_fit$cov.pars[2]
  
  Rp11 = apply(xtr_ind,1,function(x) power_exp_kernel(x_train[x[1]],x_train[x[2]],gpsc,gpls,1))
  Rp22 = apply(xts_ind,1,function(x) power_exp_kernel(x_test[x[1]],x_test[x[2]],gpsc,gpls,1))
  Rp12 = apply(xtrs_ind,1,function(x) power_exp_kernel(x_train[x[1]],x_test[x[2]],gpsc,gpls,1))
}

Rp11 = matrix(Rp11, nrow = n_train, ncol = n_train, byrow = FALSE)
Rp11 = Rp11 + diag(gps2,n_train)

Rp22 = matrix(Rp22, nrow = (n_test+n_train), ncol = (n_test+n_train), byrow = FALSE)

Rp12 = matrix(Rp12, nrow = n_train, ncol = (n_test+n_train), byrow = FALSE)

mp1 = geor_fit$trend.matrix*geor_fit$beta
mp2 = rep(geor_fit$beta,n_test+n_train)
pdgp = predict_dist_gp(y_train,mp1,mp2,Rp11,Rp22,Rp12)
fhat = pdgp$mp
fshat = sqrt(diag(pdgp$Rp))
fhat_lb = fhat - 1.96*fshat
fhat_ub = fhat + 1.96*fshat



plot_mean2d_map_viridis(xgrid[xh,], fhat,xcols = c(1,2), 
                        viridis_opt = "viridis",
                        scale_limit = c(0,40), title = "ERA5",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(0,90),
                        lon_bnds = c(-180,0)
)

sqrt(mean((tas_vec[xh] - fhat)^2))
sqrt(mean((tas_vec[xh] - fitp$mmean)^2))

plot_mean2d_map_viridis(xgrid[xh,], fhat - tas_vec[xh],xcols = c(1,2), 
                        viridis_opt = "inferno",
                        scale_limit = c(-5,5), title = "ERA5",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(0,90),
                        lon_bnds = c(-180,0)
)
