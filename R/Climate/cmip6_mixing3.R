#------------------------------------------------
# Climate CMIP6 mixing - 3 simulators
# https://stackoverflow.com/questions/59628368/failed-installation-of-package-ncdf4-in-r-studio
#------------------------------------------------
setwd("/home/johnyannotty/Documents/openbt/src")

# Load the R wrapper functions to the OpenBT library.
source("/home/johnyannotty/Documents/openbt/src/openbt.R")
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
# Mixing 3 simulators
#------------------------------------------------
# Read in simulator data
dirname = "/home/johnyannotty/NOAA_DATA/CMIP6_Interpolations/"
era5dir = "/home/johnyannotty/NOAA_DATA/ERA5/era5_avg_mon_tas/"
sim1 = nc_open(paste0(dirname,"bilinear_tas_Amon_ACCESS-CM2_historical_r1i1p1f1_gn_2014.nc"))
sim2 = nc_open(paste0(dirname,"bilinear_tas_Amon_BCC-CSM2-MR_historical_r1i1p1f1_gn_2014.nc"))
sim3 = nc_open(paste0(dirname,"bilinear_tas_Amon_MIROC-ES2L_historical_r1i1p1f2_gn_2014.nc"))

era5 = nc_open(paste0(era5dir,"data_1990-2023.nc"))

era5_lon = ncvar_get(era5,"longitude")
era5_lat = ncvar_get(era5,"latitude")
era5_time_hrs = ncvar_get(era5,"time")

tunits = ncatt_get(era5,"time","units")

# Convert time units
tustr = strsplit(tunits$value, " ")
tdstr = strsplit(unlist(tustr)[3], "-")
tmonth = as.integer(unlist(tdstr)[2])
tday = as.integer(unlist(tdstr)[3])
tyear = as.integer(unlist(tdstr)[1])
era5_time = chron(era5_time_hrs/24,origin=c(tmonth, tday, tyear))

sim1_tas = ncvar_get(sim1,"tas")
sim2_tas = ncvar_get(sim2,"tas")
sim3_tas = ncvar_get(sim3,"tas")
dim(sim1_tas)

xgrid = cbind(lon = rep(era5_lon,each = length(era5_lat)),
              lat = rep(era5_lat,length(era5_lon)))

region = "north_america"
if(region == "north_america"){
  max_lon = 310; min_lon = 235
  max_lat = 70; min_lat = 30
}else if(region == "east_hemisphere"){
  max_lon = 135; min_lon = 5
  max_lat = 75; min_lat = 0
}

h = which(xgrid[,1] >= min_lon & xgrid[,1] < max_lon & xgrid[,2] >= min_lat & xgrid[,2] < max_lat)


#psim1 = plot_pred_2d_gg2(xgrid,as.vector(t(sim1_tas[,,12])),title = "True System", 
#                        scale_vals = c(-50,0,40))

psim1 = plot_pred_2d_gg2(xgrid[h,],as.vector(t(sim1_tas[,,6]))[h],title = "True System", 
                         scale_vals = c(-5,0,35))
psim2 = plot_pred_2d_gg2(xgrid[h,],as.vector(t(sim2_tas[,,6]))[h],title = "True System", 
                         scale_vals = c(-5,5,35))
psim3 = plot_pred_2d_gg2(xgrid[h,],as.vector(t(sim3_tas[,,6]))[h],title = "True System", 
                         scale_vals = c(-5,5,37))
rm(psim1);rm(psim2);rm(psim3)


# randomly split into test and train
nlon = length(era5_lon)
nlat = length(era5_lat)
dec14_ind = which(era5_time == "12/01/14")

era5_t2m = ncvar_get(era5,"t2m",start=c(1,1,1,dec14_ind),
                     count = c(nlon,nlat,1,1))
dim(era5_t2m)
pera5 = plot_pred_2d_gg2(xgrid,as.vector(t(era5_t2m))-272.15,title = "True System", 
                         scale_vals = c(-50,0,40))
rm(pera5)

x_train_ind = sample(h, size = 200)#size = 0.05*length(h))
x_train = xgrid[x_train_ind,]
f_train = cbind(sim1=as.vector(t(sim1_tas[,,12])),
                sim2=as.vector(t(sim2_tas[,,12])),
                sim3=as.vector(t(sim3_tas[,,12])))[x_train_ind,]
y_train = as.vector(t(era5_t2m))[x_train_ind]-272.15

test_grid = TRUE
stepsz = 7 # Take every 5 points
if(test_grid){
  x_test_lon = era5_lon[which(era5_lon<max_lon & era5_lon>=min_lon)]
  x_test_lat = era5_lat[which(era5_lat<max_lat & era5_lat>=min_lat)]
  hlon = seq(1,length(x_test_lon),by=stepsz)
  hlat = seq(1,length(x_test_lat),by=stepsz)
  x_test_ind = which(xgrid[,"lon"] %in% x_test_lon[hlon] & 
                       xgrid[,"lat"] %in% x_test_lat[hlat])
}else{
  x_test_ind = which(xgrid[,1] <= max(x_train[,1]) & xgrid[,1] >= min(x_train[,1]) &
                       xgrid[,2] <= max(x_train[,2]) & xgrid[,2] >= min(x_train[,2])
  )
  x_test_ind = setdiff(x_test_ind,x_train_ind)
  x_test_ind = sample(x_test_ind, size = 5000)
}

x_test = xgrid[x_test_ind,]
f_test = cbind(sim1=as.vector(t(sim1_tas[,,12])),
               sim2=as.vector(t(sim2_tas[,,12])),
               sim3=as.vector(t(sim3_tas[,,12])))[x_test_ind,]
y_test = as.vector(t(era5_t2m))[x_test_ind]-272.15

# EDA on residuals from each simulator
sim1_resid = sim1_tas[,,12] - (era5_t2m - 272.15) 
sim2_resid = sim2_tas[,,12] - (era5_t2m - 272.15)
sim3_resid = sim3_tas[,,12] - (era5_t2m - 272.15)


rsim1 = plot_residuals_hm_gg2(xgrid[h,],as.vector(t(sim1_resid))[h],title = "Simulator 1 Residuals", 
                              scale_vals = c(-20,0,20),scale_colors = c("darkblue","gray95","darkred"))
rsim2 = plot_residuals_hm_gg2(xgrid[h,],as.vector(t(sim2_resid))[h], title="Simulator 2 Residual", 
                              scale_colors = c("darkblue","gray95","darkred"),
                              scale_vals = c(-20.0,0,20))
rsim3 = plot_residuals_hm_gg2(xgrid[h,],as.vector(t(sim3_resid))[h], title="Simulator 3 Residual", 
                              scale_colors = c("darkblue","gray95","darkred"),
                              scale_vals = c(-20.0,0,20))

rm(xgrid)
rm(sim1_tas)
rm(sim2_tas)
rm(sim3_tas)
rm(era5_t2m)

rm(sim1_resid)
rm(sim2_resid)
rm(sim3_resid)
rm(rsim1)
rm(rsim2)
rm(rsim3)

nc_close(sim1)
nc_close(sim2)
nc_close(sim3)
nc_close(era5)

# Format the data in n_train x 2 matrix
nu = 10
rho = 1
sig2_hat = max(apply(apply(f_train, 2, function(x) (x-y_train)^2),2,min))
lam = rho*sig2_hat*(nu+2)/nu
q0 = 4
fit=openbt(x_train,y_train,f_train,pbd=c(1.0,0),ntree = 30,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="physics_model",
           ndpost = 8000, nskip = 1000, nadapt = 2000, adaptevery = 200, printevery = 500,
           power = 2.0, base = 0.95, minnumbot = 3, overallsd = sqrt(sig2_hat), k = 1.5, overallnu = nu,
           summarystats = FALSE, selectp = FALSE, rpath = TRUE, q = 4.0, rshp1 = 2, rshp2 = 10,
           stepwpert = 0.1, probchv = 0.1,maxd = 4)

scp = scanpost.openbtmixing(fit,2)
scp[[1]][[112]]

#Get mixed mean function
fitp = predict.openbt(fit,x.test = x_test, f.test = f_test,tc=4, q.lower = 0.025, q.upper = 0.975)
fitw = mixingwts.openbt(fit, x.test = x_test, numwts = ncol(f_test), tc = 4)

resid = fitp$mmean - y_test
resid_lab = bquote(hat(r)*"(x)"~"= "*hat("f")["\u2020"]*"(x) - f"["\u2020"]*"(x)")
r1 = plot_residuals_hm_gg2(x_test, resid,xcols = c(1,2), title=resid_lab, scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-20.0,0,20)) #c(-3.5,0,3.5)
r1 = r1 + labs(fill = bquote(hat(r)*"(x)"))
r1 = r1 + theme(axis.text=element_text(size=12),axis.title=element_text(size=13), plot.title = element_text(size = 17),
                legend.title = element_text(size = 16))

w1_lab = bquote(hat(W)[1]*"(x)")
w1x_lab = bquote(W[1]*"(x)")
w1 = plot_wts_2d_gg2(x_test,fitw$wmean,wnum = 1,xcols = c(1,2), title=w1_lab, scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.1,0.6,1.12))
w1 = w1 + labs(fill = w1x_lab)
w1 = w1 + theme(axis.text=element_text(size=12),axis.title=element_text(size=13), plot.title = element_text(size = 17),
                legend.title = element_text(size = 16))


#w2_lab = bquote(hat(W)[2]*"(x"[1]*", x"[2]*")")
w2_lab = bquote(hat(W)[2]*"(x)")
w2x_lab = bquote(W[2]*"(x)")
w2 = plot_wts_2d_gg2(x_test,fitw$wmean,wnum = 2,xcols = c(1,2), title=w2_lab, scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.1,0.6,1.12))
w2 = w2 + labs(fill = w2x_lab)
w2 = w2 + theme(axis.text=element_text(size=12),axis.title=element_text(size=13), plot.title = element_text(size = 17),
                legend.title = element_text(size = 16))

# Third weight
w3_lab = bquote(hat(W)[3]*"(x)")
w3x_lab = bquote(W[3]*"(x)")
w3 = plot_wts_2d_gg2(x_test,fitw$wmean,wnum = 3,xcols = c(1,2), title=w2_lab, scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.1,0.6,1.12))
w3 = w3 + labs(fill = w3x_lab)
w3 = w3 + theme(axis.text=element_text(size=12),axis.title=element_text(size=13), plot.title = element_text(size = 17),
                legend.title = element_text(size = 16))


# Get wsum
wsum = 0*fitw$wdraws[[1]]
for(i in 1:ncol(f_test)){
  wsum = wsum + fitw$wdraws[[i]]
}
wsum_mean = apply(wsum,2,mean)
wsum_lb = apply(wsum,2,quantile, 0.025)
wsum_ub = apply(wsum,2,quantile, 0.975)

sqrt(mean(resid^2))
median(abs(resid))
hist(wsum_mean)
hist(fitp$sdraws[,1])


filedir = '/home/johnyannotty/NOAA_DATA/cmip6_mixing_results/north_america_dec2014'
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

saveRDS(fit_data, paste0(filedir,"/acc_bcc_miroc_10_30_23.rds"))
saveRDS(fitp$sdraws[,1], paste0(filedir,"/acc_bcc_miroc_sig1_10_30_23.rds"))

rm(fit)
rm(fitp)
rm(fitw)
rm(wsum)

