#------------------------------------------------
# CMIP6 Data Plots
#------------------------------------------------
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Samplers/sampler_functions.R")
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
# Control
#------------------------------------------------
# Read in simulator data
dirname = "/home/johnyannotty/NOAA_DATA/CMIP6_Interpolations/"
era5dir = "/home/johnyannotty/NOAA_DATA/ERA5/era5_avg_mon_tas/"
sim_list = c("bilinear_tas_Amon_ACCESS-CM2_historical_r1i1p1f1_gn_2014.nc",
  "bilinear_tas_Amon_BCC-CSM2-MR_historical_r1i1p1f1_gn_2014.nc",
  "bilinear_tas_Amon_MIROC-ES2L_historical_r1i1p1f2_gn_2014.nc",
  "bilinear_tas_Amon_CMCC-CM2-SR5_historical_r1i1p1f1_gn_2014.nc"
  #"bilinear_tas_Amon_CESM2_historical_r1i1p1f1_gn_2014.nc"
)

K = length(sim_list)
sim_tas = list()
for(i in 1:length(sim_list)){
  temp_sim = nc_open(paste0(dirname,sim_list[i]))
  temp_tas = ncvar_get(temp_sim, "tas")
  sim_tas[[i]] = temp_tas
  #assign(paste0("sim",i), temp_sim)
  #assign(paste0("sim",i,"_tas"), temp_tas)
  rm(temp_tas)
  rm(temp_sim)
}

# Era5 data
era5 = nc_open(paste0(era5dir,"data_1990-2023.nc"))
era5_lon = ncvar_get(era5,"longitude")
era5_lat = ncvar_get(era5,"latitude")
era5_time_hrs = ncvar_get(era5,"time")


# Locatation and time selection
max_lon = 310; min_lon = 235
max_lat = 70; min_lat = 30
#mon_list = c("06","12")
mon_list = c("12")
yr_list = c("14")

# Convert time units
tunits = ncatt_get(era5,"time","units")

tustr = strsplit(tunits$value, " ")
tdstr = strsplit(unlist(tustr)[3], "-")
tmonth = as.integer(unlist(tdstr)[2])
tday = as.integer(unlist(tdstr)[3])
tyear = as.integer(unlist(tdstr)[1])
era5_time = chron(era5_time_hrs/24,origin=c(tmonth, tday, tyear))


# Create lon and lat grid
xgrid = cbind(lon = rep(era5_lon,each = length(era5_lat)),
              lat = rep(era5_lat,length(era5_lon)))

# Get grid subset
h = which(xgrid[,1] >= min_lon & xgrid[,1] < max_lon & xgrid[,2] >= min_lat & xgrid[,2] < max_lat)

# Split into
nlon = length(era5_lon)
nlat = length(era5_lat)
mon_yr_list = expand.grid(mon = mon_list, yr = yr_list)
era5_mon_yr = which(era5_time %in% paste0(mon_yr_list$mon,"/01/",mon_yr_list$yr))

# Long and lat indexes
stepsz = 1 # used for test set
x_test_lon = era5_lon[which(era5_lon<max_lon & era5_lon>=min_lon)]
x_test_lat = era5_lat[which(era5_lat<max_lat & era5_lat>=min_lat)]
hlon = seq(1,length(x_test_lon),by=stepsz)
hlat = seq(1,length(x_test_lat),by=stepsz)

# Define x indexes
x_test_ind = which(xgrid[,"lon"] %in% x_test_lon[hlon] & 
                     xgrid[,"lat"] %in% x_test_lat[hlat])

n_test = length(x_test_ind)

x_train = matrix(0,nrow = 0, ncol = 3)
x_test = matrix(0,nrow = 0, ncol = 3)
y_train = c()
y_test = c()
f_train = matrix(0,nrow = 0, ncol = length(sim_list))
f_test = matrix(0,nrow = 0, ncol = length(sim_list))

for(i in 1:length(era5_mon_yr)){
  mon = as.numeric(paste(mon_yr_list[i,"mon"]))
  yr = as.numeric(paste(mon_yr_list[i,"yr"]))
  era5_t2m = ncvar_get(era5,"t2m",start=c(1,1,1,era5_mon_yr[i]),
                       count = c(nlon,nlat,1,1))
  
  x_test = rbind(x_test,cbind(xgrid[x_test_ind,],rep(i,n_test)))
  f_test_temp = matrix(0,nrow = n_test, ncol = K)
  for(j in 1:K){
    f_test_temp[,j] = as.vector(t(sim_tas[[j]][,,mon]))[x_test_ind]
  }
  f_test = rbind(f_test, f_test_temp)
  y_test = c(y_test,as.vector(t(era5_t2m))[x_test_ind]-273.15)
}

# Remove the time index if needed
if(length(era5_mon_yr) == 1){
  x_train = x_train[,1:2]
  x_test = x_test[,1:2]
}

# Residual Plots
resacc = f_test[,1] - y_test
resbcc = f_test[,2] - y_test
resmir = f_test[,3] - y_test
rescmcc = f_test[,4] - y_test

racc = plot_residuals_hm_gg2(x_test,resacc,xcols = c(1,2), title="Access", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-25,0,25)) #c(-3.5,0,3.5)
racc = racc + labs(fill = bquote(hat(r)*"(x)"), x = "", y = "Latitude") +
  theme(plot.title = element_text(size = 16), axis.title = element_text(size = 14),
        plot.margin=unit(c(0.2,0,0,0.2), "cm"))


rbcc = plot_residuals_hm_gg2(x_test,resbcc,xcols = c(1,2), title="BCC", 
                             scale_colors = c("darkblue","gray95","darkred"),
                             scale_vals = c(-25,0,25)) #c(-3.5,0,3.5)
rbcc = rbcc + labs(fill = bquote(hat(r)*"(x)"), x = "", y = "") +
  theme(plot.title = element_text(size = 16), axis.title = element_text(size = 14),
        plot.margin=unit(c(0.2,0,0,0.2), "cm"))


rmir = plot_residuals_hm_gg2(x_test,resmir,xcols = c(1,2), title="Miroc", 
                             scale_colors = c("darkblue","gray95","darkred"),
                             scale_vals = c(-25,0,25)) #c(-3.5,0,3.5)
rmir = rmir + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(size = 16), axis.title = element_text(size = 14),
        plot.margin=unit(c(-0.2,0,0,0.2), "cm"))


rcmcc = plot_residuals_hm_gg2(x_test,rescmcc,xcols = c(1,2), title="CMCC", 
                             scale_colors = c("darkblue","gray95","darkred"),
                             scale_vals = c(-25,0,25)) #c(-3.5,0,3.5)
rcmcc = rcmcc + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "") +
  theme(plot.title = element_text(size = 16), axis.title = element_text(size = 14),
        plot.margin=unit(c(-0.2,0,0,0.2), "cm"))

# 2X2 Plot
r_leg = g_legend(racc + theme(legend.position = "bottom", legend.key.size = unit(1.7,'cm'),
                              legend.key.height = unit(0.7,"cm"), 
                              legend.title = element_text(size = 14)))
grid.arrange(arrangeGrob(racc+theme(legend.position = "none", legend.key.size = unit(1.0,'cm')),
                         rbcc+theme(legend.position = "none", legend.key.size = unit(1.0,'cm')),
                         rmir+theme(legend.position = "none", legend.key.size = unit(1.0,'cm')),
                         rcmcc+theme(legend.position = "none", legend.key.size = unit(1.0,'cm')),
                         nrow = 2), nrow=2, heights = c(12,2),r_leg)
