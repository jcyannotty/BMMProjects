#------------------------------------------------
# CMIP6 Emulation Data Laoder
#------------------------------------------------
#setwd("/home/yannotty.1/openbt/src")
# Load the R wrapper functions to the OpenBT library.
#source("/home/yannotty.1/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Samplers/sampler_functions.R")
library(ncdf4)
library(chron)
library(dplyr)
library(akima)

#------------------------------------------------
# Control
#------------------------------------------------
# Read in simulator data
dirname = "/home/johnyannotty/NOAA_DATA/CMIP6_Interpolations/"
era5dir = "/home/johnyannotty/NOAA_DATA/ERA5/era5_avg_mon_tas/"
era5evdir = "/home/johnyannotty/NOAA_DATA/ERA5_Elevations/"
era5elevname = "bilinear_copernicus_altitude.nc"
sim_list = c("tas_Amon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412.nc",
             "tas_Amon_BCC-CSM2-MR_historical_r1i1p1f1_gn_185001-201412.nc",
             "tas_Amon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc",
             "tas_Amon_CMCC-CM2-SR5_historical_r1i1p1f1_gn_185001-201412.nc",
             "tas_Amon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc",
             "tas_Amon_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_185001-201412.nc",
             "tas_Amon_CanESM5_historical_r7i1p2f1_gn_185001-201412.nc",
             "tas_Amon_KIOST-ESM_historical_r1i1p1f1_gr1_185001-201412.nc"
)
K = length(sim_list)

#----------------------------------------------------------------
# Pick a simulator, area of the globe, and time period
#----------------------------------------------------------------
cmfiledir = "/home/johnyannotty/NOAA_DATA/CMIP6/historical/"
sim_num = 3
sim = nc_open(paste0(cmfiledir,sim_list[sim_num]))

lon = ncvar_get(sim,"lon")
lat = ncvar_get(sim,"lat")
lon_merc = ifelse(lon>180,lon-360,lon)

# sim_time_days = ncvar_get(sim,"time")
# 
# # Convert time units
# tunits = ncatt_get(sim,"time","units")
# 
# tustr = strsplit(tunits$value, " ")
# tdstr = strsplit(unlist(tustr)[3], "-")
# tmonth = as.integer(unlist(tdstr)[2])
# tday = as.integer(unlist(tdstr)[3])
# tyear = as.integer(unlist(tdstr)[1])
# sim_time = chron(sim_time_days,origin=c(tmonth, tday, tyear))
# sim_time = format(sim_time,"%m_%Y")

# Pick month and year
date_index = expand.grid(mon = seq(1,12, by = 1),yr = seq(1850,2014, by = 1)) 
mon_list = c(08)
yr_list = c(2014)
mon_yr_list = expand.grid(mon = mon_list, yr = yr_list)
sim_mon_yr = c()
for(i in nrow(mon_yr_list)){
  sim_mon_yr[i] = which(date_index$mon == mon_yr_list$mon[i] & 
                          date_index$yr == mon_yr_list$yr[i])
}

nlon = length(lon)
nlat = length(lat)
tas = ncvar_get(sim,"tas",start=c(1,1,sim_mon_yr),count = c(nlon,nlat,1))

tas_vec = as.vector(t(tas)) - 273.15
xgrid = cbind(lon = rep(lon_merc,each = nlat),lat = rep(lat,nlon))

# Convert long to (-180,180) scale
xs_test = xgrid
xs_test[,"lon"] = ifelse(xs_test[,"lon"] > 180,xs_test[,"lon"]-360,xs_test[,"lon"])

# Setup the test grid
stepsz = 2
max_lon = 180; min_lon = -180; max_lat = 90; min_lat = -90
x_train_lon = lon_merc[which(lon_merc<=max_lon & lon_merc>=min_lon)]
x_train_lat = lat[which(lat<=max_lat & lat>=min_lat)]
hlon = seq(1,length(x_train_lon),by=stepsz)
hlat = seq(1,length(x_train_lat),by=stepsz)
x_train_grid = expand.grid(lon = x_train_lon[hlon], lat = x_train_lat[hlat])

# Data generation
eo_split = TRUE
xh = which(xgrid[,1] <= max_lon & xgrid[,1] >= min_lon & 
             xgrid[,2] <= max_lat & xgrid[,2] >= min_lat)
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

#------------------------------------------------
# Map to elevation
#------------------------------------------------
era5dir = "/home/johnyannotty/NOAA_DATA/ERA5/era5_avg_mon_tas/"
era5evdir = "/home/johnyannotty/NOAA_DATA/ERA5_Elevations/"
era5elevname = "bilinear_copernicus_altitude.nc"
add_elev = TRUE
if(add_elev){
  era5_elev = nc_open(paste0(era5evdir,era5elevname))
  ev = ncvar_get(era5_elev,"elev")
  ev_lon = ncvar_get(era5_elev,"lon")
  ev_lat = ncvar_get(era5_elev,"lat")
  
  xx = cbind(lon = rep(ev_lon,each = length(ev_lat)),
             lat = rep(ev_lat,length(ev_lon)))
  
  yy = as.vector(t(ev))
  xxe = cbind(xx, ev = yy)
  
  lmodd = which(1:nlon %% 2 ==1)
  leven = which(1:nlon %% 2 ==0)
  
  by_odd_xx = cbind(lon = rep(lon_merc[lmodd],each = nlat),lat = rep(lat,nlon/2))
  by_even_xx = cbind(lon = rep(lon_merc[leven],each = nlat),lat = rep(lat,nlon/2))
  
  byev_odd = bilinear(ev_lon,ev_lat,ev,by_odd_xx[,1],by_odd_xx[,2])
  byev_even = bilinear(ev_lon,ev_lat,ev,by_even_xx[,1],by_even_xx[,2])
  
  byxx = rbind(by_odd_xx,by_even_xx)
  byev = c(byev_odd$z,byev_even$z)
  byall = data.frame(byxx,ev = byev)
  
  x_test = data.frame(x_test) %>% left_join(data.frame(byall))
  x_test$ev = ifelse(x_test$ev == -999, 0,x_test$ev)
  
  x_train = data.frame(x_train) %>% left_join(data.frame(byall))
  x_train$ev = ifelse(x_train$ev == -999, 0,x_train$ev)
}

usa = map_data("world",region = "USA")
world = map_data("world")
states = map_data("state")


plot_mean2d_map_viridis(byxx,byev,xcols = c(1,2),
                        viridis_opt = "viridis",
                        scale_limit = c(0,1000), title = "ERA5",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(-90,90),
                        lon_bnds = c(-180,180)
)

plot_mean2d_map_viridis(xxe[,1:2],xxe[,3],xcols = c(1,2),
                        viridis_opt = "viridis",
                        scale_limit = c(0,1000), title = "ERA5",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(-90,90),
                        lon_bnds = c(-180,180)
)


#------------------------------------------------
# Write Data
#------------------------------------------------
# Save data
filedir = '/home/johnyannotty/Documents/CMIP6_mixing/Data/'
out_data = list(
  x_train = x_train,
  y_train = y_train,
  x_test = x_test,
  y_test = y_test,
  time_pd = mon_yr_list,
  lat_bnds = c(min_lat,max_lat),
  lon_bnds = c(min_lon,max_lon)
)

#dt = gsub("/","_",format(Sys.time(), "%D_%H:%M"))
dt = gsub("/","_",format(Sys.time(), "%D"))
sn = paste(unlist(sapply(sim_list[sim_num],
                         function(x) 
                           strsplit(strsplit(x,"_Amon_")[[1]][2],"_historical")[[1]][1])
), collapse = "_")

sn = gsub("ACCESS-CM2","ACC",sn)
sn = gsub("BCC-CSM2-MR","BCC",sn)
sn = gsub("MIROC-ES2L","MIROC",sn)
sn = gsub("CMCC-CM2-SR5","CMCC",sn)
sn = gsub("CNRM-CM6-1-HR","CNRM",sn)
sn = gsub("KIOST-ESM","KIOST",sn)

desc = "World_ev_Aug2014"
ffold = "World/Emulation/"
xx = system(paste0("ls ",filedir,ffold),intern = TRUE)
#fname = paste0(sn,"_",desc,"_",dt,".rds")
fname = paste0(sn,"_",desc,"_",dt,"_n",sum(n_train),".rds")

if(fname %in% xx){
  cat("------------------------ \n \tCheck Name \n------------------------")  
}else{
  saveRDS(out_data, paste0(filedir,ffold,fname))
  cat("------------------------ \n \tSaved Data \n------------------------")
}
