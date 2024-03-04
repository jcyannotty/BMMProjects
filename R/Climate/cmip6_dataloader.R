#------------------------------------------------
# CMIP6 Data Laoder
#------------------------------------------------
#setwd("/home/yannotty.1/openbt/src")
# Load the R wrapper functions to the OpenBT library.
#source("/home/yannotty.1/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Samplers/sampler_functions.R")
library(ncdf4)
library(chron)
library(dplyr) 

#------------------------------------------------
# Control
#------------------------------------------------
# Read in simulator data
dirname = "/home/johnyannotty/NOAA_DATA/CMIP6_Interpolations/"
era5dir = "/home/johnyannotty/NOAA_DATA/ERA5/era5_avg_mon_tas/"
era5evdir = "/home/johnyannotty/NOAA_DATA/ERA5_Elevations/"
era5elevname = "bilinear_copernicus_altitude.nc"
sim_list = c("bilinear_tas_Amon_ACCESS-CM2_historical_r1i1p1f1_gn_2014.nc",
             "bilinear_tas_Amon_BCC-CSM2-MR_historical_r1i1p1f1_gn_2014.nc",
             "bilinear_tas_Amon_MIROC-ES2L_historical_r1i1p1f2_gn_2014.nc",
             "bilinear_tas_Amon_CMCC-CM2-SR5_historical_r1i1p1f1_gn_2014.nc",
             "bilinear_tas_Amon_CESM2_historical_r1i1p1f1_gn_2014.nc",
             "bilinear_tas_Amon_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_2014.nc",
             "bilinear_tas_Amon_CanESM5_historical_r7i1p2f1_gn_2014.nc",
             "bilinear_tas_Amon_KIOST-ESM_historical_r1i1p1f1_gr1_2014.nc"
             )
K = length(sim_list)

# Era5 data
era5 = nc_open(paste0(era5dir,"data_1990-2023.nc"))
era5_lon = ncvar_get(era5,"longitude")
era5_lat = ncvar_get(era5,"latitude")
era5_lon_merc = ifelse(era5_lon>180,era5_lon-360,era5_lon)
era5_time_hrs = ncvar_get(era5,"time")

# Data sizes
max_design = FALSE
test_grid = TRUE
test_equal_train = FALSE
stepsz = 2 # used for test set

# Locatation and time selection
min_lon = -180; max_lon = 180
min_lat = -90; max_lat = 90
#mon_list = c("01","02","03","04","05","06","07","08","09","10","11","12")
mon_list = c("04","08","12")
yr_list = c("14")
n_train = rep(15000,length(yr_list)*length(mon_list))

# Additional regions (lon_min, lon_max, lat_min, lat_max)
add_data = FALSE
nadd_train = 200
add_data_coords = c(min_lon, 255, min_lat, 55)

# Convert time units
tunits = ncatt_get(era5,"time","units")

tustr = strsplit(tunits$value, " ")
tdstr = strsplit(unlist(tustr)[3], "-")
tmonth = as.integer(unlist(tdstr)[2])
tday = as.integer(unlist(tdstr)[3])
tyear = as.integer(unlist(tdstr)[1])
era5_time = chron(era5_time_hrs/24,origin=c(tmonth, tday, tyear))


# Create lon and lat grid
xgrid = cbind(lon = rep(era5_lon_merc,each = length(era5_lat)),
              lat = rep(era5_lat,length(era5_lon)))

# Get grid subset
h = which(xgrid[,1] >= min_lon & xgrid[,1] < max_lon & xgrid[,2] >= min_lat & xgrid[,2] < max_lat)

# Split into
nlon = length(era5_lon)
nlat = length(era5_lat)
mon_yr_list = expand.grid(mon = mon_list, yr = yr_list)
era5_mon_yr = which(era5_time %in% paste0(mon_yr_list$mon,"/01/",mon_yr_list$yr))

# Setup the test grid
#x_test_lon = era5_lon[which(era5_lon<max_lon & era5_lon>=min_lon)]
x_test_lon = era5_lon_merc[which(era5_lon_merc<max_lon & era5_lon_merc>=min_lon)]
x_test_lat = era5_lat[which(era5_lat<max_lat & era5_lat>=min_lat)]
hlon = seq(1,length(x_test_lon),by=stepsz)
hlat = seq(1,length(x_test_lat),by=stepsz)

x_test_ind = which(xgrid[,"lon"] %in% x_test_lon[hlon] & 
                     xgrid[,"lat"] %in% x_test_lat[hlat])
n_test = length(x_test_ind)


# Build training and testing sets
x_train_ind = c()
for(i in 1:length(n_train)){
  if(max_design){
    xlist = list()
    #xlist[[1]] = sort(unique(xgrid[h,"lon"]))
    #xlist[[2]] = sort(unique(xgrid[h,"lat"]))
    xlist[[1]] = sort(unique(xgrid[x_test_ind,"lon"]))
    xlist[[2]] = sort(unique(xgrid[x_test_ind,"lat"]))
    x1diff = mean(diff(xlist[[1]]))
    x2diff = mean(diff(xlist[[1]]))
    
    #alist = c(min_lon, min_lat)
    #blist = c(max_lon, max_lat)
    Nlist = c(length(xlist[[1]]),length(xlist[[2]]))
    lon_lat_train = max_distance_design_unif_grid(Nlist,n_train[i],10)
    rownames(lon_lat_train) = NULL
    #x_train_ind_temp = cbind(xlist[[1]][lon_lat_train[,1],],xlist[[2]][lon_lat_train[,2],])
    #x_train_ind_temp = apply(lon_lat_train,1,function(x) which(x[1] == xgrid[,1] & x[2] == xgrid[,2]))
    
    # The most annoying way to get R to not repeat a few points
    x_train_ind_temp = lon_lat_train*0
    x_train_ind_temp[,1] = round((lon_lat_train[,1]-1)*x1diff + min(xlist[[1]]),2) 
    x_train_ind_temp[,2] = round((lon_lat_train[,2]-1)*x2diff + min(xlist[[2]]),2)
    
    if(nrow(unique(x_train_ind_temp)) < nrow(x_train_ind_temp)){stop("Uniqueness Error")}
    if(min(x_train_ind_temp[,1]) < min(xlist[[1]])){stop("Min Error - Lon")}
    if(max(x_train_ind_temp[,1]) > max(xlist[[1]])){stop("Max Error - Lon")}
    if(min(x_train_ind_temp[,2]) < min(xlist[[2]])){stop("Min Error - Lat")}
    if(max(x_train_ind_temp[,2]) > max(xlist[[2]])){stop("Max Error - Lat")}
    
    # Get the train indicies to be mapped to simulator output
    df1 = data.frame(xgrid)
    colnames(df1) = c("lon", "lat")
    df1[,"ind"] = 1:nrow(df1)
    
    df2 = data.frame(x_train_ind_temp)
    colnames(df2) = c("lon", "lat")
    
    df2 = df2 %>% left_join(df1, by = c("lon","lat"))
    
    x_train_ind = c(x_train_ind,df2$ind)
  }else{
    # Random design
    x_train_ind_temp = sample(x_test_ind, size = n_train[i])
    x_train_ind = c(x_train_ind,x_train_ind_temp)
  }
}

# Add to specific regions
if(add_data){
  adc = add_data_coords
  ha = which(xgrid[,1] >= adc[1] & xgrid[,1] < adc[2] & 
               xgrid[,2] >= adc[3] & xgrid[,2] < adc[4])
  ha = setdiff(ha, x_train_ind)
  xlist = list()
  xlist[[1]] = unique(xgrid[ha,"lon"])
  xlist[[2]] = sort(unique(xgrid[ha,"lat"]))
  alist = c(adc[1], adc[3])
  blist = c(adc[2], adc[4])
  lon_lat_train = max_distance_design(alist, blist, 0, nadd_train, 50, xgrid = xlist)
  x_train_ind = c(x_train_ind,
    apply(lon_lat_train,1,function(x) which(x[1] == xgrid[,1] & x[2] == xgrid[,2] ))
  )
}else{
  nadd_train = 0
}

if(test_equal_train){
  n_train = n_test
  x_train_ind = x_test_ind
}

# Get train and test sets
x_train = matrix(0,nrow = 0, ncol = 3)
x_test = matrix(0,nrow = 0, ncol = 3)
y_train = c()
y_test = c()
f_train = matrix(0,nrow = 0, ncol = length(sim_list))
f_test = matrix(0,nrow = 0, ncol = length(sim_list))

nstart = c(0,cumsum(n_train))+1
for(i in 1:length(era5_mon_yr)){
  mon = as.numeric(paste(mon_yr_list[i,"mon"]))
  yr = as.numeric(paste(mon_yr_list[i,"yr"]))
  era5_t2m = ncvar_get(era5,"t2m",start=c(1,1,1,era5_mon_yr[i]),
                       count = c(nlon,nlat,1,1))
  
  xtind = x_train_ind[nstart[i]:(nstart[i+1]-1)]
  x_train = rbind(x_train,cbind(xgrid[xtind,],rep(i,length(xtind))))
  y_train = c(y_train,as.vector(t(era5_t2m))[xtind]-273.15)
  f_train_temp = matrix(0,nrow = length(xtind), ncol = K)
  
  sim_tas = list()
  for(j in 1:length(sim_list)){
    simname = gsub("2014.nc",paste0("20",yr,".nc"),sim_list[j])
    temp_sim = nc_open(paste0(dirname,simname))
    temp_tas = ncvar_get(temp_sim, "tas")
    sim_tas[[j]] = temp_tas
    rm(temp_tas)
    rm(temp_sim)
  }
  
  for(j in 1:K){
    f_train_temp[,j] = as.vector(t(sim_tas[[j]][,,mon]))[xtind]
  }
  f_train = rbind(f_train, f_train_temp)
  
  x_test = rbind(x_test,cbind(xgrid[x_test_ind,],rep(i,n_test)))
  f_test_temp = matrix(0,nrow = n_test, ncol = K)
  for(j in 1:K){
    f_test_temp[,j] = as.vector(t(sim_tas[[j]][,,mon]))[x_test_ind]
  }
  f_test = rbind(f_test, f_test_temp)
  y_test = c(y_test,as.vector(t(era5_t2m))[x_test_ind]-273.15)
  cat("Progress: ", round(i/length(era5_mon_yr),4),"\r")
}


# Remove the time index if needed
if(length(era5_mon_yr) == 1){
  x_train = x_train[,1:2]
  x_test = x_test[,1:2]
}

# Sanity check
head(f_train)
head(f_test)
head(x_train)
tail(x_train)
head(x_test)
tail(x_test)

# Add Elevations
add_elev = TRUE
if(add_elev){
  era5_elev = nc_open(paste0(era5evdir,era5elevname))
  ev = ncvar_get(era5_elev,"elev")
  ev_lon = ncvar_get(era5_elev,"lon")
  ev_lat = ncvar_get(era5_elev,"lat")
  
  
  xx = cbind(lon = rep(ev_lon,each = length(ev_lat)),
                lat = rep(ev_lat,length(ev_lon)))
  
  #xx = expand.grid(lat = ev_lat,lon = ev_lon)
  #xx = cbind(lon = xx[,2],lat =  xx[,1])
  yy = as.vector(t(ev))
  xxe = cbind(xx, ev = yy)
  
  x_test = data.frame(x_test) %>% left_join(data.frame(xxe))
  x_test$ev = ifelse(x_test$ev == -999, 0,x_test$ev)
  
  x_train = data.frame(x_train) %>% left_join(data.frame(xxe))
  x_train$ev = ifelse(x_train$ev == -999, 0,x_train$ev)
}

# Dim Check
if(
  length(y_test) == length(era5_mon_yr)*n_test &
  nrow(x_test) == length(era5_mon_yr)*n_test & 
  nrow(f_test) == length(era5_mon_yr)*n_test &
  length(y_train) == (sum(n_train)+nadd_train) &
  nrow(x_train) == (sum(n_train)+nadd_train) &
  nrow(f_train) == (sum(n_train)+nadd_train) & 
  length(era5_mon_yr) == sum(unique(x_test[,3]) == unique(x_train[,3])) &
  sum(is.na(x_train$ev)) == 0 &
  sum(is.na(x_test$ev)) == 0
){
  cat("------------------------ \n \tPass \n------------------------")
}


# Remove large objects
rm(xgrid)
rm(sim_tas)
rm(era5_t2m)
nc_close(era5)

# Save data
filedir = '/home/johnyannotty/Documents/CMIP6_mixing/Data/'
out_data = list(
  x_train = x_train,
  y_train = y_train,
  f_train = f_train,
  x_test = x_test,
  y_test = y_test,
  f_test = f_test,
  time_pd = mon_yr_list,
  sims = sim_list,
  lat_bnds = c(min_lat,max_lat),
  lon_bnds = c(min_lon,max_lon)
)

#dt = gsub("/","_",format(Sys.time(), "%D_%H:%M"))
dt = gsub("/","_",format(Sys.time(), "%D"))
sn = paste(unlist(sapply(sim_list,
                         function(x) 
                           strsplit(strsplit(x,"_Amon_")[[1]][2],"_historical")[[1]][1])
            ), collapse = "_")

sn = gsub("ACCESS-CM2","ACC",sn)
sn = gsub("BCC-CSM2-MR","BCC",sn)
sn = gsub("MIROC-ES2L","MIROC",sn)
sn = gsub("CMCC-CM2-SR5","CMCC",sn)
sn = gsub("CNRM-CM6-1-HR","CNRM",sn)
sn = gsub("KIOST-ESM","KIOST",sn)

desc = "W_3M_2014"
ffold = "World/"
xx = system(paste0("ls ",filedir,ffold),intern = TRUE)
#fname = paste0(sn,"_",desc,"_",dt,".rds")
fname = paste0(sn,"_",desc,"_",dt,"_n",sum(n_train),".rds")

if(fname %in% xx){
  cat("------------------------ \n \tCheck Name \n------------------------")  
}else{
  saveRDS(out_data, paste0(filedir,ffold,fname))
  cat("------------------------ \n \tSaved Data \n------------------------")
}


#------------------------------------------------
# Write csv's - used for python
#------------------------------------------------
#csvfold = paste0(sn,"_",desc,"_",dt,"_n",sum(n_train),"_csvs")
csvfold = "ACC_BCC_MIROC_CMCC_CESM2_CNRM_CanESM5_KIOST_W_3M_2014_01_26_24_n45000_csvs/"
write.csv(out_data$f_train, paste0(filedir,ffold,csvfold,"f_train.csv"), row.names = FALSE)
write.csv(out_data$y_train, paste0(filedir,ffold,csvfold,"y_train.csv"), row.names = FALSE)
write.csv(out_data$x_train, paste0(filedir,ffold,csvfold,"x_train.csv"), row.names = FALSE)
write.csv(out_data$f_test, paste0(filedir,ffold,csvfold,"f_test.csv"), row.names = FALSE)
write.csv(out_data$y_test, paste0(filedir,ffold,csvfold,"y_test.csv"), row.names = FALSE)
write.csv(out_data$x_test, paste0(filedir,ffold,csvfold,"x_test.csv"), row.names = FALSE)




