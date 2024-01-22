#------------------------------------------------
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
# Elevation
#------------------------------------------------
era5evdir = "/home/johnyannotty/NOAA_DATA/ERA5_Elevations/"
era5elevname = "era5_elevations_12_31_23.nc"

era5_elev = nc_open(paste0(era5evdir,era5elevname))
ev = ncvar_get(era5_elev,"elev")
ev_lon = ncvar_get(era5_elev,"lon")
ev_lat = ncvar_get(era5_elev,"lat")

dim(ev)

xgrid = cbind(lon = rep(ev_lon,each = length(ev_lat)),
              lat = rep(ev_lat,length(ev_lon)))

pera5 = plot_pred_2d_gg2(xgrid,as.vector(ev),title = "Elevation", 
                         scale_vals = c(-150,100,1000))

min(ev)
hist(ev)


era5evdir = "/home/johnyannotty/NOAA_DATA/ERA5_Elevations/"
era5elevname = "copernicus_altitude.nc"

era5_elev = nc_open(paste0(era5evdir,era5elevname))
str(era5_elev)
ev = ncvar_get(era5_elev,"ASurf")
ev_lon = ncvar_get(era5_elev,"lon")
ev_lat = ncvar_get(era5_elev,"lat")

length(ev_lat)
length(ev_lon)
dim(ev)
head(ev_lon)
head(ev_lat)

xgrid = cbind(lon = rep(ev_lon,each = length(ev_lat)),
              lat = rep(ev_lat,length(ev_lon)))

pera5 = plot_pred_2d_gg2(xgrid,as.vector(t(ev)),title = "Elevation", 
                         scale_vals = c(-150,100,5000))

min(ev)
hist(ev)





era5evdir = "/home/johnyannotty/NOAA_DATA/ERA5_Elevations/"
era5elevname = "bilinear_copernicus_altitude.nc"

era5_elev = nc_open(paste0(era5evdir,era5elevname))
str(era5_elev)
ev = ncvar_get(era5_elev,"elev")
ev_lon = ncvar_get(era5_elev,"lon")
ev_lat = ncvar_get(era5_elev,"lat")

length(ev_lat)
length(ev_lon)
dim(ev)
head(ev_lon)
head(ev_lat)

xgrid = cbind(lon = rep(ev_lon,each = length(ev_lat)),
              lat = rep(ev_lat,length(ev_lon)))

pera5 = plot_pred_2d_gg2(xgrid,as.vector(ev),title = "Elevation", 
                         scale_vals = c(-150,100,5000))

min(ev)
hist(ev)
#nc_close(era5_elev)
