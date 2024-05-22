#------------------------------------------------
# CMIP6 Run for Unity
# sudo openconnect --user=yannotty.1  vpn.asc.ohio-state.edu
# ssh yannotty.1@unity.asc.ohio-state.edu
#------------------------------------------------
filedir = '/home/johnyannotty/Documents/CMIP6_mixing/Results/World/world_miroc_emulation_aug2014_ev_m300_03_18/'
datadir = "/home/johnyannotty/Documents/CMIP6_mixing/Data/World/"
ms = readRDS(paste0(datadir,"ACC_BCC_MIROC_CMCC_CESM2_CNRM_CanESM5_KIOST_W_3M_2014_01_26_24_n45000.rds"))
dname = "MIROC_World_ev_Aug2014_03_10_24.rds"

setwd("/home/johnyannotty/Documents/openbt/src")
source("/home/johnyannotty/Documents/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
library(rBayesianOptimization)

# Save openbt results
#openbt_file = "World_2014_8models.obt"
openbt_file = "cmip6_save.obt"
fit = openbt.load(fname = paste0(filedir,openbt_file))

#Get mixed mean function
batchsz = 10000
nt = 3 #length(unique(ms$x_test[,3]))
N = nrow(ms$x_test)/nt
x_test = ms$x_test[1:N,c(1:2,4)]

#Get mixed mean function
batchsz = 10000
nb = ceiling(N/batchsz)
for(i in 25:nb){
  i = 24
  indst = (i-1)*batchsz + 1
  indend = min(i*batchsz,N)
  fitp=predict.openbtmixing(fit,x.test = as.matrix(x_test[indst:indend,]), 
                            f.test = as.matrix(rep(1,nrow(x_test[indst:indend,]))),
                            tc=4, q.lower = 0.025, q.upper = 0.975,
                            ptype = "mean")
  
  # Save results
  fit_data = list(
    pred_mean = fitp$mmean,
    pred_ub = fitp$m.upper,
    pred_lb = fitp$m.lower,
#    x_train = out_data$x_train,
#    y_train = out_data$y_train,
    m = fit$m,
    k = fit$k,
    shp1 = fit$rshp1,
    shp2 = fit$rshp2,
    minnodesz = fit$minnumbot,
    q = 4,
    base = fit$base,
    power = fit$power,
    maxd = fit$maxd,
    nu = fit$overallnu,
    lam = fit$overalllambda,
    sims = ms$sims,
    lat_bnds = ms$lat_bnds,
    lon_bnds = ms$lon_bnds
  )
  
  # Save results
  dt = gsub("/","_",format(Sys.time(), "%D_%H:%M"))
  sn = strsplit(dname,"_")[[1]][1]
  
  if(i<10){
    bnstr = paste0("0",i)
  }else{
    bnstr = paste0(i)
  }
  saveRDS(fit_data,paste0(filedir,sn,"_batch",bnstr,".rds"))
}

fitp=predict.openbtmixing(fit,x.test = x_test[1:10,], 
                          f.test = ms$f_test[1:10,],
                          tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "sigma")
saveRDS(fitp$sdraws[,1], paste0(filedir,sn,"_sdraws",".rds"))


hist(unlist(fitp$sdraws[,1]))
plot(unlist(fitp$sdraws[,1]))


