#------------------------------------------------
# CMIP6 Run for Unity
# sudo openconnect --user=yannotty.1  vpn.asc.ohio-state.edu
# ssh yannotty.1@unity.asc.ohio-state.edu
#------------------------------------------------
#filedir = '/home/johnyannotty/Documents/CMIP6_mixing/Results/World/world_miroc_emulation_aug2014_ev_m300_0315/'
filedir = "/home/johnyannotty/Documents/CMIP6_mixing/Results/World/world_3m_2014_n45k_0315_pme/"
datadir = "/home/johnyannotty/Documents/CMIP6_mixing/Data/World/"
ms = readRDS(paste0(datadir,"ACC_BCC_MIROC_CMCC_CESM2_CNRM_CanESM5_KIOST_W_3M_2014_01_26_24_n45000.rds"))
dname = "MIROC_World_ev_Aug2014_03_10_24.rds"

setwd("/home/johnyannotty/Documents/openbt/src")
source("/home/johnyannotty/Documents/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")

# Save openbt results
#openbt_file = "World_2014_8models.obt"
openbt_file = "world8_0316.obt"
fitb = openbt.load(fname = paste0(filedir,openbt_file))

#Get mixed mean function
nt = 3 #length(unique(ms$x_test[,3]))
N = nrow(ms$x_test)/nt
x_test = ms$x_test[1:N,c(1:2,4)]

#Get mixed mean function
Nr = 40
h = sample(1:N,Nr)
fitp=predict.openbtmixing(fitb,x.test = as.matrix(x_test[h,]), 
                            f.test = as.matrix(rep(1,Nr)),
                            tc=4, q.lower = 0.025, q.upper = 0.975,
                            ptype = "mean")

Nlim = 5000

apply(fitp$mdraws[1:Nlim,],2,function(x) length(unique(x)))  

fitp$mmean ==  apply(fitp$mdraws + fitb$ymean,2,mean)  
fitp$m.5 ==  apply(fitp$mdraws + fitb$ymean,2,median)  
fitp$msd ==  apply(fitp$mdraws + fitb$ymean,2,sd)  
fitp$m.lower ==  apply(fitp$mdraws + fitb$ymean,2,quantile, 0.025)  
fitp$m.upper ==  apply(fitp$mdraws + fitb$ymean,2,quantile, 0.975)  



plot(fitp$mdraws[,10], type = "l")
which(fitp$mdraws[,1] == fitp$mdraws[4462,1])



fit$pred_mean[h] ==  apply(fitp$mdraws[1:Nlim,] + fitb$ymean,2,mean)  
fit$pred_lb[h] ==  apply(fitp$mdraws[1:Nlim,] + fitb$ymean,2,quantile, 0.025)  
fit$pred_ub[h] ==  apply(fitp$mdraws[1:Nlim,] + fitb$ymean,2,quantile, 0.975)  

#----------------------------------------------------------
# Test BMM - old results with new results
#----------------------------------------------------------
#Get mixed mean function
nt = 3 #length(unique(ms$x_test[,3]))
N = nrow(ms$x_test)/nt


#Get mixed mean function
# Helper to order result names
order_fnames = function(fls){
  fls_split = sapply(fls,function(x) strsplit(x,"batch"))
  fls_new = fls
  for(i in 1:length(fls_split)){
    if(length(fls_split[[i]])>1){
      if(nchar(fls_split[[i]][2])< 7){
        fls_new[i] = paste0(fls_split[[i]][1], "batch0",fls_split[[i]][2])
      }  
    }
  }
  fls = fls[order(fls_new)]
  return(fls)
}

filedir = "/home/johnyannotty/Documents/CMIP6_mixing/"
datadir = "Data/World/"
resdir = "Results/World/"

#ffold = "world_3m_2014_n45k_0315/"
ffold = "world_3m_2014_n45k_sparsemax_tmp13_0317_final/"
dataname = "ACC_BCC_MIROC_CMCC_CESM2_CNRM_CanESM5_KIOST_W_3M_2014_01_26_24_n45000.rds"
fls = system(paste0("ls ",filedir,resdir,ffold),intern = TRUE)
fls = order_fnames(fls)
ms = readRDS(paste0(filedir,datadir,dataname))
batch = 0
for(i in 1:length(fls)){
  if(grepl("batch",fls[i])){
    batch = batch + 1
    temp_fit = readRDS(paste0(filedir,resdir,ffold,fls[i]))
    if(batch == 1){
      fit1 = temp_fit
    }else{
      fit1 = list(
        pred_mean = c(fit1$pred_mean,temp_fit$pred_mean),
        pred_ub = c(fit1$pred_ub,temp_fit$pred_ub),
        pred_lb = c(fit1$pred_lb,temp_fit$pred_lb),
        wts_mean = rbind(fit1$wts_mean,temp_fit$wts_mean),
        wts_ub = rbind(fit1$wts_ub,temp_fit$wts_ub),
        wts_lb = rbind(fit1$wts_lb,temp_fit$wts_lb),
        wsum_mean = c(fit1$wsum_mean,temp_fit$wsum_mean),
        wsum_lb = c(fit1$wsum_lb,temp_fit$wsum_lb),
        wsum_ub = c(fit1$wsum_ub,temp_fit$wsum_ub),
        proj_mean = c(fit1$proj_mean,temp_fit$proj_mean),
        proj_ub = c(fit1$proj_ub,temp_fit$proj_ub),
        proj_lb = c(fit1$proj_lb,temp_fit$proj_lb),
        pwts_mean = rbind(fit1$pwts_mean,temp_fit$pwts_mean),
        pwts_ub = rbind(fit1$pwts_ub,temp_fit$pwts_ub),
        pwts_lb = rbind(fit1$pwts_lb,temp_fit$pwts_lb),
        delta_mean = c(fit1$delta_mean,temp_fit$delta_mean),
        delta_lb = c(fit1$delta_lb,temp_fit$delta_lb),
        delta_ub = c(fit1$delta_ub,temp_fit$delta_ub)
      )
    }
  }else if(grepl("sdraws",fls[i])){
    sfit1 = readRDS(paste0(filedir,resdir,ffold,fls[i]))
  }
}

# Number of models
K = ncol(ms$f_test)
N = nrow(ms$f_test)
Nr = 500
h = sample(1:N,Nr)
fitp=predict.openbtmixing(fitb,x.test = as.matrix(ms$x_test[h,]), 
                          f.test = ms$f_test[h,],
                          tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "mean_and_proj", temperature = temp_fit$tmp, proj_type = "euclidean")


sum(fitp$mmean ==  apply(fitp$mdraws,2,mean))  
sum(fitp$m.5 ==  apply(fitp$mdraws,2,median))
sum(fitp$msd ==  apply(fitp$mdraws,2,sd))  
sum(fitp$m.lower ==  apply(fitp$mdraws,2,quantile, 0.025)) 
sum(fitp$m.upper ==  apply(fitp$mdraws,2,quantile, 0.975))


for(j in 1:ncol(ms$f_test)){
  print(max(abs(fit1$wts_mean[h,j] - fitp$wmean[,j])))
  print(max(abs(fit1$wts_lb[h,j] - fitp$w.lower[,j])))
  print(max(abs(fit1$wts_ub[h,j] - fitp$w.upper[,j])))
}

for(j in 1:ncol(ms$f_test)){
  print(max(abs(fit1$pwts_mean[h,j] - fitp$pwmean[,j])))
  print(max(abs(fit1$pwts_lb[h,j] - fitp$pw.lower[,j])))
  print(max(abs(fit1$pwts_ub[h,j] - fitp$pw.upper[,j])))
}


max(abs(fitp$mmean - rowSums(fit1$wts_mean[h,]*ms$f_test[h,])))
max(abs(fitp$pmmean - rowSums(fit1$pwts_mean[h,]*ms$f_test[h,])))

fitp$wdraws[[2]][1:20,1]
fitp$wdraws[[2]][51:70,1]
fitp$wdraws[[2]][1:20,1]
fitp$wdraws[[2]][2220:2240,1]

hist(fitp$mdraws[,1])


#----------------------------------------------------------
# Compare fit1 and fit2
#----------------------------------------------------------
#ffold = "world_3m_2014_n45k_0315/"
ffold = "world_3m_2014_n45k_sparsemax_tmp13_0317_pme/"
dataname = "ACC_BCC_MIROC_CMCC_CESM2_CNRM_CanESM5_KIOST_W_3M_2014_01_26_24_n45000.rds"
fls = system(paste0("ls ",filedir,resdir,ffold),intern = TRUE)
fls = order_fnames(fls)
ms = readRDS(paste0(filedir,datadir,dataname))
batch = 0
for(i in 1:length(fls)){
  if(grepl("batch",fls[i])){
    batch = batch + 1
    temp_fit = readRDS(paste0(filedir,resdir,ffold,fls[i]))
    if(batch == 1){
      fit2 = temp_fit
    }else{
      fit2 = list(
        pred_mean = c(fit2$pred_mean,temp_fit$pred_mean),
        pred_ub = c(fit2$pred_ub,temp_fit$pred_ub),
        pred_lb = c(fit2$pred_lb,temp_fit$pred_lb),
        wts_mean = rbind(fit2$wts_mean,temp_fit$wts_mean),
        wts_ub = rbind(fit2$wts_ub,temp_fit$wts_ub),
        wts_lb = rbind(fit2$wts_lb,temp_fit$wts_lb),
        wsum_mean = c(fit2$wsum_mean,temp_fit$wsum_mean),
        wsum_lb = c(fit2$wsum_lb,temp_fit$wsum_lb),
        wsum_ub = c(fit2$wsum_ub,temp_fit$wsum_ub),
        proj_mean = c(fit2$proj_mean,temp_fit$proj_mean),
        proj_ub = c(fit2$proj_ub,temp_fit$proj_ub),
        proj_lb = c(fit2$proj_lb,temp_fit$proj_lb),
        pwts_mean = rbind(fit2$pwts_mean,temp_fit$pwts_mean),
        pwts_ub = rbind(fit2$pwts_ub,temp_fit$pwts_ub),
        pwts_lb = rbind(fit2$pwts_lb,temp_fit$pwts_lb),
        delta_mean = c(fit2$delta_mean,temp_fit$delta_mean),
        delta_lb = c(fit2$delta_lb,temp_fit$delta_lb),
        delta_ub = c(fit2$delta_ub,temp_fit$delta_ub)
      )
    }
  }else if(grepl("sdraws",fls[i])){
    sfit2 = readRDS(paste0(filedir,resdir,ffold,fls[i]))
  }
}


max(abs(fitp$mmean - rowSums(fit2$wts_mean[h,]*ms$f_test[h,])))
max(abs(fitp$pmmean - rowSums(fit2$pwts_mean[h,]*ms$f_test[h,])))

max(abs(fit1$pred_mean - rowSums(fit2$wts_mean*ms$f_test)))
max(abs(fit1$proj_mean - rowSums(fit2$pwts_mean*ms$f_test)))


