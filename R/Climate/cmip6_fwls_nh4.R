#------------------------------------------------
# 2D Taylor Feature weighted linear stacking
# Frequentist version of model mixing
#------------------------------------------------
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/BMM/fwls.R")
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")


filedir = '/home/johnyannotty/Documents/CMIP6_mixing/'
resdir = "Results/North_Hemisphere/FWLS/"
datadir = "Data/North_Hemisphere/"
dataname = "ACC_BCC_CESM2_CNRM_NH_6M_2014_03_05_24_n30000.rds"

ms = readRDS(paste0(filedir,datadir,dataname))

xs_train = ms$x_train
xs_train[,"lon"] = ifelse(xs_train[,"lon"] > 180,xs_train[,"lon"]-360,xs_train[,"lon"])
xs_train[,1] = ifelse(xs_train[,1] < 0,xs_train[,1]-0.25,xs_train[,1])

xs_test = ms$x_test
xs_test[,"lon"] = ifelse(xs_test[,"lon"] > 180,xs_test[,"lon"]-360,xs_test[,"lon"])
xs_test[,1] = ifelse(xs_test[,1] < 0,xs_test[,1]-0.25,xs_test[,1])

#------------------------------------------------
# Fit the model
#------------------------------------------------
K = ncol(ms$f_train)
basis = rep("linear",K)
g_train = fwls_construct_basis(xs_train, K, basis)
#g_test = fwls_construct_basis(ms$x_test, K, basis)

fit_fwls = fwls_cv(ms$y_train, ms$f_train, g_train, lambda = seq(0.01,2,by = 0.01))
#fit_fwls = fwls(ms$y_train, ms$f_train, g_train, lambda = 3)

batchsz = 2000
nbatches = ceiling(nrow(ms$f_test)/batchsz)
pred_mean = pred_lb = pred_ub = NULL
wts_mean = wts_lb = wts_ub = matrix(0, ncol = K, nrow = 0)
bn = 1
fit_fwls$lambda

# f_test = ms$f_test
# x_test = ms$x_test
# y_test = ms$y_test
# rm(ms)
for(i in 1:nbatches){
  sind = (i-1)*batchsz + 1
  eind = min(i*batchsz,nrow(ms$f_test))
  g_test_batch = fwls_construct_basis(xs_test[sind:eind,], K, basis)
  
  fitp_fwls = fwls_predict(fit_fwls,ms$f_test[sind:eind,],g_test_batch)
  pred_mean = c(pred_mean,fitp_fwls$fx)
  pred_lb = c(pred_lb,fitp_fwls$fx_lb)
  pred_ub = c(pred_ub,fitp_fwls$fx_ub)
  wts_mean = rbind(wts_mean,fitp_fwls$wx)
  wts_lb = rbind(wts_lb,fitp_fwls$wx_lb)
  wts_ub = rbind(wts_ub,fitp_fwls$wx_ub)
  
  cat("Progress: ", round(i/nbatches,4), "\r")
  
  if(i%%20 == 0 || i == nbatches){
    fit_data = list(
      pred_mean = pred_mean,
      pred_lb = pred_lb,
      pred_ub = pred_ub,
      wts_mean = wts_mean,
      wts_lb = wts_lb,
      wts_ub = wts_ub,
      sims = ms$sims,
      lat_bnds = ms$lat_bnds,
      lon_bnds = ms$lon_bnds,
      beta = fit_fwls$beta,
      sigma2 = fit_fwls$sigma2,
      lambda = fit_fwls$lambda,
      basis = basis
    )
    
    # Save results
    if(bn<10){bns = paste0("0",bn)}else{bns = paste(bn)} 
    resname = paste0("FWLS_linear/batch_",bns,".rds")
    bn = bn + 1
    saveRDS(fit_data,paste0(filedir,resdir,resname))
    pred_mean = pred_lb = pred_ub = NULL
    wts_mean = wts_lb = wts_ub = matrix(0, ncol = K, nrow = 0)
    rm(fit_data)
  }
} 

