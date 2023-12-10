#------------------------------------------------
# CMIP6 Feature weighted linear stacking
# Frequentist version of model mixing
#------------------------------------------------
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/BMM/fwls.R")
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")


filedir = "/home/johnyannotty/Documents/CMIP6_mixing/"
datadir = "Data/"
resdir = "Results/SWUSA/FWLS/"
ms = readRDS(paste0(filedir,datadir,"CESM2_CNRM_SWUSA_GRID_EV_Dec_2014_12_03_23_n3000.rds"))
#ms = readRDS(paste0(filedir,datadir,"CESM2_CNRM_NorthAmerica_EV_Dec_2014_12_05_23_n4000.rds"))

#------------------------------------------------
# Fit the model
#------------------------------------------------
basis = c("quad","quad","quad")
K = ncol(ms$f_train)
g_train = fwls_construct_basis(ms$x_train, K, basis)
g_test = fwls_construct_basis(ms$x_test, K, basis)

fit_fwls = fwls(ms$y_train, ms$f_train, g_train, lambda = 3)
fitp_fwls = fwls_predict(fit_fwls,ms$f_test,g_test)

# FWLS Residuals
resid = fitp_fwls$fx - ms$y_test
sqrt(mean(resid^2))

rb = plot_residuals_hm_gg2(ms$x_test,resid,xcols = c(1,2), title="BMM Mean Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-20,0,10)) #c(-3.5,0,3.5)
rb = rb + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")

# Simulator 1 residuals
resid1 = ms$f_test[,1] - ms$y_test
r1 = plot_residuals_hm_gg2(ms$x_test,resid1,xcols = c(1,2), title="Sim1 Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-20,0,10)) #c(-3.5,0,3.5)
r1 = r1 + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")

# Simulator 2 residuals
resid2 = ms$f_test[,2] - ms$y_test
r2 = plot_residuals_hm_gg2(ms$x_test,resid2,xcols = c(1,2), title="Sim1 Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-20,0,10)) #c(-3.5,0,3.5)
r2 = r2 + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")


max(fitp_fwls$wx);min(fitp_fwls$wx)

w1 = plot_wts_2d_gg2(ms$x_test,fitp_fwls$wx,wnum = 1,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.1,0,2.5), title = "W1")



w2 = plot_wts_2d_gg2(ms$x_test,fitp_fwls$wx,wnum = 2,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.5,0,1.5), title = "W2")

# Prediciton
p1 = plot_pred_2d_gg2(ms$x_test, fitp_fwls$fx,title = "BMM", 
                      scale_vals = c(-22,0,22) #scale_vals = c(-50,0,40)
) + labs(x = "Longitude", y = "Latitude", 
         fill = "Values")

pdag = plot_pred_2d_gg2(ms$x_test, ms$y_test,title = "ERA5", 
                        scale_vals = c(-22,0,22) ,#scale_vals = c(-50,0,40)
) + labs(x = "Longitude", y = "Latitude", 
         fill = "Values")


#-------------------------------------------------
# Save results
#-------------------------------------------------
fit_data = list(
  pred_mean = fitp_fwls$fx,
  pred_ub = fitp_fwls$fx_ub,
  pred_lb = fitp_fwls$fx_lb,
  wts_mean = fitp_fwls$wx,
  wts_ub = fitp_fwls$wx_ub,
  wts_lb = fitp_fwls$wx_lb,
  x_train = ms$x_train,
  y_train = ms$y_train,
  sims = ms$sims,
  lat_bnds = ms$lat_bnds,
  lon_bnds = ms$lon_bnds,
  beta = fit_fwls$beta,
  beta_cov = fit_fwls$beta_cov,
  sigma2 = fit_fwls$sigma2,
  lambda = fit_fwls$lambda,
  basis = basis
)

# Save results
dt = gsub("/","_",format(Sys.time(), "%D"))
sn = paste(unlist(sapply(ms$sims,
                         function(x) strsplit(strsplit(x,"_Amon_")[[1]][2],"_historical")[[1]][1])
), collapse = "_")
sn = gsub("ACCESS-CM2","ACC",sn)
sn = gsub("BCC-CSM2-MR","BCC",sn)
sn = gsub("MIROC-ES2L","MIROC",sn)
sn = gsub("CMCC-CM2-SR5","CMCC",sn)
sn = gsub("CNRM-CM6-1-HR","CNRM",sn)
desc = "NorthAmerica_Dec_2014"

resname = paste0(sn,"_FWLS_quad_",desc,"_",dt,".rds")
saveRDS(fit_data,paste0(filedir,resdir,resname))


#------------------------------------------------
# Cross Validation
#------------------------------------------------
