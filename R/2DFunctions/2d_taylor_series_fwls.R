#------------------------------------------------
# 2D Taylor Feature weighted linear stacking
# Frequentist version of model mixing
#------------------------------------------------
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/BMM/fwls.R")
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")


filedir = '/home/johnyannotty/Documents/Dissertation/results/2d_functions/'
ms = readRDS(paste0(filedir,"ms_taylor_sincos.rds"))
dfit = readRDS(paste0(filedir,"sincos_res_dpath_09_15_23.rds"))
rfit = readRDS(paste0(filedir,"sincos_res_rpath_09_15_23.rds"))

#------------------------------------------------
# Fit the model
#------------------------------------------------
basis = c("linear","linear")
K = ncol(ms$f_train)
g_train = fwls_construct_basis(ms$x_train, K, basis)
g_test = fwls_construct_basis(ms$x_test, K, basis)

fit_fwls = fwls_cv(ms$y_train, ms$f_train, g_train, lambda = seq(0.1,20,by = 0.1))
#fit_fwls = fwls(ms$y_train, ms$f_train, g_train, lambda = 3)
fitp_fwls = fwls_predict(fit_fwls,ms$f_test,g_test)

# FWLS Residuals
fit_fwls$lambda
resid = fitp_fwls$fx - ms$f0_test
sqrt(mean(resid^2))

rb = plot_residuals_hm_gg2(ms$x_test,resid,xcols = c(1,2), title="FWLS Mean Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-3,0,3)) #c(-3.5,0,3.5)

# Simulator 1 residuals
resid1 = ms$f_test[,1] - ms$f0_test
r1 = plot_residuals_hm_gg2(ms$x_test,resid1,xcols = c(1,2), title="Sim1 Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-20,0,10)) #c(-3.5,0,3.5)
r1 = r1 + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")

# Simulator 2 residuals
resid2 = ms$f_test[,2] - ms$f0_test
r2 = plot_residuals_hm_gg2(ms$x_test,resid2,xcols = c(1,2), title="Sim1 Residuals", 
                           scale_colors = c("darkblue","gray95","darkred"),
                           scale_vals = c(-20,0,10)) #c(-3.5,0,3.5)
r2 = r2 + labs(fill = bquote(hat(r)*"(x)"), x = "Longitude", y = "Latitude")


max(fitp_fwls$wx);min(fitp_fwls$wx)

w1 = plot_wts_2d_gg2(ms$x_test,fitp_fwls$wx,wnum = 1,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.1,0.8,1.6), title = "W1")

w2 = plot_wts_2d_gg2(ms$x_test,fitp_fwls$wx,wnum = 2,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(-0.1,0.8,1.6), title = "W2")

# Prediciton
p1 = plot_pred_2d_gg2(ms$x_test, fitp_fwls$fx,title = "BMM", 
                      scale_vals = c(-2.6,0,2.6) #scale_vals = c(-50,0,40)
) + labs(x = "Longitude", y = "Latitude", 
         fill = "Values")

pdag = plot_pred_2d_gg2(ms$x_test, ms$f0_test,title = "ERA5", 
                        scale_vals = c(-2.6,0,2.6) ,#scale_vals = c(-50,0,40)
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
resname = "sincos_fwls_linear_02_06_24.rds"
saveRDS(fit_data,paste0(filedir,resname))


#------------------------------------------------
# Cross Validation
#------------------------------------------------
