#----------------------------------------------------------
# Stacking methods for specifc heat example
#----------------------------------------------------------
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/BMM/fwls.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/BMM/hierarchical_stacking_stan.R")

library(loo)
library(rstan)

#----------------------------------------------------------
#----------------------------------------------------------
# Load Data 
filedir = "/home/johnyannotty/Documents/Dissertation/results/EFT/spheat_model_fit/"
standir = "/home/johnyannotty/Documents/BMMProjects/R/EFT/"
ms = readRDS(paste0(filedir,"ms_n20_01_23_24.rds"))
sgf = readRDS(paste0(filedir,"fsg_fit_n4_01_23_24.rds"))
lgf = readRDS(paste0(filedir,"flg_fit_n4_01_23_24.rds"))

f_train = cbind(fsg = ms$fsg_train, flg = ms$flg_train)
f_test = cbind(fsg = ms$fsg_test, flg = ms$flg_test)

n_train = nrow(f_train)
n_test = nrow(f_test)

sig2_hat1 = median(sgf$post[,3])
sig2_hat2 = median(lgf$post[,3])

# compute the log densities
lpd1 = log(dnorm(ms$y_train,f_train[,1],sd = sqrt(sig2_hat1)))
lpd2 = log(dnorm(ms$y_train,f_train[,2],sd = sqrt(sig2_hat2)))


#----------------------------------------------------------
# FWLS
#----------------------------------------------------------
basis = c("cubic","cubic")
K = ncol(f_train)
g_train = fwls_construct_basis(ms$g_train, K, basis)
g_test = fwls_construct_basis(ms$g_test, K, basis)

fit_fwls = fwls_cv(ms$y_train, f_train, g_train, lambda = seq(0.1,20,by = 0.1))
#fit_fwls = fwls(ms$y_train, ms$f_train, g_train, lambda = 3)
fitp_fwls = fwls_predict(fit_fwls,f_test,g_test)

# FWLS Residuals
fit_fwls$lambda
resid = fitp_fwls$fx - ms$f0_test
sqrt(mean(resid^2))

plot(ms$g_test, ms$f0_test, type = 'l')
lines(ms$g_test,fitp_fwls$fx, col = "purple")
lines(ms$g_test,fitp_fwls$fx_lb, col = "orange")
lines(ms$g_test,fitp_fwls$fx_ub, col = "orange")

plot(ms$g_test, fitp_fwls$wx[,1], col = "red", ylim = c(-2,2))
lines(ms$g_test, fitp_fwls$wx[,2], col = "blue")

fit_data = list(
  pred_mean = fitp_fwls$fx,
  pred_ub = fitp_fwls$fx_ub,
  pred_lb = fitp_fwls$fx_lb,
  wts_mean = fitp_fwls$wx,
  wts_ub = fitp_fwls$wx_ub,
  wts_lb = fitp_fwls$wx_lb,
  x_train = ms$x_train,
  y_train = ms$y_train,
  beta = fit_fwls$beta,
  beta_cov = fit_fwls$beta_cov,
  sigma2 = fit_fwls$sigma2,
  lambda = fit_fwls$lambda,
  basis = basis
)

# Save results
resname = "spheat_fwls_quad_02_08_24.rds"
saveRDS(fit_data,paste0(filedir,resname))


#----------------------------------------------------------
# HS
#----------------------------------------------------------
# Define the basis function
xtrain_basis = as.matrix(cbind(ms$g_train, ms$g_train^2))
xtest_basis = as.matrix(cbind(ms$g_test, ms$g_test^2))
head(xtrain_basis)
head(xtest_basis)

# Fit the HS Model using stan
xt = normx(xtrain_basis)
xts = normx(xtest_basis)
lpd_point = cbind(lpd1,lpd2)
stan_data = list(X = xt, N=n_train, d=ncol(xt), d_discrete=0,
                 lpd_point = lpd_point, K=ncol(lpd_point),tau_mu = 0.5,
                 tau_sigma = 1, tau_discrete = 0.5 , tau_con = 1)
fiths = stan(paste0(standir,"spheat_hs.stan"), data = stan_data, cores = 4, iter = 5000,
             warmup = 1000, chains = 1, pars = c("mu","mu_0","sigma","beta_con","sigma_con"),
             include = TRUE, algorithm = "NUTS", control = list(max_treedepth = 4))

# Posteriors
mu = extract(fiths)$mu
mu0 = extract(fiths)$mu_0
sigma = extract(fiths)$sigma
beta_con = extract(fiths)$beta_con # N x K-1 x p
sigma_con = extract(fiths)$sigma_con

dim(beta_con)

hist(mu0)
hist(mu)
hist(beta_con[,1,1]) 
plot(beta_con[,1,2], type = 'l')

# Predictions
res = predict_hs_stan(fiths,xts,f_test,tau_mu=0.5, tau_dis=0.5, tau_sig=1, qlower = 0.025, 
                      qupper = 0.975)

plot(ms$g_test, ms$f0_test, type = "l")
lines(ms$g_test, res$pred_mean, col = "purple")
lines(ms$g_test, res$pred_lb, col = "orange")
lines(ms$g_test, res$pred_ub, col = "orange")

plot(ms$g_test, res$wts_mean[,1], col = "red", ylim = c(-2,2))
lines(ms$g_test, res$wts_mean[,2], col = "blue")


# Save Results
fit_data = list(
  pred_mean = res$pred_mean,
  pred_ub = res$pred_ub,
  pred_lb = res$pred_lb,
  wts_mean = res$wts_mean,
  wts_ub = res$wts_ub,
  wts_lb = res$wts_lb,
  hs_beta = beta_con,
  hs_mu = mu,
  hs_mu0 = mu0,
  hs_sig_beta = sigma_con,
  hs_sig = sigma,
  lpd = lpd_point
)

# Save results
resname = "spheat_hs_quad_02_08_24.rds"
saveRDS(fit_data,paste0(filedir,resname))

#----------------------------------------------------------
# CS
#----------------------------------------------------------
csw = stacking_weights(lpd_point)
cs_pred = f_test[,1]*csw[1] + f_test[,2]*csw[2]
plot(ms$g_test, ms$f0_test, type = "l")
lines(ms$g_test, cs_pred, col = "purple")

fit_data = list(
  pred_mean = cs_pred,
  wts_mean = csw,
  lpd = lpd_point
)

resname = "spheat_cs_02_08_24.rds"
saveRDS(fit_data,paste0(filedir,resname))
