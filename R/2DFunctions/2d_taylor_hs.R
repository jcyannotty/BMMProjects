#------------------------------------------------
# 2D Taylor Hierarhcical stacking
#------------------------------------------------
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Samplers/sampler_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Samplers/density_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/BMM/hierarchical_stacking_stan.R")

library(plotly)
library(viridis)
library(latex2exp)
library(rstan)
library(loo)

# Set directory
standir = '/home/johnyannotty/Documents/BMMProjects/R/2DFunctions/'
filedir = '/home/johnyannotty/Documents/Dissertation/results/2d_functions/'
ms = readRDS(paste0(filedir,"ms_taylor_sincos.rds"))

#-------------------------------------------------
# Fit each individual model (Closed form available)
# Ml: Y ~ N(fl(x),sig^2), sig2 ~ Lam*nu/Chi_nu^2
#-------------------------------------------------
nu = 10; lam = 2
s2 = rscinvchi2(10000,nu,lam)
hist(sqrt(s2), main = "Prior draws of sigma")
nupost = lampost = 0

n_train = nrow(ms$x_train)
K = ncol(ms$f_test)
r2matrix = matrix(0,nrow = n_train, ncol = K)
y_loopd = matrix(0,nrow = n_train, ncol = K)

for(l in 1:K){
  # Data info
  r2matrix[,l] = (ms$y_train - ms$f_train[,l])^2
  ssr = sum(r2matrix[,l])
  
  # Posterior info (stored for informative purposes)
  nupost[l] = n_train + nu
  lampost[l] = (nu*lam + ssr)/nupost[l]
  
  # LOO PD
  nuloo = nupost[l] - 1
  lamloo = (nu*lam + ssr-r2matrix[,l])/nuloo
  y_loopd[,l] = dtscaled(ms$y_train, df=nuloo, mean=ms$f_train[,l], scale = sqrt(lamloo)) 
}


#-------------------------------------------------
# Define the basis function
#-------------------------------------------------
xtrain_basis = as.matrix(cbind(ms$x_train, ms$x_train^2))
xtest_basis = as.matrix(cbind(ms$x_test, ms$x_test^2))
head(xtrain_basis)
head(xtest_basis)

#-------------------------------------------------
# Fit the HS Model using stan
#-------------------------------------------------
lpd_point = log(y_loopd)
xt = normx(xtrain_basis)
xts = normx(xtest_basis)
stan_data = list(X = xt, N=n_train, d=ncol(xt), d_discrete=0,
                 lpd_point = lpd_point, K=ncol(lpd_point),tau_mu = 0.5,
                 tau_sigma = 1, tau_discrete = 0.5 , tau_con = 1)
fiths = stan(paste0(standir,"2d_taylor_hs.stan"), data = stan_data, cores = 4, iter = 5000,
             warmup = 1000, chains = 1, pars = c("mu","mu_0","sigma","beta_con","sigma_con"),
             include = TRUE, algorithm = "NUTS", control = list(max_treedepth = 4))


#w_fit = extract(fiths, pars=w)$w

mu = extract(fiths)$mu
mu0 = extract(fiths)$mu_0
sigma = extract(fiths)$sigma
beta_con = extract(fiths)$beta_con # N x K-1 x p
sigma_con = extract(fiths)$sigma_con

dim(beta_con)

hist(mu)
hist(beta_con[,1,1]) 
hist(beta_con[,1,2])
hist(beta_con[,1,3])
plot(beta_con[,1,3], type = "l")

#w_fit = extract(fiths)$w
#res = fitted_hs_stan(w_fit,ms$f_train)
res = predict_hs_stan(fiths,xts,ms$f_test,tau_mu=0.5, tau_dis=0.5, tau_sig=1, qlower = 0.025, 
                      qupper = 0.975)


# FWLS Residuals
resid = res$pred_mean - ms$y_test
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


max(res$wts_mean);min(res$wts_mean)

w1 = plot_wts_2d_gg2(ms$x_test,res$wts_mean,wnum = 1,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(0,0.5,1), title = "W1")



w2 = plot_wts_2d_gg2(ms$x_test,res$wts_mean,wnum = 2,xcols = c(1,2), 
                     scale_colors = c("black","red2","yellow"),
                     scale_vals = c(0,0.5,1), title = "W2")

# Prediciton
p1 = plot_pred_2d_gg2(ms$x_test, res$pred_mean,title = "BMM", 
                      scale_vals = c(-30,0,22) #scale_vals = c(-50,0,40)
) + labs(x = "Longitude", y = "Latitude", 
         fill = "Values")

pdag = plot_pred_2d_gg2(ms$x_test, ms$y_test,title = "ERA5", 
                        scale_vals = c(-30,0,22) ,#scale_vals = c(-50,0,40)
) + labs(x = "Longitude", y = "Latitude", 
         fill = "Values")

#-------------------------------------------------
# Save results
#-------------------------------------------------
fit_data = list(
  pred_mean = res$pred_mean,
  pred_ub = res$pred_ub,
  pred_lb = res$pred_lb,
  wts_mean = res$wts_mean,
  wts_ub = res$wts_ub,
  wts_lb = res$wts_lb,
  x_train = ms$x_train,
  y_train = ms$y_train,
  sims = ms$sims,
  lat_bnds = ms$lat_bnds,
  lon_bnds = ms$lon_bnds,
  hs_beta = beta_con,
  hs_mu = mu,
  hs_mu0 = mu0,
  hs_sig_beta = sigma_con,
  hs_sig = sigma
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

resname = paste0(sn,"_HS_quad_",desc,"_",dt,".rds")
saveRDS(fit_data,paste0(filedir,resdir,resname))
