#------------------------------------------------
# Generating data from Trees at 1 am
#------------------------------------------------
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/TreeModels/tree_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/TreeModels/rpath_cvg.R")
source("/home/johnyannotty/Documents/RandomPathBART/rpath_prior_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Plotting/computer_expt_plots.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/kernels.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/gp_utils.R")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")

setwd("/home/johnyannotty/Documents/openbt/temp/openbt-current-rpath")
setwd("/home/johnyannotty/Documents/")

library(geoR)

#------------------------------------------------
# Conditional on 1 trees and gamma
#------------------------------------------------
# Set x
n_train = 200
n_test = 300
xmin = -1; xmax = 1
ymin = -5; ymax = 5
kappa = 1
#s = 0.05
tau2 = ((ymax - ymin)/(2*kappa))^2
x_train = seq(xmin+0.05,xmax-0.05,length = n_train) 
x_test = seq(xmin+0.01,xmax-0.01,length = n_test) 
x_all = c(x_train, x_test)

alpha = 0.95
beta = 2.0
#tau2 = 1
gam = 0.5
xbnds = matrix(c(xmin,xmax),nrow = 1)
hgrid = seq(0, 1.5, length = 30)

xcuts = create_cuts(xmin,xmax,200)

tree=init_tree()
bout1 = birth(tree$nid1,1,0.15)

tree = append_tree(tree, bout1)
bout2 = birth(tree$nid2,1,-0.4)

tree = append_tree(tree, bout2)
bout3 = birth(tree$nid3,1,0.55)

tree = append_tree(tree, bout3)
bout4 = birth(tree$nid4,1,-0.65)

tree = append_tree(tree, bout4)
bout8 = birth(tree$nid8,1,-0.80)

tree = append_tree(tree, bout8)
bout7 = birth(tree$nid7,1,0.8)

tree = append_tree(tree, bout7)
bout15 = birth(tree$nid15,1,0.9)

tree = append_tree(tree, bout15)
names(tree)

bnv = getbots(tree)


# Compute the covariance function
# igrid = expand.grid(i1 = 1:n_train, i2 = 1:n_train)
# phix = phix_tree(tree,x_train,xbnds,gam)
# rpc = apply(igrid,1,function(x) rpath_cov(phix[x[1],],phix[x[2],],tau2)) 
# rpc = matrix(rpc, nrow = n_train, ncol = n_train, byrow = FALSE)
# diag(rpc) = tau2 
# dim(rpc)

# Compute the covariance function
n_all = n_train + n_test
igrid_all = expand.grid(i1 = 1:n_all, i2 = 1:n_all)
phix_all = phix_tree(tree,x_all,xbnds,gam)
phix_test = phix_all[(n_train+1):(n_train+n_test),]

rpc = apply(igrid_all,1,function(x) rpath_cov(phix_all[x[1],],phix_all[x[2],],tau2)) 
rpc = matrix(rpc, nrow = n_all, ncol = n_all, byrow = FALSE)
diag(rpc) = tau2 

plot_phix(x_test, phix_all[(n_train+1):(n_train+n_test),])
plot_kernel_viridis(rpc[(n_train+1):n_all,(n_train+1):n_all]/tau2)


#------------------------------------------------
# Generate Data
#------------------------------------------------
sig = 0.005 
R11 = rpc[1:n_train,1:n_train]  + diag(sig^2,n_train)
R22 = rpc[(n_train+1):(n_train+n_test),(n_train+1):(n_train+n_test)]
R12 = rpc[1:n_train,(n_train+1):(n_train+n_test)]

plot_kernel_viridis(R22/tau2)

m1 = rep(0,n_train)
m2 = rep(0,n_test)

set.seed(311)
y_train = sample_gp(m1,R11)
plot(x_train, y_train)

# True Predictive Distribution 
pdgp = predict_dist_gp(y_train,m1,m2,R11,R22,R12)
f0_test = pdgp$mp

plot(x_test, f0_test, type = "l")
points(x_train, y_train)


#------------------------------------------------
# Conditional Variogram Time
#------------------------------------------------
# Empirical vg
ugrid = seq(0.02,1.5, length = 40)
vgyhat = variog(coords = cbind(x_train,rep(0,n_train)), data = y_train, uvec = ugrid)

# Compute the variogram
hgrid = seq(min(ugrid), max(ugrid), length = 120)
vg_means = matrix(0, nrow = 0, ncol = length(hgrid))
vg_params = matrix(0, nrow = 0, ncol = 2)

gam_list = c(0.1,0.5)
for(j in 1:length(gam_list)){
  k = 1; vg_gam = gam_list[j]
  N = 5000
  vg_tau2 = (diff(range(y_train))/(2*k))^2
  svg = svg_tree(tree,vg_gam,vg_tau2,xbnds,hgrid,N,0)
  
  vg_means = rbind(vg_means,svg$vgfhat)
  vg_params = rbind(vg_params,c(k,vg_gam))
  
  klist = c(1.1,1.25,1.5,1.6,1.7,1.75,2.0,2.25)
  for(k0 in klist){
    vg_means = rbind(vg_means,svg$vgfhat*(1/k0^2))
    vg_params = rbind(vg_params, c(k = k0, gam = vg_gam))
  }
  print("Round over....")
}

plot(vgyhat$u,vgyhat$v, ylim = c(0,50))
lines(hgrid, vg_means[27,])
abline(h = var(y_train))

vg_params

#-------------------------------------------------
# Batch Fit
#-------------------------------------------------
nu = 30
q0 = 4

# Train a bart model
fit=train.openbtmixing(x_train,y_train,as.matrix(rep(1,n_train)),pbd=c(1.0,0),ntree = 1,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="physics_model",
                       ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
                       power = 1.0, base = 0.95, minnumbot = 5, overallsd = 0.2, k = 2.0, overallnu = nu,
                       summarystats = FALSE, rpath = TRUE, q = q0, rshp1 = 15, rshp2 = 135,
                       stepwpert = 0.1, probchv = 0.1, batchsize = 10000)


#Get mixed mean function
fitp=predict.openbtmixing(fit,x.test = x_test, f.test = as.matrix(rep(1,n_test)), q.lower = 0.025, q.upper = 0.975,
                          ptype = "mean_and_sigma", tc=4)

# Sigma
hist(fitp$sdraws[,1])
plot(fitp$sdraws[,1], type = 'l')

# compare vs true curve
sqrt(mean((fitp$mmean - f0_test)^2))

plot_mean1d(x_test, pmean = fitp$mmean, plb = fitp$m.lower, pub = fitp$m.upper,
            amean = c(f0_test), colors = c("black","purple2"),
            apts_x = x_train, apts_y = y_train)


#------------------------------------------------
# Save Results
#------------------------------------------------
fit_data = list(
  pred_mean = fitp$mmean,
  pred_ub = fitp$m.upper,
  pred_lb = fitp$m.lower,
  m = fit$m,
  k = fit$k,
  shp1 = fit$rshp1,
  shp2 = fit$rshp2,
  minnodesz = fit$minnumbot,
  q = q0,
  base = fit$base,
  power = fit$power,
  maxd = fit$maxd,
  nu = fit$overallnu,
  lam = fit$overalllambda,
  sigma = fitp$smean[1],
  x_train = x_train,
  y_train = y_train
)

# Kernel
out = list(tree = tree, gam = gam, phix = phix_test,rpc_kernel = R22/tau2,
           m = m,alpha = alpha, beta = beta, tau2 = tau2, 
           x = x_train, y_train = y_train, x_test = x_test, f0_test = f0_test,
           xmin = xmin, xmax = xmax, ymin_start = ymin, ymax_start = ymax)


# Variogram out
vg_out = list(y_train = y_train, x_train = x_train, h = hgrid, vge = vgyhat$v,
              vgu = vgyhat$u,vg_means = vg_means, vg_params = vg_params)


# Save results
filedir = "/home/johnyannotty/Documents/Dissertation/results/rpath_bart/variogram_gen/"
saveRDS(fit_data,paste0(filedir,"pred/rpath_res2_m1_k1.rds"))
saveRDS(fitp$sdraws[,1], paste0(filedir,"pred/rpath_sdraws2_m1_k1.rds"))

saveRDS(out, paste0(filedir,"vg/tree4_info.rds"))
saveRDS(vg_out, paste0(filedir,"vg/vg_results4_n200.rds"))


#================================================
# Now generate data from 1000 tree kernel 
#================================================
filedir = '/home/johnyannotty/Documents/Dissertation/results/rpath_bart/variogram_exs_1d/cov/'
kern_data = readRDS(paste0(filedir,"rpc_list3_02_21_24.rds"))
m = length(kern_data$tree_list)
n = nrow(kern_data$rpc_list[[1]])

n_train = 200
n_test = 200

rpc = kern_data$rpc_list[[length(kern_data$rpc_list)]]
train_ind = 2*(0:(n_train-1)) + 1
test_ind = 2*(1:n_train)
x_train = kern_data$x[train_ind]
x_test = kern_data$x[test_ind]

#------------------------------------------------
# Generate Data
#------------------------------------------------
sig = 0.005 
R11 = rpc[train_ind,train_ind] + diag(sig^2,n_train)
R22 = rpc[test_ind,test_ind]
R12 = rpc[train_ind,test_ind]

plot_kernel_viridis(R22, vopt = "inferno",k_lim = c(-0.01,0.5,1.01))

m1 = rep(0,n_train)
m2 = rep(0,n_test)

set.seed(31)
y_train = sample_gp(m1,R11)
plot(x_train, y_train)

# True Predictive Distribution 
pdgp = predict_dist_gp(y_train,m1,m2,R11,R22,R12)
f0_test = pdgp$mp

plot(x_test, f0_test, type = "l")
points(x_train, y_train)


#------------------------------------------------
# Variogram Time
#------------------------------------------------
ugrid = seq(0.02,1.0, length = 40)
h_grid = ugrid
vgyhat = variog(coords = cbind(x_train,rep(0,n_train)), data = y_train ,uvec = ugrid)

plot(vgyhat$u,vgyhat$v)

xmin = -1; xmax = 1
xbnds = matrix(c(xmin,xmax), nrow = 1, byrow = FALSE)

param_grid = expand.grid(k = c(1,1.1,1.2), a1 = c(2,5,10), a2 = c(10,25,30), pwr = c(1,1.25,1.5,2))
vg_means = matrix(0, nrow = 0, ncol = length(h_grid))

# Set outfile directory
for(i in 1:nrow(param_grid)){
  vg = variogram.openbtmixing(xbnds,h_grid,10000,1,
                              k=param_grid[i,"k"],
                              0.95,
                              power = param_grid[i,"pwr"],
                              a1 = param_grid[i,"a1"],
                              a2 = param_grid[i,"a2"],
                              4, 
                              type = 'b',
                              ymin = min(y_train),
                              ymax = max(y_train),
                              sigma2 = 0.001,
                              ncut = 200,
                              gam = NULL)

  vg_means = rbind(vg_means, vg$vmean/2)
  #vg_params[[length(vg_params)+1]] = vg$params
}

# Semi-variogram
plot(h_grid,vg_means[12,], type = "l", ylim = c(0,1))
points(vgyhat$u,vgyhat$v)
abline(h = var(y_train), col = "grey")


#-------------------------------------------------
# Batch Fit
#-------------------------------------------------
nu = 20
q0 = 4

# Train a bart model
fit=train.openbtmixing(x_train,y_train,as.matrix(rep(1,n_train)),pbd=c(1.0,0),ntree =20,ntreeh=1,numcut=300,tc=4,model="mixbart",modelname="physics_model",
                       ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
                       power = 1.25, base = 0.95, minnumbot = 1, overallsd = sqrt(0.001), k = 1.1, overallnu = nu,
                       summarystats = FALSE, rpath = TRUE, q = q0, rshp1 = 2, rshp2 = 30,
                       stepwpert = 0.1, probchv = 0.1, batchsize = 10000)


#Get mixed mean function
fitp=predict.openbtmixing(fit,x.test = x_test, f.test = as.matrix(rep(1,n_test)),tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "mean_and_sigma")


g = gammapost.openbtmixing(fit)
hist(g[,1])
plot(g[,1], type = 'l')
apply(g,2,mean)


# Sigma
hist(fitp$sdraws[,1])
plot(fitp$sdraws[,1], type = 'l')

# compare vs true curve
sqrt(mean((fitp$mmean - f0_test)^2))

plot_mean1d(x_test, pmean = fitp$mmean, plb = fitp$m.lower, pub = fitp$m.upper,
            amean = c(f0_test), colors = c("black","purple2"),
            apts_x = x_train, apts_y = y_train)

plot_mean1d(x_test, pmean = fitp$mmean, plb = fitp$m.lower, pub = fitp$m.upper,
            colors = c("black","purple2"),
            apts_x = x_train, apts_y = y_train)


#------------------------------------------------
# Save Results
#------------------------------------------------
fit_data = list(
  pred_mean = fitp$mmean,
  pred_ub = fitp$m.upper,
  pred_lb = fitp$m.lower,
  m = fit$m,
  k = fit$k,
  shp1 = fit$rshp1,
  shp2 = fit$rshp2,
  minnodesz = fit$minnumbot,
  q = q0,
  base = fit$base,
  power = fit$power,
  maxd = fit$maxd,
  nu = fit$overallnu,
  lam = fit$overalllambda,
  sigma = fitp$smean[1],
  x_train = x_train,
  y_train = y_train,
  x_test = x_test,
  f0_test = f0_test
)

# Kernel
out = list(tree_list = tree_list, gam_list = gam_list, phix_list = phix_list,
           rpc_kernel = R22/((ymax - ymin)/(2*kappa))^2, 
           m = m,alpha = alpha, beta = beta,a1 = a1, a2 = a2,k=kappa,tau = tau, 
           x = x_train, y_train = y_train, x_test = x_test,f0_test = f0_test,
           xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)


# Variogram out
vg_out = list(y_train = y_train, x_train = x_train, h = h_grid, vge = vgyhat$v,
              vg_means = vg_means, vg_params = param_grid)

# Save results
filedir = "/home/johnyannotty/Documents/Dissertation/results/rpath_bart/variogram_exs_1d/"
saveRDS(fit_data,paste0(filedir,"pred/rpath_rpc3d_03_30.rds"))
saveRDS(fitp$sdraws[,1], paste0(filedir,"pred/rpath_rpc3d_sdraws_03_30.rds"))
saveRDS(vg_out, paste0(filedir,"vg/vg_rpc3_n200.rds"))


#------------------------------------------------
# Prediction with gp
#------------------------------------------------
#geodata = list(data = y_train, coords = cbind(x_train,runif(n_train,-0.001,0.001)))
geodata = list(data = y_train, coords = cbind(x_train,rep(0,n_train)))
geor_fit = likfit(geodata, cov.model = "gaussian", ini.cov.pars = c(0.5,0.1),fix.nugget = FALSE,
                  nugget = 0.01) 

xtr_ind = expand.grid(1:n_train,1:n_train)
xts_ind = expand.grid(1:n_test,1:n_test)
xtrs_ind = expand.grid(1:n_train,1:n_test)

if(geor_fit$cov.model == "gaussian"){
  gps2 = geor_fit$tausq
  gpsc = sqrt(geor_fit$cov.pars[1])
  gpls = geor_fit$cov.pars[2]^2/2
  
  Rp11 = apply(xtr_ind,1,function(x) sqr_exp_kernel(x_train[x[1]],x_train[x[2]],gpsc,gpls))
  Rp22 = apply(xts_ind,1,function(x) sqr_exp_kernel(x_test[x[1]],x_test[x[2]],gpsc,gpls))
  Rp12 = apply(xtrs_ind,1,function(x) sqr_exp_kernel(x_train[x[1]],x_test[x[2]],gpsc,gpls))
}else{
  gps2 = geor_fit$tausq
  gpsc = sqrt(geor_fit$cov.pars[1])
  gpls = geor_fit$cov.pars[2]
  
  Rp11 = apply(xtr_ind,1,function(x) power_exp_kernel(x_train[x[1]],x_train[x[2]],gpsc,gpls,1))
  Rp22 = apply(xts_ind,1,function(x) power_exp_kernel(x_test[x[1]],x_test[x[2]],gpsc,gpls,1))
  Rp12 = apply(xtrs_ind,1,function(x) power_exp_kernel(x_train[x[1]],x_test[x[2]],gpsc,gpls,1))
}

Rp11 = matrix(Rp11, nrow = n_train, ncol = n_train, byrow = FALSE)
Rp11 = Rp11 + diag(gps2,n_train)

Rp22 = matrix(Rp22, nrow = n_test, ncol = n_test, byrow = FALSE)

Rp12 = matrix(Rp12, nrow = n_train, ncol = n_test, byrow = FALSE)

mp1 = geor_fit$trend.matrix*geor_fit$beta
mp2 = rep(geor_fit$beta,n_test)
pdgp = predict_dist_gp(y_train,mp1,mp2,Rp11,Rp22,Rp12)
fhat = pdgp$mp
fshat = sqrt(diag(pdgp$Rp))
fhat_lb = fhat - 1.96*fshat
fhat_ub = fhat + 1.96*fshat

sqrt(mean((fhat - f0_test)^2))

plot_mean1d(x_test, pmean = fhat, plb = fhat_lb, pub = fhat_ub,
            amean = c(f0_test), colors = c("black","purple2"),
            apts_x = x_train, apts_y = y_train)


# GP Results
gp_data = list(
  fhat = fhat,
  fshat = fshat,
  fhat_lb = fhat_lb,
  fhat_ub = fhat_ub,
  sc = gpsc,
  ls = gpls,
  beta = geor_fit$beta,
  model = geor_fit$cov.model,
  sigma = sqrt(gps2)
)

filedir = "/home/johnyannotty/Documents/Dissertation/results/rpath_bart/variogram_gen/"
saveRDS(gp_data,paste0(filedir,"pred/gp_sqr_exp_rpc3.rds"))

