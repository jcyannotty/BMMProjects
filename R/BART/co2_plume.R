#------------------------------------------------
# CO2 Plume Data
#------------------------------------------------
setwd("/home/johnyannotty/Documents/openbt/src")

source("/home/johnyannotty/Documents/openbt/src/openbt.R")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/kernels.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Plotting/computer_expt_plots.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/gp_utils.R")

#------------------------------------------------
filedir = "/home/johnyannotty/Documents/Dissertation/data_sets/"
load(paste0(filedir,"co2plume.dat_.gz"))
load(paste0(filedir,"co2holdout.dat_.gz"))
head(co2plume)
plot(co2plume)

x_train = cbind(st = co2plume$stack_inerts, tm = co2plume$time)
y_train = co2plume$co2

x_test = cbind(st = co2holdout$stack_inerts, tm = co2holdout$time)
y_test = co2holdout$co2

n_train = length(y_train)
n_test = length(y_test)

#------------------------------------------------
# VG Optimized
#------------------------------------------------
vg_loss = function(pvec){
  #k = pvec[1]*(ubnds[1] - lbnds[1]) + lbnds[1]
  #a1 = pvec[2]*(ubnds[2] - lbnds[2]) + lbnds[2]
  #a2 = pvec[3]*(ubnds[3] - lbnds[3]) + lbnds[3]
  #pwr = pvec[4]*(ubnds[4] - lbnds[4]) + lbnds[4]
  #lam = pvec[4]*(ubnds[4] - lbnds[4]) + lbnds[4]
  #sig2 = nu*(lam)/(nu + 2)
  
  k = (exp(pvec[1])/(exp(pvec[1]) + 1))*(ubnds[1] - lbnds[1]) + lbnds[1]
  a1 = (exp(pvec[2])/(exp(pvec[2]) + 1))*(ubnds[2] - lbnds[2]) + lbnds[2]
  a2 = (exp(pvec[3])/(exp(pvec[3]) + 1))*(ubnds[3] - lbnds[3]) + lbnds[3]
  lam = (exp(pvec[4])/(exp(pvec[4]) + 1))*(ubnds[4] - lbnds[4]) + lbnds[4]
  #a1 = exp(pvec[2]) + lbnds[2]
  #a2 = exp(pvec[3]) + lbnds[3]
  #sig2 = nu*(lam)/(nu + 2)
  sig2 = lam
  vg = variogram.openbtmixing(xbnds,hgrid,N,1,
                              k=k,
                              0.95,
                              power = pwr,
                              a1 = a1,
                              a2 = a2,
                              4,
                              ncut = ncut,
                              beta = 0,
                              sigma2 = sig2,
                              maxd = 999,
                              type = "b",
                              ymin = ymin,
                              ymax = ymax
  )
  score = mean((vg$vmean/2-emp_vg)^2) 
  return(score)
}

# Global options
nu = 3
pwr = 1
xbnds = t(apply(x_train,2,range))
N = 10000
ncut = 200
ymin = min(y_train)
ymax = max(y_train)
lbnds = c(1,2,5,0.1)
ubnds = c(5,10,70,2)

# Empirical
hgrid = seq(0.01,0.8, by = 0.02)
vgyhat = variog(coords = x_train, data = y_train, uvec = hgrid)
emp_vg = vgyhat$v
plot(vgyhat$u,vgyhat$v)
abline(h = var(y_train))

# Test the loss
#pvec0 = c(0.75,0.1,0.5,0.5)
#pvec0 = rep(0.1,4)
pvec0 = c(log(0.20)/log(0.8),log(0.1)/log(0.9),log(0.8)/log(0.2),1)
vg_loss(pvec0)

# Run Optim
#vgoptim = optim(pvec0, vg_loss, upper = rep(1,4), method = "L-BFGS-B",
#                lower = rep(0,4), control = list(maxit = 200))
vgoptim = optim(pvec0, vg_loss,method = "SANN",control = list(maxit = 500))
vgoptim$par
vgoptim$value
res = 0 
#for(j in 1:4){res[j] = vgoptim$par[j]*(ubnds[j]-lbnds[j])+lbnds[j]}

pvec = vgoptim$par
for(j in 1:4){res[j] = (exp(pvec[j])/(exp(pvec[j]) + 1))*(ubnds[j] - lbnds[j]) + lbnds[j]}

#res[2] = exp(pvec[2]) + lbnds[2]
#res[3] = exp(pvec[3]) + lbnds[3]
#res[4] = (exp(pvec[4])/(exp(pvec[4]) + 1))*(ubnds[4] - lbnds[4]) + lbnds[4]


sig2 = nu*res[4]/(nu+2)
vg = variogram.openbtmixing(xbnds,hgrid,10000,1,
                            k=2.0,#res[1],
                            0.95,
                            power = 1,#pwr,
                            a1 = 2,#res[2],
                            a2 = 60,#res[3],
                            4,
                            ncut = ncut,
                            beta = 0,
                            sigma2 = sig2,
                            maxd = 10,
                            type = "b",
                            ymin = ymin,
                            ymax = ymax
)

plot(hgrid,vg$vmean/2, type = 'l', ylim = c(0,2000))
points(hgrid,emp_vg)
abline(h = var(y_train))

#res1 = res # hgrid = seq(0.02,0.4, by = 0.02)
#res2 = res # hgrid = seq(0.02,0.2, by = 0.02)

#-------------------------------------------------
# Batch Fit
#-------------------------------------------------
# Train a bart model
q0 = 4
res[1] = 3.0; res[2] = 2; res[3] = 60; pwr = 1
nu = 1
fit=train.openbtmixing(x_train,y_train-mean(y_train),as.matrix(rep(1,n_train)),pbd=c(1.0,0),ntree = 75,ntreeh=1,
                       numcut=300,tc=4,model="mixbart",modelname="co2_plume",
                       ndpost = 10000, nskip = 2000, nadapt = 5000, adaptevery = 500, printevery = 500,
                       power = pwr, base = 0.95, minnumbot = 2, overallsd = sqrt(sig2), k = res[1], overallnu = nu,
                       summarystats = FALSE, rpath = TRUE, q = q0, rshp1 = res[2], rshp2 = res[3],
                       stepwpert = 0.1, probchv = 0.1, batchsize = 5000, maxd = 10)

#Get mixed mean function
fitp=predict.openbtmixing(fit,x.test = x_test, f.test = as.matrix(rep(1,n_test)),tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "mean_and_sigma", proj_type = "softmax", temperature = 0.2)


# Sigma
hist(unlist(fitp$sdraws[,1]))
plot(unlist(fitp$sdraws[,1]))

# Gamma
#gp = gammapost.openbtmixing(fit)
#hist(gp[,4])

# on Testing data
rdf = data.frame(x_test,y_test - fitp$mmean-mean(y_train))
rownames(rdf) = NULL
colnames(rdf)[3] = 'r'

rdf %>% ggplot() + geom_point(aes(x = st, y = tm, color = r), size = 3) + 
  scale_color_gradientn(limits = c(-10,10), colors = c("navy","grey95","darkred")) + 
  theme_bw() + theme(axis.line = element_line(color = "grey70"),
                     panel.border = element_blank(),plot.title = element_text(hjust = 0.5),
                     legend.position = "right", legend.text = element_text(size = 10)) 

sqrt(mean((fitp$mmean + mean(y_train) - y_test)^2))

# Trace plots of selected predictions/fits
hist(unlist(fitp$mdraws[,1]))
plot(unlist(fitp$mdraws[,188]), type = 'l')

# Coverage
mean(((y_test - mean(y_train)) > fitp$m.lower) & ((y_test - mean(y_train)) < fitp$m.upper))


# Save results
fit_data = list(
  pred_mean = fitp$mmean + mean(y_train),
  pred_ub = fitp$m.upper + mean(y_train),
  pred_lb = fitp$m.lower + mean(y_train),
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
  lam = fit$overalllambda
)

# Save results
filedir = "/home/johnyannotty/Documents/RandomPathBART/BART_Results/"
saveRDS(fit_data,paste0(filedir,"co2_res5_02_07_24.rds"))
saveRDS(fitp$sdraws[,1], paste0(filedir,"co2_sdraws5_02_07_24.rds"))


#------------------------------------------------
# On training data
#------------------------------------------------
fitp=predict.openbtmixing(fit,x.test = x_train, f.test = as.matrix(rep(1,n_train)),tc=4, q.lower = 0.025, q.upper = 0.975,
                          ptype = "mean_and_sigma", proj_type = "softmax", temperature = 0.2)

sqrt(mean((fitp$mmean + mean(y_train) - y_train)^2))

df = data.frame(x_train,fitp$mmean+mean(y_train))
rownames(df) = NULL
colnames(df)[3] = 'fit'

df %>% ggplot() + geom_point(aes(x = st, y = tm, color = fit),size = 3) + 
  scale_color_viridis(option = "viridis", limits = c(0,220)) + 
  theme_bw() + theme(axis.line = element_line(color = "grey70"),
                     panel.border = element_blank(),plot.title = element_text(hjust = 0.5),
                     legend.position = "right", legend.text = element_text(size = 10)) 

co2plume %>% ggplot() + geom_point(aes(x = stack_inerts, y = time, color = co2),size = 3) + 
  scale_color_viridis(option = "viridis", limits = c(0,220)) + 
  theme_bw() + theme(axis.line = element_line(color = "grey70"),
                     panel.border = element_blank(),plot.title = element_text(hjust = 0.5),
                     legend.position = "right", legend.text = element_text(size = 10)) 

rdf = data.frame(x_train,y_train - fitp$mmean-mean(y_train))
rownames(rdf) = NULL
colnames(rdf)[3] = 'r'

rdf %>% ggplot() + geom_point(aes(x = st, y = tm, color = r), size = 3) + 
  scale_color_viridis(option = "turbo", limits = c(-15,15)) + 
  theme_bw() + theme(axis.line = element_line(color = "grey70"),
                     panel.border = element_blank(),plot.title = element_text(hjust = 0.5),
                     legend.position = "right", legend.text = element_text(size = 10)) 


#------------------------------------------------
# VG Batch Fit
#------------------------------------------------
# Empirical Variogram
n_train = length(y_train)
hgrid = seq(0.02,1, by = 0.025)
vgyhat = variog(coords = x_train, data = y_train ,uvec = hgrid)
plot(vgyhat$u,vgyhat$v)
abline(h = var(y_train))

# x bounds and h grid
xbnds = t(apply(x_train,2,range))
param_list = list(k = c(1.0,1.5,2.0,3.0), a1 = c(2,5,10), a2 = c(10,20,40,60),power = c(0.5,1,2))
param_grid = expand.grid(param_list)
ncut = 200
vg_out = matrix(0,nrow = nrow(param_grid), ncol = length(hgrid))
sig2 = 0.25
for(j in 1:nrow(param_grid)){
  vg = variogram.openbtmixing(xbnds,hgrid,10000,1,
                              k=param_grid[j,"k"],
                              0.95,
                              power = param_grid[j,"power"],
                              a1 = param_grid[j,"a1"],
                              a2 = param_grid[j,"a2"],
                              4,
                              ncut = ncut,
                              beta = 0,
                              sigma2 = sig2,
                              maxd = 10,
                              type = "b",
                              ymin = min(y_train),
                              ymax = max(y_train)
  )
  cat("Progress: ", round(j/nrow(param_grid),4)*100)
  vg_out[j,] = vg$vmean
}

# Semi-variogram
plot(hgrid,vg_out[112,]/2, type = "l", ylim = c(0,700))
points(vgyhat$u,vgyhat$v)
abline(h = var(y_train), col = "grey")
#points(hgrid,sapply(h-grid,function(h) (1-power_exp_kernel(0,h,1,0.3,1))), pch = 3, col = 'red')

out = list(vg = vg_out, xbnds = xbnds, h = hgrid, pgrid = param_grid,
           sigma2 = sig2,vgyhat = 2*vgyhat$v, vgyhat_uvec = vgyhat$u)
saveRDS(out, paste0(filedir,"Variograms/vgb_co2_02_08_24.rds"))
