#------------------------------------------------
# CMIP6 Data Pre-Analysis
#------------------------------------------------
setwd("/home/johnyannotty/Documents/openbt/src")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/kernels.R")

library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)
library(plotly)
library(viridis)
library(latex2exp)
library(geoR)

filedir = "/home/johnyannotty/Documents/CMIP6_mixing/"
datadir = "Data/World/"
resdir = "Results/World/"
ffold = "nh_6month_2014_n30k/"
#dataname = "CESM2_CNRM_NWH_12M_2012-2014_01_26_24_n72000.rds" #ACC_BCC_CESM2_CNRM_NH_6M14_01_06_24_n30000.rds"
dataname = "ACC_BCC_MIROC_CMCC_CESM2_CNRM_CanESM5_KIOST_W_3M_2014_01_26_24_n45000.rds"

ms = readRDS(paste0(filedir,datadir,dataname))
K = ncol(ms$f_train)

#-------------------------------------------------
# Theoretical Variogram
#-------------------------------------------------
# Theoretical
xbnds = t(apply(ms$x_test,2,range))
#xbnds = matrix(c(xmin,xmax),ncol = 2, byrow = TRUE)
# drop time
xbnds = xbnds[-3,]
hgrid = seq(0.5,80,by = 1.5)

# Grid Search
param_list = list(k = c(0.5,0.8,1,2), a1 = c(2,5,10), a2 = c(5,20,40),
                  power = c(0.5,1))
param_grid = expand.grid(param_list)
vg_matrix = matrix(0,nrow = nrow(param_grid), ncol = length(hgrid))
for(i in 1:nrow(param_grid)){
  vg = variogram.openbtmixing(xbnds,hgrid,10000,1,
                              k=param_grid[i,"k"],
                              0.95,
                              power = param_grid[i,"power"],
                              a1 = param_grid[i,"a1"],
                              a2 = param_grid[i,"a2"],
                              4,
                              maxd = 5)
  vg_matrix[i,] = vg$vmean
  cat("Progress: ", round(i/nrow(vg_matrix),4)*100)
}

plot(h_grid,vg$vmean, type = "l", ylim = c(0,0.3))

# Save results
out = list(vg = vg_matrix, xbnds = xbnds, h = hgrid, pgrid = param_grid,
           dataref = paste0(filedir,datadir,dataname))
saveRDS(out, paste0(filedir,"Variograms/wts_nwhem_n30000.rds"))


#-------------------------------------------------
# Precision Weights for Variogram
#-------------------------------------------------
# Precision Weight calculation
pw_se = apply(ms$f_train,2,function(x) (ms$y_train-x)^2)
pw_denom = rowSums(1/pw_se)
pw = (1/pw_se)/pw_denom
#rowSums(pw)

# Get nearest neightbors
#knn = 2
#d = abs(outer(ms$x_train,ms$x_train, "-"))
#d = outer(x_train[,1],x_train[,1], "-")^2 + outer(x_train[,2],x_train[,2], "-")^2
#nn = t(matrix(apply(d,1,function(x) order(x)[2:(knn+1)]),nrow = knn))
#pw_nn = 0*pw

#for(i in 1:nrow(ms$x_train)){
  #ind = c(i,nn[i,])
  #pw_se_nn = apply(pw_se[ind,],2,sum)
  #pw_nn[i,] = (1/pw_se_nn)/sum(1/pw_se_nn)
#}

rowSums(pw_nn)
plot(pw[,2])
plot(pw_nn[,1])

xdata = ms$x_train[,c(1,2,4)]
vghat = c() 
vghat_nn = c()
#h_grid = seq(0.9,max(x_train)-min(x_train),by = 1)
ugrid = seq(0.5,80,by = 1)
K = ncol(ms$f_train)
for(j in 1:K){
  geor_vghat = variog(coords = xdata, data = pw[,j],uvec = ugrid)
  vgh = 2*geor_vghat$v
  vghat = cbind(vghat,vgh)
  
  #geor_vghat = variog(coords = xdata, data = pw_nn[,j],uvec = ugrid)
  #vgh = 2*geor_vghat$v
  #vghat_nn = cbind(vghat_nn,vgh)
}

colnames(vghat) = paste0("vw",1:K)
colnames(vghat_nn) = paste0("vw",1:K)


# Plot the empirical variogram for the mean weight
vghat_mean = rowMeans(vghat)
plot(geor_vghat$u,vghat_mean, ylim = c(0,0.2))

plot(geor_vghat$u,vghat[,3], ylim = c(0,0.2))
abline(h = 2*var(pw[,3]))

# Save results
#out = list(vghat = vghat,vghat_nn = vghat_nn,hgrid = geor_vghat$u, 
#           what = pw, what_nn = pw_nn, x = x_train,y = y_train)
#saveRDS(out, paste0(filedir,"pwtrig3_emp_variogram_precwts_11_15_23.rds"))



# For climate models
ind = 4
plot(ms$x_train[,1],ms$x_train[,2], 
     col = ifelse(pw[,ind]>0.75,"black",
                  ifelse(pw[,ind]>0.5,"darkred",
                        ifelse(pw[,ind]>0.25,"orange2","yellow"))), pch = 16,cex = 0.5)

#-------------------------------------------------
# Variogram in terms of y
#-------------------------------------------------
# In terms of y
h = sample(1:nrow(ms$x_train), size = 10000)
xdata = ms$x_train[h,c(1,2,4)]
ugrid = seq(0.5,80,by = 1.5)

#geor_vghat = variog(coords = xdata, data = ms$y_train - rowMeans(ms$f_train), uvec = ugrid)
#geor_vghat = variog(coords = xdata, data = ms$y_train, uvec = ugrid)

vgf_list = list()
for(j in 1:K){
  vg_tmp = variog(coords = xdata, data = ms$f_train[h,j], uvec = ugrid)
  vgf_list[[j]] = vg_tmp
}

geor_vghaty = variog(coords = xdata, data = ms$y_train[h], uvec = ugrid)

plot(geor_vghaty$u, geor_vghaty$v*2, ylim = c(0,1700), type = 'l', col='black')
for(j in 1:K){
  lines(geor_vghaty$u, vgf_list[[j]]$v*2, col=j+1)
}

evg_f = matrix(0, nrow = length(geor_vghaty$u), ncol = K)  
for(j in 1:K){
  evg_f[,j] = 2*vgf_list[[j]]$v
}
colnames(evg_f) = paste0("f",1:K)

empvg_out = list(u = geor_vghaty$u, evg_y = 2*geor_vghaty$v, 
     evg_f = evg_f,
     dataref = paste0(filedir,datadir,dataname))

saveRDS(empvg_out, paste0(filedir,"Variograms/empvg_world_n45000.rds"))

#-------------------------------------------------
# Variogram in terms of fx
hgrid = seq(0.5,80,by = 1.5)
p = 3
#sqr_exp_kernel(0,rep(10/p,p),1,c(70,20))

sc2 = 265
nug = 40
plot(geor_vghaty$u, geor_vghaty$v*2, ylim = c(0,1000), type = 'l', col='black')
lines(geor_vghaty$u, vgf_list[[8]]$v*2, col='red')
Rxvec = sapply(hgrid,function(x) sqr_exp_kernel(0,rep(x/sqrt(p)),sqrt(sc2),c(700,700,700)))
points(hgrid,2*(sc2-Rxvec + nug))

plot(Rxvec)

plot(hgrid,200-Rxvec)

# x bounds
xbnds = t(apply(ms$x_test,2,range))
xbnds = xbnds[-3,]
#h_grid = seq(0.1,80,by = 0.1)
xgrid = cbind(runif(10000,xbnds[1,1],xbnds[1,2]),
              runif(10000,xbnds[2,1],xbnds[2,2]),
              runif(10000,xbnds[3,1],xbnds[3,2])
              )

# Grid Search
param_list = list(k = c(0.5,0.8,1), a1 = c(2,10), a2 = c(20,40),
                  power = c(0.5,1))
param_grid = expand.grid(param_list)
vg_matrix = matrix(0,nrow = nrow(param_grid), ncol = length(hgrid))

for(i in 1:nrow(param_grid)){
  vg = variogram.openbtmixing(xbnds[1:3,],hgrid,10000,1,
                              k=param_grid[i,"k"],
                              0.95,
                              power = param_grid[i,"power"],
                              a1 = param_grid[i,"a1"],
                              a2 = param_grid[i,"a2"],
                              4,
                              ncut = 1000,
                              beta = 0.25,
                              sigma2 = 1,
                              maxd = 5,
                              xgrid = xgrid,
                              rgrid = Rxvec,
                              fmean = 7.35,
                              rscale = 320,
                              type = "y"
  )
  vg_matrix[i,] = vg$vmean
  cat("Progress: ", round(i/nrow(vg_matrix),4)*100)
}


vg = variogram.openbtmixing(xbnds[1:3,],hgrid,10000,1,
                            k=1.0,
                            0.95,
                            power = 0.5,
                            a1 = 2,
                            a2 = 40,
                            4,
                            ncut = 1000,
                            beta = 0.25,
                            sigma2 = 1,
                            maxd = 5,
                            xgrid = xgrid,
                            rgrid = Rxvec,
                            fmean = 7.35,
                            rscale = 220+55,
                            type = "y"
)

plot(hgrid,vg$vmean, type = "l", ylim = c(100,1000))
lines(geor_vghaty$u, geor_vghaty$v*2,col='blue')
lines(geor_vghat1$u, geor_vghat1$v*2, col='red')
points(hgrid,2*(200-Rxvec) + 110)
abline(h = 2*var(ms$f_train[,1]))


#-------------------------------------------------
# Variogram in terms of fx
hgrid = seq(0.5,80,by = 1.5)
p = 3
#sqr_exp_kernel(0,rep(10/p,p),1,c(70,20))
sc0 = 300
nug0 = 35
lam0 = 400
plot(geor_vghaty$u, geor_vghaty$v*2, ylim = c(0,1700), type = 'l', col='black')
lines(geor_vghat2$u, geor_vghat1$v*2, col='red')
Rxvec = sapply(hgrid,function(x) sqr_exp_kernel(0,rep(x/sqrt(p)),sqrt(sc0),c(lam0,lam0,lam0)))
points(hgrid, 2*(sc0-Rxvec + nug0))

# x bounds
xbnds = t(apply(ms$x_test,2,range))
xbnds = xbnds[-3,]
#h_grid = seq(0.1,80,by = 0.1)
xgrid = cbind(runif(10000,xbnds[1,1],xbnds[1,2]),
              runif(10000,xbnds[2,1],xbnds[2,2]),
              runif(10000,xbnds[3,1],xbnds[3,2])
)

# Grid Search
fmeans = apply(ms$f_train,2,mean)
# North Hemi
#scale_list = c(150,150,150,150)
#nug_list = c(70,70,65,75)
#lc_list = c(450,450,480,480)

# South Hemi
#scale_list = c(300,300,300,300)
#nug_list = c(35,25,30,30)
#lc_list = c(400,400,400,400)

# NW Hemi
#scale_list = c(170,170)
#nug_list = c(50,60)
#lc_list = c(600,600)

# World
scale_list = c(240,230,225,240,240,265,265,265)
nug_list = c(60,50,60,70,60,65,60,40)
lc_list = c(500,500,700,620,600,700,700,700)

Rxvec_list = list()
for(j in 1:K){
  Rxvec_list[[j]] = sapply(hgrid,function(x) sqr_exp_kernel(0,rep(x/sqrt(p)),sqrt(scale_list[j]),rep(lc_list[j],p)))
}

param_list = list(k = c(1.5,1.7,1.8), a1 = c(2), a2 = c(40),power = c(0.5))
param_grid = expand.grid(param_list)
vg_out = matrix(0,nrow = nrow(param_grid), ncol = length(hgrid))
vg_matrix = matrix(0,nrow = K, ncol = length(hgrid))
for(j in 1:nrow(param_grid)){
  for(i in 1:ncol(ms$f_train)){
    vg = variogram.openbtmixing(xbnds[1:3,],hgrid,10000,1,
                                k=param_grid[j,"k"],
                                0.95,
                                power = param_grid[j,"power"],
                                a1 = param_grid[j,"a1"],
                                a2 = param_grid[j,"a2"],
                                4,
                                ncut = 1000,
                                beta = 1/K,
                                sigma2 = 0.55,
                                maxd = 5,
                                xgrid = xgrid,
                                rgrid = Rxvec_list[[i]],
                                fmean = fmeans[i],
                                rscale = scale_list[i]+nug_list[i],
                                type = "y"
    )
    vg_matrix[i,] = vg$vmean
  }
  cat("Progress: ", round(j/nrow(param_grid),4)*100)
  vg_out[j,] = apply(vg_matrix,2,sum)
}

plot(hgrid, apply(vg_matrix,2,sum), type = 'l', ylim = c(100,1000))
lines(geor_vghaty$u, geor_vghaty$v*2,col='blue')

plot(hgrid,vg_out[1,], type = "l", ylim = c(100,1000))
lines(geor_vghaty$u, geor_vghaty$v*2,col='blue')

# Save results
out = list(vg = vg_out, xbnds = xbnds, h = hgrid, pgrid = param_grid,sigma2 = 0.55,
           fmean = fmeans, Rxvec_list = Rxvec_list,
           scale2_list = scale_list, len_scale_list = lc_list,
           nug_list = nug_list, sigma2 =1, kernels = "sqr_exponential",
           dataref = paste0(filedir,datadir,dataname))
saveRDS(out, paste0(filedir,"Variograms/vgy_world_2_n45000.rds"))


#-------------------------------------------------
# Scratch....
library(GPfit)
apply(ms$f_train,2,mean)

n_train = nrow(ms$x_train)
h = sample(n_train, size = 300)
h = sort(h)

hp = sample(n_train, size = 2000)
hp = sort(hp)

x_scale = apply(ms$x_train,2,function(x) (x - min(x))/(max(x) - min(x)))

gp1 = GP_fit(x_scale[h,], ms$f_train[h,1])
gp1$beta
gp1$sig2
gp1$correlation_param

gp1_pred = predict(gp1, x_scale[hp,])
mean((gp1_pred$Y_hat - ms$f_train[hp,1])^2)
max((gp1_pred$Y_hat - ms$f_train[hp,1])^2)

x_scale = apply(ms$x_train,2,function(x) (x - min(x))/(max(x) - min(x)))
hmat = expand.grid(h,h)
Rx = apply(hmat,1,function(x) sqr_exp_kernel(x_scale[x[1],],x_scale[x[2],],1,1))
Rx = matrix(Rx, nrow = length(h))

