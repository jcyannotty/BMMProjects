#------------------------------------------------
# EFT Interpolators
#------------------------------------------------
source("/home/johnyannotty/Documents/SciData/NuclearPhysics/R/heat_2d_ising.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/LinearRegression/bayesian_regression.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/kernels.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/gp_utils.R")
#source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/gp_regression.R")

#------------------------------------------------
# Set the data
gvec = seq(0.01,5,length=200)
gset1 = c(0.02,0.1,0.15,0.25)
gset2 = c(3,3.5,4,4.5,5)
gset = c(gset1,gset2)

gset01 = gset/5
gvec01 = gvec/5

#fsg = heat_c2x_sg(gset1)
#flg = heat_c2x_lg(gset2)
#fg = c(fsg,flg)

L = 8
fsg = heat_sg_exp(gset1,10,L)
flg = heat_lg_exp(gset2,10,L)
fg = c(fsg,flg)

xvec = 0.5*log(gvec+1)
dz2 = sapply(xvec,function(x) d2_logz(x,L,0.05))

fsg_grid = heat_sg_exp(gvec,10,L)
flg_grid = heat_lg_exp(gvec,10,L)

#------------------------------------------------
# Polynomial - order 6
#------------------------------------------------
# Set train and test basis
pbasis = cbind(gset^0,gset,gset^2,gset^3,gset^4,gset^5,gset^6)
pbasis = cbind(gset01^0,gset01^4)#)gset^2,gset^3,gset^4,gset^5,gset^6)
colnames(pbasis) = paste0("g",0:(ncol(pbasis)-1))

ptest = cbind(gvec^0,gvec,gvec^2,gvec^3,gvec^4,gvec^5,gvec^6)
ptest = cbind(gvec01^0,gvec01^4)#,gvec^2,gvec^3,gvec^4,gvec^5,gvec^6)
colnames(ptest) = paste0("g",0:(ncol(ptest)-1))

# Set the priors and fit
mu_vec = c(mean(fg),rep(0,ncol(pbasis)-1))
V = diag(5,ncol(pbasis))
fpoly_post = bayes_reg(pbasis,fg,mu_vec,V,30,0.05,5000,1000)
fpoly_mean = predict_bayes(fpoly_post,ptest)

# Trace plot....
plot(fpoly_post$sig2,type = "l")

# Get predictions
pm = apply(fpoly_mean$post_dist,2,mean)
pm_lb = apply(fpoly_mean$post_dist,2,quantile, 0.025)
pm_ub = apply(fpoly_mean$post_dist,2,quantile, 0.975)

plot(gvec,dz2/L^2, type = 'l')
lines(gvec,fsg_grid,col = "red",lty = "dashed")
lines(gvec,flg_grid,col = "blue",lty = "dashed")
lines(gvec,pm,col = "green",lty = "dashed")
lines(gvec,pm_lb,col = "green",lty = "dotted")
lines(gvec,pm_ub,col = "green",lty = "dotted")

#lmdata = data.frame(fg = fg, pbasis[,-1])
#fpoly = lm(fg~., data = lmdata)
#summary(fpoly)
#fgtest = predict(fpoly,data.frame(ptest))


#------------------------------------------------
# Polynomial - order 4
#------------------------------------------------
# Set train and test basis
pbasis = cbind(gset^0,gset,gset^2,gset^3,gset^4)
colnames(pbasis) = paste0("g",0:(ncol(pbasis)-1))

ptest = cbind(gvec^0,gvec,gvec^2,gvec^3,gvec^4)
colnames(ptest) = paste0("g",0:(ncol(ptest)-1) )

# Set the priors and fit
mu_vec = c(mean(fg),rep(0,ncol(pbasis)-1))
V = diag(5,ncol(pbasis))
fpoly_post = bayes_reg(pbasis,fg,mu_vec,V,30,0.05,5000,1000)
fpoly_mean = predict_bayes(fpoly_post,ptest)

# Trace plot....
plot(fpoly_post$sig2,type = "l")

# Get predictions
pm = apply(fpoly_mean$post_dist,2,mean)
pm_lb = apply(fpoly_mean$post_dist,2,quantile, 0.025)
pm_ub = apply(fpoly_mean$post_dist,2,quantile, 0.975)

plot(gvec,dz2/L^2, type = 'l')
lines(gvec,heat_c2x_sg(gvec),col = "red",lty = "dashed")
lines(gvec,heat_c2x_lg(gvec),col = "blue",lty = "dashed")
lines(gvec,pm,col = "green",lty = "dashed")
lines(gvec,pm_lb,col = "green",lty = "dotted")
lines(gvec,pm_ub,col = "green",lty = "dotted")

#lmdata = data.frame(fg = fg, pbasis[,-1])
#fpoly = lm(fg~., data = lmdata)
#summary(fpoly)

#------------------------------------------------
# GP With different Kernels
#------------------------------------------------
# Train kernel
lam = 1
sc = 0.3
gvec_gp = c(0,0.005,seq(0.021,5.2,length=70))
gset_gp = gset[-c(2,3)]
fg_gp = fg[-c(2,3)]

g1grid = expand.grid(gset_gp,gset_gp)
R11 = apply(g1grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc,lam))
R11 = matrix(R11, nrow = length(gset_gp), ncol = length(gset_gp))
dim(R11)
#Rx = Rx + diag(sig2,n)

g2grid = expand.grid(gvec,gvec)
R22 = apply(g2grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc,lam))
R22 = matrix(R22, nrow = length(gvec), ncol = length(gvec))
dim(R22)

g12grid = expand.grid(gset_gp,gvec)
R12 = apply(g12grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc,lam))
R12 = matrix(R12, nrow = length(gset_gp), ncol = length(gvec), byrow = FALSE)
dim(R12)

m1 = rep(mean(fg_gp),length(fg_gp))
m2 = rep(mean(fg_gp),length(gvec))
pdgp = predict_dist_gp(fg_gp,m1,m2,R11,R22,R12)
pdgp_lb = pdgp$mp - 2*diag(pdgp$Rp)^0.5
pdgp_ub = pdgp$mp + 2*diag(pdgp$Rp)^0.5

# Plot
plot(gvec,dz2/L^2, type = 'l')
points(gset_gp,fg_gp)
lines(gvec,heat_c2x_sg(gvec),col = "red",lty = "dashed")
lines(gvec,heat_c2x_lg(gvec),col = "blue",lty = "dashed")
lines(gvec,pdgp$mp,col = "green",lty = "dashed")
lines(gvec,pdgp_lb,col = "green",lty = "dotted")
lines(gvec,pdgp_ub,col = "green",lty = "dotted")

#------------------------------------------------
# GP with trend and different Kernels
#------------------------------------------------
# Get the BLUP 
gset01_gp = gset_gp/5
#X = model.matrix(~I(((gset01_gp*5-2.5)/2.5)^1))
#Xtilde = model.matrix(~I(((gvec01*5-2.5)/2.5)^1))
X = model.matrix(~I(sin(2*gset01_gp*pi))+I(cos(2*gset01_gp*pi)))
Xtilde = model.matrix(~I(sin(2*gvec01*pi))+I(cos(2*gvec01*pi)))

R11_inv = solve(R11)
B = solve(t(X)%*%R11_inv%*%X)%*%t(X)%*%R11_inv
M = Xtilde%*%B + t(R12)%*%R11_inv%*%(diag(1,length(gset01_gp)) - X%*%B)

betahat = solve(t(X)%*%R11_inv%*%X)%*%t(X)%*%R11_inv%*%fg_gp
fhat = M%*%fg_gp 
fhatsd = sqrt(diag(M%*%R11%*%t(M)))

fhat_ub = fhat + 2*fhatsd
fhat_lb = fhat - 2*fhatsd

plot(gvec,dz2/L^2, type = 'l', col="black")
points(gset_gp,fg_gp)
lines(gvec,fsg_grid,col = "red",lty = "dashed")
lines(gvec,flg_grid,col = "blue",lty = "dashed")
lines(gvec,pdgp$mp,col = "green",lty = "dashed")
lines(gvec,pdgp_lb,col = "green",lty = "dotted")
lines(gvec,pdgp_ub,col = "green",lty = "dotted")
lines(gvec,fhat,col = "orange",lty = "dashed")
lines(gvec,fhat_lb,col = "orange",lty = "dotted")
lines(gvec,fhat_ub,col = "orange",lty = "dotted")
lines(gvec,Xtilde%*%betahat,col = "purple",lty = "dashed")
points(gvec,predict(ss,gvec01)$y)


#------------------------------------------------
# Two different GPs
#------------------------------------------------
# GPs with squared exp
gp_modelset = list() 
lam = c(0.4,1,4)
sc = c(0.7,0.4,0.6)
#gvec_gp = c(0,0.005,seq(0.021,5.2,length=70))
gvec_gp = seq(0,5.2,length = 100)
gset_gp = gset[-c(2,3)]
fg_gp = fg[-c(2,3)]

for(j in 1:length(sc)){
  g1grid = expand.grid(gset_gp,gset_gp)
  R11 = apply(g1grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc[j],lam[j]))
  R11 = matrix(R11, nrow = length(gset_gp), ncol = length(gset_gp))
  dim(R11)
  #Rx = Rx + diag(sig2,n)
  
  g2grid = expand.grid(gvec_gp,gvec_gp)
  R22 = apply(g2grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc[j],lam[j]))
  R22 = matrix(R22, nrow = length(gvec_gp), ncol = length(gvec_gp))
  dim(R22)
  
  g12grid = expand.grid(gset_gp,gvec_gp)
  R12 = apply(g12grid,1,function(x) sqr_exp_kernel(x[1],x[2],sc[j],lam[j]))
  R12 = matrix(R12, nrow = length(gset_gp), ncol = length(gvec_gp), byrow = FALSE)
  dim(R12)
  
  m1 = rep(mean(fg_gp),length(fg_gp))
  m2 = rep(mean(fg_gp),length(gvec_gp))
  pdgp = predict_dist_gp(fg_gp,m1,m2,R11,R22,R12)
  pdgp_lb = pdgp$mp - 2*diag(pdgp$Rp)^0.5
  pdgp_ub = pdgp$mp + 2*diag(pdgp$Rp)^0.5
  
  gp_modelset[[j]] = data.frame(m = pdgp$mp, lb = pdgp_lb, ub = pdgp_ub)
}

plot(gvec,dz2/L^2, type = 'l', col="black")
points(gset_gp,fg_gp)
lines(gvec,heat_c2x_sg(gvec),col = "red",lty = "dashed")
lines(gvec,heat_c2x_lg(gvec),col = "blue",lty = "dashed")
lines(gvec_gp,gp_modelset[[1]]$m,col = "green",lty = "dashed")
lines(gvec_gp,gp_modelset[[1]]$lb,col = "green",lty = "dotted")
lines(gvec_gp,gp_modelset[[1]]$ub,col = "green",lty = "dotted")
lines(gvec_gp,gp_modelset[[2]]$m,col = "orange",lty = "dashed")
lines(gvec_gp,gp_modelset[[2]]$lb,col = "orange",lty = "dotted")
lines(gvec_gp,gp_modelset[[2]]$ub,col = "orange",lty = "dotted")
lines(gvec_gp,gp_modelset[[3]]$m,col = "cyan",lty = "dashed")
lines(gvec_gp,gp_modelset[[3]]$lb,col = "cyan",lty = "dotted")
lines(gvec_gp,gp_modelset[[3]]$ub,col = "cyan",lty = "dotted")


# Get SVD and PC
f_grid = cbind(f1=gp_modelset[[1]]$m,f2=gp_modelset[[2]]$m,f3=gp_modelset[[3]]$m)
svd_out = svd(f_grid-matrix(colMeans(f_grid),nrow = nrow(f_grid), ncol = 3, byrow = TRUE))
ev = svd_out$u%*%diag(svd_out$d)

plot(gvec,dz2/L^2, type = 'l', col="black", ylim = c(-10,10))
points(gset_gp,fg_gp)
lines(gvec,heat_c2x_sg(gvec),col = "red",lty = "dashed")
lines(gvec,heat_c2x_lg(gvec),col = "blue",lty = "dashed")
lines(gvec_gp,ev[,1],col = "green",lty = "dashed")
lines(gvec_gp,ev[,2],col = "orange",lty = "dashed")
lines(gvec_gp,ev[,3],col = "cyan",lty = "dashed")


#------------------------------------------------
# Different kernels
#------------------------------------------------
# GPs with squared exp
gset1 = c(0.02,0.2,0.4)
gset2 = c(2,3.5,4,4.5)
gset = c(gset1,gset2)

L = 8
fsg = heat_sg_exp(gset1,10,L)
flg = heat_lg_exp(gset2,10,L)
fg = c(fsg,flg)


gp_modelset = list() 
p1 = rep(0.4,5)
p2 = c(2,1,1.5,2,1.5)
p3 = c(NA,1.5,1.3,0.5,1)
p4 = c(NA,NA,NA,NA,0)
klist = c("se","m","pe","rq","w")
gvec_gp = seq(0,5.2,length = 100)
gset_gp = gset
fg_gp = fg
#gset_gp = gset[-c(2,3)]
#fg_gp = fg[-c(2,3)]

for(j in 1:length(klist)){
  # Get Grids
  g1grid = expand.grid(gset_gp,gset_gp)
  g2grid = expand.grid(gvec_gp,gvec_gp)
  g12grid = expand.grid(gset_gp,gvec_gp)
  
  # Kernel
  if(klist[j] == "se"){
    R11 = apply(g1grid,1,function(x) sqr_exp_kernel(x[1],x[2],p1[j],p2[j]))
    R22 = apply(g2grid,1,function(x) sqr_exp_kernel(x[1],x[2],p1[j],p2[j]))
    R12 = apply(g12grid,1,function(x) sqr_exp_kernel(x[1],x[2],p1[j],p2[j]))
  }else if(klist[j] == "m"){
    R11 = apply(g1grid,1,function(x) matern_kernel(x[1],x[2],p1[j],p2[j],p3[j]))
    R22 = apply(g2grid,1,function(x) matern_kernel(x[1],x[2],p1[j],p2[j],p3[j]))
    R12 = apply(g12grid,1,function(x) matern_kernel(x[1],x[2],p1[j],p2[j],p3[j]))
  }else if(klist[j] == "pe"){
    R11 = apply(g1grid,1,function(x) power_exp_kernel(x[1],x[2],p1[j],p2[j],p3[j]))
    R22 = apply(g2grid,1,function(x) power_exp_kernel(x[1],x[2],p1[j],p2[j],p3[j]))
    R12 = apply(g12grid,1,function(x) power_exp_kernel(x[1],x[2],p1[j],p2[j],p3[j]))
  }else if(klist[j] == "rq"){
    R11 = apply(g1grid,1,function(x) rational_quad_kernel(x[1],x[2],p1[j],p2[j],p3[j]))
    R22 = apply(g2grid,1,function(x) rational_quad_kernel(x[1],x[2],p1[j],p2[j],p3[j]))
    R12 = apply(g12grid,1,function(x) rational_quad_kernel(x[1],x[2],p1[j],p2[j],p3[j]))
  }else if(klist[j] == "w"){
    R11 = apply(g1grid,1,function(x) wendland_kernel(x[1],x[2],p1[j],p2[j],p3[j],p4[j]))
    R22 = apply(g2grid,1,function(x) wendland_kerenl(x[1],x[2],p1[j],p2[j],p3[j],p4[j]))
    R12 = apply(g12grid,1,function(x) wendland_kerenl(x[1],x[2],p1[j],p2[j],p3[j],p4[j]))
  }
  
  # Reshape
  R11 = matrix(R11, nrow = length(gset_gp), ncol = length(gset_gp))
  R22 = matrix(R22, nrow = length(gvec_gp), ncol = length(gvec_gp))
  R12 = matrix(R12, nrow = length(gset_gp), ncol = length(gvec_gp), byrow = FALSE)
  
  # Fit and Predict
  m1 = rep(mean(fg_gp),length(fg_gp))
  m2 = rep(mean(fg_gp),length(gvec_gp))
  pdgp = predict_dist_gp(fg_gp,m1,m2,R11,R22,R12)
  pdgp_lb = pdgp$mp - 2*diag(pdgp$Rp)^0.5
  pdgp_ub = pdgp$mp + 2*diag(pdgp$Rp)^0.5
  
  gp_modelset[[j]] = data.frame(m = pdgp$mp, lb = pdgp_lb, ub = pdgp_ub)
}

par(mfrow = c(3,2))
for(j in 1:length(klist)){
  plot(gvec,dz2/L^2, type = 'l', col="black", main = paste(klist[j],j))
  points(gset_gp,fg_gp)
  lines(gvec,fsg_grid,col = "red",lty = "dashed")
  lines(gvec,flg_grid,col = "blue",lty = "dashed")
  lines(gvec_gp,gp_modelset[[j]]$m,col = "orange",lty = "dashed")
  lines(gvec_gp,gp_modelset[[j]]$lb,col = "orange",lty = "dotted")
  lines(gvec_gp,gp_modelset[[j]]$ub,col = "orange",lty = "dotted")
}
#lines(gvec_gp,gp_modelset[[2]]$m,col = "orange",lty = "dashed")
#lines(gvec_gp,gp_modelset[[2]]$lb,col = "orange",lty = "dotted")
#lines(gvec_gp,gp_modelset[[2]]$ub,col = "orange",lty = "dotted")
#lines(gvec_gp,gp_modelset[[3]]$m,col = "cyan",lty = "dashed")
#lines(gvec_gp,gp_modelset[[3]]$lb,col = "cyan",lty = "dotted")
#lines(gvec_gp,gp_modelset[[3]]$ub,col = "cyan",lty = "dotted")


# Get SVD and PC
f_grid = cbind(f1=gp_modelset[[1]]$m,f2=gp_modelset[[2]]$m,f3=gp_modelset[[3]]$m)
svd_out = svd(f_grid-matrix(colMeans(f_grid),nrow = nrow(f_grid), ncol = 3, byrow = TRUE))
ev = svd_out$u%*%diag(svd_out$d)

plot(gvec,dz2/L^2, type = 'l', col="black", ylim = c(-10,10))
points(gset_gp,fg_gp)
lines(gvec,heat_c2x_sg(gvec),col = "red",lty = "dashed")
lines(gvec,heat_c2x_lg(gvec),col = "blue",lty = "dashed")
lines(gvec_gp,ev[,1],col = "green",lty = "dashed")
lines(gvec_gp,ev[,2],col = "orange",lty = "dashed")
lines(gvec_gp,ev[,3],col = "cyan",lty = "dashed")

