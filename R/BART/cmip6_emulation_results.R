#------------------------------------------------
# Climate CMIP6 results loader - from Unity
#------------------------------------------------
setwd("/home/johnyannotty/Documents/openbt/src")
source("/home/johnyannotty/Documents/openbt/src/openbt_mixing.R")
source("/home/johnyannotty/Documents/openbt/R/eft_mixing_helper_functions.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Plotting/computer_expt_plots.R")
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/GaussianProcesses/kernels.R")


library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)
library(plotly)
library(viridis)
library(geoR)
library(latex2exp)

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


#------------------------------------------------
# Variogram
#------------------------------------------------
# Empirical Variogram
filedir = "/home/johnyannotty/Documents/CMIP6_mixing/"
datadir = "Data/World/"
emudir = "Data/World/Emulation/"
resdir = "Results/World/"
vgdir = "Variograms/"

dataname = "ACC_BCC_MIROC_CMCC_CESM2_CNRM_CanESM5_KIOST_W_3M_2014_01_26_24_n45000.rds"
emudataname = "MIROC_World_ev_Aug2014_03_10_24.rds"
ffold = "world_miroc_emulation_aug2014_ev_m300/"

ms = readRDS(paste0(filedir,datadir,dataname))
ems = readRDS(paste0(filedir,emudir,emudataname))

svgdata = readRDS(paste0(filedir,vgdir,"svgb_miroc_aug2014.rds"))
vg_out = svgdata$svgb


# Empirical Variogram
x_train = rbind(ems$x_train,ems$x_test)
y_train = c(ems$y_train,ems$y_test)
hgrid = seq(0.5,400,by = 5)
vgyhat = variog(coords = x_train, data = y_train ,uvec = hgrid)
plot(vgyhat$u,vgyhat$v, ylim = c(0,500))
abline(h = var(y_train))

# x bounds
xbnds = t(apply(x_train,2,range))
xbnds = xbnds[1:2,]

param_list = list(k = c(1.5,1.7,2.0), a1 = c(2,5), a2 = c(20,40,65),power = c(1.0,1.5,1.75))
param_list = list(k = c(1.7), a1 = c(2), a2 = c(65),power = c(1.75))
param_grid = expand.grid(param_list)
vg_out = matrix(0,nrow = nrow(param_grid), ncol = length(hgrid))
sig2 = 0.05
for(j in 1:nrow(param_grid)){
  vg = variogram.openbtmixing(xbnds,hgrid,10000,1,
                              k=param_grid[j,"k"],
                              0.95,
                              power = param_grid[j,"power"],
                              a1 = param_grid[j,"a1"],
                              a2 = param_grid[j,"a2"],
                              4,
                              ncut = 500,
                              beta = 0,
                              sigma2 = sig2,
                              maxd = 999,
                              type = "b",
                              ymin = min(y_train),
                              ymax = max(y_train)
  )
  cat("Progress: ", round(j/nrow(param_grid),4)*100)
  vg_out[j,] = vg$vmean/2
}

# Semi-variogram
plot(hgrid,vg_out[1,], type = "l", ylim = c(0,800))
lines(hgrid,vg_out[15,])
points(vgyhat$u,vgyhat$v)
abline(h = var(y_train), col = "grey")

#vgm1 = likfit(coords = x_train[,1:2], data = y_train, cov.model = "exponential", ini.cov.pars = c(sqrt(380),40))
vgm1 = variofit(vario = vgyhat, cov.model = "gaussian", ini.cov.pars = c(380,40), nugget = 0.1)
vgm1

# Other svg models 
exp_cov = list()
exp_params = list(s2 = c(600,400), ls = c(100,70), pwr = 1) #400, 70 and 600, 100
exp_cov[[1]] = sapply(hgrid, function(x) power_exp_kernel(x,0,sqrt(exp_params$s2[1]),
                                                          exp_params$ls[1],exp_params$pwr))
exp_cov[[2]] = sapply(hgrid, function(x) power_exp_kernel(x,0,sqrt(exp_params$s2[2]),
                                                          exp_params$ls[2],exp_params$pwr))

sqrexp_params = list(s2 = 391, ls = 1800) # 380, 1200 and
sqrexp_cov = sapply(hgrid, function(x) sqr_exp_kernel(x,0,sqrt(sqrexp_params$s2),sqrexp_params$ls))

mat_params = list(s2 = 440, ls = 50, nu = 3/2) #380, 40
mat_cov = sapply(hgrid, function(x) matern_kernel(x,0,sqrt(mat_params$s2),mat_params$ls,mat_params$nu))

svg = cbind(exp1 = exp_params$s2[1] - exp_cov[[1]], exp2 = exp_params$s2[2] - exp_cov[[2]], 
            sqrexp = sqrexp_params$s2 - sqrexp_cov,
            mat = mat_params$s2 - mat_cov)

pb = plot_mean1d(hgrid,pmean = t(vg_out[c(49,50),]), apts_x = hgrid, apts_y = vgyhat$v, apts_col = 'grey50',
                 apts_sz = 2.0, x_lab = "h", y_lab = TeX("$\\nu(h)$"), apts_alpha = 0.5,
                 colors = c("red","blue")) + 
  geom_hline(yintercept = var(y_train), col = "grey40", alpha = 0.4, linewidth = 1) + 
  coord_cartesian(xlim = c(0,250), ylim = c(0,600)) + 
  theme(legend.position = "none", axis.text = element_text(size = 13), 
        axis.title.y = element_blank(), axis.title.x = element_text(size = 13)) +
  labs(title = "Random Path")

pe = plot_mean1d(hgrid,pmean = svg[,c(1,2)], apts_x = hgrid, apts_y = vgyhat$v, apts_col = 'grey50',
                 apts_sz = 1.8, x_lab = "h", y_lab = TeX("$\\nu(h)$"),apts_alpha = 0.5,colors = c("red","blue")) + 
  geom_hline(yintercept = var(y_train), col = "grey40", alpha = 0.4, linewidth = 1) + 
  coord_cartesian(xlim = c(0,250), ylim = c(0,600)) + 
  theme(legend.position = "none", axis.text = element_text(size = 13), 
        axis.title.y = element_blank(), axis.title.x = element_text(size = 13)) + 
  labs(title = "Exponential")

ps = plot_mean1d(hgrid,pmean = svg[,3], apts_x = hgrid, apts_y = vgyhat$v, apts_col = 'grey50',
                 apts_sz = 1.8, x_lab = "h", y_lab = TeX("$\\nu(h)$"), apts_alpha = 0.5,colors = "green3") + 
  geom_hline(yintercept = var(y_train), col = "grey40", alpha = 0.4, linewidth = 1) + 
  coord_cartesian(xlim = c(0,250), ylim = c(0,600)) + 
  theme(legend.position = "none", axis.text = element_text(size = 13), 
        axis.title.y = element_blank(), axis.title.x = element_text(size = 13)) + 
  labs(title = "Squared Exponential")

pm = plot_mean1d(hgrid,pmean = svg[,4], apts_x = hgrid, apts_y = vgyhat$v, apts_col = 'grey50',
                 apts_sz = 1.8, x_lab = "h", y_lab = TeX("$\\nu(h)$"), apts_alpha = 0.5,colors = "green3") + 
  geom_hline(yintercept = var(y_train), col = "grey40", alpha = 0.4, linewidth = 1) + 
  coord_cartesian(xlim = c(0,250), ylim = c(0,600)) + 
  theme(legend.position = "none", axis.text = element_text(size = 13), 
        axis.title.y = element_blank(), axis.title.x = element_text(size = 13)) + 
  labs(title = "Matern")


grid.arrange(pb,pe,pm,nrow = 1, left = grid.text(TeX("$\\nu(h)$")))

svgb_miroc = list(param_grid = param_grid, hgrid = hgrid, dataname = emudataname,
                  svgy = vgyhat$v, svgb = vg_out/2, label = "semivargiogram",
                  exp_params = exp_params, exp_cov = exp_cov,
                  sqrexp_params = sqrexp_params, exp_cov = sqrexp_cov,
                  mat_params = mat_params, mat_cov = mat_cov)

#saveRDS(svgb_miroc, paste0(filedir,vgdir,"svgb_miroc_aug2014.rds"))

#------------------------------------------------
# Results loader
#------------------------------------------------
fls = system(paste0("ls ",filedir,resdir,ffold),intern = TRUE)
fls = order_fnames(fls)
batch = 0
for(i in 1:length(fls)){
  if(grepl("batch",fls[i])){
    batch = batch + 1
    temp_fit = readRDS(paste0(filedir,resdir,ffold,fls[i]))
    if(batch == 1){
      fit = temp_fit
    }else{
      fit = list(
        pred_mean = c(fit$pred_mean,temp_fit$pred_mean),
        pred_ub = c(fit$pred_ub,temp_fit$pred_ub),
        pred_lb = c(fit$pred_lb,temp_fit$pred_lb),
        wts_mean = rbind(fit$wts_mean,temp_fit$wts_mean),
        wts_ub = rbind(fit$wts_ub,temp_fit$wts_ub),
        wts_lb = rbind(fit$wts_lb,temp_fit$wts_lb),
        wsum_mean = c(fit$wsum_mean,temp_fit$wsum_mean),
        wsum_lb = c(fit$wsum_lb,temp_fit$wsum_lb),
        wsum_ub = c(fit$wsum_ub,temp_fit$wsum_ub),
        proj_mean = c(fit$proj_mean,temp_fit$proj_mean),
        proj_ub = c(fit$proj_ub,temp_fit$proj_ub),
        proj_lb = c(fit$proj_lb,temp_fit$proj_lb),
        pwts_mean = rbind(fit$pwts_mean,temp_fit$pwts_mean),
        pwts_ub = rbind(fit$pwts_ub,temp_fit$pwts_ub),
        pwts_lb = rbind(fit$pwts_lb,temp_fit$pwts_lb),
        delta_mean = c(fit$delta_mean,temp_fit$delta_mean),
        delta_lb = c(fit$delta_lb,temp_fit$delta_lb),
        delta_ub = c(fit$delta_lb,temp_fit$delta_ub)
      )
    }
  }else if(grepl("sdraws",fls[i])){
    sfit = readRDS(paste0(filedir,resdir,ffold,fls[i]))
  }
}

#------------------------------------------------
# Results
#------------------------------------------------
# Prediciton
usa = map_data("world",region = "USA")
world = map_data("world")
states = map_data("state")

# Convert long to (-180,180) scale
ntp = 3
tp_len = nrow(ms$f_test)/ntp
for(i in 1:ntp){
  assign(paste0("h",i),(tp_len*(i-1)+1):(tp_len*i))
}
hlist = list(h1,h2,h3)

xs_test = ms$x_test
xs_test[,"lon"] = ifelse(xs_test[,"lon"] > 180,xs_test[,"lon"]-360,xs_test[,"lon"])
xs_test[,1] = ifelse(xs_test[,1] < 0,xs_test[,1]-0.25,xs_test[,1])

pb = plot_mean2d_map_viridis(xs_test[h1,], fit$pred_mean,xcols = c(1,2), 
                        viridis_opt = "viridis",
                        scale_limit = c(-65,42), title = "Random Path Emulation",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(-90,90),
                        lon_bnds = c(-180,180)
)

pb = pb + theme(axis.text = element_text(size = 13),axis.title = element_text(size = 13),
                    legend.position = "bottom", legend.key.width = unit(2.1,'cm'),
                    axis.title.y = element_blank()) + 
  labs(x = "Longitude", fill = "Celsius")

pleg = g_legend(pb)


pdag = plot_mean2d_map_viridis(xs_test[h2,], ms$f_test[h2,3],xcols = c(1,2), 
                        viridis_opt = "viridis",
                        scale_limit = c(-65,42), title = "Bilinear Interpolation",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(-90,90),
                        lon_bnds = c(-180,180)
)

pdag = pdag + theme(axis.text = element_text(size = 13),axis.title = element_text(size = 13),
           legend.position = "none", legend.key.width = unit(2.2,'cm'),
           axis.title.y = element_blank()) + 
  labs(x = "Longitude", fill = "Celsius")

resid = ms$f_test[h2,3] - fit$pred_mean
rmin = -6; rmax = 6
rscale = turbo(5)
rscale[3] = 'white' #c("darkblue","gray95","darkred")
rp = plot_mean2d_map_gradient(xs_test[h2,1:2],ifelse(resid>5,5,ifelse(resid < (-5),-5,resid)), xcols = c(1,2), 
                               title = "Mean Residuals", 
                               scale_colors = rscale,
                               scale_vals = c(-5,-2,0,2,5),
                               maps_list = list(data.frame(world),data.frame(states)),
                               maps_cols = c("grey30","grey30"),
                               lat_bnds = ms$lat_bnds,
                               lon_bnds = c(min(xs_test[,'lon']),max(xs_test[,'lon']))
)

rp = rp + theme(axis.text = element_text(size = 13),axis.title = element_text(size = 13),
                legend.position = "bottom", legend.key.width = unit(1.5,'cm'),
                axis.title.y = element_blank()) + 
  labs(x = "Longitude", fill = "Celsius")

rleg = g_legend(rp)  

grid.arrange(arrangeGrob(arrangeGrob(pdag,pb+theme(legend.position = "none"),nrow = 1),pleg, 
                         nrow = 2, heights = c(1,0.1)),
             arrangeGrob(rp+theme(legend.position = "none"), rleg, nrow = 2, heights = c(1,0.1)),
             nrow = 1, widths = c(1,0.5), left = "Latitude")



plot_mean2d_map_viridis(xs_test[h1,], fit$pred_ub - fit$pred_lb,xcols = c(1,2), 
                        viridis_opt = "inferno",
                        scale_limit = c(0,15), title = "ERA5",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(-90,90),
                        lon_bnds = c(-180,180)
)



plot_mean2d_map_viridis(xs_test[h1,], cp_ind,xcols = c(1,2), 
                        viridis_opt = "inferno",
                        scale_limit = c(-1,1), title = "ERA5",
                        maps_list = list(data.frame(world),data.frame(states)),
                        maps_cols = c("grey30","grey30"),
                        lat_bnds = c(-90,90),
                        lon_bnds = c(-180,180)
)


cp_ind = ifelse(fit$pred_ub >= ms$f_test[h2,3] & fit$pred_lb <= ms$f_test[h2,3],1,0)
cp = mean(fit$pred_ub >= ms$f_test[h2,3] & fit$pred_lb <= ms$f_test[h2,3])
sqrt(mean(resid^2))


# Error Standard Deviation
par(mfrow = c(1,2))
hist(unlist(sfit), xlab = TeX("Posterior Draws of $\\sigma$"), 
     main = TeX("Posterior Draws of $\\sigma$"))
hist(fit$pred_ub - fit$pred_lb, xlab = "Credible Interval Width",xlim = c(0,15), 
     main = "Credible Interval Width",font.main = 1)
#plot(unlist(sfit), type = 'l')

# Plot smooth pahs
xlon = unique(ms$x_test[,1])[40] #77
xh = h2[which(ms$x_test[h2,1] == xlon)]
xhm = xh - length(h2)

# Indivudal plots
plot(ms$x_test[xh,2], fit$pred_mean[xhm], type = 'l')
lines(ms$x_test[xh,2], fit$pred_lb[xhm], type = 'l', col = "red")
lines(ms$x_test[xh,2], fit$pred_ub[xhm], type = 'l', col = "red")
points(ms$x_test[xh,2],ms$f_test[xh,3], cex = 0.3)

ms$x_test[which.max(resid),1]


# First Draw at .....
xlon = unique(ms$x_test[,1])[1] #77
xh = h2[which(ms$x_test[h2,1] == xlon)]
xhm = xh - length(h2)

pd1 = plot_mean1d(ms$x_test[xh,2], fit$pred_mean[xhm], plb = fit$pred_lb[xhm],
                  pub = fit$pred_ub[xhm], amean = ms$f_test[xh,3], y_lim = c(-60,40), 
                  line_type = c("solid","solid"), colors = c("grey30","red"))
pd1 = pd1 + theme(axis.text = element_text(size = 13),axis.title = element_text(size = 13),
                legend.position = "bottom", legend.key.width = unit(1.5,'cm'),
                axis.title.y = element_blank()) + 
  labs(x = "Latitude", fill = "Celsius")


# Second Draw at .....
xlon = 20 #unique(ms$x_test[,1])[1] #77
xh = h2[which(ms$x_test[h2,1] == xlon)]
xhm = xh - length(h2)

pd2 = plot_mean1d(ms$x_test[xh,2], fit$pred_mean[xhm], plb = fit$pred_lb[xhm],
                  pub = fit$pred_ub[xhm], amean = ms$f_test[xh,3], y_lim = c(-60,40),
                  line_type = c("solid","solid"), colors = c("grey30","red"))
pd2 = pd2 + theme(axis.text = element_text(size = 13),axis.title = element_text(size = 13),
                  legend.position = "bottom", legend.key.width = unit(1.5,'cm'),
                  axis.title.y = element_blank()) + 
  labs(x = "Latitude", fill = "Celsius")


# Third Draw at .....
xlon = 77.5 #unique(ms$x_test[,1])[1] #77
xh = h2[which(ms$x_test[h2,1] == xlon)]
xhm = xh - length(h2)

pd3 = plot_mean1d(ms$x_test[xh,2], fit$pred_mean[xhm], plb = fit$pred_lb[xhm],
                  pub = fit$pred_ub[xhm], amean = ms$f_test[xh,3], y_lim = c(-60,40), 
                  line_type = c("solid","solid"), colors = c("grey30","red"))
pd3 = pd3 + theme(axis.text = element_text(size = 13),axis.title = element_text(size = 13),
                  legend.position = "bottom", legend.key.width = unit(1.5,'cm'),
                  axis.title.y = element_blank()) + 
  labs(x = "Latitude", fill = "Celsius")


grid.arrange(pd1 + theme(legend.position = "none") + labs(title = TeX("$0^\\circ$ Longitude")),
             pd2 + theme(legend.position = "none") + labs(title = TeX("$20^\\circ$ Longitude")),
             pd3 + theme(legend.position = "none") + labs(title = TeX("$77.5^\\circ$ Longitude")),
             nrow = 1, left = textGrob("Temperature", gp=gpar(fontsize=18),rot = 90))




