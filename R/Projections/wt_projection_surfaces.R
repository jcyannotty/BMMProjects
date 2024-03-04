#------------------------------------------------
# Projection Geometry
#------------------------------------------------
source("/home/johnyannotty/Documents/BayesToolBox/bayestb/Optimization/optimization_helpers.R")
source("/home/johnyannotty/Documents/openbt/R/polynomials.R")

library(lattice)
n = 15
w = seq(-2,2,length = n)
w1 = rep(w, each = n)
w2 = rep(w, n)
wgrid = cbind(w1,w2)

simplex_l2(c(0,2))

p1_l2 = apply(wgrid,1,function(x) simplex_l2(x)[1]) 
wireframe(p1_l2 ~ w1 * w2,scales=list(arrows=FALSE, axis=list(text=list(cex=2))),
          xlab = list(TeX("$w_1$"),rot = 30),
          ylab = list(TeX("$w_2$"),rot = -30),zlab = TeX("$\\alpha_1$"),
          main = "Euclidean Weight 1")


p1_sm1 = softmax(wgrid)
wireframe(p1_sm1[,1] ~ w1 * w2,scales=list(arrows=FALSE, axis=list(text=list(cex=2))),
          xlab = list(TeX("$w_1$"),rot = 30),
          ylab = list(TeX("$w_2$"),rot = -30),zlab = TeX("$\\alpha_1$"),
          main = "Softmax (T=1) Weight 1")

p1_sm05 = softmax(wgrid/0.5)
wireframe(p1_sm05[,1] ~ w1 * w2,scales=list(arrows=FALSE, axis=list(text=list(cex=2))),
          xlab = list(TeX("$w_1$"),rot = 30),
          ylab = list(TeX("$w_2$"),rot = -30),zlab = TeX("$\\alpha_1$"),
          main = "Softmax (T=0.5) Weight 1")

p1_sm025 = softmax(wgrid/0.25)
wireframe(p1_sm025[,1] ~ w1 * w2,scales=list(arrows=FALSE, axis=list(text=list(cex=2))),
          xlab = list(TeX("$w_1$"),rot = 30),
          ylab = list(TeX("$w_2$"),rot = -30),zlab = TeX("$\\alpha_1$"),
          main = "Softmax (T=0.25) Weight 1")


data
