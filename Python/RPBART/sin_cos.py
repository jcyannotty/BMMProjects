import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

os.getcwd()
os.chdir("/home/johnyannotty/Documents/openbt/")
#sys.path.append("/home/johnyannotty/Documents/openbt/openbtmixing_pypkg/")
#sys.path.append("/home/johnyannotty/Documents/openbt/src/")

from openbtmixing_pypkg import OpenbtRpath

## Functions for design points 
# n1 = number of bins in the x1 dimension
# n2 = number of bins in the x2 dimension
# n = n1*n2 is the total training size
def grid_2d_design(n1,n2, xmin = [-1,-1], xmax = [1,1]):
  # Generate n uniform rvs
  n = n1*n2
  ux = np.random.uniform(0,1,n)
  uy = np.random.uniform(0,1,n)

  # Dimensions for each rectangle
  x1_len = (xmax[0] - xmin[0])/n1
  x2_len = (xmax[1] - xmin[1])/n2
  xgrid = [[x, y] for x in range(n1) for y in range(n2)]
  xgrid = np.array(xgrid).transpose()

  # Get points
  x1 = ux*x1_len + x1_len*xgrid[0] + xmin[0]
  x2 = uy*x2_len + x2_len*xgrid[1] + xmin[1]

  # Join data
  xdata = np.array([x1,x2]).transpose()
  return xdata

n_train = 100 
n_test = 25
xoff = 0.10
xmin = -np.pi
xmax = np.pi
s = 0.1

# Generate data
x_train = grid_2d_design(n1=10,n2=10, xmin = [xmin,xmin], xmax = [xmax,xmax])

# Get test data
x1_test = np.outer(np.linspace(-np.pi, np.pi, n_test), np.ones(n_test))
x2_test = x1_test.copy().transpose()
f0_test = (np.sin(x1_test) + np.cos(x2_test))
x_test = np.array([x1_test.reshape(x1_test.size,),x2_test.reshape(x2_test.size,)]).transpose()


# Generate True function
np.random.seed(99)
f0_train = np.sin(x_train[:,0]) + np.cos(x_train[:,1])
y_train = f0_train + np.random.normal(0,s,n_train)

rpb = OpenbtRpath(local_openbt_path = os.getcwd() + "/src/")

nu = 40
rho = 1
sig2_hat = 0.1
lam = sig2_hat#*(nu+2)/nu
q0 = 4

rpb.set_prior(
    k=1,
    ntree=20,
    nu=nu,
    sighat=np.sqrt(sig2_hat),
    rpath = True,
    a1 = 2,
    a2 = 5,
    q = 4,
    power = 1.0,
    ymin = np.min(y_train),
    ymax = np.max(y_train))


fit = rpb.train(
        x_train=x_train,
        y_train=y_train,
        ndpost=5000,
        nadapt=5000,
        nskip=2000,
        adaptevery=500,
        minnumbot=3,
        tc = 4,
        numcut = 300)


fitp=rpb.predict(x_test, ci=0.95)

# rpb.sp_train
# rpb.sp_pred

fitp["sigma"]["mean"][0]

np.sqrt(np.mean((fitp["pred"]["mean"].reshape(625,1) - f0_test.reshape(625,1))**2)) 
pmean = fitp["pred"]["mean"]


cmap = plt.get_cmap('turbo')
fig, ax = plt.subplots(1,2, figsize = (12,5))

pcm1 = ax[0].pcolormesh(f0_test.transpose(),cmap = cmap, vmin = -2.0, vmax = 2.0,shading='gouraud')
ax[0].set_title("True System", size = 16)
ax[0].set(xlabel = "$x_1$", ylabel = "$x_2$")
ax[0].xaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[0].xaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[0].yaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[0].yaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
fig.colorbar(pcm1,ax = ax[0])


# Predicted mean
pcm2 = ax[1].pcolormesh(pmean.reshape(x1_test.shape).transpose(),cmap = cmap, vmin = -2.0, vmax = 2.0,shading='gouraud')
ax[1].set_title("Posterior Mean Prediction", size = 16)
ax[1].set(xlabel = "$x_1$", ylabel = "$x_2$")
ax[1].xaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[1].xaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[1].yaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[1].yaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))

fig.colorbar(pcm2,ax = ax[1])
fig.suptitle("Figure 8: True versus Predicted System", size = 18)
plt.show()




#pd.DataFrame(x_train).to_csv("/home/johnyannotty/Downloads/tempx", header = False, index = False)
#pd.DataFrame(y_train).to_csv("/home/johnyannotty/Downloads/tempy", header = False, index = False)


from scipy.stats import spearmanr
