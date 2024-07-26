import os
import numpy as np
import pandas as pd

os.getcwd()
os.chdir("/home/johnyannotty/Documents/openbt/")
#sys.path.append("/home/johnyannotty/Documents/openbt/openbtmixing_pypkg/")
#sys.path.append("/home/johnyannotty/Documents/openbt/src/")


from openbtmixing_pypkg import Openbtmix

filedir = "/home/johnyannotty/Documents/Dissertation/results/2d_functions/"

f_train = pd.read_csv(filedir + "sincos_4k_ftrain_rpath_05_07_24.txt").values
x_train = pd.read_csv(filedir + "sincos_4k_xtrain_rpath_05_07_24.txt").values
y_train = pd.read_csv(filedir + "sincos_4k_ytrain_rpath_05_07_24.txt").values

f_test = pd.read_csv(filedir + "sincos_4k_ftest_rpath_05_07_24.txt").values
x_test = pd.read_csv(filedir + "sincos_4k_xtest_rpath_05_07_24.txt").values
f0_test = pd.read_csv(filedir + "sincos_4k_f0test_rpath_05_07_24.txt").values

f_test.shape

mix = Openbtmix(local_openbt_path = os.getcwd() + "/src/")

nu = 20
rho = 1
minse = [np.min((f_train[:,i].reshape(100,1)-y_train)**2)  for i in range(f_train.shape[1])]
sig2_hat = np.max(minse)
lam = sig2_hat*(nu+2)/nu
q0 = 4

mix.set_prior(
    k=1,
    ntree=10,
    nu=20,
    sighat=np.sqrt(sig2_hat),
    inform_prior=False,
    rpath = True,
    a1 = 2,
    a2 = 10,
    q = 4,
    power = 1.0)

mix.local_openbt_path
fit = mix.train(
        x_train=x_train,
        y_train=y_train,
        f_train=f_train,
        ndpost=10000,
        nadapt=2000,
        nskip=2000,
        adaptevery=500,
        minnumbot=3,
        tc = 4,
        numcut = 500)

fitp=mix.predict(x_test, f_test, ci=0.95)

fitp.keys()
np.sqrt(np.mean((fitp["pred"]["mean"].reshape(625,1) - f0_test)**2)) 

fitw=mix.predict_weights(x_test, ci=0.95)


import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
wmean = fitw["wts"]["mean"]

# The posterior mean weight functions
cmap_hot = plt.get_cmap('hot')
w1 = wmean.transpose()[0]
w2 = wmean.transpose()[1]
w3 = wmean.transpose()[1]
w4 = wmean.transpose()[1]

w1_mean = wmean.transpose()[0]
w1_mean = w1_mean.reshape((25,25)).transpose()

w2_mean = wmean.transpose()[1]
w2_mean = w2_mean.reshape((25,25)).transpose()

w3_mean = wmean.transpose()[2]
w3_mean = w3_mean.reshape((25,25)).transpose()

w4_mean = wmean.transpose()[3]
w4_mean = w4_mean.reshape((25,25)).transpose()


# Posterior Mean resiudals
n_test = 625
cmap_rb = plt.get_cmap("RdBu")
fig, ax = plt.subplots(1,4, figsize = (24,8))

pcm0 = ax[0].pcolormesh(w1_mean,cmap = cmap_hot, vmin = -0.05, vmax = 1.05)
ax[0].set_title("Posterior Mean of $w_1(x)$\n", size = 22)
#ax[1].set(xlabel = "$x_1$", ylabel = "$x_2")
ax[0].set_xlabel("$x_1$",size=18)
ax[0].set_ylabel("$x_2$",size=18)
ax[0].xaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[0].xaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[0].yaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[0].yaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[0].tick_params(labelsize = 13)
fig.colorbar(pcm0,ax = ax[0],location="bottom")

pcm1 = ax[1].pcolormesh(w2_mean,cmap = cmap_hot, vmin = -0.05, vmax = 1.05)
ax[1].set_title("Posterior Mean of $w_1(x)$\n", size = 22)
#ax[1].set(xlabel = "$x_1$", ylabel = "$x_2")
ax[1].set_xlabel("$x_1$",size=18)
ax[1].set_ylabel("$x_2$",size=18)
ax[1].xaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[1].xaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[1].yaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[1].yaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[1].tick_params(labelsize = 13)
fig.colorbar(pcm1,ax = ax[1],location="bottom")

pcm2 = ax[2].pcolormesh(w3_mean,cmap = cmap_hot, vmin = -0.05, vmax = 1.05)
ax[2].set_title("Posterior Mean of $w_2(x)$\n", size = 22)
ax[2].set_xlabel("$x_1$",size=18)
ax[2].set_ylabel("$x_2$",size=18)
ax[2].xaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[2].xaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[2].yaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[2].yaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[2].tick_params(labelsize = 13)
fig.colorbar(pcm2,ax = ax[2],location="bottom")

pcm3 = ax[3].pcolormesh(w4_mean,cmap = cmap_hot, vmin = -0.05, vmax = 1.05)
ax[3].set_title("Posterior Mean of $w_1(x)$\n", size = 22)
#ax[1].set(xlabel = "$x_1$", ylabel = "$x_2")
ax[3].set_xlabel("$x_1$",size=18)
ax[3].set_ylabel("$x_2$",size=18)
ax[3].xaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[3].xaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[3].yaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[3].yaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[3].tick_params(labelsize = 13)
fig.colorbar(pcm3,ax = ax[3],location="bottom")
plt.show()