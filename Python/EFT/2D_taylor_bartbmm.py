import numpy as np
from scipy.special import factorial

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


os.getcwd()
os.chdir("/home/johnyannotty/Documents/openbt/")

from openbtmixing_pypkg import Openbtmix


class polynomal_model():
    '''
        Polynomial models class. Used to define a function of the form

        .. math::
            f(x) = c(x-a)^p + b
    '''

    def __init__(self, a=0, b=0, c=1, p=1):
        '''
            Parameters:
            -----------
            :param float a: center parameter.
            :param float b: shift parameter.
            :param float c: scale parameter.
            :param float p: power parameter.

        '''
        self.a = a
        self.b = b
        self.c = c
        self.p = p

    def evaluate(self, x):
        '''
            Evaluate the polynomial at a grid of x's. The standard deviation
            output is set to 1 by default.

            Parameters:
            -----------
            :param np.ndarray x: design matrix.

            Returns:
            --------
            :returns: mean and standard deviation of the
                model at the x grid points.
            :rtype: np.ndarray, np.ndarray
            :return values: mean predictions.
            :return values: standard deviation of the predictions.

        '''

        if isinstance(x, list):
            x = np.array(x)
        m = self.c * (x - self.a)**self.p + self.b
        if len(m.shape) == 1:
            m = m.reshape(m.shape[0], 1)
        s = np.array([1] * x.shape[0]).reshape(m.shape[0], 1)

        return m, s


# Taylor Expansions
class sin_exp():
    '''
        Sine Taylor series expansion model class.
    '''

    def __init__(self, k, x0):
        '''
        Parameters:
        -----------
        :param int k: the degree of the expansion.
        :param float x0: the center of the expansion.

        '''

        self.k = k
        self.x0 = x0

    def evaluate(self, x):
        '''
        Evaluate the Taylor Series at a grid of x's. The standard deviation
        output is set to 1 by default.

        Parameters:
        -----------
        :param np.ndarray x: design matrix.

        Returns:
        --------
        :returns: mean and standard deviation of the model
            at the x grid points.
        :rtype: np.ndarray, np.ndarray
        :return values: mean predictions.
        :return values: standard deviation of the predictions.

        '''

        # Check type of x
        if isinstance(x, list):
            x = np.array(x)

        # Get degree list for polynomial expansion
        deg = np.linspace(0, self.k, self.k + 1)

        # Get every 4th term (derivative repaets every 4 terms)
        h0 = deg[deg % 4 == 0]
        h1 = deg[deg % 4 == 1]
        h2 = deg[deg % 4 == 2]
        h3 = deg[deg % 4 == 3]

        # Compute taylor series:
        xc = x - self.x0
        xc0 = xc.repeat(h0.shape[0]).reshape(x.shape[0], h0.shape[0])
        ts = np.sum(
            np.sin(
                self.x0) *
            np.power(
                xc0,
                h0) /
            factorial(h0),
            axis=1)

        xc1 = xc.repeat(h1.shape[0]).reshape(x.shape[0], h1.shape[0])
        ts = ts + np.sum(np.cos(self.x0) * np.power(xc1,
                         h1) / factorial(h1), axis=1)

        xc2 = xc.repeat(h2.shape[0]).reshape(x.shape[0], h2.shape[0])
        ts = ts + np.sum(-np.sin(self.x0) *
                         np.power(xc2, h2) / factorial(h2), axis=1)

        xc3 = xc.repeat(h3.shape[0]).reshape(x.shape[0], h3.shape[0])
        ts = ts + np.sum(-np.cos(self.x0) *
                         np.power(xc3, h3) / factorial(h3), axis=1)

        if len(ts.shape) == 1:
            ts = ts.reshape(ts.shape[0], 1)

        s = np.array([1] * x.shape[0]).reshape(ts.shape[0], 1)

        return ts, s


class cos_exp():
    '''
        Cosine Taylor series expansion model class.
    '''

    def __init__(self, k, x0):
        '''
        Parameters:
        -----------
        :param int k: the degree of the expansion.
        :param float x0: the center of the expansion.

        '''
        self.k = k
        self.x0 = x0

    def evaluate(self, x):
        '''
        Evaluate the Taylor series at a grid of x's. The standard deviation
        output is set to 1 by default.

        Parameters:
        -----------
        :param np.ndarray x: design matrix.

        Returns:
        --------
        :returns: mean and standard deviation of the model at the
            x grid points.
        :rtype: np.ndarray, np.ndarray
        :return values: mean predictions.
        :return values: standard deviation of the predictions.

        '''

        # Check type of x
        if isinstance(x, list):
            x = np.array(x)

        # Get degree list for polynomial expansion
        deg = np.linspace(0, self.k, self.k + 1)

        # Get every 4th term (derivative repaets every 4 terms)
        h0 = deg[deg % 4 == 0]
        h1 = deg[deg % 4 == 1]
        h2 = deg[deg % 4 == 2]
        h3 = deg[deg % 4 == 3]

        # Compute taylor series:
        xc = x - self.x0
        xc0 = xc.repeat(h0.shape[0]).reshape(x.shape[0], h0.shape[0])
        ts = np.sum(
            np.cos(
                self.x0) *
            np.power(
                xc0,
                h0) /
            factorial(h0),
            axis=1)

        xc1 = xc.repeat(h1.shape[0]).reshape(x.shape[0], h1.shape[0])
        ts = ts + np.sum(-np.sin(self.x0) *
                         np.power(xc1, h1) / factorial(h1), axis=1)

        xc2 = xc.repeat(h2.shape[0]).reshape(x.shape[0], h2.shape[0])
        ts = ts + np.sum(-np.cos(self.x0) *
                         np.power(xc2, h2) / factorial(h2), axis=1)

        xc3 = xc.repeat(h3.shape[0]).reshape(x.shape[0], h3.shape[0])
        ts = ts + np.sum(np.sin(self.x0) * np.power(xc3,
                         h3) / factorial(h3), axis=1)

        if len(ts.shape) == 1:
            ts = ts.reshape(ts.shape[0], 1)

        s = np.array([1] * x.shape[0]).reshape(ts.shape[0], 1)

        return ts, s


class sin_cos_exp():
    '''
    Taylor series expansion of

    .. math::
            f(x) = \\sin(x_1) + \\cos(x_2)

    '''

    def __init__(self, ks, kc, xs, xc):
        '''
        Parameters:
        -----------
        :param int ks: the degree of the sine expansion.
        :param int kc: the degree of the cosine expansion.
        :param float xs: the center of the sine expansion.
        :param float xc: the center of the cosine expansion.

        '''

        self.ks = ks
        self.xs = xs

        self.kc = kc
        self.xc = xc

    def evaluate(self, x):
        '''
        Evaluate the model at a grid of x's. The standard deviation
        output is set to 1 by default.

        Parameters:
        -----------
        :param np.ndarray x: design matrix.

        Returns:
        --------
        :returns: mean and standard deviation of the model at the x
            grid points.
        :rtype: np.ndarray, np.ndarray
        :return values: mean predictions.
        :return values: standard deviation of the predictions.

        '''

        # Check type of x
        if isinstance(x, list):
            x = np.array(x)

        # Sine Part
        # Get degree list for polynomial expansion
        deg = np.linspace(0, self.ks, self.ks + 1)

        # Get every 4th term (derivative repaets every 4 terms)
        h0 = deg[deg % 4 == 0]
        h1 = deg[deg % 4 == 1]
        h2 = deg[deg % 4 == 2]
        h3 = deg[deg % 4 == 3]

        # Compute taylor series:
        xctr = x.transpose()[0] - self.xs
        xc0 = xctr.repeat(h0.shape[0]).reshape(x.shape[0], h0.shape[0])
        tss = np.sum(
            np.sin(
                self.xs) *
            np.power(
                xc0,
                h0) /
            factorial(h0),
            axis=1)

        xc1 = xctr.repeat(h1.shape[0]).reshape(x.shape[0], h1.shape[0])
        tss = tss + np.sum(np.cos(self.xs) * np.power(xc1,
                           h1) / factorial(h1), axis=1)

        xc2 = xctr.repeat(h2.shape[0]).reshape(x.shape[0], h2.shape[0])
        tss = tss + np.sum(-np.sin(self.xs) *
                           np.power(xc2, h2) / factorial(h2), axis=1)

        xc3 = xctr.repeat(h3.shape[0]).reshape(x.shape[0], h3.shape[0])
        tss = tss + np.sum(-np.cos(self.xs) *
                           np.power(xc3, h3) / factorial(h3), axis=1)

        # Cosine Part
        # Get degree list for polynomial expansion
        deg = np.linspace(0, self.kc, self.kc + 1)

        # Get every 4th term (derivative repaets every 4 terms)
        h0 = deg[deg % 4 == 0]
        h1 = deg[deg % 4 == 1]
        h2 = deg[deg % 4 == 2]
        h3 = deg[deg % 4 == 3]

        # Compute taylor series:
        xctr = x.transpose()[1] - self.xc
        xc0 = xctr.repeat(h0.shape[0]).reshape(x.shape[0], h0.shape[0])
        tsc = np.sum(
            np.cos(
                self.xc) *
            np.power(
                xc0,
                h0) /
            factorial(h0),
            axis=1)

        xc1 = xctr.repeat(h1.shape[0]).reshape(x.shape[0], h1.shape[0])
        tsc = tsc + np.sum(-np.sin(self.xc) *
                           np.power(xc1, h1) / factorial(h1), axis=1)

        xc2 = xctr.repeat(h2.shape[0]).reshape(x.shape[0], h2.shape[0])
        tsc = tsc + np.sum(-np.cos(self.xc) *
                           np.power(xc2, h2) / factorial(h2), axis=1)

        xc3 = xctr.repeat(h3.shape[0]).reshape(x.shape[0], h3.shape[0])
        tsc = tsc + np.sum(np.sin(self.xc) * np.power(xc3,
                           h3) / factorial(h3), axis=1)

        # Add the sine and cosine parts
        ts = tss + tsc

        if len(ts.shape) == 1:
            ts = ts.reshape(ts.shape[0], 1)

        s = np.array([1] * x.shape[0]).reshape(ts.shape[0], 1)

        return ts, s


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


# Generate Data
# Pick the grid dimensions
nx1 = 10; nx2 = 8
n_train = nx1*nx2
#x_train = grid_2d_design(nx1,nx2,[-np.pi,-np.pi],[np.pi,np.pi])

# Get test data
n_test = 30
x1_test = np.outer(np.linspace(-3, 3, n_test), np.ones(n_test))
x2_test = x1_test.copy().transpose()
x_test = np.array([x1_test.reshape(x1_test.size,),
                    x2_test.reshape(x1_test.size,)]).transpose()


# Generate observational data with true standard deviation of 0.1
#f0_train = np.sin(x_train.transpose()[0]) + np.cos(x_train.transpose()[1])
#y_train = f0_train + np.random.normal(0,0.1,n_train)


y_train = np.loadtxt("/home/johnyannotty/Documents/Taweret/test/bart_bmm_test_data/2d_y_train.txt").reshape(80, 1)
x_train = np.loadtxt("/home/johnyannotty/Documents/Taweret/test/bart_bmm_test_data/2d_x_train.txt").reshape(80, 2)
x_train = x_train.reshape(2, 80).transpose()
f0_train = np.sin(x_train.transpose()[0]) + np.cos(x_train.transpose()[1])


# Define the model set
f1 = sin_cos_exp(7,10,np.pi,np.pi) # 7th order sin(x1) + 10th order cos(x2)
f2 = sin_cos_exp(13,6,-np.pi,-np.pi) # 13th order sin(x1) + 6th order cos(x2)
model_dict = {'model1':f1, 'model2':f2}

f_train = np.concatenate([f1.evaluate(x_train)[0],f2.evaluate(x_train)[0]], axis = 1)
f_test = np.concatenate([f1.evaluate(x_test)[0],f2.evaluate(x_test)[0]], axis = 1)
f0_test = (np.sin(x1_test) + np.cos(x2_test))

# Fit the BMM Model
# Initialize the Trees class instance
mix = Openbtmix(local_openbt_path = os.getcwd() + "/src/")

nu = 10
sig2_hat = 0.01**2/(7/5)

mix.set_prior(
    k=2.5,
    ntree=30,
    nu=nu,
    sighat=np.sqrt(sig2_hat),
    inform_prior=False,
    rpath = False,
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
pmean = fitp["pred"]["mean"]

np.sqrt(np.mean((fitp["pred"]["mean"].reshape(n_test**2,1) - f0_test.reshape(n_test**2,1))**2)) 

fitw=mix.predict_weights(x_test, ci=0.95)


# Define color map
cmap = plt.get_cmap('viridis')

# Heat map comparing the surfaces
fig, ax = plt.subplots(1,2, figsize = (12,5))

pcm1 = ax[0].pcolormesh(f0_test.transpose(),cmap = cmap, vmin = -2.5, vmax = 2.5,shading='gouraud')
ax[0].set_title("True System", size = 16)
ax[0].set(xlabel = "$x_1$", ylabel = "$x_2$")
ax[0].xaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[0].xaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[0].yaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[0].yaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
fig.colorbar(pcm1,ax = ax[0])


# Predicted mean
pcm2 = ax[1].pcolormesh(pmean.reshape(x1_test.shape).transpose(),cmap = cmap, vmin = -2.5, vmax = 2.5,shading='gouraud')
ax[1].set_title("Posterior Mean Prediction", size = 16)
ax[1].set(xlabel = "$x_1$", ylabel = "$x_2$")
ax[1].xaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[1].xaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[1].yaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[1].yaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))

fig.colorbar(pcm2,ax = ax[1])
fig.suptitle("Figure 8: True versus Predicted System", size = 18)
plt.show()


# Create Figure 10: The posterior mean weight functions
cmap_hot = plt.get_cmap('hot')
wmean = fitw["wts"]["mean"]
w1 = wmean.transpose()[0]
w2 = wmean.transpose()[1]

w1_mean = wmean.transpose()[0]
w1_mean = w1_mean.reshape(x1_test.shape).transpose()

w2_mean = wmean.transpose()[1]
w2_mean = w2_mean.reshape(x1_test.shape).transpose()

w_sum = w1_mean + w2_mean 

fig, ax = plt.subplots(1,2, figsize = (12,5))
pcm0 = ax[0].pcolormesh(w1_mean,cmap = cmap_hot, vmin = -0.05, vmax = 1.05,shading='gouraud')
ax[0].set_title("Posterior Mean of $w_1(x)$", size = 14)
ax[0].set(xlabel = "$x_1$", ylabel = "$x_2")
ax[0].xaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[0].xaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[0].yaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[0].yaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
fig.colorbar(pcm0,ax = ax[0])

pcm1 = ax[1].pcolormesh(w2_mean,cmap = cmap_hot, vmin = -0.05, vmax = 1.05,shading='gouraud')
ax[1].set_title("Posterior Mean of $w_2(x)$", size = 14)
ax[1].set(xlabel = "$x_1$", ylabel = "$x_2$")
ax[1].xaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[1].xaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
ax[1].yaxis.set_major_locator(ticker.FixedLocator(np.round(np.linspace(0, n_test, 6),3)))
ax[1].yaxis.set_major_formatter(ticker.FixedFormatter(np.round(np.linspace(-np.pi, np.pi, 6),3)))
fig.colorbar(pcm1,ax = ax[1])
fig.suptitle("Figure 10: Posterior Mean Weight Functions", size = 16)
plt.show()


pm = np.loadtxt("/home/johnyannotty/Documents/Taweret/test/bart_bmm_test_data/2d_pmean.txt")
perr = np.mean(np.abs(pmean - pm))

wm = np.loadtxt("/home/johnyannotty/Documents/Taweret/test/bart_bmm_test_data/2d_wmean.txt")
werr = np.mean(np.abs(wmean - wm))