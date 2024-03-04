"""
Name: cmip6_fwls.py
Desc: Apply nn feature weighted linear stacking to climate models from cmip6
"""

import sys
import pandas as pd
import numpy as np
import torch
import pickle
import matplotlib.pyplot as plt
import importlib


#import nns
#importlib.reload(nns)

sys.path.append("/home/johnyannotty/Documents/nnstacking/nnstacking")
from nns import NNS

#-----------------------------------------------------
# Fix up the nns python package
#-----------------------------------------------------
class FWLSNN(NNS):      
    def predict_new(self,x_test, f_test):
        with torch.no_grad():
            self._check_dims(x_test, np.empty((1,1)))
            f_test = torch.from_numpy(f_test)

            self.neural_net.eval()
            nnx = np.array(x_test, dtype='f4')
            nnx = torch.from_numpy(nnx)

            nnpred = np.array(f_test, dtype='f4')
            nnpred = torch.from_numpy(nnpred)

            if self.gpu:
                nnx = nnx.cuda()
                nnpred = nnpred.cuda()

            output = self._ensemblize(nnx, nnpred)

            return output.data.cpu().numpy()

#-----------------------------------------------------
# Data load
#-----------------------------------------------------
filedir = "/home/johnyannotty/Documents/CMIP6_mixing/"
datadir = "Data/World/"
resdir = "Results/World/"
dataname = "ACC_BCC_MIROC_CMCC_CESM2_CNRM_CanESM5_KIOST_W_3M_2014_01_26_24_n45000" 
csvdir = dataname+"_csvs/"

f_train = pd.read_csv(filedir+datadir+csvdir+"f_train.csv").values
y_train = pd.read_csv(filedir+datadir+csvdir+"y_train.csv").values
x_train = pd.read_csv(filedir+datadir+csvdir+"x_train.csv").values
f_test = pd.read_csv(filedir+datadir+csvdir+"f_test.csv").values
y_test = pd.read_csv(filedir+datadir+csvdir+"y_test.csv").values
x_test = pd.read_csv(filedir+datadir+csvdir+"x_test.csv").values

n_train = x_train.shape[0]
n_test = x_test.shape[0]
K = f_train.shape[1]
sn = csvdir.split("_")[0:K]

f_train = f_train.reshape(n_train,1,K)
f_test = f_test.reshape(n_test,1,K)

#-----------------------------------------------------
# Model Fit
#-----------------------------------------------------
method = "UNNS"
nlayers = 3
hsz = 100
delta = False

nnsmodel = FWLSNN(ensemble_method = method, nworkers = 4, num_layers = nlayers,
              hidden_size = hsz, es = True, es_validation_set = 5000, nepoch = 2, gpu = False,
              ensemble_addition = delta, estimators = sn)

nnsmodel.fit(x_train, y_train, predictions = f_train, x_val = x_test[0:10,:], y_val = y_test[0:10,:], f_val = f_test[0:10,:,:])



#nnsmodel.__dict__.keys()
#nnsmodel.__dict__["best_loss_val"]
#nnsmodel.__dict__["loss_history_train"][0:5]
#nnsmodel.__dict__["loss_history_validation"][0:5]

nns_pred = nnsmodel.predict_new(x_test, f_test)
nns_wts = nnsmodel.get_weights(x_test)
np.sqrt(np.mean((nns_pred - y_test)**2))

nit = len(nnsmodel.__dict__["loss_history_validation"])
nnsmodel.__dict__["loss_history_train"][(nit-5):]
nnsmodel.__dict__["loss_history_validation"][(nit-5):]


#-----------------------------------------------------
# Save Data
#-----------------------------------------------------
out_dict = {}
out_dict["loss_history_train"] = nnsmodel.__dict__["loss_history_train"]
out_dict["loss_history_val"] = nnsmodel.__dict__["loss_history_validation"]
out_dict["pred"] = nns_pred
out_dict["wts"] = nns_wts
out_dict["sn"] = sn
out_dict["hp"] = {"nlayers":nlayers,"hsz":hsz, "method":method, "delta":delta}

outfile = filedir + resdir + dataname + "_nnfwls_tp1.pickle"
with open(outfile, 'wb') as pickle_file:
    pickle.dump(out_dict,pickle_file)



#A = np.array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
#A.reshape(4,1,3)[1][0][1]
