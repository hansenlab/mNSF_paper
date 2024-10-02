


dir_mNSF_functions='/users/ywang/Hansen_projects/scRNA/mNSFH_2023_06_14/mNSF-main'
dir_output='/dcs04/hansen/data/ywang/ST/data_10X_ST/mouse_Sagittal_spaceRanger1_1_0/out/'

########################################################################
########################################################################
import sys
sys.path.append(dir_mNSF_functions)

from scanpy import read_h5ad


import mNSF

from mNSF import process_multiSample

from mNSF.NSF import preprocess
from mNSF.NSF import misc
#from mNSF.NSF import visualize
#from mNSF import training_multiSample
from mNSF import training_multiSample
from mNSF import process_multiSample
from mNSF.NSF import visualize

#from tensorflow.data import Dataset

from os import path
#import pandas
import os
import numpy as np
import tensorflow as tf
import pandas as pd
import sys 
import pickle




########################################################################

os.chdir('/dcs04/legacy-dcs01-hansen/hansen_lab1/ywang/ST/April8_2022_NSF/nsf-paper-main_9_ori')
pth = "simulations/ggblocks_lr"
dpth = path.join(pth,"data")
plt_pth = path.join(pth,"results")
mpth = path.join(pth,"models")
pp = path.join(mpth,"pp_ss")#list_fit[0].generate_pickle_path("constant",base=mpth)
########################################################################
L=4



########################################################################3
################### step 0  Data loading
########################################################################
nsample=2

list_D=list()
list_X=list()
for ksample in range(0,nsample):
	ad=read_h5ad(path.join(dpth,"ggblocks_lr_"+str(ksample+1)+"_SSlayer.h5ad"))
	D,_ = preprocess.anndata_to_train_val(ad, layer="counts", train_frac=1.0,
                                      flip_yaxis=False)
	D["X"] = D["X"]
	D["Z"] = D["X"]
	list_D.append(D)
	list_X.append(D["X"])


list_Dtrain=process_multiSample.get_listDtrain(list_D)
list_sampleID=process_multiSample.get_listSampleID(list_D)



########################################################################3
################### step 1 initialize model
########################################################################

list_fit = process_multiSample.ini_multiSample(list_D,L)


########################################################################
################### step 2 fit model
########################################################################


list_fit=training_multiSample.train_model_mNSF(list_fit,pp,list_Dtrain,list_D)


# save the fitted model
process_multiSample.save_object(list_fit, path.join(plt_pth,'sim_SS_list_fit.pkl') )




########################################################################3
################### step 3 save and plot results
########################################################################
## save the loadings
#loadings=visualize.get_loadings(list_fit[0])
#DF = pd.DataFrame(loadings) 
#DF.to_csv(("loadings.csv"))
inpf12=process_multiSample.interpret_npf_v3(list_fit,list_X,S=100,lda_mode=False)


W = inpf12["loadings"]
#Wdf=pd.DataFrame(W*inpf12["totals1"][:,None], index=ad.var.index, columns=range(1,L+1))

Wdf=pd.DataFrame(W*inpf12["totalsW"][:,None],  columns=range(1,L+1))
Wdf.to_csv(path.join(plt_pth,"sim_SS_loadings_spde.csv"))



## save the factors
#inpf12 = process_multiSample.interpret_npf_v3(list_fit,list_X,S=100,lda_mode=False)
Factors = inpf12["factors"][:,0:L]

for k in range(0,nsample):
	indices=list_sampleID[k]
	indices=indices.astype(int)
	Factors_df = pd.DataFrame(Factors[indices,:]) 
	Factors_df.to_csv(path.join(plt_pth,"sim_SS_factors_sample"+str(k+1)+".csv"))






## plot the factors
hmkw = {"figsize":(7,.9),"bgcol":"white","subplot_space":0.1,"marker":"s","s":10}
for ksample in range(0,nsample):
	indices=list_sampleID[k]
	indices=indices.astype(int)
	fig,axes=visualize.multiheatmap(list_D[ksample]["X"],Factors[indices,:], (1,L), cmap="Blues", **hmkw)
	fig.savefig(path.join(plt_pth,"sim_SS_factors_sample"+str(ksample+1)+".png"),bbox_inches='tight')






	

