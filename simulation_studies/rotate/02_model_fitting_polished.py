

import pandas as pd
import os

os.chdir('/users/ywang/Hansen_projects/mNSF/revision/simulation')
dtp = "float32"
#from mNSF import process_multiSample
import mNSF
import numpy as np
import matplotlib.pyplot as plt
#from os import path
from pandas import get_dummies
from anndata import AnnData
from scanpy import pp
from os import path
from mNSF.NSF import preprocess,misc,visualize
#from utils import misc,preprocess,visualize
#from mNSF.NSF import misc,preprocess,visualize
from mNSF import process_multiSample
from scanpy import read_h5ad
from mNSF import training_multiSample
from scanpy import pp

########################################################################
mpth = path.join("models")
misc.mkdir_p(mpth)
pp = path.join(mpth,"pp",str(2))#list_fit[0].generate_pickle_path("constant",base=mpth)
misc.mkdir_p(pp)


########################################################################
L=4

########################################################################3
################### step 0  Data loading
########################################################################
nsample=2

list_D=list()
list_X=list()
for ksample in range(0,nsample):
	ad=read_h5ad(path.join("ggblocks_lr_"+str(ksample+1)+".h5ad"))
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
process_multiSample.save_object(list_fit, path.join('sim_rotated_list_fit.pkl') )




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
Wdf.to_csv(path.join("sim_rotated_loadings_spde.csv"))



## save the factors
#inpf12 = process_multiSample.interpret_npf_v3(list_fit,list_X,S=100,lda_mode=False)
Factors = inpf12["factors"][:,0:L]

for k in range(0,nsample):
	indices=list_sampleID[k]
	indices=indices.astype(int)
	Factors_df = pd.DataFrame(Factors[indices,:]) 
	Factors_df.to_csv(path.join("sim_rotated_factors_sample"+str(k+1)+".csv"))


## plot the factors
hmkw = {"figsize":(7,.9),"bgcol":"white","subplot_space":0.1,"marker":"s","s":10}
for ksample in range(0,nsample):
	indices=list_sampleID[k]
	indices=indices.astype(int)
	fig,axes=visualize.multiheatmap(list_D[ksample]["X"],Factors[indices,:], (1,L), cmap="Blues", **hmkw)
	fig.savefig(path.join("sim_rotated_factors_sample"+str(ksample+1)+".png"),bbox_inches='tight')






	

