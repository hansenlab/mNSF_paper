

########################################################################
dir_mNSF_functions='/users/ywang/Hansen_projects/scRNA/mNSF_2023_10_20/'
#dir_mNSF_functions='/users/ywang/Hansen_projects/scRNA/mNSF_2023_10_20/mNSF-main'
dir_output="//dcs04/hansen/data/ywang/revision/SlideseqCerebellum/"

########################################################################
########################################################################
import sys
sys.path.append(dir_mNSF_functions)

#from scanpy import read_h5ad

import random
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

from scanpy import pp

sys.path.append(dir_output)
os.chdir(dir_output)



########################################################################
########################################################################
L=10


nsample=11

mpth = path.join("models")
misc.mkdir_p(mpth)
pp = path.join(mpth,"pp",str(2))#list_fit[0].generate_pickle_path("constant",base=mpth)
misc.mkdir_p(pp)


########################################################################3
################### step 0  Data loading
########################################################################

list_D=list()
list_X=list()


Y=pd.read_csv(path.join('/dcs04/hansen/data/ywang/ST/data_10X_ST//mouse_Sagittal/put/Y_features_sele_sample'+ str(2+1) +'_v2_460genes.csv')) #_500genes
X=pd.read_csv(path.join('/dcs04/hansen/data/ywang/ST/data_10X_ST/mouse_Sagittal_spaceRanger1_1_0/out/X_sample'+ str(2+1) +'.csv'))
D=process_multiSample.get_D(X,Y)
D["sz"]=D["sz"]/10
list_D.append(D)
list_X.append(D["X"])


for k in range(0,10):
	print(k)
	X=pd.read_csv(path.join("Xsub_xenium_sub"+str(k+1)+".csv")).iloc[:,1:3]
	Y=pd.read_csv(path.join("Ysub_xenium_sub"+str(k+1)+".csv")).iloc[:,1:461]
	D=process_multiSample.get_D(X,Y)
	D["sz"]=D["sz"]/10
	list_D.append(D)
	list_X.append(D["X"])



list_Dtrain=process_multiSample.get_listDtrain(list_D)
list_sampleID=process_multiSample.get_listSampleID(list_D)


# inducing points, 70% of total spots for each sample
ksample = 0
random.seed(111)
ninduced=round(list_D[ksample]['X'].shape[0] * 0.35)
random.seed(222)
print(ninduced)
D=list_D[ksample]
rd_ = random.sample(range(0, D['X'].shape[0]), ninduced)
list_D[ksample]["Z"]=D['X'][rd_ ,:]


for ksample in range(1,11):
	random.seed(111)
	ninduced=round(list_D[ksample]['X'].shape[0] * 0.35)
	random.seed(222)
	print(ninduced)
	D=list_D[ksample]
	rd_ = random.sample(range(0, D['X'].shape[0]), ninduced)
	list_D[ksample]["Z"]=D['X'][rd_ ,:]



########################################################################3
################### step 1 initialize model
########################################################################
#lik="nb"

list_fit=process_multiSample.ini_multiSample(list_D,L,"nb")

process_multiSample.save_object(list_fit, 'list_fit_nb_2samples_szMean_L10_saveMem_ini.pkl') 

########################################################################
################### step 2 fit model
########################################################################


list_fit=training_multiSample.train_model_mNSF(list_fit,pp,list_Dtrain,list_D, S=1)

#list_fit=training_multiSample.train_model_mNSF(list_fit[0:2],pp,list_Dtrain[0:2],list_D[0:2],num_epochs=1, S=1)
# save the fitted model
process_multiSample.save_object(list_fit, 'list_fit_nb_2samples_szMean_L10_saveMem.pkl') 



########################################################################
with open( 'list_fit_nb_2samples_szMean_L10_saveMem.pkl', 'rb') as inp:
              list_fit = pickle.load(inp)

  
########################################################################
################### step 3 save and plot results
########################################################################
inpf12=process_multiSample.interpret_npf_v3(list_fit,list_X,S=1,lda_mode=False)


W = inpf12["loadings"]
#Wdf=pd.DataFrame(W*inpf12["totals1"


Wdf=pd.DataFrame(W*inpf12["totalsW"][:,None],  columns=range(1,L+1))
Wdf.to_csv(path.join("loadings_L10_saveMem.csv"))


## save the factors
#inpf12 = process_multiSample.interpret_npf_v3(list_fit,list_X,S=100,lda_mode=False)
Factors = inpf12["factors"][:,0:L]

for k in range(0,nsample):
	indices=list_sampleID[k]
	indices=indices.astype(int)
	Factors_df = pd.DataFrame(Factors[indices,:]) 
	Factors_df.to_csv(path.join(dir_output,"factors_nb_szMean_sample_s"+str(k+1)+"_L10_saveMem.csv"))


#




	

