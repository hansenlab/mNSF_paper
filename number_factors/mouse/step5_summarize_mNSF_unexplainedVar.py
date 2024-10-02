

dir_mNSF_functions='/users/ywang/Hansen_projects/scRNA/mNSF_2023_10_20/'
dir_output='/dcs04/hansen/data/ywang/ST/data_10X_ST/mouse_Sagittal_spaceRanger1_1_0/out/'

########################################################################
########################################################################
import sys
sys.path.append(dir_mNSF_functions)


import mNSF

from mNSF import process_multiSample

from mNSF.NSF import preprocess
from mNSF.NSF import misc
from mNSF.NSF import visualize
#from mNSF import training_multiSample
from mNSF import training_multiSample
from mNSF import process_multiSample
from scanpy import read_h5ad
#from tensorflow.data import Dataset
from os import path
#import pandas
import os
import numpy as np
import tensorflow as tf
import pandas as pd
import sys 
import pickle



sys.path.append(dir_output)
os.chdir(dir_output)


########################################################################
########################################################################
L=16
nsample=4

dpth='data'

pth="/dcs04/hansen/data/ywang/ST/data_10X_ST/mouse_Sagittal_spaceRanger1_1_0/out/"
mpth = path.join(pth,"models_")
misc.mkdir_p(mpth)
pp = path.join(mpth,"pp")#list_fit[0].generate_pickle_path("constant",base=mpth)
misc.mkdir_p(pp)


########################################################################3
################### step 0  Data loading
########################################################################

list_D=list()
list_X=list()


import random
for ksample in range(0,nsample):
	Y=pd.read_csv(path.join(dpth,'/dcs04/hansen/data/ywang/ST/data_10X_ST//mouse_Sagittal/put/Y_features_sele_sample'+ str(ksample+1) +'_v2_500genes.csv')) #_500genes
	#Y=pd.read_csv(path.join(dpth,'/dcs04/hansen/data/ywang/ST/data_10X_ST/mouse_Sagittal_spaceRanger1_1_0/out/Y_sample'+ str(ksample+1) +'.csv'))
	X=pd.read_csv(path.join(dpth,'/dcs04/hansen/data/ywang/ST/data_10X_ST/mouse_Sagittal_spaceRanger1_1_0/out/X_sample'+ str(ksample+1) +'.csv'))
	D=process_multiSample.get_D(X,Y)
	D["sz"]=D["sz"]/10
	list_D.append(D)
	list_X.append(D["X"])




########################################################################3
################### step 0  Data loading
########################################################################
list_Dtrain=process_multiSample.get_listDtrain(list_D) #,nbatch=1
list_sampleID=process_multiSample.get_listSampleID(list_D)



import random
for ksample in range(0,nsample):
	random.seed(222)
	ninduced=round(list_D[ksample]['X'].shape[0] * 0.1)
	print(ninduced)
	D=list_D[ksample]
	list_D[ksample]["Z"]=D['X'][random.sample(range(0, D['X'].shape[0]-1), ninduced) ,:]
	
	


########################################################################
# get deviance for each sample
for L in [8,16,24,32,40]:# finished running
	print(L)
	with open( 'list_fit_500selectedFeatures_dev_interpolated_10percent_szMean_L'+str(L)+'.pkl', 'rb') as inp:
              list_fit = pickle.load(inp)
	vec_dev = np.zeros(4)
	for ksample in range(0,nsample):
		dev_mnsf=visualize.gof(list_fit[ksample],list_D[ksample],Dval=None,S=10,plot=False)
		vec_dev[ksample]=dev_mnsf['tr']['mean']
	vec_dev_df = pd.DataFrame(vec_dev) 
	vec_dev_df.to_csv(path.join("vec_dev_500selectedFeatures_dev_interpolated_10percent_szMean_L"+str(L)+".csv"))




