

########################################################################
dir_mNSF_functions='/users/ywang/Hansen_projects/scRNA/mNSF_2023_10_20/'
#dir_mNSF_functions='/users/ywang/Hansen_projects/scRNA/mNSF_2023_10_20/mNSF-main'
dir_output="/dcs04/hansen/data/ywang/ST/DLPFC/PASTE_out_keepAll_scTransform/"

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


nsample=12

mpth = path.join("models")
misc.mkdir_p(mpth)
pp = path.join(mpth,"pp",str(2))#list_fit[0].generate_pickle_path("constant",base=mpth)
misc.mkdir_p(pp)


########################################################################3
################### step 0  Data loading
########################################################################

list_D=list()
list_X=list()

for ksample in range(0,nsample):
	Y=pd.read_csv(path.join('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data//Y_features_sele_sample'+str(ksample+1)+'_500genes.csv'))
	X=pd.read_csv(path.join('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data///X_allSpots_sample'+str(ksample+1)+'.csv'))
	D=process_multiSample.get_D(X,Y)
	list_D.append(D)
	list_X.append(D["X"])
	

list_Dtrain=process_multiSample.get_listDtrain(list_D)
list_sampleID=process_multiSample.get_listSampleID(list_D)


# inducing points, 70% of total spots for each sample
for ksample in range(0,nsample):
	random.seed(111)
	ninduced=round(list_D[ksample]['X'].shape[0] * 0.35)
	random.seed(222)
	print(ninduced)
	D=list_D[ksample]
	rd_ = random.sample(range(0, D['X'].shape[0]), ninduced)
	list_D[ksample]["Z"]=D['X'][rd_ ,:]



########################################################################
# get deviance for each sample

#{'tr': {'mean': 1.0590545, 'argmax': 2, 'max': 4.113069, 'med': 1.3276968}}

########################################################################
for disp in ['0001', '001', '01', 1, 10, 100]:# still need to run： 001, 01
for disp in [ '10']:# still need to run： 01
	print(disp)
	with open( 'list_fit_nb_12samples_szMean_L10_fullData_disp'+str(disp)+'.pkl', 'rb') as inp:
              list_fit = pickle.load(inp)
	vec_dev = np.zeros(12)
	for ksample in range(0,nsample):
		dev_mnsf=visualize.gof(list_fit[ksample],list_D[ksample],Dval=None,S=10,plot=False)
		vec_dev[ksample]=dev_mnsf['tr']['mean']
	vec_dev_df = pd.DataFrame(vec_dev) 
	vec_dev_df.to_csv(path.join("vec_dev_500selectedFeatures_dev_interpolated_35percent_szMean_disp"+str(disp)+".csv"))


