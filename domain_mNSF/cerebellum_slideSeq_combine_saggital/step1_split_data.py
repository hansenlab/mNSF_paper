

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


nsample=2

mpth = path.join("models")
misc.mkdir_p(mpth)
pp = path.join(mpth,"pp",str(2))#list_fit[0].generate_pickle_path("constant",base=mpth)
misc.mkdir_p(pp)


########################################################################3
################### step 0  Data loading
########################################################################


Y=pd.read_csv(path.join('/dcs04/hansen/data/ywang/revision/SpatialPCA_example_data/SlideseqCerebellum/Y.csv')) #_500genes
X=pd.read_csv(path.join('/dcs04/hansen/data/ywang/revision/SpatialPCA_example_data/SlideseqCerebellum/X.csv')) #_500genes


for k in range(0,9):
	print(k)
	st = k*2500
	end = st + 2500
	Xsub=pd.DataFrame(X.iloc[st:end,:])
	Ysub=pd.DataFrame(Y.iloc[st:end,:])
	Xsub.to_csv(path.join("Xsub_xenium_sub"+str(k+1)+".csv"))
	Ysub.to_csv(path.join("Ysub_xenium_sub"+str(k+1)+".csv"))

k=9
st = k*2500
end = 25415
Xsub=pd.DataFrame(X.iloc[st:end,:])
Ysub=pd.DataFrame(Y.iloc[st:end,:])
Xsub.to_csv(path.join("Xsub_xenium_sub"+str(k+1)+".csv"))
Ysub.to_csv(path.join("Ysub_xenium_sub"+str(k+1)+".csv"))



	

