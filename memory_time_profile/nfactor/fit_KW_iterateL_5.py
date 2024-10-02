### mNSF 4
#conda activate mNSF4
#export CUDA_VISIBLE_DEVICES=0
#export LD_LIBRARY_PATH=/usr/local/cuda/lib64

### KW edit of Yi fit.py, to benchmark memory usage
# iterate through various L
# Note: in training_multiSample.py line 271 in _train_model_fixed_lr, set num_epochs = 10 for benchmarking
import sys
import os
#%%
dir_output="//dcs04/hansen/data/ywang/ST/memory_time_profile/"
sys.path.append(dir_output)
os.chdir(dir_output)

dir_data="/dcs04/hansen/data/ywang/ST/DLPFC/processed_Data/"




########################################################################
########################################################################
import sys
import mNSF
from mNSF.NSF import preprocess, misc
from mNSF import process_multiSample, training_multiSample
from scanpy import read_h5ad
import os
from os import path
import numpy as np
import tensorflow as tf
import pandas as pd
import sys 
import pickle
import time
from csv import writer


########################################################################
########################################################################
L=5  ## CHANGE THIS
nsample=2
epochs = 20
legacy = False # Use legacy optimizer if tensorflow 2.12.0 +
tries = 1
lr = 0.0001


dpth='data'
pth="iterateL_" + str(L)
misc.mkdir_p(pth)
mpth = path.join(pth,"models")
misc.mkdir_p(mpth)
pp = path.join(mpth,"pp")#list_fit[0].generate_pickle_path("constant",base=mpth)
misc.mkdir_p(pp)

output_memory_list = []
output_memory_columns = ["L", "peak_before_first_iteration", "peak_after_training", "runtime"]

memory_outfile = "memoryUsage_mouseSagittal_n3_iterL_runSeperately.csv"

if not os.path.isfile(memory_outfile):
    with open(memory_outfile, 'a') as file:
        writer_object = writer(file)
        writer_object.writerow(output_memory_columns)
        file.close()

      


########################################################################3
################### step 0  Data loading
########################################################################

list_D=list()
list_X=list()
for ksample in range(0,nsample):
	Y=pd.read_csv(path.join(dir_data,'Y_small_forMemTimeProfile.csv'))
	X=pd.read_csv(path.join(dir_data,'Y_small_forMemTimeProfile.csv'))
	D=process_multiSample.get_D(X,Y)
	list_D.append(D)
	list_X.append(X)


	
list_Dtrain=process_multiSample.get_listDtrain(list_D)
list_sampleID=process_multiSample.get_listSampleID(list_D)


########################################################################3
################### step 1 initialize model
########################################################################

list_fit = process_multiSample.ini_multiSample(list_D,L)

# measuring here - would need to build into train_model_mNSF step
# can set max iteration to 10 and then measure
########################################################################3
################### step 2 fit model
########################################################################



output_memory_list.append(L)


start = time.time()




output_memory_list.append(tf.config.experimental.get_memory_info("GPU:0")['peak']/(1024*1024*1024))

list_fit=training_multiSample.train_model_mNSF(list_fit,pp,list_Dtrain,list_D, legacy=legacy, maxtry = tries, num_epochs=epochs, tol=0, lr=lr) # set tol to 0 to prevent convergence

output_memory_list.append(tf.config.experimental.get_memory_info("GPU:0")['peak']/(1024*1024*1024))


end = time.time()
output_memory_list.append(end-start)

with open(memory_outfile, 'a') as file:
    writer_object = writer(file)
    writer_object.writerow(output_memory_list)
    file.close()







