#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/users/ywang/anaconda/envs/mNSF4/lib

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


# %%
rng = np.random.default_rng(222)
def generate_gg_blocks():
    A = np.zeros( [ 4 , 36 ] )
    A[0, [ 1 , 6 , 7 , 8 , 13 ] ] = 1
    A[1, [ 3 , 4 , 5 , 9 , 11 , 15 , 16 , 17  ] ] = 1
    A[2, [ 18 , 24 , 25 , 30 , 31 , 32 ] ] = 1
    A[3, [ 21 , 22 , 23 , 28 , 34 ] ] = 1
    return A

ncopy = 5
nside = ncopy*6
N = (nside)**2
N=900 #number of spots per sample
X = misc.make_grid(N)
X[:,1] = -X[:,1] #make the display the same
X = preprocess.rescale_spatial_coords(X)
J = 500
L = 4
A = generate_gg_blocks()
A = A.reshape((L,6,6))
A = np.kron(A,np.ones((1,ncopy,ncopy)))
Ftrue = A.reshape((L,N)).T #NxL

X1 = X
X2= -X
X2 = -X2



#%%
w = rng.choice(L,J,replace=True)
Wtrue = 10.9*get_dummies(w).to_numpy(dtype=dtp) #JxL indicator matrix
v = rng.choice(4,J,replace=True)
Vtrue = 8.9*get_dummies(v).to_numpy(dtype=dtp) #Jx2 indicator matrix
Utrue = rng.binomial(1,0.2,size=(N,4))
UVt = Utrue @ Vtrue.T
Lambda_true = 0.2+Ftrue @ Wtrue.T + UVt #NxJ


Y1 = rng.negative_binomial(10,10/(Lambda_true+10))
Y2 = rng.negative_binomial(20,10/(Lambda_true+10))



Yss = rng.choice(Y1,replace=False,axis=1,size=4)
Yss = Y1[:,(0,2,3,4)]


Yss = rng.choice(Y2,replace=False,axis=1,size=4)
Yss = Y2[:,(0,2,3,4)]



#%% Save as anndata - sample 1
ad1 = AnnData(Y1,obsm={"spatial":X1,"factors":Ftrue},varm={"loadings":Wtrue})
ad1.layers = {"counts":ad1.X.copy()} #store raw counts before normalization changes ad1.X
pp.normalize_total(ad1, inplace=True, layers=None, key_added="sizefactor")
pp.log1p(ad1)
ad1.write_h5ad(path.join("ggblocks_lr_1_noiseDif.h5ad"),compression="gzip")

DF = pd.DataFrame(Ftrue) 
DF.to_csv("ggblocks_lr_1_noiseDif_sample1.csv")



#%% Save as anndata - sample 2
ad2 = AnnData(Y2,obsm={"spatial":X2,"factors":Ftrue},varm={"loadings":Wtrue})
ad2.layers = {"counts":ad2.X.copy()} #store raw counts before normalization changes ad2.X
pp.normalize_total(ad2, inplace=True, layers=None, key_added="sizefactor")
pp.log1p(ad2)
ad2.write_h5ad(path.join("ggblocks_lr_2_noiseDif.h5ad"),compression="gzip")
DF = pd.DataFrame(Ftrue) 
DF.to_csv("ggblocks_lr_1_noiseDif_sample2.csv")





#%% make  plots - sample 1
hmkw = {"figsize":(4,.9),"bgcol":"white","subplot_space":0.1,"marker":"s","s":10}
#fig,axes=visualize.heatmap(ad1.obsm["spatial"],ad1.obsm["factors"][:,3],s=50,marker="s",figsize=(4.1,4))
fig,axes=visualize.multiheatmap(X1, ad1.obsm["factors"], (1,4), cmap="Blues", **hmkw)
fig.savefig(path.join("ggblocks_factors_sample1_noiseDif.png"),bbox_inches='tight')




#%% make  plots - sample 2
#visualize.heatmap(ad2.obsm["spatial"],ad2.obsm["factors"][:,3],s=50,marker="s",figsize=(4.1,4))
fig,axes=visualize.multiheatmap(X2, ad2.obsm["factors"], (1,4), cmap="Blues", **hmkw)
#fig,axes=visualize.heatmap(ad1.obsm["spatial"],ad2.obsm["factors"][:,3],s=50,marker="s",figsize=(4.1,4))
fig.savefig(path.join("ggblocks_factors_sample2_noiseDif.png"),bbox_inches='tight')













