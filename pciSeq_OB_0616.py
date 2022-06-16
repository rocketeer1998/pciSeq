# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 14:14:41 2022

@author: Rocketeer

"""

####################
####Cellpose
####################

import os
os.getcwd() 
os.chdir("D:/5.ISSDemo_OlfactoryBulb_0326") 
from urllib.parse import urlparse
import matplotlib.pyplot as plt
#import matplotlib as mpl
#%matplotlib inline
#mpl.rcParams['figure.dpi'] = 300
import numpy as np
#print(np.__version__)
from cellpose import utils, io
from cellpose import models, io
# DEFINE CELLPOSE MODEL
# model_type='cyto' or model_type='nuclei'
model = models.Cellpose(gpu=False, model_type='cyto') 
import skimage
from scipy import ndimage as ndi
from skimage import (
    io, color, feature, filters, measure, morphology, segmentation, util
)
from skimage.segmentation import watershed, expand_labels
import scipy
from scipy.sparse import (
    coo_matrix, save_npz, load_npz
)
import tifffile
import pandas as pd
from skimage.filters import threshold_multiotsu
from skimage.measure import label, regionprops

### Step1-Define cell_pose_segemenation_to_coo() function
def cell_pose_segemenation_to_coo(image, diam, expanded_distance):
    
    '''
    function to segement nuclei and expand segemented nuclei using cell pose. 
    cellpose is a generalist algorithm for cellular segmentation. the function returns a 
    coo object which saves the outlines of the cells. 
    
    the input for the function is the dapi stained images, the diameter of the objects 
    to segment and the distance the to expand (we generally do a 10% expansion). 
    
    the segementation startegy employs the scikit-image package and watershed
    segementation. 
    from the docs: 
    "The watershed is a classical algorithm used for segmentation, that is, for separating 
    different objects in an image. Starting from user-defined markers, the watershed algorithm
    treats pixels values as a local topography (elevation).The algorithm floods basins from the 
    markers until basins attributed to different markers meet on watershed lines. In many cases, 
    markers are chosen as local minima of the image, from which basins are flooded."
    '''
    
    # run cell pose segementation on the objects, and get masks
    masks_nuclei, flows, styles, diams = model.eval(image,diameter=diam)
    
    distance = ndi.distance_transform_edt(masks_nuclei)
    
    local_max_coords = feature.peak_local_max(distance, min_distance=7)
    
    local_max_mask = np.zeros(distance.shape, dtype=bool)
    local_max_mask[tuple(local_max_coords.T)] = True

    # find markers
    markers = measure.label(local_max_mask)

    # run watershed segementation
    segmented_cells = segmentation.watershed(-distance, markers, mask=masks_nuclei)
    seg1 = measure.label(segmented_cells)
    expanded = expand_labels(seg1, distance=expanded_distance)
    expanded_new = expanded.astype('uint32')
    coo = coo_matrix(expanded_new)
        
    return expanded_new, coo

### Step2-Watershed
filepath="D:/5.ISSDemo_OlfactoryBulb_0326/1.Raw/DAPI/DAPI-1_ch00.tif"
image = tifffile.imread(filepath)

fig, ax = plt.subplots(figsize=(35,30))
ax.imshow(image, cmap='gray')
ax.set_title('OlfactoryBulb')
plt.show()

coo_file = cell_pose_segemenation_to_coo(image, diam=46.2, expanded_distance=0.1) #46.2 pixels=15 um; 15min
coo_file[0]
coo_file[1] #label image matrix
np.savez("CellSeg.npz",coo_file[1])

####################
####pciSeq
####################

import numpy as np
import pandas as pd
import skimage.color
import matplotlib.pyplot as plt
from scipy.sparse import load_npz, coo_matrix
import copy
import loompy
import pciSeq
print(pciSeq.__version__)

### Part1-Segmentation coo_matrix file
#coo = load_npz("D:/5.ISSDemo_OlfactoryBulb_0326/CellSeg.npz")
coo = copy.deepcopy(coo_file[1])
print('The image has %d cells' % len(set(coo.data)))
# The image has 13265 cells

rgb_label_image = skimage.color.label2rgb(coo.toarray(), bg_label=0)

plt.style.use('dark_background')
plt.figure(figsize=(30, 30), dpi=300)
imgplot = plt.imshow(rgb_label_image)
plt.axis('off')
plt.show()

### Part2-Spots file
iss_spots = pd.read_csv("D:/5.ISSDemo_OlfactoryBulb_0326/4.pciSeq_0616/spots.csv")
iss_spots.head()

plt.figure(figsize=(30, 30), dpi=300)
imgplot = plt.imshow(rgb_label_image)
plt.scatter(iss_spots.x, iss_spots.y, s=0.1)
plt.axis('off')
plt.show()

### Part3-scRNA-seq data
ds = loompy.connect("D:/5.ISSDemo_OlfactoryBulb_0326/4.pciSeq_0616/l1_olfactory.loom")
mx = ds[(ds.ra.Gene == ["Slc17a7"]) | 
   (ds.ra.Gene == ["Gad1"]) |
   (ds.ra.Gene == ["Pcp4"]) |
   (ds.ra.Gene == ["Penk"]) |
   (ds.ra.Gene == ["Cck"]) |
   (ds.ra.Gene == ["Doc2g"]) |
   (ds.ra.Gene == ["Cdhr1"]) |
   (ds.ra.Gene == ["Nrsn1"]) |
   (ds.ra.Gene == ["Kcnj4"]) |
   (ds.ra.Gene == ["Camk4"]) |
   (ds.ra.Gene == ["Gsn"]) |
   (ds.ra.Gene == ["Nr2f2"]) |
   (ds.ra.Gene == ["Kctd12"]) |
   (ds.ra.Gene == ["Fabp7"]) |
   (ds.ra.Gene == ["Nrgn"]) |
   (ds.ra.Gene == ["Ptn"]), :]
mx.shape
mx = pd.DataFrame(mx,columns=ds.ca.Class)
mx = mx.rename(index={0: 'Slc17a7',
                      1:'Gad1',
                      2:'Pcp4',
                      3:'Penk',
                      4:'Cck',
                      5:'Doc2g',
                      6:'Cdhr1',
                      7:'Nrsn1',
                      8:'Kcnj4',
                      9:'Camk4',
                      10:'Gsn',
                      11:'Nr2f2',
                      12:'Kctd12',
                      13:'Fabp7',
                      14:'Nrgn',
                      15:'Ptn'})
mx.to_csv("Matrix.csv") 
ds.close()

### Main Step
opts = {'save_data': True}
cellData, geneData = pciSeq.fit(iss_spots, coo, mx, opts) #16:23-

### Result
### Part1-cellData
cellData.head()
### cell_number,(x,y)
cellData.iloc[0].Cell_Num, (int(cellData.iloc[0].X), int(cellData.iloc[0].Y)) 
### The genes and their counts for cell1
pd.DataFrame(zip(cellData.iloc[0].Genenames, cellData.iloc[0].CellGeneCount), columns=['Gene names', 'Gene counts'])
### The cell types of that cell1 and the corresponding probabilities
pd.DataFrame(zip(cellData.iloc[0].ClassName, cellData.iloc[0].Prob), columns=['Class name', 'Prob'])

### Part2-geneData
geneData.head()
