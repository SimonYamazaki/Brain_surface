#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 16:32:56 2021

@author: simonyamazaki
"""

#%% load the brain surface

import nibabel as nib
import numpy as np

fname = "/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/tutorial_data_20190918_1558/buckner_data/tutorial_subjs/004/surf/rh.pial";

surf = nib.freesurfer.read_geometry(fname)

#points are the vertices
points = surf[0];

#edges are the points of each triangle
edges = surf[1];



#%% Load the computed measures needed for plotting

import numpy as np
from pylab import get_cmap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import nibabel as nib

tv2 = np.load("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/triangle_vertices.npy")

tv2 = tv2[:10000,:,:]

GC2 = np.load("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/GC2.npy")
MC2 = np.load("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/MC2.npy")

GCF = np.load("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/GCF.npy")
MCF = np.load("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/MCF.npy")


#%% Regular surface in green


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(points[:,0], points[:,1], points[:,2],triangles=edges, color = [0.5,1,0.5])
#fig.colorbar(surf)

ax.set_xlim(points[:,0].min(),points[:,0].max())
ax.set_ylim(points[:,1].min(),points[:,1].max())
ax.set_zlim(points[:,2].min(),points[:,2].max())

ax.view_init(elev=20, azim=120)

plt.savefig('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/brain_angle1.png')  

plt.show()


#%% Regular surface in green from another angle 

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(points[:,0], points[:,1], points[:,2],triangles=edges, color = [0.5,1,0.5])
#fig.colorbar(surf)

ax.set_xlim(points[:,0].min(),points[:,0].max())
ax.set_ylim(points[:,1].min(),points[:,1].max())
ax.set_zlim(points[:,2].min(),points[:,2].max())

ax.view_init(elev=20, azim=50)

plt.savefig('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/brain_angle2.png')  

plt.show()

#%%  Regular surface in green magnified

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(points[:,0], points[:,1], points[:,2],triangles=edges, color = [0.5,1,0.5])
#fig.colorbar(surf)

p1v = abs(points[:,0].min()) + abs(points[:,0].max())
p1v = abs(points[:,1].min()) + abs(points[:,1].max())
p1v = abs(points[:,2].min()) + abs(points[:,2].max())

ax.set_xlim(points[:,0].min()+30,points[:,0].max()-30)
ax.set_ylim(points[:,1].min()+84.5,points[:,1].max()-84.5)
ax.set_zlim(points[:,2].min()+53,points[:,2].max()-53)

ax.view_init(elev=20, azim=50)

plt.savefig('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/brain_zoom.png')  


