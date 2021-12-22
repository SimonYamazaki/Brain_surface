#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:41:49 2021

@author: simonyamazaki
"""

#%%

import nibabel as nib
import numpy as np

fname = "/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/tutorial_data_20190918_1558/buckner_data/tutorial_subjs/004/surf/rh.pial";

surf = nib.freesurfer.read_geometry(fname)

points = surf[0];
edges = surf[1];

#%%


Total_area = 0
Total_vol = 0

p0 = np.array([30,0,0])

for ii,p in enumerate(edges):
    
    print(ii)
    
    p1 = points[p[0],:]
    p2 = points[p[1],:]
    p3 = points[p[2],:]
    
    v1 = p1-p2
    v2 = p1-p3
    v3 = p1-p0
    
    #volume
    p_matrix = np.concatenate((v1.reshape(3,1),v2.reshape(3,1),v3.reshape(3,1)),axis=1)
    T_vol = 1/6*np.linalg.det(p_matrix)
    Total_vol += T_vol
    
    #area
    v4 = np.cross(v1,v2)
    T_area = 0.5*np.linalg.norm(v4)
    Total_area += T_area
    
    #triangular_vertices
    tv = np.concatenate((p1.reshape(3,1),p2.reshape(3,1),p3.reshape(3,1)),axis=1).T
    tv = tv.reshape(1,3,3)
    
    if ii == 0:
        tv2 = tv
    else:
        tv2 = np.concatenate((tv2,tv),axis=0)
        
else:
    np.save("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/triangle_vertices.npy",tv2)
        

tarea = f"The total surface area: {Total_area}"
tvol = f"The total volumne: {Total_vol}"

print(tarea)
print(tvol)

with open('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/measures.txt', 'a') as f:
    f.writelines([tarea,tvol])

#%% Curvature 
    
import math 

GC_list = []
MC_list = []

#curvature for each vertex
for ii,p in enumerate(points):
    
    print(ii)
    
    #gauss curvature
    T_idx = np.argwhere(edges==ii)[:,0]
    
    Tp_idx = edges[T_idx][edges[T_idx]!=ii]
    Tp_idxU = np.unique(Tp_idx)
    Tp_idx2 = Tp_idx.reshape(int(len(Tp_idx)/2),2)
    
    #T_points = points[Tp_idx]
    #T_edges = T_points - p
    
    T_area_sum = 0
    T_angle_sum = 0
    T_angle = []
    ab_cross_n = []
    
    for n,pj in enumerate(Tp_idx2):
        
        pj1 = p - points[pj[0]]
        pj2 = p - points[pj[1]]
        
        T_angle_sum += np.arccos(np.dot(pj1,pj2) / (np.linalg.norm(pj1)*np.linalg.norm(pj2)))
        T_angle.append(np.arccos(np.dot(pj1,pj2) / (np.linalg.norm(pj1)*np.linalg.norm(pj2))))
        
        ab_cross = np.cross(pj1,pj2)
        ab_cross_n.append(ab_cross / np.linalg.norm(ab_cross))
        
        T_area_sum += 1/2*np.linalg.norm(ab_cross)
    
    GC = (2*math.pi - T_angle_sum) / ( 1/3*T_area_sum )
    
    GC_list.append(GC)
    
    
    #mean curvature 
    edge_sum = np.zeros(3)
    AWPN = np.zeros(3)
    
    #for the edge k in the 1 ring around p
    for jj,k in enumerate(Tp_idxU):
        
        AWPN += T_angle[jj]*ab_cross_n[jj]
        
        #Tp_idx_with_common_edge
        ce = np.argwhere(Tp_idx2==k)
        p1_idx = Tp_idx2[ce[0,0],1-ce[0,1]]
        p2_idx = Tp_idx2[ce[1,0],1-ce[1,1]]
        
        idx3 = [p1_idx, k, p2_idx]
        
        T_points2 = points[idx3]
        T_edges2 = T_points2 - p
        
        outer_edge1 = T_points2[0]-T_points2[1]
        outer_edge2 = T_points2[2]-T_points2[1]
        
        alpha = np.arccos(np.dot(outer_edge1,T_edges2[0]) / (np.linalg.norm(outer_edge1)*np.linalg.norm(T_edges2[0])))
        beta = np.arccos(np.dot(outer_edge2,T_edges2[2]) / (np.linalg.norm(outer_edge2)*np.linalg.norm(T_edges2[2])))
        
        edge_sum +=  (np.cos(alpha)/np.sin(alpha) + np.cos(beta)/np.sin(beta)) * (p - T_points2[1])

    T_1ring_area = T_area_sum*4
    
    AWPN = AWPN / np.linalg.norm(AWPN)
    
    MCN = edge_sum / T_1ring_area
    
    k1 = MCN[0] / AWPN[0]
    k2 = MCN[1] / AWPN[1]
    k3 = MCN[2] / AWPN[2]
    
    kappa = np.array([k1,k2,k3]).mean()
    
    MC_list.append(kappa)
        
else:
    MC2 = np.array(MC_list)
    GC2 = np.array(GC_list)
    np.save("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/GC2.npy",GC2)
    np.save("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/MC2.npy",MC2)
    
    MCF = []
    GCF = []

    for e in edges:
        MCF.append(MC2[e].mean())
        GCF.append(GC2[e].mean())
    
    MCF = np.array(MCF)
    GCF = np.array(GCF)
    
    np.save("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/GCF.npy",GCF)
    np.save("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/MCF.npy",MCF)
    
    np.savetxt("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/GCF.csv", GCF, delimiter=",")
    np.savetxt("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/MCF.csv", MCF, delimiter=",")

with open('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/measures.txt', 'a') as f:
    f.writelines([f"Mean GC: {GC2.mean()}",f"Mean MC: {MC2.mean()}"])



#%% plot the triangulated brain 

import numpy as np
from pylab import get_cmap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import nibabel as nib

fname = "/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/tutorial_data_20190918_1558/buckner_data/tutorial_subjs/004/surf/rh.pial";
surf = nib.freesurfer.read_geometry(fname,read_metadata=True)
points = surf[0];
edges = surf[1];

tv2 = np.load("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/triangle_vertices.npy")
print("was here")

tv2 = tv2[:10000,:,:]

GC2 = np.load("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/GC2.npy")
MC2 = np.load("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/MC2.npy")

GCF = np.load("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/GCF.npy")
MCF = np.load("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/MCF.npy")


#%% Regular surface 


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(points[:,0], points[:,1], points[:,2],triangles=edges, color = [0.5,1,0.5])
#fig.colorbar(surf)

ax.set_xlim(points[:,0].min(),points[:,0].max())
ax.set_ylim(points[:,1].min(),points[:,1].max())
ax.set_zlim(points[:,2].min(),points[:,2].max())

ax.view_init(elev=20, azim=120)

plt.savefig('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/brain_angle1.png')  

#plt.show()


#%% Regular surface from another angle 

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(points[:,0], points[:,1], points[:,2],triangles=edges, color = [0.5,1,0.5])
#fig.colorbar(surf)

ax.set_xlim(points[:,0].min(),points[:,0].max())
ax.set_ylim(points[:,1].min(),points[:,1].max())
ax.set_zlim(points[:,2].min(),points[:,2].max())

ax.view_init(elev=20, azim=50)

plt.savefig('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/brain_angle2.png')  

#plt.show()

#%% zoomed 

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


#%% Colored surface

import numpy as np
from pylab import get_cmap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

cmap = get_cmap('plasma')

#cmap = get_cmap('Spectral')

MCplot = MCF[:tv2.shape[0]]
face_color_index = MCplot/MCplot.max()

# 0.5 centered colors
GCplot = GCF[:tv2.shape[0]]
maxC = np.array([abs(GCplot.min()), abs(GCplot.max())]).max()
GCplot = GCplot+maxC
face_color_index = GCplot/GCplot.max()

collection = Poly3DCollection(tv2, facecolors=cmap(face_color_index))

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.add_collection(collection)

ax.set_xlim(points[:,0].min(),points[:,0].max())
ax.set_ylim(points[:,1].min(),points[:,1].max())
ax.set_zlim(points[:,2].min(),points[:,2].max())

#plt.savefig('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/brain_colors.png')  
plt.show()


#%% Colored surface zoom binary, GC

import numpy as np
from pylab import get_cmap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

cmap = get_cmap('plasma')

#cmap = get_cmap('Spectral')

#binary
GCplot = GCF[:tv2.shape[0]]
GCplot[GCplot<0] = 0
GCplot[GCplot>0] = 1
face_color_index = GCplot

collection = Poly3DCollection(tv2, facecolors=cmap(face_color_index))

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.add_collection(collection)

ax.set_xlim(20,30)
ax.set_ylim(-85,-75)
ax.set_zlim(-15,-5)

#plt.savefig('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/brain_colors.png')  
plt.show()


#%% Colored surface zoom binary, MC

import numpy as np
from pylab import get_cmap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

cmap = get_cmap('plasma')

#cmap = get_cmap('Spectral')

#binary
GCplot = MCF[:tv2.shape[0]]
GCplot[GCplot<0] = 0
GCplot[GCplot>0] = 1
face_color_index = GCplot

collection = Poly3DCollection(tv2, facecolors=cmap(face_color_index))

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.add_collection(collection)

ax.set_xlim(20,30)
ax.set_ylim(-85,-75)
ax.set_zlim(-15,-5)

#plt.savefig('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/brain_colors.png')  
plt.show()


#%%
"""
from pylab import get_cmap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection='3d')

cmap = get_cmap('Spectral')

face_color_index = np.linspace(0, 1, num=edges.shape[0])
colors = cmap(face_color_index)

ax.plot_trisurf(points[:,0], points[:,1], points[:,2], triangles = edges )
ax.set_fc(colors)

plt.show()

"""


#%% FROM STACK OVERFLOW


"""
# Setup is the same

import itertools
import matplotlib.pyplot as plt
from pylab import get_cmap
from matplotlib.tri import Triangulation, LinearTriInterpolator
import numpy as np
from scipy import stats
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def simplex(n_vals):
    base = np.linspace(0, 1, n_vals, endpoint=False)
    coords = np.asarray(list(itertools.product(base, repeat=3)))
    return coords[np.isclose(coords.sum(axis=-1), 1.0)]

sim = simplex(20)
pdf = stats.dirichlet([1.1, 1.5, 1.3]).pdf(sim.T)

# For shorter notation we define x, y and z:

x = sim[:, 0]
y = sim[:, 1]
z = sim[:, 2]

# Creating a triangulation object and using it to extract the actual triangles. 
# Note if it is necessary that no patch will be vertical (i.e. along the z direction)

tri = Triangulation(x, y)

triangle_vertices = np.array([np.array([[x[T[0]], y[T[0]], z[T[0]]],
                                        [x[T[1]], y[T[1]], z[T[1]]], 
                                        [x[T[2]], y[T[2]], z[T[2]]]]) for T in tri.triangles])

# Finding coordinate for the midpoints of each triangle. 
# This will be used to extract the color

midpoints = np.average(triangle_vertices, axis = 1)
midx = midpoints[:, 0]
midy = midpoints[:, 1]

# Interpolating the pdf and using it with the selected cmap to produce the color RGB vector for each face. 
# Some roundoff and normalization are needed

face_color_function = LinearTriInterpolator(tri, pdf)
face_color_index = face_color_function(midx, midy)
face_color_index[face_color_index < 0] = 0
face_color_index /= np.max(pdf)

cmap = get_cmap('Spectral')

# Creating the patches and plotting

face_color_index = face_color_index._data
face_color_index = np.linspace(0,1,len(face_color_index))

collection = Poly3DCollection(triangle_vertices, facecolors=cmap(face_color_index))


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.add_collection(collection)
plt.show()

"""


#%% Colored surface zoom
"""
import numpy as np
from pylab import get_cmap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

cmap = get_cmap('plasma')

#cmap = get_cmap('Spectral')

MCplot = MCF[:tv2.shape[0]]
face_color_index = MCplot/MCplot.max()

GCplot = GCF[:tv2.shape[0]]
maxC = np.array([abs(GCplot.min()), abs(GCplot.max())]).max()
GCplot = GCplot+maxC

face_color_index = GCplot/GCplot.max()

face_color_index = np.log(face_color_index)
maxC = np.array([abs(face_color_index.min()), abs(face_color_index.max())]).max()
face_color_index = face_color_index+maxC


collection = Poly3DCollection(tv2, facecolors=cmap(face_color_index))

print("was here1")

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.add_collection(collection)
print("was here2")

#ax.color_bar(face_color_index)

ax.set_xlim(0,10)
ax.set_ylim(-55,-45)
ax.set_zlim(-35,-25)

"""
ax.set_xlim(20,30)
ax.set_ylim(-85,-75)
ax.set_zlim(-15,-5)
"""

print("was here3")

#plt.savefig('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/brain_colors.png')  
plt.show()

"""