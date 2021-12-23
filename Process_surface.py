#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:41:49 2021

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


#%% compute surface area and volume


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
    
    #reshape the surface data used for a matplotlib library
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

#write the results to a file 
with open('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/measures.txt', 'a') as f:
    f.writelines([tarea,tvol])

#%% Compute curvature s
    
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
    #write the curvature at each vertex to a numpy file
    MC2 = np.array(MC_list)
    GC2 = np.array(GC_list)
    np.save("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/GC2.npy",GC2)
    np.save("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/MC2.npy",MC2)
    
    #the average curvature of each triangle
    #average of the 3 vertices in each triangle
    MCF = []
    GCF = []

    for e in edges:
        MCF.append(MC2[e].mean())
        GCF.append(GC2[e].mean())
    
    MCF = np.array(MCF)
    GCF = np.array(GCF)
    
    #write the average curvature of each triangle to numpy and csv file
    np.save("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/GCF.npy",GCF)
    np.save("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/MCF.npy",MCF)
    
    np.savetxt("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/GCF.csv", GCF, delimiter=",")
    np.savetxt("/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/MCF.csv", MCF, delimiter=",")


#write the average of positive and negative curvatures to file 
GC2_pos = GC2[GC2>0].mean()
GC2_neg = GC2[GC2<0].mean()

MC2_pos = MC2[MC2>0].mean()
MC2_neg = MC2[MC2<0].mean()

with open('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/measures.txt', 'a') as f:
    f.writelines([f"GC2_pos: {GC2_pos}",f"GC2_neg: {GC2_neg)}",f"MC2_pos: {MC2_pos}",f"MC2_neg: {MC2_neg)}"])


#write the average curvature of all vertices to file
with open('/Users/simonyamazaki/Documents/1_M/diff_geo1/Project/measures.txt', 'a') as f:
    f.writelines([f"Mean GC: {GC2.mean()}",f"Mean MC: {MC2.mean()}"])



#%% Load the computed measures 

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


#%% Colored surface
#not included in report
"""
import numpy as np
from pylab import get_cmap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

cmap = get_cmap('plasma')


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
"""

#%% Colored surface zoom binary, GC
#not included in report
"""
import numpy as np
from pylab import get_cmap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

cmap = get_cmap('plasma')

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
"""

#%% Colored surface zoom binary, MC
#not included in report

"""
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
"""

