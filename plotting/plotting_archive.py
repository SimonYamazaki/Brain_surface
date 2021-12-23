#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 16:34:21 2021

@author: simonyamazaki
"""



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


