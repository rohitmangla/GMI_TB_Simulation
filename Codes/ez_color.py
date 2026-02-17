def cmap(n1=None,n2=None,ncolor=None,r=None):
  '''
  This is a function to creat the colormap. modified from Chuntao's code
  01/04/2016
  '''
  import matplotlib.colors as color
  import xarray as xr
  import numpy as np
  import matplotlib.pyplot as plt
  if n1==None:n1=0
  if n2==None:n2=60
  if ncolor==None:ncolor=59
  if r==None:r=0
  colors=xr.open_dataset('D:/Datasets/GPM_PF/ez_color.nc')
  a=np.array([[colors['R'][n1:n2]/255.],[colors['G'][n1:n2]/255.],[colors['B'][n1:n2]/255.]])
  a=a.transpose().reshape((n2-n1),3)
  if r==0:
     cm=color.ListedColormap(a,N=ncolor)
# reverse the color map
  if r==1:
     a=a[::-1]
     cm=color.ListedColormap(a,N=ncolor)
  #if r==2:
  #   colors=read_libs.read_rgb('D:/Datasets/GPM_PF/topo_15lev.rgb')
#     a=np.array([[colors['R']],[colors['G']],[colors['B']]])
#     a=a.transpose().reshape(16,3)
     cm=color.ListedColormap(a,N=ncolor)
  return cm

'''
cm=cmap(r=2)
# test the color map
import numpy as np
import matplotlib.pyplot as plt
x = np.random.random(50)
y = np.random.random(50)
c = np.random.random(50)  # color of points
s = 500 * np.random.random(50)  # size of points

fig, ax = plt.subplots()
im = ax.scatter(x, y, c=c, s=s, cmap=cm)

fig.colorbar(im, ax=ax)

# set the color limits - not necessary here, but good to know how.
im.set_clim(0.0, 1.0)
plt.show()
'''
