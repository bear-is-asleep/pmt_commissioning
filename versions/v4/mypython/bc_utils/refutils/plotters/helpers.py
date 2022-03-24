import matplotlib.pyplot as plt
from random import choice
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from datetime import date
import sys
import seaborn as sns
import matplotlib
from scipy import optimize
from mpl_toolkits.mplot3d import Axes3D
from bc_utils.refutils import pic
from numpy import pi
import glob
from PIL import Image
#matplotlib.rcParams['axes.unicode_minus'] = False

xls = 16 #axis size
tls = 20 #title size
lls = 16 #legend size
tbs=10 #Text box size
#Small displacement for text
small = 5
matplotlib.rcParams['xtick.labelsize'] = 13
matplotlib.rcParams['ytick.labelsize'] = 13 #Set to same as axis size 
matplotlib.rcParams['axes.labelsize'] = 'large'
matplotlib.rcParams['legend.fontsize'] = 'large'
matplotlib.rcParams['axes.titlesize'] = 'x-large'
matplotlib.rcParams['figure.figsize'] = [9,7]


def use_science_style():
  if 'science' in plt.style.available:
    plt.style.use(['science','no-latex'])
  else:
    os.system('pip install SciencePlots -q')
    plt.style.reload_library()
    plt.style.use(['science','sns'])

#Organization
def make_plot_dir():
    day = date.today().strftime("%d_%m_%Y")
    isDir = os.path.isdir("Plots/Plots_"+day)
    if isDir == False:
        os.system("mkdir -p Plots_" +day)
        os.system("mv -n Plots_" +day+"/ Plots/")

def save_plot(fname,fig=None,ftype='.jpg',dpi=500):
    day = date.today().strftime("%d_%m_%Y")
    if fig == None:
      plt.savefig(f'{fname}{ftype}',bbox_inches = "tight",dpi=dpi)
    else:
      fig.savefig(f'{fname}{ftype}',bbox_inches = "tight",dpi=dpi)
    #print(os.getcwd())
    os.system("mv " + fname + "* Plots/Plots_" +day+"/")

def plot_stuff():
  matplotlib.rcParams['xtick.labelsize'] = 13
  matplotlib.rcParams['ytick.labelsize'] = 13 #Set to same as axis size 
  matplotlib.rcParams['axes.labelsize'] = 'large'
  matplotlib.rcParams['legend.fontsize'] = 'large'
  matplotlib.rcParams['axes.titlesize'] = 'x-large'
  matplotlib.rcParams['figure.figsize'] = [9,7]
  make_plot_dir()
  use_science_style()

def make_plane_meshes(x,yl=-200,yu=200,zl=0,zu=500,yiter=400,ziter=500,xiter=200,n_ref=1):
  ys = np.linspace(yl,yu,2) #Get all y refs.
  zs = np.linspace(zl,zu,2) #Get all z refs.
  yy,zz = np.meshgrid(ys,zs) #Mesh grid
  xx = np.full(yy.shape,x) #Populate mesh grid
  xs = np.linspace(-n_ref*xiter,+n_ref*xiter,2*n_ref+1) #Get all x refs.
  ys = np.linspace(-n_ref*yiter,yiter*n_ref,2*n_ref+1) #Get all y refs.
  zs = np.linspace(-ziter*n_ref,+ziter*n_ref,2*n_ref+1) #Get all z refs.
  #print(xx,yy,zz,xs,ys,zs)
  grids = []
  for xal in xs:
    for yal in ys:
      for zal in zs:
        grids.append([xx+xal,yy+yal,zz+zal,0]) #Add mesh grid to list
  return grids

def plot_ref(xkey,ykey,zkey,ckey,plane_grid,n_ref=1,x=None,y=None,z=None,
coord_ref_df=pd.DataFrame(),elev=20,azim=40,axis='on',title='',cmap='viridis_r',plot_planes=True):
  if coord_ref_df.empty:
    if x != None and y != None and z != None:
      coord_ref_df = pic.single_plane_reflection(x,y,z,x_refs=n_ref,y_refs=n_ref,z_refs=n_ref,initialize=True)
    else: 
      raise Exception('We need some coordinates for this lame dataframe')
      
  fig = plt.figure(figsize=(9,7))
  ax = Axes3D(fig,auto_add_to_figure=False)
  fig.add_axes(ax)

  x = coord_ref_df.loc[:,xkey]
  y = coord_ref_df.loc[:,ykey]
  z = coord_ref_df.loc[:,zkey]
  c = coord_ref_df.loc[:,ckey]

  if n_ref > 0:
    s = 200/n_ref
  else: 
    s= 200
  im = ax.scatter3D(x,z,y,c=c,s=s,cmap=cmap,edgecolors='black',alpha=0.5)
  if plot_planes:
    for grid in plane_grid:
      plane_ref = ax.plot_surface(grid[0],grid[2],grid[1],alpha=0.3)
  #ax.axis('off')
  ax.view_init(elev,azim)
  if axis != 'off':
    fig.colorbar(im)
    ax.set(xlabel='x',ylabel='z',zlabel='y')
    ax.set_title(title,fontsize=tls)
  else:
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
  ax.axis(axis)
  return ax

def make_ref_gif(cwd,n_images,n_ref,a=200,b=1,c=100):
  xs = np.linspace(0,2*b,n_images)
  ts = -a*(xs-b)**2 + a*b**2
  #Filepaths
  day = date.today().strftime("%d_%m_%Y")
  fp_in = 'ref*.png'
  fp_out = f'n_ref__{n_ref}.gif'

  #print(fp_in)
  os.chdir(cwd)
  print(os.getcwd())

  imgs = glob.glob(fp_in)
  #print(imgs,sortstrings_numerically(imgs))
  imgs = pic.sortstrings_numerically(imgs)
  frames = [Image.open(image) for image in imgs]
  #imgs = natsorted(fp_in, key=lambda y: y.lower())#sort alphanumeric in ascending order
  img = frames[0]
  img.save(fp=fp_out,format='GIF',append_images=frames,save_all=True,duration=100,loop=0)

  #subprocess.call("convert -delay 150 -loop 5 " + fp_in + " " + fp_out, shell=True)
  #os.system("start output.gif")
  
    
  

def make_images(params,xiter=0,yiter=2,ziter=2.5,eleviter=.02,azimiter=.01,n_ref_iter=0,
  n_images=100,cmap='viridis_r',plot_planes=True):
  #params = [0,0,0,0,0,1] #initialize params x,y,z,elev,azim,ref
  i = 0
  n_ref_float = params[-1]
  elev = params[3]
  azim = params[4]
  while i < n_images: #make n_images images
    y = 200*np.sin((i*yiter+params[1])*pi/180)
    z = 250*np.sin((i*ziter+params[2])*pi/180)+250
    x = 200*np.cos((i*xiter+params[0]*pi/180))+100
    if x > 200: #we don't like stragglers
      x=0
    #elev = 90*np.sin(eleviter*i*pi/180)
    #azim = 90*np.sin(azimiter*i*pi/180)
    elev = elev + eleviter
    azim = azim + azimiter
    n_ref_float = n_ref_float+i*n_ref_iter
    n_ref = int(n_ref_float) #n_ref

    coord_ref_df = pic.single_plane_reflection(x,y,z,x_refs=n_ref,y_refs=n_ref,z_refs=n_ref,initialize=True)
    if plot_planes:
      grids = make_plane_meshes(x,n_ref=n_ref)
    else:
      grids = None
    ax = plot_ref('x','y','z','tot_ref',grids,n_ref=n_ref,
    coord_ref_df=coord_ref_df,axis='off',azim=azim,elev=elev,title=i,cmap=cmap,
    plot_planes=plot_planes)
    save_plot(f'ref{i}',ftype='.png',dpi=200)
    plt.close()

    #iterate x,y,z
    i += 1