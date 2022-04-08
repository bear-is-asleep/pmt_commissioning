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
from bc_utils.utils import pic
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

def remove_outliers(arr,max_dev=4):
  #Remove outliers based on standard deviations of data to keep
  median,std = np.median(arr),np.std(arr)
  zero_based = abs(arr-median)
  return arr[zero_based < max_dev * std]

def max_bin_height(ax,bins):
  #Get max bin heigth
  max_bin = 0
  for bar, b0, b1 in zip(ax.containers[0], bins[:-1], bins[1:]):
    if bar.get_height() > max_bin:
      max_bin = bar.get_height()
  return max_bin

def make_lines():
  plt.axhline(y=62.5,linestyle='--',linewidth=3)
  plt.axhline(y=-62.5,linestyle='--',linewidth=3)
  plt.axvline(x=125,linestyle='--',linewidth=3)
  plt.axvline(x=250,linestyle='--',linewidth=3)
  plt.axvline(x=375,linestyle='--',linewidth=3)

def convert_p_str(parameters):
  s = ''
  last_key = list(parameters)[-1]
  for key in parameters.keys():
    if last_key == key:
      s+= f"{key} = {parameters[key][0]} {parameters[key][1]}"
    else:
      s+= f"{key} = {parameters[key][0]} {parameters[key][1]}\n"
  return s

def kde_1dhist(df,titles,xlabels,ylabels,keys,show,save,save_name,fit_gaus,trim_outliers):
  #Make kde plots
  for (i,title) in enumerate(titles):
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot()
    arr = df[keys[i]].to_numpy()
    w = 0.2
    n = math.ceil((arr.max() - arr.min())/w) #set bin width
    counts,_ = np.histogram(arr,n)
    plt.hist(arr,fill=True,bins=n,label='Data')
    if trim_outliers[i]:
      arr = remove_outliers(arr,max_dev=0.06)
      w = 0.2
      n = math.ceil((arr.max() - arr.min())/w) #set bin width
      counts,_ = np.histogram(arr,n)
      plt.hist(arr,fill=True,bins=n,label='Fit Data')
    if fit_gaus[i]: #Fit curve to gaussian and plot results
      """x,popt,pcov = fit_gaussian(arr)
      perr = np.sqrt(np.diag(pcov))
      params = {'$A$':[f'{popt[0]:.2f}$\pm${perr[0]:.1e}',''],
                '$\mu$':[f'{popt[1]:.2f}$\pm${perr[1]:.1e}',''],
                '$\sigma$':[f'{popt[2]:.2f}$\pm${perr[2]:.1e}','']}
      """
      mu, std = norm.fit(arr)
      params = {'$\mu$':[f'{mu:.2f}',''],
                '$\sigma$':[f'{std:.2f}','']}
      s = convert_p_str(params)
      ax.text(0.75, 0.85,s,transform=ax.transAxes,bbox=dict(facecolor='Pink',alpha=0.5),horizontalalignment='left',fontsize=lls)
      # Plot the PDF.
      xmin, xmax = plt.xlim()
      x = np.linspace(xmin, xmax, int(1e5))
      p = norm.pdf(x, mu, std)
      scale = max(counts)/max(p)
      plt.plot(x, scale*p, 'k', linewidth=1,label='Gaussian Fit') 
      #ax.plot(x,gaussian(x,popt[0],popt[1],popt[2]),label='Gaussian fit')

    plt.xlabel(xlabels[i],fontsize=xls)
    plt.ylabel(ylabels[i],fontsize=xls)
    plt.title(titles[i],fontsize=tls)
    ax.legend(fontsize=lls)
    if i == 0 or i == 2:
      plt.xlim(-2, 5)
      ax.axvline(x=0,linestyle='-.',color='red')
    if i == 1:
      plt.xlim(-0.1, 4)
      ax.axvline(x=1,linestyle='-.',color='red')
    if save[i]:
      save_plot(save_name[i])
    if show[i]:
      plt.show()
    else:
      plt.close(fig)

def plot_TPC(tpc,label,label_title,df,coating=2,cmap='viridis',return_plot=False):
  #If coating is 2, plot both, coating 0 for coated, coating 1 for uncoated
  if df.shape[0] != 120:
    print('df needs to contain only one event, or combined events')
    return None


  #plt.figure()
  #sns.distplot(df[label].to_numpy())
  #plt.show()
  #Plot 2d hist with colorscale as label
  fig = plt.figure(figsize=(10,8))
  ax = fig.add_subplot()
  ax.set_facecolor('cyan')
  make_lines()



  data_points = []
  for _,line in df.iterrows():
    skip = False #Skips text for PMTs that are filtered
    det_type = int(line['ophit_opdet_type'])
    det_ch = str(int(line['ophit_ch']))
    x = line['ophit_opdet_x']
    y = line['ophit_opdet_y']
    z = line['ophit_opdet_z']
    c = line[label]
    data_line = [x,y,z,c]
    if tpc == 0 and x < 0:
      #Apply coating cut
      if coating == 2:
        if det_type == 0 or det_type == 1:
          data_points.append(data_line)
          skip = True #Keeps good PMTs
      elif coating == 1 or coating == 0:
        if det_type == coating:
          data_points.append(data_line)
          skip = True #Keeps good PMTs
      if z > 250 and skip:
        if det_ch == 166 or det_ch == 160:
          ax.text(z-2*small,y-2*small,det_ch,fontsize=tbs)
        else:
          ax.text(z-2*small,y+small,det_ch,fontsize=tbs)
      if z < 250 and skip:
        if det_ch == 150 or det_ch == 156:
          ax.text(z+0.15*small,y-2*small,det_ch,fontsize=tbs)
        else:
          ax.text(z+0.15*small,y+small,det_ch,fontsize=tbs)
    
    if tpc == 1 and x > 0: #186 included (for some reason)
      #Apply coating cut
      if coating == 2:
        if det_type == 0 or det_type == 1:
          data_points.append(data_line)
          skip = True #Keeps good PMTs
      elif coating == 1 or coating == 0:
        if det_type == coating:
          data_points.append(data_line)
          skip = True #Keeps good PMTs
      if z > 250 and skip: 
        if det_ch == 166 or det_ch == 160:
          ax.text(z-2*small,y-2*small,det_ch,fontsize=tbs)
        else:
          ax.text(z-2*small,y+small,det_ch,fontsize=tbs)
      if z < 250 and skip:
        if det_ch == 150 or det_ch == 156:
          ax.text(z+0.15*small,y-2*small,det_ch,fontsize=tbs)
        else:
          ax.text(z+0.15*small,y+small,det_ch,fontsize=tbs)

  data_points = np.asarray(data_points)
  #print(data_points[:,3])
  sc = ax.scatter(data_points[:,2],data_points[:,1],c=data_points[:,3],cmap=cmap,s=80,alpha=0.7)
  ax.margins(x=0.05)
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cbar = fig.colorbar(sc,cax=cax,ax=ax)
  #cbar.set_label(f'{label_title}',rotation=270,fontsize=xls)
  ax.set_xlabel('Z [cm]',fontsize=xls)
  ax.set_ylabel('Y [cm]',fontsize=xls)
  ax.set_title(f'{label_title} TPC{tpc}',fontsize = tls)
  if return_plot:
    return fig,ax,sc
  else:
    return fig,ax

#Plot muon tracks
def plot_tracks(df,x1_key0,y1_key0,x2_key0,y2_key0,x1_key1,y1_key1,x2_key1,y2_key1,ax='None',fig='None'):
  #Plot muon tracks with coordinate x as horizantel axis, y as vertical
  #Use keys

  #Set figure 
  if ax == 'None' and fig == 'None':
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot()
    plot_text = True
  else: plot_text = False
  #TPC0
  x1_0 = df[x1_key0].to_numpy()
  y1_0 = df[y1_key0].to_numpy()
  x2_0 = df[x2_key0].to_numpy()
  y2_0 = df[y2_key0].to_numpy()
  #TPC1
  x1_1 = df[x1_key1].to_numpy()
  y1_1 = df[y1_key1].to_numpy()
  x2_1 = df[x2_key1].to_numpy()
  y2_1 = df[y2_key1].to_numpy()
  
  #Keep track of event ID
  events = df['event'].to_numpy()


  #Create labels for plot
  if 'x' in x1_key0:
    xlabel = 'x [cm]'
  if 'y' in x1_key0:
    xlabel = 'y [cm]'
  if 'z' in x1_key0:
    xlabel = 'z [cm]'
  if 'x' in y1_key0:
    ylabel = 'x [cm]'
  if 'y' in y1_key0:
    ylabel = 'y [cm]'
  if 'z' in y1_key0:
    ylabel = 'z [cm]'
  title = f'$\mu$ track trajectory'

  #print(events.shape,x1_1.shape)
  events0 = []
  events1 = []
  for i in range(events.shape[0]):
    #Small offset for text
    smallx = choice([-1*small,small,0])
    smally = choice([-1*small,small,0])
    #In TPC0
    if x1_1[i] == -999 and x2_1[i] == -999 and y1_1[i] == -999 and y2_1[i] == -999:
     #and events[i] not in events0:
      #Slope and intercept for text
      events0.append(events[i])
      xs = [x1_0[i],x2_0[i]]
      ys = [y1_0[i],y2_0[i]]
      xtext,ytext = (xs[1]+xs[0])/2+smallx,(ys[1]+ys[0])/2+smally
      if plot_text:
        ax.plot(xs,ys)
        ax.set_xlabel(xlabel,fontsize=xls)
        ax.set_ylabel(ylabel,fontsize=xls)
        ax.set_title(title,fontsize=tls)
        ax.text(xtext,ytext,f'{int(events[i])}',fontsize=tbs)
      else: 
        ax.plot(xs,ys,ls='--',c='r',label=r'$\mu$ tracks')
        #ax.legend(fontsize=lls)
        #plt.text(xtext,ytext,f'{int(events[i])}',fontsize=tbs)
    #In TPC1
    if x1_0[i] == -999 and x2_0[i] == -999 and y1_0[i] == -999 and y2_0[i] == -999:
      # and events[i] not in events1:
      #Slope and intercept for text
      events1.append(events[i])
      xs = [x1_1[i],x2_1[i]]
      ys = [y1_1[i],y2_1[i]]
      xtext,ytext = (xs[1]+xs[0])/2+smallx,(ys[1]+ys[0])/2+smally
      if plot_text:
        ax.plot(xs,ys)
        ax.set_xlabel(xlabel,fontsize=xls)
        ax.set_ylabel(ylabel,fontsize=xls)
        ax.set_title(title,fontsize=tls)
        ax.text(xtext,ytext,f'{int(events[i])}',fontsize=tbs)
      else: 
        ax.plot(xs,ys,ls='--',c='r',label=r'$\mu$ tracks')
        #ax.legend(fontsize=lls)
        #plt.text(xtext,ytext,f'{int(events[i])}',fontsize=tbs)

  #save_plot('tracks')
  #plt.show()

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
coord_ref_df=pd.DataFrame(),elev=20,azim=40,axis='on',title=''):
  if coord_ref_df.empty:
    if x != None and y != None and z != None:
      coord_ref_df = pic.single_plane_reflection(x,y,z,x_refs=n_ref,y_refs=n_ref,z_refs=n_ref,initialize=True)
    else: 
      raise Exception('We need some coordinates for this lame dataframe')
  if plane_grid:
    #print('Making Reflection inception')
    ttt = 0
  else:
    raise Exception('We need a plane grid')
      
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
  im = ax.scatter3D(x,z,y,c=c,s=s,cmap='viridis_r',edgecolors='black',alpha=0.6)
  for grid in plane_grid:
    plane_ref = ax.plot_surface(grid[0],grid[2],grid[1],alpha=0.3)
  #ax.axis('off')
  ax.view_init(elev,azim)
  if axis != 'off':
    fig.colorbar(im)
  else:
    ax.set(xlabel='x',ylabel='z',zlabel='y')
    ax.set_title(title,fontsize=tls+15)
  ax.axis(axis)
  return ax
