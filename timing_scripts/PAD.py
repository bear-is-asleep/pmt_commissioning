import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.pyplot as plt
import sys
sys.path.append('/sbnd/app/users/brindenc/mypython') #My utils path
from bc_utils.pmtutils import pic as pmtpic
from bc_utils.pmtutils import plotters as pmtplotters
from bc_utils.hitutils import pic as hitpic
from bc_utils.hitutils import plotters as hitplotters
from bc_utils.utils import pic,plotters
import pandas as pd
import matplotlib.image as mpimg

#Define style
#plt.style.use('science')

#Define global params
global sc
global cax

#Define extra functions
def findvmax(df,coating,df_label):
  #Decide max colorscale
  #print(coating,df_label,df.head())
  if coating == -1:
    vmax =  df.loc[:,df_label].values.max()
  elif coating == 1 or coating == 0 or coating == 2 or coating == 3:
    vmax = df[df.loc[:,'ophit_opdet_type'] == coating].loc[:,df_label].values.max()
  elif coating == 2.5: #All x-arapucas
    max1 = df[df.loc[:,'ophit_opdet_type'] == 2].loc[:,df_label].values.max()
    max2 = df[df.loc[:,'ophit_opdet_type'] == 3].loc[:,df_label].values.max()
    vmax = max([max1,max2]) #Max between two values
  elif coating == 0.5: #All pmts
    max1 = df[df.loc[:,'ophit_opdet_type'] == 0].loc[:,df_label].values.max()
    max2 = df[df.loc[:,'ophit_opdet_type'] == 1].loc[:,df_label].values.max()
    vmax = max([max1,max2]) #Max between two values
  return vmax

def load_dfs(which_sample,tpc,ind=0,loadg4=True,loadmuon=True):
  #Make dataframes
  scalarPE_df = pd.read_pickle(f'data/scalarPE{tpc}_{which_sample}_df.pkl')
  vectorPE_df = pd.read_pickle(f'data/vectorPE{tpc}_{which_sample}_df.pkl')
  indeces = scalarPE_df.index.drop_duplicates()

  index = indeces[int(argList[0])]

  scalarPE_df = scalarPE_df.sort_index()
  vectorPE_df = vectorPE_df.sort_index()

  #Keep single event
  scalardf = scalarPE_df.loc[index]
  vectordf = vectorPE_df.loc[index]

  dfs = [scalardf,vectordf]

  #Load dataframes for plotting
  if loadmuon:
    muon_df = pd.read_pickle(f'data/muon_{which_sample}_df.pkl')
    muon_df = muon_df.sort_index()
    #Clean data - print trajectory through tpc
    muon_plot_df = muon_df.drop(['nmuontrks','muontrk_t0'],axis=1)
    muon_plot_df = pmtpic.get_muon_tracks(muon_plot_df)
    muon_plot_df = muon_plot_df.set_index(['run','subrun','event'])
    muon_plot_df = muon_plot_df.loc[index]
    #Print statements
    pic.print_stars()
    print(f'Making plot for Run {index[0]} Subrun {index[1]} Event {index[2]} TPC{tpc}')
    if tpc == 0:
      x1 = muon_plot_df.loc[:,'muontrk_x1_0'].values
      x2 = muon_plot_df.loc[:,'muontrk_x2_0'].values
      y1 = muon_plot_df.loc[:,'muontrk_y1_0'].values
      y2 = muon_plot_df.loc[:,'muontrk_y2_0'].values
      z1 = muon_plot_df.loc[:,'muontrk_z1_0'].values
      z2 = muon_plot_df.loc[:,'muontrk_z2_0'].values
    if tpc == 1:
      x1 = muon_plot_df.loc[:,'muontrk_x1_1'].values
      x2 = muon_plot_df.loc[:,'muontrk_x2_1'].values
      y1 = muon_plot_df.loc[:,'muontrk_y1_1'].values
      y2 = muon_plot_df.loc[:,'muontrk_y2_1'].values
      z1 = muon_plot_df.loc[:,'muontrk_z1_1'].values
      z2 = muon_plot_df.loc[:,'muontrk_z2_1'].values

    for i in range(len(x1)):
        print(f'Muon {i} moves from [{x2[i]:.1f},{y2[i]:.1f},{z2[i]:.1f}] to [{x1[i]:.1f},{y1[i]:.1f},{z1[i]:.1f}]')

    pic.print_stars()
    dfs.append(muon_df)
  if loadg4:
    g4_df = pd.read_pickle(f'data/g4_{which_sample}_df.pkl')
    g4_df = g4_df.sort_index()
    dfs.append(g4_df)

  return dfs,index

#Pass event to display
argList = sys.argv[1:]

#Handle people who don't know how to run this
if len(argList)<2:
  raise Exception(f'You should specify an event for the first input.\nSecondly you should specify the tpc (0 or 1)\nThe last value should specify sample type (ew or fb)')
else: 
  tpc = int(argList[1])
  ind = int(argList[0])
  which_sample = argList[2]
if which_sample == 'gun':
  #Make dataframes
  scalarPE_df = pd.read_pickle(f'../top_bottom/data/scalarPE_df.pkl')
  vectorPE_df = pd.read_pickle(f'../top_bottom/data/vectorPE_df.pkl')
  indeces = scalarPE_df.index.drop_duplicates()

  index = indeces[int(argList[0])]

  #Load dataframes for plotting
  muon_df = pd.read_pickle(f'../top_bottom/data/muon_df_tpc0.pkl')
else:
  #Make dataframes
  dfs,index = load_dfs(which_sample,tpc,ind=ind)
  scalardf = dfs[0]
  vectordf = dfs[1]
  muon_df = dfs[2]
  g4_df = dfs[3]
  # scalarPE_df = pd.read_pickle(f'data/scalarPE{tpc}_{which_sample}_df.pkl')
  # vectorPE_df = pd.read_pickle(f'data/vectorPE{tpc}_{which_sample}_df.pkl')
  # indeces = scalarPE_df.index.drop_duplicates()

  # index = indeces[int(argList[0])]

  # #Load dataframes for plotting
  # muon_df = pd.read_pickle(f'data/muon_{which_sample}_df.pkl')
  # g4_df = pd.read_pickle(f'data/g4_{which_sample}_df.pkl')

# scalarPE_df = scalarPE_df.sort_index()
# vectorPE_df = vectorPE_df.sort_index()

# muon_t0 = muon_df.loc[index,'muontrk_t0']

# #Clean data
# muon_plot_df = muon_df.drop(['nmuontrks','muontrk_t0'],axis=1)
# muon_plot_df = pmtpic.get_muon_tracks(muon_plot_df)
# muon_plot_df = muon_plot_df.set_index(['run','subrun','event'])

# #Keep single event
# scalardf = scalarPE_df.loc[index]
# vectordf = vectorPE_df.loc[index]
# #print(muon_plot_df,index)
# muon_plot_df = muon_plot_df.loc[index]

# #Print statements
# pic.print_stars()
# print(f'Making plot for Run {index[0]} Subrun {index[1]} Event {index[2]} TPC{tpc}')
# if tpc == 0:
#   x1 = muon_plot_df.loc[:,'muontrk_x1_0'].values
#   x2 = muon_plot_df.loc[:,'muontrk_x2_0'].values
#   y1 = muon_plot_df.loc[:,'muontrk_y1_0'].values
#   y2 = muon_plot_df.loc[:,'muontrk_y2_0'].values
#   z1 = muon_plot_df.loc[:,'muontrk_z1_0'].values
#   z2 = muon_plot_df.loc[:,'muontrk_z2_0'].values
# if tpc == 1:
#   x1 = muon_plot_df.loc[:,'muontrk_x1_1'].values
#   x2 = muon_plot_df.loc[:,'muontrk_x2_1'].values
#   y1 = muon_plot_df.loc[:,'muontrk_y1_1'].values
#   y2 = muon_plot_df.loc[:,'muontrk_y2_1'].values
#   z1 = muon_plot_df.loc[:,'muontrk_z1_1'].values
#   z2 = muon_plot_df.loc[:,'muontrk_z2_1'].values

# for i in range(len(x1)):
#     print(f'Muon {i} moves from [{x2[i]:.1f},{y2[i]:.1f},{z2[i]:.1f}] to [{x1[i]:.1f},{y1[i]:.1f},{z1[i]:.1f}]')

# pic.print_stars()

#This will probably go in a plotter function
trights = vectordf.loc[index,'tright'].drop_duplicates().values #Use this for slider value
tleft = scalardf.loc[:,'tleft'].drop_duplicates().values[0] #First time PE is seen
#Constants/initialize
cdict = {'PE':'summed_PE',
          'Channels':'ophit_opch',
          'Coatings':'ophit_opdet_type'} #Dictionary for coloring and labeling points
tdict = {'All':-1,
          'X-ARAPUCA':2.5,
          'PMT':0.5,
          'Coated PMT':0,
          'Uncoated PMT':1,
          'VIS X-ARAPUCA':2,
          'VUV X-ARAPUCA':3} #Dictionary for coating type
tpcdict = {'TPC0 - West APA':0,
           'TPC1 - East APA':1} #Dictionary for TPC
df_label='summed_PE'
coating=-1 #All coatings

#PAD
textcolor = 'white'
labelcolor = 'black'
facecolor = 'black'
figcolor = 'gray'
cmap = 'hot'
markboxes = False
#Lines
plotg4 = True
plotmuon = True
show_legend = True
thigh=1000
tlow=-1000
#Slider/boxes
fillcolor = 'teal'
boxfacecolor = 'lightgoldenrodyellow'
handlesize = 50



#Initial params
tright = trights[3] #Initialize right bound
title = f'PE for t$\in$[{tleft:.3f},{tright:.3f}] $\mu$s '
init_cs = vectordf[vectordf.loc[:,'tright'] == tright]#Initial colors, specify key later?
df = pd.concat([init_cs,scalardf],axis=1)
df = df.loc[:,~df.columns.duplicated()] #Remove duplicate columns
if df_label == 'summed_PE':
  vmax = findvmax(df,coating,'tot_PE')
else:
  vmax = findvmax(df,coating,df_label)

#Make initial plot
fig,ax,sc,cax = pmtplotters.interactive_TPC(tpc,df_label,title,df,coating=coating,cmap=cmap,
  return_plot=True,normalize=False,facecolor=facecolor,ax=None,fig=None,vmax=vmax,text_label=df_label,
  textcolor=textcolor,markboxes=markboxes,figcolor=figcolor)
#pmtplotters.plot_tracks(muon_df,'muontrk_z1_0','muontrk_y1_0','muontrk_z2_0','muontrk_y2_0',
#            'muontrk_z1_1','muontrk_y1_1','muontrk_z2_1','muontrk_y2_1',indeces=[index],ax=ax,fig=fig,tpc=tpc)
pmtplotters.plot_g4_muon(g4_df,muon_df,index,thigh,tlow,x='z',y='y',show_vtx=False,remove_other=True,fig=fig,ax=ax,
    show_legend=show_legend,small=0,checktpc=True,tpc=tpc,plotg4=plotg4,plotmuon=plotmuon) #plotting package

#Vertical slider
axtright = fig.add_axes([0.05,0.9,0.9,0.02],facecolor=boxfacecolor)
tright_slider = Slider(
ax=axtright,
label=r'$t_{max}$',
valmin=trights[0],
valmax=trights[-1],
valstep=trights, #Set valstep to constant binwidth
valinit=tright,
orientation='horizontal',
initcolor=None,
color=fillcolor
#handle_style={'facecolor'=fillcolor}
)

def update(val):
  global ax
  global fig
  global sc
  global cax
  global df
  global title
  global df_label
  global coating
  global tright
  global tpc
  global vectordf
  global scalardf
  del ax.texts[:] #Remove all text from figure
  tright = round(val,5)
  title = f'PE for t$\in$[{tleft:.3f},{tright:.3f}] $\mu$s '
  #ax.clear()
  cs = vectordf[round(vectordf.loc[:,'tright'],5) == tright]#Initial colors, specify key later?
  #print(cs.head(),val)
  df = pd.concat([cs,scalardf],axis=1)
  df = df.loc[:,~df.columns.duplicated()] #Remove duplicate columns
  sc.remove()
  cax.remove()
  if df_label == 'summed_PE':
    vmax = findvmax(df,coating,'tot_PE')
  else:
    vmax = findvmax(df,coating,df_label)
  fignew,axnew,scnew,caxnew = pmtplotters.interactive_TPC(tpc,df_label,title,df,coating=coating,cmap=cmap,
    return_plot=True,normalize=False,facecolor=facecolor,ax=ax,fig=fig,vmax=vmax,text_label=df_label,
    textcolor=textcolor,markboxes=markboxes,figcolor=figcolor)
  fig = fignew
  ax = axnew
  sc = scnew
  cax = caxnew
  
  
#Register slider value
tright_slider.on_changed(update)
# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
# resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
# button = Button(resetax, 'Reset', hovercolor='0.975')

# def reset(evt):
#   tright_slider.reset()
#   tright_slider.on_changed(update)

# button.on_clicked(reset)

#Select colors/text
radiocolor=boxfacecolor
axradio1 = fig.add_axes([0.825,0.1,0.1,0.15],facecolor=radiocolor)
radio1 = RadioButtons(axradio1,('PE','Channels','Coatings'))
def colorbuttons(label): #Color button update function
  global ax
  global fig
  global sc
  global cax
  global df
  global title
  global df_label
  global coating
  global tright
  global tpc
  global vectordf
  global scalardf
  df_label = cdict[label] #Pass this key into interactive tpc

  del ax.texts[:] #Remove all text from figure
  sc.remove() #Remove scatter poitns
  cax.remove() #Remove coloscale
  if df_label == 'summed_PE':
    vmax = findvmax(df,coating,'tot_PE')
  else:
    vmax = findvmax(df,coating,df_label)
  fignew,axnew,scnew,caxnew = pmtplotters.interactive_TPC(tpc,df_label,title,df,coating=coating,cmap=cmap,
    return_plot=True,normalize=False,facecolor=facecolor,ax=ax,fig=fig,vmax=vmax,text_label=df_label,
    textcolor=textcolor,markboxes=markboxes,figcolor=figcolor)
  fig = fignew
  ax = axnew
  sc = scnew
  cax = caxnew
  plt.draw()
radio1.on_clicked(colorbuttons)

axradio2 = fig.add_axes([0.825,0.3,0.1,0.15],facecolor=radiocolor)
radio2 = RadioButtons(axradio2,('All','X-ARAPUCA','PMT','Coated PMT','Uncoated PMT','VIS X-ARAPUCA','VUV X-ARAPUCA'))
def coatingbuttons(label): #Color button update function
  global ax
  global fig
  global sc
  global cax
  global df
  global title
  global df_label
  global coating
  global tright
  global tpc
  global vectordf
  global scalardf
  coating = tdict[label] #Pass this key into interactive tpc
  #Decide max colorscale
  if df_label == 'summed_PE':
    vmax = findvmax(df,coating,'tot_PE')
  else:
    vmax = findvmax(df,coating,df_label)
    

  del ax.texts[:] #Remove all text from figure
  sc.remove() #Remove scatter poitns
  cax.remove() #Remove coloscale
  fignew,axnew,scnew,caxnew = pmtplotters.interactive_TPC(tpc,df_label,title,df,coating=coating,cmap=cmap,
    return_plot=True,normalize=False,facecolor=facecolor,ax=ax,fig=fig,vmax=vmax,text_label=df_label,
    textcolor=textcolor,markboxes=markboxes,figcolor=figcolor)
  fig = fignew
  ax = axnew
  sc = scnew
  cax = caxnew
  plt.draw()
radio2.on_clicked(coatingbuttons)

axradio3 = fig.add_axes([0.825,0.5,0.1,0.15],facecolor=radiocolor)
radio3 = RadioButtons(axradio3,('TPC0 - West APA','TPC1 - East APA'))
def tpcbuttons(label): #Color button update function
  global ax
  global fig
  global sc
  global cax
  global df
  global title
  global df_label
  global coating
  global tright
  global tpc
  global vectordf
  global scalardf
  tpc = tpcdict[label] #Pass this key into interactive tpc
  dfs,index = load_dfs(which_sample,tpc,ind=ind,loadmuon=True,loadg4=False)
  scalardf = dfs[0]
  vectordf = dfs[1]
  muon_df = dfs[2]

  cs = vectordf[vectordf.loc[:,'tright'] == tright]#Initial colors, specify key later?
  df = pd.concat([cs,scalardf],axis=1)
  df = df.loc[:,~df.columns.duplicated()] #Remove duplicate columns

  title = f'PE for t$\in$[{tleft:.3f},{tright:.3f}] $\mu$s '
  #Decide max colorscale
  if df_label == 'summed_PE':
    vmax = findvmax(df,coating,'tot_PE')
  else:
    vmax = findvmax(df,coating,df_label)
  del ax.texts[:] #Remove all text from figure
  sc.remove() #Remove scatter poitns
  cax.remove() #Remove coloscale
  del ax.lines[:] #Remove tracks
  if markboxes:
    pmtplotters.make_lines(ax=ax) #Add boundary lines again
  fignew,axnew,scnew,caxnew = pmtplotters.interactive_TPC(tpc,df_label,title,df,coating=coating,cmap=cmap,
    return_plot=True,normalize=False,facecolor=facecolor,ax=ax,fig=fig,vmax=vmax,text_label=df_label,
    textcolor=textcolor,markboxes=markboxes,figcolor=figcolor)
  pmtplotters.plot_g4_muon(g4_df,muon_df,index,thigh,tlow,x='z',y='y',show_vtx=False,remove_other=True,fig=fignew,ax=axnew,
    show_legend=show_legend,small=0,checktpc=True,tpc=tpc,plotg4=plotg4,plotmuon=plotmuon) #plotting package
  fig = fignew
  ax = axnew
  sc = scnew
  cax = caxnew
  plt.draw()
radio3.on_clicked(tpcbuttons)

axsbnd = fig.add_axes([0.825,0.7,0.1,0.15])
sbnd = mpimg.imread('Images/SBND-color.jpg')
axsbnd.imshow(sbnd)
axsbnd.axis('off');



plt.show()
