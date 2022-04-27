from cProfile import label
from curses import init_color
import numpy as np
from matplotlib.widgets import Slider, Button
import matplotlib.pyplot as plt
import sys
sys.path.append('/sbnd/app/users/brindenc/mypython') #My utils path
from bc_utils.pmtutils import pic as pmtpic
from bc_utils.pmtutils import plotters as pmtplotters
from bc_utils.hitutils import pic as hitpic
from bc_utils.hitutils import plotters as hitplotters
from bc_utils.utils import pic,plotters
import pandas as pd

#Define global params
global sc
global cax

#Load dataframes for plotting
scalarPE_df = pd.read_pickle('data/scalarPE_df.pkl')
vectorPE_df = pd.read_pickle('data/vectorPE_df.pkl')
muon_df = pd.read_pickle('data/muon_df_tpc0.pkl')

#Clean data
muon_plot_df = muon_df.drop(['nmuontrks','muontrk_t0'],axis=1)
muon_plot_df = pmtpic.get_muon_tracks(muon_plot_df)
muon_plot_df = muon_plot_df.set_index(['run','subrun','event'])

indeces = scalarPE_df.index.drop_duplicates()

#Pass event to display
argList = sys.argv[1:]

#Handle people who don't know how to run this
if len(argList)<1 or int(argList[0]) > len(indeces):
  raise Exception(f'You should specify an event (choose number between 0 and {indeces.shape[0]-1})')
else: 
  index = indeces[int(argList[0])]
  print(f'Making plot for Run {index[0]} Subrun {index[1]} Event {index[2]}')

#Keep single event
scalardf = scalarPE_df.loc[index]
vectordf = vectorPE_df.loc[index]
muon_df = muon_plot_df.loc[index]


#This will probably go in a plotter function
trights = vectordf.loc[:,'tright'].drop_duplicates().values #Use this for slider value

#Constants
vmax = vectordf.loc[:,'summed_PE'].values.max() #Max value for colorbar
coating=2 #Both coatings
tpc = 0 #Initialize tpc
tleft = scalardf.loc[:,'tleft'].drop_duplicates().values[0] #First time PE is seen

#Initial params
init_tright = trights[2] #Initialize right bound
title = f'PE for t$\in$[{tleft:.3f},{init_tright:.3f}] $\mu$s '
init_cs = vectordf[vectordf.loc[:,'tright'] == init_tright]#Initial colors, specify key later?
df = pd.concat([init_cs,scalardf],axis=1)
df = df.loc[:,~df.columns.duplicated()] #Remove duplicate columns


#Make initial plot
fig,ax,sc,cax = pmtplotters.interactive_TPC(tpc,'summed_PE',title,df,coating=coating,cmap='viridis',
  return_plot=True,normalize=False,facecolor='cyan',ax=None,fig=None,vmax=vmax)
pmtplotters.plot_tracks(muon_df,'muontrk_z1_0','muontrk_y1_0','muontrk_z2_0','muontrk_y2_0',
            'muontrk_z1_1','muontrk_y1_1','muontrk_z2_1','muontrk_y2_1',ax=ax,fig=fig,indeces=[index])

#Vertical slider
axtright = fig.add_axes([0.1,0.25,0.0225,0.63])
tright_slider = Slider(
ax=axtright,
label=r'$t_{max}$',
valmin=trights[0],
valmax=trights[-1],
valstep=trights, #Set valstep to constant binwidth
valinit=init_tright,
orientation='vertical',
initcolor=None
)

def update(val):
  global sc
  global cax
  val = round(val,5)
  title = f'PE for t$\in$[{tleft:.3f},{val:.3f}] $\mu$s '
  #ax.clear()
  cs = vectordf[round(vectordf.loc[:,'tright'],5) == val]#Initial colors, specify key later?
  #print(cs.head(),val)
  dfnew = pd.concat([cs,scalardf],axis=1)
  dfnew = dfnew.loc[:,~dfnew.columns.duplicated()] #Remove duplicate columns
  sc.remove()
  cax.remove()
  _,_,scnew,caxnew = pmtplotters.interactive_TPC(tpc,'summed_PE',title,dfnew,coating=coating,cmap='viridis',
  return_plot=True,normalize=False,facecolor='cyan',ax=ax,fig=fig,label_PMTs=False,vmax=vmax)
  #Keep updating scatterplot
  sc = scnew
  cax = caxnew
  #plt.show()
  
#Register slider value
tright_slider.on_changed(update)
# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')

def reset(evt):
  tright_slider.reset()
  tright_slider.on_changed(update)
button.on_clicked(reset)

plt.show()

