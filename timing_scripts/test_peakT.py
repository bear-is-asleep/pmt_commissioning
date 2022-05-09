import sys
import numpy as np
import pandas as pd
sys.path.append('/sbnd/app/users/brindenc/mypython') #My utils path
from bc_utils.pmtutils import pic as pmtpic
from bc_utils.pmtutils import plotters as pmtplotters
from bc_utils.utils import pic,plotters
from datetime import date
from time import time
import matplotlib.pyplot as plt


#Date info
day = date.today().strftime("%Y_%m_%d")

#Load dataframes
op_ew_df = pd.read_pickle('data/op_ew_df.pkl')
op_fb_df = pd.read_pickle('data/op_fb_df.pkl')
pmts = pd.read_pickle('data/PMT_info.pkl')

#Constants
tpc=1
bw = 0.002 #check two digitizer readouts
digitize_rate = 0.002 #Digitize every 2 ns
number_pmts=10 #Number of pmts that need to be greater than the threshold
thresholdPE=100 #PE for PMT to pass threshold
t_width = 0.05 #+- this value for time window width
max_plots = 5 #Specify max number of plots to make
learn_rate = 20
minthreshold = 10

#Plot labels
xlabel = 'peakT [us]'
ylabel = 'PE'

#Get indeces
indeces = op_ew_df.index.drop_duplicates()

cnt = 0
pic.print_stars()
for index in indeces:
  if cnt >= max_plots: break
  start = time()
  #Run info
  run = index[0]
  subrun = index[1]
  event = index[2]

  #PeakT info and getting corresponding pe info
  peakT = pmtpic.get_peakT(op_ew_df,pmts,tpc,index,bw,number_pmts=number_pmts,thresholdPE=thresholdPE)
  temp_threshold = thresholdPE #Set temp threshold to handle events that do not pass, we will decrease it until we get the right peakT
  while peakT is None:
    if temp_threshold <= minthreshold:
      peakT = -9999 #Set dumby value for pmts that don't see the proper amount of light
      break
    temp_threshold = temp_threshold - learn_rate
    peakT = pmtpic.get_peakT(op_ew_df,pmts,tpc,index,bw,number_pmts=number_pmts,thresholdPE=temp_threshold) #Get peakT from function
  op_df = op_ew_df.loc[index]
  op_df = op_df[op_df.loc[:,'ophit_peakT']>peakT-t_width]
  op_df = op_df[op_df.loc[:,'ophit_peakT']<peakT+t_width]
  xs = op_df.loc[:,'ophit_peakT'].values
  ys = op_df.loc[:,'ophit_pe'].values
  cs = op_df.loc[:,'ophit_opdet'].values

  #Set bounds of plot
  leftbound = peakT-10*digitize_rate
  rightbound = peakT+20*digitize_rate

  #Make plot
  fig = plt.figure(figsize=(9,4))
  ax = fig.add_subplot()
  ax.scatter(xs,ys,c=cs)
  ax.set_xlabel(xlabel,fontsize=15)
  ax.set_ylabel(ylabel,fontsize=15)
  ax.axhline(temp_threshold,c='red',ls='--')
  ax.axvline(peakT,c='black',ls='-.')
  ax.axvline(leftbound,c='green',ls='-')
  ax.axvline(rightbound,c='green',ls='-')
  ax.set_title(f'PeakT Location and Threshold\nRun {run} Subrun {subrun} Event {event}')
  _,y_max = ax.get_ylim()
  #ax.text(peakT+15*digitize_rate,0.75*y_max,f'peakT = {peakT:.3f} us')

  #Set parameters text box
  parameters = {'peak T (us) ':f'{peakT:.3f}',
                #'Time slice (us) ':f'[{leftbound:.3f},{rightbound:.3f}]',
                'Threshold PE ':f'{temp_threshold}'
  }
  stattext = plotters.convert_p_str(parameters)
  props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
  ax.text(0.5, 0.9, stattext, transform=ax.transAxes, fontsize=15,
        verticalalignment='top', bbox=props)


  plotters.save_plot(f'peakT__{run}_{subrun}_{event}',folder_name='peakT_v1')

  end = time()
  print(f'peakT = {peakT:.3f} us\nRun {run} Subrun {subrun} Event {event} time: {end-start:.2f}s')

  cnt+=1

