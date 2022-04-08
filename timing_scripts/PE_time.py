import sys
import os
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('/sbnd/app/users/brindenc/mypython') #My utils path
from bc_utils.pmtutils import pic as pmtpic
from bc_utils.pmtutils import plotters as pmtplotters
from bc_utils.utils import pic,plotters
import pandas as pd

#Load dataframes
op_df_tpc0 = pd.read_pickle('data/op_df_tpc0.pkl')
op_df_tpc1 = pd.read_pickle('data/op_df_tpc1.pkl')
#crt_df_tpc0 = pd.read_pickle('data/crt_df_tpc0.pkl')
#crt_df_tpc1 = pd.read_pickle('data/crt_df_tpc1.pkl')
#ctrk_df_tpc0 = pd.read_pickle('data/ctrk_df_tpc0.pkl')
#ctrk_df_tpc1 = pd.read_pickle('data/ctrk_df_tpc1.pkl')
#muon_df_tpc0 = pd.read_pickle('data/muon_df_tpc0.pkl')
#muon_df_tpc1 = pd.read_pickle('data/muon_df_tpc1.pkl')

#Get ready to plot!
plotters.plot_stuff()

#Iterable params
bws = [0.01] #binwidth: [nleft,nright]
pmts = [0,1,6,166] #Test coating, no coating and specific channels
tpcs = [0,1] 
left = 0 #left bound value
right = 1 #right bound value

cnter = -1 #This is to keep track of bw label for plot
for bw in bws:
  cnter += 1
  for pmt in pmts:
    for tpc in tpcs:
      for i in range(10): #iterate through tb
        if tpc == 0:
          df = op_df_tpc0.copy()
        elif tpc == 1:
          df = op_df_tpc1.copy()
        else:
          df = op_df_tpc0.copy()
        index = df.index.drop_duplicates()[i]
        bincenters,yvals,title = pmtplotters.get_xy_bins(df,'ophit_peakT','ophit_pe',index,bw,pmt=pmt,tpc=tpc)
        ax,fig = pmtplotters.make_bar_scatter_plot(bincenters,yvals,bw,title=title,
                left=left,right=right)
        ax.grid()
        ax.axvline(x=0,ls='--',c='k')
        ax.axvline(x=0.25,ls='--',c='k')
        plotters.save_plot(f'pe_time_pmt{pmt}_tpc{tpc}_bw{cnter}_event{i}')
        plt.close()

"""
#Iterable params
bws = [0.05,0.2,0.5] #binwidth: [nleft,nright]
pmts = [0,1,166,167] #Test coating, no coating and specific channels
tpcs = [0,1] 
left = -1 #left bound value
right = 8 #right bound value

cnter = -1 #This is to keep track of bw label for plot
for bw in bws:
  cnter += 1
  for pmt in pmts:
    for tpc in tpcs:
      for i in range(2): #iterate through tb
        if tpc == 0:
          df = op_df_tpc0.copy()
        elif tpc == 1:
          df = op_df_tpc1.copy()
        else:
          df = op_df_tpc0.copy()
        index = df.index.drop_duplicates()[i]
        bincenters,yvals,title = pmtplotters.get_xy_bins(df,'ophit_peakT','ophit_pe',index,bw,pmt=pmt,tpc=tpc)
        ax,fig = pmtplotters.make_bar_scatter_plot(bincenters,yvals,bw,title=title,
                left=left,right=right)
        ax.grid()
        plotters.save_plot(f'pe_time_pmt{pmt}_tpc{tpc}_bw{cnter}_event{i}')
        plt.close()
"""






