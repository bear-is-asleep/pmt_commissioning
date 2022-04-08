from locale import normalize
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('/sbnd/app/users/brindenc/mypython') #My utils path
from bc_utils.pmtutils import pic as pmtpic
from bc_utils.pmtutils import plotters as pmtplotters
from bc_utils.hitutils import pic as hitpic
from bc_utils.hitutils import plotters as hitplotters
from bc_utils.utils import pic,plotters
import pandas as pd

plotters.plot_stuff()

#Load and clean data
PE_tpc0_tright0_df = pd.read_pickle('data/summedPE_tpc0_tright0df.pkl')
PE_tpc0_tright1_df = pd.read_pickle('data/summedPE_tpc0_tright1df.pkl')
PE_tpc0_tright2_df = pd.read_pickle('data/summedPE_tpc0_tright2df.pkl')
PE_tpc0_tright3_df = pd.read_pickle('data/summedPE_tpc0_tright3df.pkl')
PE_tpc0_tright4_df = pd.read_pickle('data/summedPE_tpc0_tright4df.pkl')
PE_tpc0_tright5_df = pd.read_pickle('data/summedPE_tpc0_tright5df.pkl')
muon_tpc0_df = pd.read_pickle('data/muon_df_tpc0.pkl')
muon_plot_tpc0_df = muon_tpc0_df.drop(['nmuontrks','muontrk_t0'],axis=1)
muon_plot_tpc0_df = pmtpic.get_muon_tracks(muon_plot_tpc0_df)
muon_plot_tpc0_df = muon_plot_tpc0_df.set_index(['run','subrun','event'])
hit_tpc0_df = pd.read_pickle('data/hit_df_tpc0.pkl')
hit_bw = 20


#Parameters
indeces = PE_tpc0_tright0_df.index.drop_duplicates()[:10]
tpcs = [0]
trights = [0.25,0.5,1,2,10,1000] #Time in us to sum over this window of time

#We save dataframes based on binwidth, I don't see a good reason to fix it, it works..
PE_tpc0_dfs = [PE_tpc0_tright0_df,PE_tpc0_tright1_df,
PE_tpc0_tright2_df,PE_tpc0_tright3_df,PE_tpc0_tright4_df,
PE_tpc0_tright5_df] #Enumerate over multiple binwidths

for tind,PE_tpc0_df in enumerate(PE_tpc0_dfs):
  for i,index in enumerate(indeces):
    for tpc in tpcs:
      fig,ax = pmtplotters.plot_TPC(tpc,'summed_PE',rf'PE for $t \in $[0,{trights[tind]}] $\mu$s'+f'\nSubrun {index[1]}, Event {index[2]}', 
      PE_tpc0_df.loc[index],normalize=True)
      ax,fig = pmtplotters.plot_tracks(muon_plot_tpc0_df.loc[index],'muontrk_z1_0','muontrk_y1_0','muontrk_z2_0','muontrk_y2_0',
            'muontrk_z1_1','muontrk_y1_1','muontrk_z2_1','muontrk_y2_1',ax=ax,fig=fig,indeces=[index])
      plotters.save_plot(f'summedPE_tpc{tpc}_trk{i}_tright{tind}')
      plt.close()
      if tind == 0:
        bincenters,yvals,title = hitpic.get_xy_bins(hit_tpc0_df,'hit_wire','hit_charge',index,hit_bw,tpc=0,plane=2)
        ax,fig = hitplotters.make_bar_scatter_plot(bincenters,yvals,hit_bw,title='')
        plotters.save_plot(f'Qwire_tpc{tpc}_trk{i}')
        plt.close()



