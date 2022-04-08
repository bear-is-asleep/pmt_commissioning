from matplotlib.style import library
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from datetime import date
import sys
sys.path.append('/sbnd/app/users/brindenc/mypython')
from numpy import pi,exp,sqrt
import uproot
from bc_utils.utils import pic,plotters
from time import time
import matplotlib
import awkward as ak
from sklearn.linear_model import LinearRegression
import seaborn as sn

plotters.plot_stuff()
plt.style.use(['science','no-latex']) #Turn this off if you don't have science style plots
cwd = '/sbnd/app/users/brindenc/analyze_sbnd/PMT/v09_43_00/ew_filtered/read_model/' #Use this if u cant tie ur shoez
cwd = os.getcwd()+'/'

#Get Dataframes
pmt_hits_df_trunc = pd.read_pickle(cwd+'pmt_hits_df_trunc.pkl')
pmt_hits_df = pd.read_pickle(cwd+'pmt_hits_df.pkl')
mtrks = pd.read_pickle(cwd+'muon_tracks.pkl').drop_duplicates()


#TPC plots
pic.print_stars()
print('HEY WE\'RE MAKING THESE PLOTS HERE!!')
for tpc in range(2):
  fig,ax,sc = plotters.plot_TPC(tpc,'pred_err','$ (N_{pred}-N_{obs})/N_{obs}$ cmax=3',pmt_hits_df_trunc,return_plot=True)
  sc.set_clim(vmax=3)
  plotters.save_plot(f'pred_err_TPC{tpc}_cmax3')
  plt.close()
  fig,ax,sc = plotters.plot_TPC(tpc,'pred_err','$ (N_{pred}-N_{obs})/N_{obs}$ cmax=1',pmt_hits_df_trunc,return_plot=True)
  sc.set_clim(vmax=1)
  plotters.save_plot(f'pred_err_TPC{tpc}_cmax1')
  plt.close()
  fig,ax = plotters.plot_TPC(tpc,'ophit_pred','$N_{pred}$',pmt_hits_df_trunc)
  plotters.save_plot(f'ophit_pred_TPC{tpc}')
  plt.close()
  fig1,ax1 = plotters.plot_TPC(tpc,'ophit_obs','$N_{obs}$',pmt_hits_df_trunc)
  plotters.save_plot(f'ophit_obs_TPC{tpc}')
  plt.close()
  fig,ax = plotters.plot_TPC(tpc,'pred_err','$ (N_{pred}-N_{obs})/N_{obs}$',pmt_hits_df_trunc)
  plotters.save_plot(f'pred_err_TPC{tpc}')
  plt.close()
  fig,ax = plotters.plot_TPC(tpc,'ophit_pred_noatt','$N_{pred} \ R_I(x)=1$',pmt_hits_df_trunc)
  plotters.save_plot(f'ophit_pred_TPC{tpc}_noatt')
  plt.close()
  fig,ax = plotters.plot_TPC(tpc,'pred_err_noatt','$ (N_{pred}-N_{obs})/N_{obs} \ R_I(x)=1$',pmt_hits_df_trunc)
  plotters.save_plot(f'pred_err_TPC{tpc}_noatt')
  plt.close()
  fig,ax = plotters.plot_TPC(tpc,'ophit_pred_nodisp','$N_{pred} \ d(x)=1$',pmt_hits_df_trunc)
  plotters.save_plot(f'ophit_pred_TPC{tpc}_nodisp')
  plt.close()
  fig,ax = plotters.plot_TPC(tpc,'pred_err_nodisp','$ (N_{pred}-N_{obs})/N_{obs} \ d(x)=1$',pmt_hits_df_trunc)
  plotters.save_plot(f'pred_err_TPC{tpc}_nodisp')
  plt.close()
  fig,ax = plotters.plot_TPC(tpc,'ophit_pred_noref','$N_{pred}$ (No reflections)',pmt_hits_df_trunc)
  plotters.save_plot(f'ophit_pred_TPC{tpc}_noref')
  plt.close()
  fig,ax = plotters.plot_TPC(tpc,'pred_err_noref','$ (N_{pred}-N_{obs})/N_{obs}$ (No reflections)',pmt_hits_df_trunc)
  plotters.save_plot(f'pred_err_TPC{tpc}_noref')
  plt.close()
  

#Plot 15 tracks
mtrks_trunc = mtrks.iloc[:15]
plotters.plot_tracks(mtrks_trunc,'muontrk_x1_0','muontrk_y1_0','muontrk_x2_0','muontrk_y2_0','muontrk_x1_1','muontrk_y1_1','muontrk_x2_1','muontrk_y2_1');
plotters.save_plot('xy_trks');
plt.close()
plotters.plot_tracks(mtrks_trunc,'muontrk_x1_0','muontrk_z1_0','muontrk_x2_0','muontrk_z2_0','muontrk_x1_1','muontrk_z1_1','muontrk_x2_1','muontrk_z2_1');
plotters.save_plot('xz_trks');
plt.close()
plotters.plot_tracks(mtrks_trunc,'muontrk_z1_0','muontrk_y1_0','muontrk_z2_0','muontrk_y2_0','muontrk_z1_1','muontrk_y1_1','muontrk_z2_1','muontrk_y2_1');
plotters.save_plot('zy_trks');
plt.close()

#Save params to csv
#fit_params = pd.read_pickle(cwd+'fit_params.pkl')
#fit_params.to_csv(cwd+'PMT_fit_params.csv')



