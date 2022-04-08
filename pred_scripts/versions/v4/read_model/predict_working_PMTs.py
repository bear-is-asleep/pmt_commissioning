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

#Load DFs
pmt_hits_df_trunc = pd.read_pickle(cwd+'pmt_hits_df_trunc.pkl')
pmt_hits_df = pd.read_pickle(cwd+'pmt_hits_df.pkl')
mtrks = pd.read_pickle(cwd+'muon_tracks.pkl').drop_duplicates()
fit_params = pd.read_pickle(cwd+'fit_params.pkl')

#Set bad channels
bad_chs = [145,117,144,116,226,305]
reduce = 0.4 #Reduce observed hits by this amount
bad_df = pic.fake_bad_pmt(pmt_hits_df_trunc,bad_chs,reduce=reduce)

#Check functionality
good_df = pic.check_functionality(pmt_hits_df_trunc)
bad_df = pic.check_functionality(bad_df)

#Make diagnostic plots
for tpc in range(2):
  fig,ax = plotters.plot_TPC(tpc,'functioning','Properly functioning PMTs (fake bad PMT)',bad_df)
  plotters.save_plot(f'bad_pmt_functionality{tpc}')
  plt.close()
  fig,ax = plotters.plot_TPC(tpc,'functioning','Properly functioning PMTs',good_df)
  plotters.save_plot(f'good_pmt_functionality{tpc}')
  plt.close()
  fig1,ax1 = plotters.plot_TPC(tpc,'ophit_obs','$N_{obs}$ (fake bad PMT)',bad_df)
  plotters.save_plot(f'ophit_obs_TPC{tpc}_badpmts')
  plt.close()
  fig1,ax1 = plotters.plot_TPC(tpc,'pred_err','$ (N_{pred}-N_{obs})/N_{obs}$ (fake bad PMT)',bad_df)
  plotters.save_plot(f'pred_err_TPC{tpc}_badpmts')
  plt.close()












