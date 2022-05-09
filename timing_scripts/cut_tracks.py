"""
Brinden Carlson
5/2/22

Make cuts based on tracks
"""

import sys
import os
import pandas as pd
sys.path.append('/sbnd/app/users/brindenc/mypython') #My utils path
from bc_utils.pmtutils import pic as pmtpic
from bc_utils.pmtutils import plotters as pmtplotters
from bc_utils.utils import pic,plotters
from time import time

#Constants
apa_threshold = 100 #x1 and x2 need to be greater than this to be considered 'near-apa'

#Initialize list of dataframes
op_ew_dfs = []
muon_ew_dfs = []
op_fb_dfs = []
muon_fb_dfs = []

#Load dfs - everything must happen in the for loop
for i in range(1,31):
  #Check if files exist
  ew = False 
  fb = False
  if os.path.exists(f'data/muon_df_ew__precut{i}.pkl') and os.path.exists(f'data/op_df_ew__precut{i}.pkl'):
    ew = True
  if os.path.exists(f'data/muon_df_fb__precut{i}.pkl') and os.path.exists(f'data/op_df_fb__precut{i}.pkl'):
    fb = True
  start = time()
  #East-west i know im duplicating code let me live
  if ew:
    indeces_drop = []
    muon_df_ew = pd.read_pickle(f'data/muon_df_ew__precut{i}.pkl')
    op_df_ew = pd.read_pickle(f'data/op_df_ew__precut{i}.pkl')

    for row,line in muon_df_ew.iterrows():
      if round(line['muontrk_x1'],1) == round(line['muontrk_x2'],1): #Drop rows with matching x1 and x2 (check nearest 1st digit, python has rounding issues)
        indeces_drop.append(row)
      if abs(line['muontrk_x1']) < apa_threshold or abs(line['muontrk_x2']) < apa_threshold: #Keep near-apa muons 
        indeces_drop.append(row)
      if line['nmuontrks'] > 1: #Only one track
        indeces_drop.append(row)
    
    indeces_drop = list(set(indeces_drop)) #Drop duplicate events
    op_ew_dfs.append(op_df_ew.drop(indeces_drop)) #Append good events to list
    muon_ew_dfs.append(muon_df_ew.drop(indeces_drop)) #Append good events to list

  #Front back
  if fb:
    indeces_drop = []
    muon_df_fb = pd.read_pickle(f'data/muon_df_fb__precut{i}.pkl')
    op_df_fb = pd.read_pickle(f'data/op_df_fb__precut{i}.pkl')

    for row,line in muon_df_fb.iterrows():
      if round(line['muontrk_x1'],1) == round(line['muontrk_x2'],1): #Drop rows with matching x1 and x2
        indeces_drop.append(row)
      if abs(line['muontrk_x1']) < apa_threshold or abs(line['muontrk_x2']) < apa_threshold: #Keep near-apa muons 
        indeces_drop.append(row)
      if line['nmuontrks'] > 1: #Only one track
        indeces_drop.append(row)
    
    indeces_drop = list(set(indeces_drop)) #Drop duplicate events
    op_fb_dfs.append(op_df_fb.drop(indeces_drop)) #Append good events to list
    muon_fb_dfs.append(muon_df_fb.drop(indeces_drop)) #Append good events to list

  end = time()
  print(f'File {i} elapsed(s): {end-start:.2f}')

#Concatenate all dataframes together
op_ew_df_all = pd.concat(op_ew_dfs)
muon_ew_df_all = pd.concat(muon_ew_dfs)
op_fb_df_all = pd.concat(op_fb_dfs)
muon_fb_df_all = pd.concat(muon_fb_dfs)

#Save dataframes
op_ew_df_all.to_pickle('data/op_ew_df.pkl')
muon_ew_df_all.to_pickle('data/muon_ew_df.pkl')
op_fb_df_all.to_pickle('data/op_fb_df.pkl')
muon_fb_df_all.to_pickle('data/muon_fb_df.pkl')




  



