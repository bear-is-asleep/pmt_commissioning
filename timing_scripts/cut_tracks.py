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

#Input sample name
argList = sys.argv[1:]
sample = argList[0]

#Constants
apa_threshold = 100 #x1 and x2 need to be greater than this to be considered 'near-apa'

#Initialize list of dataframes
op_dfs = []
muon_dfs = []
crt_dfs = []
g4_dfs = []

#Function to cut events
def muon_track_cuts(label,sample_type,op_dfs,muon_dfs,crt_dfs,g4_dfs):
  indeces_keep = []
  muon_df = pd.read_pickle(f'data/muon_df_{sample_type}__precut{label}.pkl')
  op_df = pd.read_pickle(f'data/op_df_{sample_type}__precut{label}.pkl')
  crt_df = pd.read_pickle(f'data/crt_df_{sample_type}__precut{label}.pkl')
  g4_df = pd.read_pickle(f'data/g4_df_{sample_type}__precut{label}.pkl')

  #Check event info
  op_inds = op_df.index.drop_duplicates().values
  crt_inds = crt_df.index.drop_duplicates().values
  g4_inds = g4_df.index.drop_duplicates().values

  for row,line in muon_df.iterrows():
    if round(line['muontrk_x1'],1) != round(line['muontrk_x2'],1): #Drop rows with matching x1 and x2 (check nearest 1st digit, python has rounding issues)
      indeces_keep.append(row)
    #if abs(line['muontrk_x1']) < apa_threshold or abs(line['muontrk_x2']) < apa_threshold: #Keep near-apa muons 
    #  indeces_drop.append(row)
    #if row in op_inds and row in crt_inds and row in g4_inds:
    #  indeces_keep.append(row)
    #if line['nmuontrks'] > 1: #Only one track
    #  indeces_drop.append(row)
  indeces_keep = list(set(indeces_keep)) #Drop duplicate events
  op_dfs.append(op_df.loc[indeces_keep]) #Append good events to list
  muon_dfs.append(muon_df.loc[indeces_keep]) #Append good events to list
  #crt_dfs.append(crt_df.loc[indeces_keep]) #Append good events to list
  g4_dfs.append(g4_df.loc[indeces_keep]) #Append good events to list

#Load dfs - everything must happen in the for loop
for i in range(1,6):
  #Check if files exist
  if not (os.path.exists(f'data/muon_df_{sample}__precut{i}.pkl') and 
  os.path.exists(f'data/op_df_{sample}__precut{i}.pkl') and 
  os.path.exists(f'data/crt_df_{sample}__precut{i}.pkl') and 
  os.path.exists(f'data/g4_df_{sample}__precut{i}.pkl')):
    continue #Keep going if even one of these files doesn't exist. Future versions will not have g4 stuff
  start = time()
  #Append dfs to list for existing files
  muon_track_cuts(i,sample,op_dfs,muon_dfs,crt_dfs,g4_dfs)

  end = time()
  print(f'File {i} elapsed(s): {end-start:.2f}')

#Concatenate all dataframes together
op_df_all = pd.concat(op_dfs)
muon_df_all = pd.concat(muon_dfs)
#crt_df_all = pd.concat(crt_dfs)
g4_df_all = pd.concat(g4_dfs)

#Save dataframes
op_df_all.to_pickle(f'data/op_{sample}_df.pkl')
muon_df_all.to_pickle(f'data/muon_{sample}_df.pkl')
#crt_df_all.to_pickle(f'data/crt_{sample}_df.pkl')
g4_df_all.to_pickle(f'data/g4_{sample}_df.pkl')







  



