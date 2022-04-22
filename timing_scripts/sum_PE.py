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
pmts = pd.read_pickle('data/PMT_info.pkl')

#Parameters
bw = 0.0025 #time step

channels = pmts.loc[:,'ophit_opdet'].drop_duplicates().values #Channels
#Run information
indeces_tpc0 = op_df_tpc0.index.drop_duplicates().values
indeces_tpc1 = op_df_tpc1.index.drop_duplicates().values

columns = ['run','subrun','event','ophit_ch','ophit_opdet_type',
          'op_tpc','tleft','tright','summed_PE','tot_PE','ophit_opdet_x',
          'ophit_opdet_y','ophit_opdet_z']

tleft = 0.1
trights = [0.21,0.22,0.225,0.23,0.24,0.25,1000] #Time in us to sum over this window of time

for tind,tright in enumerate(trights):
  #TPC0          
  arr_tpc0 = np.zeros((len(indeces_tpc0),len(channels),len(columns)))
  for i,index in enumerate(indeces_tpc0): #Iterate over events
    for j,ch in enumerate(channels):  #Iterate over channels
      bincenters,yvals,_ = pmtplotters.get_xy_bins(op_df_tpc0,'ophit_peakT','ophit_pe',index,bw,pmt=ch,tpc=0)  
      #Find total PE in time window
      boolean_array = np.logical_and(bincenters >= tleft, bincenters <= tright) #find indeces between left and right
      inds_in_range = np.where(boolean_array)[0]
      yvals_trunc=yvals[inds_in_range]
      arr_tpc0[i][j][0] = index[0] #run
      arr_tpc0[i][j][1] = index[1] #subrun
      arr_tpc0[i][j][2] = index[2] #event
      arr_tpc0[i][j][3] = ch #channel/pmt
      #Get coating type
      arr_tpc0[i][j][4] = pmts[pmts.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_type'].drop_duplicates().values[0]
      #Get PMT tpc
      arr_tpc0[i][j][5] = op_df_tpc0[op_df_tpc0.loc[:,'ophit_opdet'] == ch].loc[:,'op_tpc'].drop_duplicates().values[0]
      arr_tpc0[i][j][6] = tleft #left bound
      arr_tpc0[i][j][7] = tright #right bound
      arr_tpc0[i][j][8] = yvals_trunc.sum() #Bounded sum of PE
      arr_tpc0[i][j][9] = yvals.sum() #Total sum of PE
      arr_tpc0[i][j][10] = op_df_tpc0[op_df_tpc0.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_x'].drop_duplicates().values[0]
      arr_tpc0[i][j][11] = op_df_tpc0[op_df_tpc0.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_y'].drop_duplicates().values[0]
      arr_tpc0[i][j][12] = op_df_tpc0[op_df_tpc0.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_z'].drop_duplicates().values[0]


  #Convert to dataframe
  arr_tpc0 = arr_tpc0.reshape(len(indeces_tpc0)*len(channels),len(columns))
  summedPE_tpc0_df = pd.DataFrame(arr_tpc0,columns=columns)
  summedPE_tpc0_df = summedPE_tpc0_df.set_index(['run','subrun','event'])
  #Save pkl files
  summedPE_tpc0_df.to_pickle(f'data/summedPE_tpc0_tright{tind}df.pkl')


