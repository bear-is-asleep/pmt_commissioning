#Make dataframe to sum up bins
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('/sbnd/app/users/brindenc/mypython') #My utils path
from bc_utils.pmtutils import pic as pmtpic
from bc_utils.pmtutils import plotters as pmtplotters
from bc_utils.utils import pic,plotters
import pandas as pd
from time import time

#Load dataframes
op_df = pd.read_pickle('data/op_df_tpc0.pkl')
#op_df_tpc1 = pd.read_pickle('data/op_df_tpc1.pkl')
pmts = pd.read_pickle('data/PMT_info.pkl')

#Parameters
bw = 0.002 #time step 0.002 is minimum
thresholdPE = 5 #Min threshold for light to be seen in an event

channels = pmts.loc[:,'ophit_opdet'].drop_duplicates().values #Channels
#Run information
indeces = op_df.index.drop_duplicates().values
indeces = indeces[:10] #Grab first 10 events for now
#indeces_tpc1 = op_df_tpc1.index.drop_duplicates().values

columns = ['run','subrun','event','ophit_ch','ophit_opdet_type',
          'op_tpc','tleft','tot_PE','ophit_opdet_x',
          'ophit_opdet_y','ophit_opdet_z']
vector_columns = ['summed_PE','tright','run','subrun','event','ophit_ch']

window_size=0.02 #Time slice window to look at (100 us), must be evenly divisible by bw
trights_size = int(window_size/bw) #Length of trights array


#TPC0          
scalar_arr = np.zeros((len(indeces),len(channels),len(columns))) #Array with single value for each element
vector_arr = np.zeros((len(indeces),len(channels),trights_size,len(vector_columns))) #Array with vector for each element
for i,index in enumerate(indeces): #Iterate over events
  event_start = time()
  #Get first time that PMT sees light
  bincenters,allPE,_ = pmtplotters.get_xy_bins(op_df,'ophit_peakT','ophit_pe',index,bw,pmt=2,tpc=0)
  peakT = bincenters[np.argmax(allPE>thresholdPE)]-3*bw/2 #Window befor first time where PMT sees light
  trights = np.arange(peakT,peakT+window_size,bw) #Time in us to sum over this window of time
  for j,ch in enumerate(channels):  #Iterate over channels
    channel_start = time()
    bincenters,yvals,_ = pmtplotters.get_xy_bins(op_df,'ophit_peakT','ophit_pe',index,bw,pmt=ch,tpc=0)  
    #firstT = bincenters[0]-bw/2 #First populated time bin edge
    yvals_trunc = np.zeros(len(trights)) #Initialize y_vals, these will store PE for each tright
    for tind,tright in enumerate(trights):
      tright_start = time()
      if tright > peakT:
        #Find total PE in time window
        boolean_array = np.logical_and(bincenters >= peakT, bincenters <= tright) #find indeces between left and right
        inds_in_range = np.where(boolean_array)[0]
        yvals_trunc[tind] = yvals[inds_in_range].sum()
    tright_end = time()
    #print(f'Bin sum time: {tright_end-tright_start} s')
    #Vector arr
    vector_arr[i][j][:,0] = yvals_trunc #PE in region for channel
    vector_arr[i][j][:,1] = trights #right bound values
    vector_arr[i][j][:,2] = np.full(len(trights),index[0]) #run
    vector_arr[i][j][:,3] = np.full(len(trights),index[1]) #subrun
    vector_arr[i][j][:,4] = np.full(len(trights),index[2]) #event 
    vector_arr[i][j][:,5] = np.full(len(trights),ch) #ch
    #Scalar arr
    scalar_arr[i][j][0] = index[0] #run
    scalar_arr[i][j][1] = index[1] #subrun
    scalar_arr[i][j][2] = index[2] #event
    scalar_arr[i][j][3] = ch #channel/pmt
    #Get coating type
    scalar_arr[i][j][4] = pmts[pmts.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_type'].drop_duplicates().values[0]
    #Get PMT tpc
    scalar_arr[i][j][5] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'op_tpc'].drop_duplicates().values[0]
    scalar_arr[i][j][6] = peakT #left bound
    #scalar_arr[i][j][7] = trights #right bound
    #scalar_arr[i][j][8] = yvals_trunc #Bounded sum of PE
    scalar_arr[i][j][7] = yvals.sum() #Total sum of PE
    scalar_arr[i][j][8] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_x'].drop_duplicates().values[0]
    scalar_arr[i][j][9] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_y'].drop_duplicates().values[0]
    scalar_arr[i][j][10] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_z'].drop_duplicates().values[0]
    channel_end = time()
    #print(f'Channel {ch} time: {channel_end-channel_start} s')
  event_end = time()
  print(f'Run {index[0]} Subrun {index[1]} Event {index[2]} time: {event_end-event_start:.2f} s')


#Convert to dataframe
scalar_arr = scalar_arr.reshape(len(indeces)*len(channels),len(columns))
vector_arr = vector_arr.reshape(len(indeces)*len(channels)*len(trights),len(vector_columns))
scalar_df = pd.DataFrame(scalar_arr,columns=columns)
vector_df = pd.DataFrame(vector_arr,columns=vector_columns)
#summedPE_df = pd.concat([scalar_df,vector_df],axis=1) #Vertically concatenate two arrays


scalar_df = scalar_df.set_index(['run','subrun','event'])
vector_df = vector_df.set_index(['run','subrun','event'])
#Save pkl files
scalar_df.to_pickle(f'data/scalarPE_df.pkl')
vector_df.to_pickle(f'data/vectorPE_df.pkl')
#Dataframes can't handle multidimensional data of different sizes very well, so we'll just make two of these to save data





