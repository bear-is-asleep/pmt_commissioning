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

#Select which sample to look at, ew or fb
argList = sys.argv[1:]
which_sample = argList[0]
if which_sample == 'ew': #east west
  op_df = pd.read_pickle('data/op_ew_df.pkl')
  muon_df = pd.read_pickle('data/muon_ew_df.pkl')
elif which_sample == 'fb': #front back
  op_df = pd.read_pickle('data/op_fb_df.pkl')
  muon_df = pd.read_pickle('data/muon_fb_df.pkl')
else: #Error
  op_df = pd.read_pickle('data/op_fb_df.pkl')
  muon_df = pd.read_pickle('data/muon_df_fb.pkl')
  #raise Exception('You need to pass the sample test, either ew or fb')

#Run information
indeces = op_df.index.drop_duplicates().values #Set indeces to 

#Load dataframes
#pmts = pd.read_pickle('data/PMT_info.pkl')
pmts = pd.read_pickle('/sbnd/data/users/brindenc/analyze_sbnd/PDS/PMT_ARAPUCA_info.pkl')

#Parameters
bw = 0.002 #time step 0.002 is minimum
#Set window size, number of binwidths to left and right of peakT
leftshift = 2
rightshift = 2
number_pmts=5 #Number of pmts that need to be greater than the threshold (10)
thresholdPE=200 #PE for PMT to pass threshold (100)
learn_rate = 25 #Decrease threshold by this amount for none type (20)
minthreshold = 100 #minimum accepted threshold (20)
tbins = 50 #Number of tbins to use when checking t0 from muon track dataproduct

channels = pmts.loc[:,'ophit_opdet'].drop_duplicates().values #Channels
#indeces = indeces[:10] #Grab first 10 events for now
#indeces_tpc1 = op_df_tpc1.index.drop_duplicates().values

columns = ['run','subrun','event','ophit_ch','ophit_opdet_type',
          'op_tpc','tleft','tot_PE','ophit_opdet_x',
          'ophit_opdet_y','ophit_opdet_z']
vector_columns = ['summed_PE','tright','run','subrun','event','ophit_ch']

#window_size=0.04 #Time slice window to look at (40 us), must be evenly divisible by bw
trights_size = leftshift+rightshift+1 #Length of trights array


#TPC0       
muon_df0 = muon_df[muon_df.loc[:,'muontrk_tpc'] == 0].sort_index()
muon_indeces0 = muon_df0.index.drop_duplicates()
scalar_arr0 = np.zeros((len(indeces),len(channels),len(columns))) #Array with single value for each element
vector_arr0 = np.zeros((len(indeces),len(channels),trights_size,len(vector_columns))) #Array with vector for each element

#TPC1
muon_df1 = muon_df[muon_df.loc[:,'muontrk_tpc'] == 1].sort_index()
muon_indeces1 = muon_df1.index.drop_duplicates()
scalar_arr1 = np.zeros((len(indeces),len(channels),len(columns))) #Array with single value for each element
vector_arr1 = np.zeros((len(indeces),len(channels),trights_size,len(vector_columns))) #Array with vector for each element

#Temporary limit events to save time
cnt = 0
pic.print_stars()
#print(indeces,muon_df0.index.drop_duplicates().values)
for i,index in enumerate(indeces): #Iterate over events
  if i < 5: continue #temp skip
  event_start = time() #Time the event
  #Checks if type is float or array
  to_value0=True
  to_value1=True
  #Checks if there's a muon track in this TPC
  has_muon0=True
  has_muon1=False #temp set to false
  #Get muon info for this run
  if index in muon_indeces0:
    muon0 = muon_df0.loc[index] 
  else:
    has_muon0 = False
  if index in muon_indeces1:
    muon1 = muon_df1.loc[index] 
  else:
    has_muon1 = False
  #Check one TPC at a time

  #TPC 0 - is there a muon there?
  if has_muon0:
    if isinstance(muon0,pd.core.series.Series): #Can't convert float to .values
      to_value0=False
    if to_value0:
      peakT0s = muon0['muontrk_t0'].values
    else:
      peakT0s = muon0['muontrk_t0']
    
    #trights are timing bin ranges to sum over 
    trights0 = np.zeros((len(peakT0s),leftshift+rightshift+1)) #+1 for love  
    for i in range(len(peakT0s)):
      #tright0 = trights0[i] 
      tright0 = np.arange(peakT0s[i]-leftshift*bw,peakT0s[i]+rightshift*bw,bw)
      #print(len(trights0[i]),len(tright0))

      #Check for mismatch in size
      if len(tright0) == trights_size+1:
        #print('I did  it')
        trights0 = np.delete(trights0,-1,axis=1) #remove last element for all bins (not a big deal to lose 2 ns)
      if len(tright0) == trights_size-1:
        #print('WTF')
        tright0 = np.insert(tright0,0,peakT0s[i]-(leftshift-1)*bw) #add a another bin to the left, this also shouldn't do anything bad
      trights0[i] = tright0
    for k,peakT0 in enumerate(peakT0s):
      for j,ch in enumerate(channels):  #Iterate over channels
        channel_start = time()

        #TPC0 for each channel
        bincenters0,yvals0,_ = pmtpic.get_xy_bins(op_df,'ophit_peakT','ophit_pe',index,bw,pmt=ch,tpc=0,xmin=trights0[k,0],xmax=trights0[k,-1])
        print(trights0,trights0.shape)
        print(ch,bincenters0, yvals0,yvals0.shape,bincenters0.shape)
        #Vector arr - TPC0
        vector_arr0[i][j][:,0] = yvals0 #PE in region for channel
        vector_arr0[i][j][:,1] = trights0[k] #right bound values
        vector_arr0[i][j][:,2] = np.full(len(yvals0),index[0]) #run
        vector_arr0[i][j][:,3] = np.full(len(yvals0),index[1]) #subrun
        vector_arr0[i][j][:,4] = np.full(len(yvals0),index[2]) #event 
        vector_arr0[i][j][:,5] = np.full(len(yvals0),ch) #ch
        vector_arr0[i][j][:,6] = np.full(len(yvals0),k) #Track number
        #Scalar arr
        scalar_arr0[i][j][0] = index[0] #run
        scalar_arr0[i][j][1] = index[1] #subrun
        scalar_arr0[i][j][2] = index[2] #event
        scalar_arr0[i][j][3] = ch #channel/pmt
        #Get coating type
        scalar_arr0[i][j][4] = pmts[pmts.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_type'].drop_duplicates().values[0]
        #Get PMT tpc
        scalar_arr0[i][j][5] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'op_tpc'].drop_duplicates().values[0]
        scalar_arr0[i][j][6] = peakT0 #left bound
        scalar_arr0[i][j][7] = yvals0.sum() #Total sum of PE
        scalar_arr0[i][j][8] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_x'].drop_duplicates().values[0]
        scalar_arr0[i][j][9] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_y'].drop_duplicates().values[0]
        scalar_arr0[i][j][10] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_z'].drop_duplicates().values[0]
        scalar_arr0[i][j][11] = k #Track number



  if has_muon1:
    if isinstance(muon1,pd.core.series.Series): 
      to_value1=False
    if to_value1:
      peakT1s = muon1['muontrk_t0'].values
      trights1 = np.zeros(len(peakT1s,leftshift+rightshift))
      for i in range(len(peakT1s)):
        trights1[i] = np.arange(peakT1s[i]-leftshift*bw,peakT1s[i]+rightshift*bw,bw)
    else:
      peakT1s = muon1['muontrk_t0']
      trights1 = np.arange(peakT1s-leftshift*bw,peakT1s+rightshift*bw,bw)


  #if index == (1,26,586) or index == (1,27,18): continue #Temp check
  #cnt+=1
  #if cnt > 5: break #Break out of loop after x events
  
  #Get first time that PMT sees light - TPC0
  #peakT0 = pmtpic.get_peakT(op_df,pmts,0,index,bw,number_pmts=number_pmts,thresholdPE=thresholdPE) #Get peakT from function
  
  
  #temp_threshold = thresholdPE #Set temp threshold to handle events that do not pass, we will decrease it until we get the right peakT
  #while peakT0 is None:
  #  if temp_threshold <= minthreshold:
  #    peakT0 = -9999 #Set dumby value for pmts that don't see the proper amount of light
  #    break
  #  temp_threshold = temp_threshold - learn_rate
  #  peakT0 = pmtpic.get_peakT(op_df,pmts,0,index,bw,number_pmts=number_pmts,thresholdPE=temp_threshold) #Get peakT from function
  



  #trights0 = np.arange(peakT0-leftshift*bw,peakT0+rightshift*bw,bw) #Time in us to sum over this window of time

  #Get first time that PMT sees light - TPC1
  #peakT1 = pmtpic.get_peakT(op_df,pmts,1,index,bw,number_pmts=number_pmts,thresholdPE=thresholdPE) #Get peakT from function
  
  #temp_threshold = thresholdPE #Set temp threshold to handle events that do not pass, we will decrease it until we get the right peakT
  #while peakT1 is None:
  #  if temp_threshold <= minthreshold:
  #    peakT1 = -9999 #Set dumby value for pmts that don't see the proper amount of light
  #    break
  #  temp_threshold = temp_threshold - learn_rate
  #  peakT1 = pmtpic.get_peakT(op_df,pmts,1,index,bw,number_pmts=number_pmts,thresholdPE=temp_threshold) #Get peakT from function
  #trights1 = np.arange(peakT1-leftshift*bw,peakT1+rightshift*bw,bw) #Time in us to sum over this window of time

  #Error handling for cases where trights(0,1) > trights_size. We will remove the last element
  #if len(trights0) == trights_size+1:
  #  trights0 = np.delete(trights0,-1) #remove last element
  #if len(trights1) == trights_size+1:
  #  trights1 = np.delete(trights1,-1) #remove last element
  
        
      #firstT = bincenters[0]-bw/2 #First populated time bin edge
      #yvals_trunc0 = np.zeros(trights_size) #Initialize y_vals, these will store PE for each tright
      # for tind,tright in enumerate(trights0):
      #   tright_start = time()
      #   if tright > peakT0:
      #     #Find total PE in time window
      #     boolean_array = np.logical_and(bincenters0 >= peakT0, bincenters0 <= tright) #find indeces between left and right
      #     inds_in_range = np.where(boolean_array)[0]
      #     yvals_trunc0[tind] = yvals0[inds_in_range].sum()

    #TPC1 for each channel
    # bincenters1,yvals1,_ = pmtplotters.get_xy_bins(op_df,'ophit_peakT','ophit_pe',index,bw,pmt=ch,tpc=1)  
    # #firstT = bincenters[0]-bw/2 #First populated time bin edge
    # yvals_trunc1 = np.zeros(trights_size) #Initialize y_vals, these will store PE for each tright
    # for tind,tright in enumerate(trights1):
    #   tright_start = time()
    #   if tright > peakT1:
    #     #Find total PE in time window
    #     boolean_array = np.logical_and(bincenters1 >= peakT1, bincenters1 <= tright) #find indeces between left and right
    #     inds_in_range = np.where(boolean_array)[0]
    #     yvals_trunc1[tind] = yvals1[inds_in_range].sum()
    # tright_end = time()
    #print(f'Bin sum time: {tright_end-tright_start} s')

    #Vector arr - TPC0
    # vector_arr0[i][j][:,0] = yvals_trunc0 #PE in region for channel
    # vector_arr0[i][j][:,1] = trights0 #right bound values
    # vector_arr0[i][j][:,2] = np.full(len(trights0),index[0]) #run
    # vector_arr0[i][j][:,3] = np.full(len(trights0),index[1]) #subrun
    # vector_arr0[i][j][:,4] = np.full(len(trights0),index[2]) #event 
    # vector_arr0[i][j][:,5] = np.full(len(trights0),ch) #ch
    # #Scalar arr
    # scalar_arr0[i][j][0] = index[0] #run
    # scalar_arr0[i][j][1] = index[1] #subrun
    # scalar_arr0[i][j][2] = index[2] #event
    # scalar_arr0[i][j][3] = ch #channel/pmt
    # #Get coating type
    # scalar_arr0[i][j][4] = pmts[pmts.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_type'].drop_duplicates().values[0]
    # #Get PMT tpc
    # scalar_arr0[i][j][5] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'op_tpc'].drop_duplicates().values[0]
    # scalar_arr0[i][j][6] = peakT0 #left bound
    # scalar_arr0[i][j][7] = yvals0.sum() #Total sum of PE
    # scalar_arr0[i][j][8] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_x'].drop_duplicates().values[0]
    # scalar_arr0[i][j][9] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_y'].drop_duplicates().values[0]
    # scalar_arr0[i][j][10] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_z'].drop_duplicates().values[0]

    # #Vector arr - TPC1
    # vector_arr1[i][j][:,0] = yvals_trunc1 #PE in region for channel
    # vector_arr1[i][j][:,1] = trights1 #right bound values
    # vector_arr1[i][j][:,2] = np.full(len(trights1),index[0]) #run
    # vector_arr1[i][j][:,3] = np.full(len(trights1),index[1]) #subrun
    # vector_arr1[i][j][:,4] = np.full(len(trights1),index[2]) #event 
    # vector_arr1[i][j][:,5] = np.full(len(trights1),ch) #ch
    # #Scalar arr
    # scalar_arr1[i][j][0] = index[0] #run
    # scalar_arr1[i][j][1] = index[1] #subrun
    # scalar_arr1[i][j][2] = index[2] #event
    # scalar_arr1[i][j][3] = ch #channel/pmt
    # #Get coating type
    # scalar_arr1[i][j][4] = pmts[pmts.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_type'].drop_duplicates().values[0]
    # #Get PMT tpc
    # scalar_arr1[i][j][5] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'op_tpc'].drop_duplicates().values[0]
    # scalar_arr1[i][j][6] = peakT1 #left bound
    # scalar_arr1[i][j][7] = yvals1.sum() #Total sum of PE
    # scalar_arr1[i][j][8] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_x'].drop_duplicates().values[0]
    # scalar_arr1[i][j][9] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_y'].drop_duplicates().values[0]
    # scalar_arr1[i][j][10] = op_df[op_df.loc[:,'ophit_opdet'] == ch].loc[:,'ophit_opdet_z'].drop_duplicates().values[0]
    channel_end = time()
    #print(f'Channel {ch} time: {channel_end-channel_start} s')
  event_end = time()
  print(f'Run {index[0]} Subrun {index[1]} Event {index[2]} time: {event_end-event_start:.2f} s')


#Convert to dataframe - TPC0
scalar_arr0 = scalar_arr0.reshape(len(indeces)*len(channels),len(columns))
vector_arr0 = vector_arr0.reshape(len(indeces)*len(channels)*len(trights0[0]),len(vector_columns))
scalar_df0 = pd.DataFrame(scalar_arr0,columns=columns)
vector_df0 = pd.DataFrame(vector_arr0,columns=vector_columns)

#Convert to dataframe - TPC1
scalar_arr1 = scalar_arr1.reshape(len(indeces)*len(channels),len(columns))
vector_arr1 = vector_arr1.reshape(len(indeces)*len(channels)*len(trights1),len(vector_columns))
scalar_df1 = pd.DataFrame(scalar_arr1,columns=columns)
vector_df1 = pd.DataFrame(vector_arr1,columns=vector_columns)

#TPC0
scalar_df0 = scalar_df0.set_index(['run','subrun','event'])
vector_df0 = vector_df0.set_index(['run','subrun','event'])
#Save pkl files
scalar_df0.to_pickle(f'data/scalarPE0_{which_sample}_df.pkl')
vector_df0.to_pickle(f'data/vectorPE0_{which_sample}_df.pkl')

#TPC1
scalar_df1 = scalar_df1.set_index(['run','subrun','event'])
vector_df1 = vector_df1.set_index(['run','subrun','event'])
#Save pkl files
scalar_df1.to_pickle(f'data/scalarPE1_{which_sample}_df.pkl')
vector_df1.to_pickle(f'data/vectorPE1_{which_sample}_df.pkl')







