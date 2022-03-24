import numpy as np
from random import choice
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import math
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from datetime import date
import sys
import seaborn as sns
import matplotlib
from scipy import optimize
from numpy import sqrt,exp,arctan,cos,sin,pi,arccos
import uproot

#Constants
hc = 1240 #eV nm
r_e = 2.818e-15 #m
alpha = 1/137 #alpha
m_e = 511e3 #eV/c^2
n_A = 1.3954/6.6335209e-23 #n_A/cm^3
m_u = 105.658 #MeV

def print_stars():
  print('\n******************************\n')

#Gaussian fit
def gaussian(x,a, mean, stddev):
  #Ex: popt, pcov = optimize.curve_fit(gaussian, x, data)
  return a*np.exp(-((x - mean) / 4 / stddev)**2)

def fit_gaussian(data):
  #Fits gaussian to 1d data and returns parameters
  steps=len(data)
  x = np.linspace(min(data),max(data),steps)
  popt, pcov = optimize.curve_fit(gaussian, x, data)
  return x,popt,pcov

def sum_combine(arr,every_other=2,col_options=[]):
  """
  DESC: Somes int=every_other for every other row
  arr: 2D array
  every_other: How many rows to combine, i.e. every_other=2 combines every other 
               line
  col_options: Enter column options 1D array same number of columns as arr.
               Set column to 0 for sum of column from i to i+every_other 
  """
  
  if arr.ndim == 1:
    arr = np.reshape(arr,(arr.size,1))
    #print (arr.shape)

  result = [] #Store result here
  extra_rows=[] #Store extra rows here
  if not col_options.any():
    col_options = np.zeros(arr.shape[1]) #Set colunn default
  for (i,line) in enumerate(arr): #Iterate rows
    if np.mod(i+1,every_other) == 0: #Check if it's row to sum over and store
      row = [] #Store row here
      for (j,ele) in enumerate(line):
        val = 0
        temp = []
        if col_options[j] == 0:
          for k in range(i-every_other+1,i+1):
            val+=arr[k,j]
          row.append(val)
      
        if col_options[j] == 1:
          for k in range(i-every_other+1,i+1):
            temp.append(arr[k,j])
          val = np.mean(temp)
          row.append(val)
        if col_options[j] == 2:
          for k in range(i-every_other+1,i+1):
            temp.append(arr[k,j])
          val = np.median(temp)

          row.append(val)
      result.append(row)
      extra_rows = []
    elif i == np.shape(arr)[0]-1:
     #print(i)
      extra_rows.append(line)
      return np.asarray(np.vstack((result,extra_rows))) #Returns extra rows
    else:
      extra_rows.append(line)
  return np.asarray(result)


def truncate_df(df,keys_to_sum = {''},keys_to_median = {''}):
  #Truncate df for average values across all events and runs, len(df) = number of PMTs
  cols = df.shape[1]
  col_options = np.ones(cols)
  
  #set sum_combine to 0 for sum, 1 for mean, 2 for med
  for i,key in enumerate(df.keys()):
    if key in keys_to_sum:
      col_options[i]=0
    if key in keys_to_median:
      col_options[i]=2
  #Extract all PMT channel information and store in single df
  all_trunc = []
  for i in range(int(df['ophit_ch'].min()),int(df['ophit_ch'].max())+1):
    df_temp = df.loc[df['ophit_ch']==i]
    if df_temp.empty:
      continue
    np_temp = df_temp.to_numpy()
    np_trunc = sum_combine(np_temp,every_other=np_temp.shape[0],col_options=col_options)
    all_trunc.append(np_trunc)
  all_trunc = np.squeeze(all_trunc)
  return  pd.DataFrame(data=all_trunc,columns=df.keys())

#Get total POT
def get_pot(rootname,treename):
    #Input nueback dataframe
    POT_tree = uproot.open(f'{rootname}:{treename}/POT;1')
    pot = POT_tree.arrays('potbnb',library='np')
    return pot['potbnb'].sum()

#Return g4 dataframe
def get_g4_df(rootname,treename,offset=0):
    #Input root and tree name to return dataframe with g4 information
    #Offset run creates indexing for additional backgrounds
    g4_tree = uproot.open(f'{rootname}:{treename}/Event;1')
    with uproot.open(f'{rootname}:{treename}/Event;1') as file:
        keys = file.keys()
    g4_keys = keys[0:28]
    g4 = g4_tree.arrays(g4_keys,library='pd')
    runs = list(g4.loc[:,'run'])
    runs = [run+offset for run in runs] #Add offset to run value for additional background
    subruns = list(g4.loc[:,'subrun'])
    events = list(g4.loc[:,'event'])
    if offset != 0:
      arrays = [runs,subruns,events]
      tuples = list(zip(*arrays)) #Make multiindex
      index = pd.MultiIndex.from_tuples(tuples, names=['run','subrun','event'])
      g4 = g4.drop(['run','subrun','event'],axis=1)
      columns = g4.columns
      g4_data = g4.values
      g4 = pd.DataFrame(g4_data,index=index,columns=columns)
    else:
      g4 = g4.set_index(['run','subrun','event'])
    indeces = list(g4.index.values) 
    return g4,indeces

#Return genie dataframe
def get_genie_df(rootname,treename):
    #Input root and tree name to return dataframe with genie information
    genie_tree = uproot.open(f'{rootname}:{treename}/Event;1')
    with uproot.open(f'{rootname}:{treename}/Event;1') as file:
        keys = file.keys()
    genie_keys = [key for key in keys if 'genie' in key]
    genie_keys.append('run')
    genie_keys.append('subrun')
    genie_keys.append('event')
    genie = genie_tree.arrays(genie_keys,library='pd')
    runs = list(genie.loc[:,'run'])
    subruns = list(genie.loc[:,'subrun'])
    events = list(genie.loc[:,'event'])
    genie = genie.set_index(['run','subrun','event'])
    indeces = list(genie.index.values) 
    return genie,indeces

#Return dataframe keys
def get_keys(rootname,treename):
  with uproot.open(f'{rootname}:{treename}/Event;1') as file:
    keys = file.keys()
  return keys


def no_particle(g4,pdg=2212):
    #Input g4 dataframe, return indeces that don't contain the particle
    g4 = g4.sort_index()
    indeces = []
    index_g4 = np.unique(g4.index.values)
    for index in index_g4:
        if pdg not in list(g4.loc[index,'pdg']): #If there are no particles
            indeces.append(index)
    return indeces

def calc_T(df):
    #df is g4 df, see g4 keys
    #We're going to be including T and theta_e
    E = df['Eng']
    m = df['Mass']
    df['T'] = sqrt(E**2-m**2)

def number_events(df):
    #Return number of events in a dataframe
    return df.index.drop_duplicates().shape[0]

def get_pot_normalized_df(df,target_pot,pot,events,seed=420):
  #Return dataframe with randomly selected events, normalized to target pot
  #Also returns total number of events
  np.random.seed(seed)
  n_keep = int(np.round(target_pot/pot*events))
  n_drop = events - n_keep
  if n_drop < 0:
    raise Exception('You need more events chief!')
  index = df.index.drop_duplicates()
  drop_indices = np.random.choice(index, n_drop, replace=False)
  return df.drop(drop_indices),n_keep

def get_eventype_count(df,pdg,pdg_key='pdg'):
  #Returns number of events with a pdg
  #Also returns indeces with pdg
  #Also returns fraction of events with pdg
  indeces = df.index.drop_duplicates()
  indeces = list(indeces)
  events = len(indeces)
  cnt = 0 #Counter for events
  has_pdg = []
  for index in indeces:
    vals = df.loc[index,pdg_key].values
    #print(vals)
    if pdg in vals:
      cnt += 1
      has_pdg.append(index)
  return cnt,has_pdg,cnt/events

def calc_thetat(df,method=0,return_key='theta_t',theta_xz_key='theta_xz',theta_yz_key='theta_yz',
px_key='Px',py_key='Py',pz_key='Pz',p_key='P'):
  #Calculate transverse angle (off z-axis) for given dataframe, appends it to dataframe
  #Method 0 is using angles
  #Method 1 is using momenta
  if return_key in df.columns:
    df.drop(return_key,axis=1)
    #print('Key already exists in dataframe silly!')
  if method == 0:
    #Get values of angles
    theta_xz = df.loc[:,theta_xz_key].values
    theta_yz = df.loc[:,theta_yz_key].values
    thetat = sqrt(theta_xz**2+theta_yz**2) #Not sure if this works, makes sense to me though
  elif method == 1:
    #Get momenta values
    px = df.loc[:,px_key].values
    py = df.loc[:,py_key].values
    pz = df.loc[:,pz_key].values
    if p_key == 'None':
      p = sqrt(px**2+py**2+pz**2)
    else:
      p = df.loc[:,p_key].values
    thetat = pi/2 - arccos(sqrt(px**2+py**2)/p)
  df.loc[:,return_key] = thetat
  return df

def calc_Etheta(df,return_key='E_theta^2',E_key='Eng',theta_t_key='theta_t'):
  #Calc E_etheta^2 for electrons
  if return_key in df.columns:
    df.drop(return_key,axis=1)
    #print('Key already exists in dataframe silly!')
  #Get values
  E = df.loc[:,E_key].values
  theta_t = df.loc[:,theta_t_key].values
  E_theta2 = E*theta_t**2
  df.loc[:,return_key] = E_theta2
  return df

def get_scat_type(df,pdg_key='pdg',return_key='scat_type'):
  #Appends scat type to dataframe
  if return_key in df.columns:
    df.drop(return_key,axis=1)
    #print('Key already exists in dataframe silly!')
  indeces = df.index.drop_duplicates()
  indeces = list(indeces)
  types = [] #Scattering types
  for index in indeces:
    temp_df = df.loc[index]
    #Apply pdg check
    pdgs = temp_df.loc[:,pdg_key].values
    if 14 in pdgs: #nu mu
      types.extend(np.full(len(pdgs),0))
    elif 12 in pdgs: #nu e
      types.extend(np.full(len(pdgs),1))
    if -14 in pdgs: #nu bar mu
      types.extend(np.full(len(pdgs),2))
    elif -12 in pdgs: #nu bar e
      types.extend(np.full(len(pdgs),3))
  df.loc[:,return_key] = types
  return df

def get_signal_background(scat,back):
  #Get signal to background ratio
  #print(len(list(scat.index.drop_duplicates())),len(list(back.index.drop_duplicates())))
  return len(list(scat.index.drop_duplicates()))/len(list(back.index.drop_duplicates()))

def make_cuts(df,pdg_key='pdg',method=0,Etheta=0.03,Etheta_key='E_theta^2'):
  #Make background cuts, returns cut index
  #Method 0: E theta^2 cut
  e_df = df[df.loc[:,pdg_key] == 11] #Filter by electron
  indeces = list(e_df.index.drop_duplicates()) #Get event indeces
  keep_indeces = [] #keep indeces if they pass the cut
  for index in indeces:
    if method == 0: #E theta^2 cut
      E_theta = e_df.loc[index,Etheta_key]
      #print(isinstance(E_theta,np.floating),type(E_theta))
      if isinstance(E_theta,np.floating):
        if e_df.loc[index,Etheta_key] < Etheta: #Less than cutoff
          keep_indeces.append(index)
      else:
        if e_df.loc[index,Etheta_key].values[0] < Etheta: #Less than cutoff
          #print(e_df.loc[index,Etheta_key].values[0] < Etheta,Etheta,e_df.loc[index,Etheta_key].values[0])
          keep_indeces.append(index)
  
  return df.loc[keep_indeces]

def get_electron_count(df,pdg_key='pdg',return_key='e_count',drop_duplicates=False):
  #Get electron count for each event
  indeces = df.index.drop_duplicates()
  indeces = list(indeces)
  cnts=[] #Count number of electrons in each event
  drop_index = [] #Drop indeces with multiple electrons
  for index in indeces:
    temp_df = df.loc[index]
    pdgs = list(temp_df.loc[:,pdg_key].values)
    es = pdgs.count(11) #Number of electrons
    if drop_duplicates and es > 1:
      drop_index.append(index)
      continue
    else: 
      cnts.extend(np.full(len(pdgs),es)) #append number of electrons in each event
  if drop_duplicates:
    df = df.drop(drop_index)
  df.loc[:,return_key] = cnts
  return df

#Return shower dataframe
def get_shw_df(rootname,treename,offset=0):
    #Input root and tree name to return dataframe with g4 information
    #Offset run creates indexing for additional backgrounds
    shw_tree = uproot.open(f'{rootname}:{treename}/Event;1')
    with uproot.open(f'{rootname}:{treename}/Event;1') as file:
        keys = file.keys()
    shw_keys = [key for key in keys if 'shw' in key]
    print(shw_keys)
    shw = shw_tree.arrays(shw_keys,library='pd')
    runs = list(shw.loc[:,'run'])
    runs = [run+offset for run in runs] #Add offset to run value for additional background
    subruns = list(shw.loc[:,'subrun'])
    events = list(shw.loc[:,'event'])
    if offset != 0:
      arrays = [runs,subruns,events]
      tuples = list(zip(*arrays)) #Make multiindex
      index = pd.MultiIndex.from_tuples(tuples, names=['run','subrun','event'])
      shw = shw.drop(['run','subrun','event'],axis=1)
      columns = shw.columns
      shw_data = shw.values
      shw = pd.DataFrame(shw_data,index=index,columns=columns)
    else:
      shw = shw.set_index(['run','subrun','event'])
    indeces = list(shw.index.values) 
    return shw,indeces
  
