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

#Return hits and calo info dataframe
def get_hits_df(rootname,treename,branchname='Event;1',offset=0):
    #Input root and tree name to return dataframe with g4 information
    #Offset run creates indexing for additional backgrounds
    hits_tree = uproot.open(f'{rootname}:{treename}/{branchname}')
    with uproot.open(f'{rootname}:{treename}/{branchname}') as file:
        keys = file.keys()
    hits_keys = [key for key in keys if 'hit_' in key]
    hits_keys.extend(['run','subrun','event'])
    hits = hits_tree.arrays(hits_keys,library='pd')
    runs = list(hits.loc[:,'run'])
    runs = [run+offset for run in runs] #Add offset to run value for additional background
    subruns = list(hits.loc[:,'subrun'])
    events = list(hits.loc[:,'event'])
    if offset != 0:
      arrays = [runs,subruns,events]
      tuples = list(zip(*arrays)) #Make multiindex
      index = pd.MultiIndex.from_tuples(tuples, names=['run','subrun','event'])
      hits = hits.drop(['run','subrun','event'],axis=1)
      columns = hits.columns
      hits_data = hits.values
      hits = pd.DataFrame(hits_data,index=index,columns=columns)
    else:
      hits = hits.set_index(['run','subrun','event'])
    indeces = list(hits.index.values) 
    return hits,indeces


