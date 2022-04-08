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
import getopt
import random
pd.options.mode.chained_assignment = None  # default='warn'

sbnd_version = 'v09_43_00'
#cwd = '/sbnd/app/users/brindenc/analyze_sbnd/PMT/'+sbnd_version+'/analysis_pipeline/' #Put your working directory here
cwd = os.getcwd()+'/'
sys.path.append(cwd)
print(f'Working in {cwd}')

#Constants/Default values
make_ophits = True
train = 0.8

QE = 0.25 # cathode at 390 nm https://www.hamamatsu.com/us/en/product/type/R5912/index.html
Tw = 0.92 #Transmissivity of PMT https://www.azom.com/article.aspx?ArticleID=4765
Rw = 0 #Reflectivity of PMT
#Columns for pmt info
columns = ['ophit_ch','ophit_pred','ophit_obs','pred_err','run','subrun',
          'event','ophit_opdet_type','distance','ophit_opdet_x','ophit_opdet_y',
          'ophit_opdet_z','opdet_tpc','muontrk_tpc','muontrk_x1', 'muontrk_y1', 'muontrk_z1', 'muontrk_x2',
          'muontrk_y2', 'muontrk_z2','mhits','mtype','mean_distance','muontrk_dx','muontrk_length'] 

#Handle args passed
argList = sys.argv[1:]

#Options
options = 'fa:'
longoptions = ''

if len(argList) > 0:
  args,values = getopt.getopt(argList,options)
  for opt,arg in args:
    curArg = opt
    curVal = arg
    #print(args,curVal,values)
    if 'a' in curArg:
      print(f'Analyzing {curVal}')
      file = curVal
      do_fit = False
    elif 'f' in curArg:
      print(f'Fitting {curVal} on {int(train*100)}% of dataset')
      file = curVal
      do_fit = True
else:
  file = 'hitdumper_ew_filtered_one.root'
  do_fit = False
  print(f'No args passed, defaulting to analyzing {file}')

#Functions

#Load trees
hitdumper_tree = uproot.open(file+':hitdumper/hitdumpertree;1')
#print('Got hitdumper tree')

#PMT information
pmts = pd.read_pickle(cwd+'PMT_info.pkl')
#print('Got pmt info')

#Muon information
muon_keys = ['run','subrun','event','muontrk_t0', 'muontrk_x1', 'muontrk_y1', 'muontrk_z1', 'muontrk_x2',
       'muontrk_y2', 'muontrk_z2', 'muontrk_theta_xz', 'muontrk_theta_yz',
       'muontrk_tpc', 'muontrk_type','nhits','nmuontrks']
muontracks = hitdumper_tree.arrays(muon_keys,library='pd')
#print('Got muon tracks')
muontracks = muontracks.drop_duplicates()
muontracks = muontracks.set_index(['run','subrun','event'])
muontracks = muontracks.sort_index()

#Op hits
if make_ophits:
  #Make ophits
  #print('Making OP hits ')
  op_hit_keys = ['run','subrun','event','ophit_opdet','ophit_opdet_type']
  op_hits = hitdumper_tree.arrays(op_hit_keys)
  op_hits = ak.to_pandas(op_hits) #If it ain't broke don't fix it
  #print('Got op hits')
  op_hits = op_hits.set_index(['run','subrun','event'])
  op_hits = op_hits.sort_index()
  op_hits.to_pickle(cwd+'op_hits.pkl')
else:
  #op_hits = pd.read_pickle(cwd + 'op_hits.pkl')
  print('Loaded op hits')

#Get test train split indeces
indeces = list(op_hits.index.drop_duplicates())

def remove_list_duplicates(lis):
  #Remove duplicates from list
  res = []
  [res.append(x) for x in lis if x not in res]
  return res


if do_fit:
  train_indeces = random.sample(indeces,int(train*len(indeces)))
  train_indeces.sort()
  test_indeces = [ind for ind in indeces if ind not in train_indeces]
  test_indeces.sort()
  pic.print_stars()
  print(f'Training on {len(train_indeces)} events, testing on {len(test_indeces)} events')

  #Train events, subruns, runs
  train_events = remove_list_duplicates([ind[2] for ind in train_indeces])
  train_events.sort()
  train_subruns = remove_list_duplicates([ind[1] for ind in train_indeces])
  train_subruns.sort()
  train_runs = remove_list_duplicates([ind[0] for ind in train_indeces])
  train_runs.sort()

  #Test events, subruns, runs
  events = remove_list_duplicates([ind[2] for ind in test_indeces])
  events.sort()
  subruns = remove_list_duplicates([ind[1] for ind in test_indeces])
  subruns.sort()
  runs = remove_list_duplicates([ind[0] for ind in test_indeces])
  runs.sort()

  print('Sorting DFs')
  #Get training info
  op_hits_train = op_hits.loc[train_indeces].sort_index()
  muontracks_train = muontracks.loc[train_indeces].sort_index()

  #Make fit dataframe
  pic.make_analysis_df(pmts,muontracks_train,op_hits_train,columns,
  train_events,train_subruns,train_runs,cwd,fit_results=True,move_dfs=True) #fit
  
  #Get test info
  op_hits = op_hits.loc[test_indeces].sort_index()
  muontracks = muontracks.loc[test_indeces].sort_index()

  #Analyze on test results
  pic.make_analysis_df(pmts,muontracks,op_hits,columns,
  events,subruns,runs,cwd,fit_results=False,move_dfs=True) #analyze
  

else:
  #Events, subruns, runs
  events = remove_list_duplicates([ind[2] for ind in indeces])
  events.sort()
  subruns = remove_list_duplicates([ind[1] for ind in indeces])
  subruns.sort()
  runs = remove_list_duplicates([ind[0] for ind in indeces])
  runs.sort()

  #Analyze run
  pic.make_analytical_df(pmts,muontracks,op_hits,events,subruns,runs,cwd,move_dfs=True) 









