#!/sbnd/data/users/brindenc/.local/bin/python3.9
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('/sbnd/app/users/brindenc/mypython') #My utils path
from bc_utils.pmtutils import pic as pmtpic
from bc_utils.pmtutils import plotters as pmtplotters
#from bc_utils.hitutils import pic as hitpic
#from bc_utils.hitutils import plotters as hitplotters
from bc_utils.utils import pic,plotters
import uproot

#You can also specify a cwd
cwd = os.getcwd()
data_dir = '/pnfs/sbnd/persistent/users/brindenc/analyze_sbnd/nue/v09_43_00/data'

#Constants/parameters
tlow = -1.5 #ms
thigh = 2.5 #ms
cutoff = 10 #Max number of plots
#dpi = 300 #image resolution

#Load files
hd_tree = uproot.open(f'{cwd}/hitdumper_tree.root:hitdumper/hitdumpertree;1')
keys = hd_tree.keys()

#g4 information
g4_keys = pic.g4_keys
g4_keys.extend(['run','subrun','event'])
g4 = hd_tree.arrays(g4_keys,library='pd')
g4 = g4.set_index(['run','subrun','event'])
g4_indeces = g4.index.drop_duplicates()

#muon information
muon_keys = [key for key in keys if 'muon' in key]
muon_keys.extend(['run','subrun','event'])
muon = hd_tree.arrays(muon_keys,library='pd')
muon = muon.set_index(['run','subrun','event'])
muon_indeces = muon.index.drop_duplicates()

#g4 first cuts for tracks
g4.loc[:,'theta_yx'] = np.arctan2(g4.loc[:,'Py'],g4.loc[:,'Px']) #calculate angle using momentum
temp = g4[g4.loc[:,'pathlen']>0] #Length it was in the detector for
temp = temp[abs(temp.loc[:,'pdg']) == 13] #Muon
temp = temp[temp.loc[:,'status'] == 1] #Primary particle
g4_cut = temp.copy() #First cut

#Treadout cuts 
g4_cut = pmtpic.find_cosmicentrance(g4_cut,method=1) #find where cosmics enter detector, important for treadout, method 1 uses start and end point for angles
g4_cut = pmtpic.get_treadout(g4_cut) #get treadout
temp = g4_cut[g4_cut.loc[:,'treadout'] < thigh*1e6] #convert to ns
temp = temp[temp.loc[:,'treadout'] > tlow*1e6] #convert to ns
g4_cut = temp.copy() #Second cut

#Make plots
for i in range(len(muon_indeces)): #iterate over all tracks
  if i > cutoff: #introduce cutoff, don't overburden me with plots
    break
  else:
    for x in ['x','z']:
      for y in ['x','y','z']: #iterate over all coordinate pairs
        if x != y: #don't plot single coordinate, that's boring
          pmtplotters.plot_g4_muon(g4_cut,muon,muon_indeces[i],thigh,tlow,x=x,y=y,show_vtx=False) #plotting package




