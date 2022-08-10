from random import sample
import sys
import os

from matplotlib.style import library
sys.path.append('/sbnd/app/users/brindenc/mypython') #My utils path
from bc_utils.pmtutils import pic as pmtpic
from bc_utils.pmtutils import plotters as pmtplotters
from bc_utils.utils import pic,plotters
import uproot
import awkward as ak
import numpy as np

#Pass file number to script, we'll be doing this one file at a time then combining after we filter the data
argList = sys.argv[1:]
file_id = argList[0] #File number, it could also be all for all data (currently not working)
sample = argList[1]

#You can also specify a cwd
cwd = os.getcwd()

#Import root files
tree = uproot.open(cwd+f'/data/{sample}_{file_id}.root:hitdumper/hitdumpertree;1')
keys = tree.keys() #they have the same keys, don't worry

#Constants
tlow = -5 #ms
thigh = 10 #ms

#print('Opened trees')

#Get op,crt,muon,g4 info
run_info_keys = ['run','subrun','event']
g4keys = ['no_primaries',
 'geant_list_size',
 'pdg',
 'status',
 'Eng',
 'EndE',
 'Mass',
 'Px',
 'Py',
 'Pz',
 'P',
 'StartPointx',
 'StartPointy',
 'StartPointz',
 'StartT',
 'EndT',
 'EndPointx',
 'EndPointy',
 'EndPointz',
 'theta_xz',
 'theta_yz',
 'pathlen',
 'NumberDaughters',
 'TrackId',
 'Mother']
g4keys.extend(run_info_keys)

opkeys = [key for key in keys if 'op' in key]
opkeys.extend(run_info_keys)
opkeys.remove('nophits')

crtkeys = [key for key in keys if 'crt' in key]
crtkeys.extend(run_info_keys)

muonkeys = [key for key in keys if 'muon' in key]
muonkeys.extend(run_info_keys)

#print('Making arrays')

op_df = tree.arrays(opkeys,library='pd')
op_df = op_df.set_index(run_info_keys)

#print('Op dataframe')

crt_df = tree.arrays(crtkeys,library='pd')
crt_df = crt_df.set_index(run_info_keys)

muon_df = tree.arrays(muonkeys,library='pd')
muon_df = muon_df.set_index(run_info_keys)

g4_df = tree.arrays(g4keys,library='pd')
g4 = g4_df.set_index(run_info_keys)

#g4 first cuts for tracks
g4.loc[:,'theta_yx'] = np.arctan2(g4.loc[:,'Py'],g4.loc[:,'Px']) #calculate angle using momentum
temp = g4[g4.loc[:,'pathlen']>0] #Length it was in the detector for
temp = temp[abs(temp.loc[:,'pdg']) == 13] #Muon
temp = temp[temp.loc[:,'status'] == 1] #Primary particle
g4_cut = temp.copy() #First cut

#Treadout cuts 
g4_cut = pmtpic.find_cosmicentrance(g4_cut) #find where cosmics enter detector, important for treadout
g4_cut = pmtpic.get_treadout(g4_cut) #get treadout
temp = g4_cut[g4_cut.loc[:,'treadout'] < thigh*1e6] #convert to ns
temp = temp[temp.loc[:,'treadout'] > tlow*1e6] #convert to ns
g4_df = temp.copy() #Second cut


#print('Made dataframes')

#Add tpc info to op info
op_df.loc[:,'op_tpc'] = op_df.loc[:,'ophit_opch'].values%2 #if the channel is odd, the tpc is 1, this works out nicely

#Employ early cuts

#Save dataframes
op_df.to_pickle(f'data/op_df_{sample}__precut{file_id}.pkl')
crt_df.to_pickle(f'data/crt_df_{sample}__precut{file_id}.pkl')
muon_df.to_pickle(f'data/muon_df_{sample}__precut{file_id}.pkl')
g4_df.to_pickle(f'data/g4_df_{sample}__precut{file_id}.pkl')




