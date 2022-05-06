import sys
import os

from matplotlib.style import library
sys.path.append('/sbnd/app/users/brindenc/mypython') #My utils path
from bc_utils.pmtutils import pic as pmtpic
from bc_utils.pmtutils import plotters as pmtplotters
from bc_utils.utils import pic,plotters
import uproot
import awkward as ak

#Pass file number to script, we'll be doing this one file at a time then combining after we filter the data
argList = sys.argv[1:]
file_id = argList[0] #File number, it could also be all for all data (currently not working)

#You can also specify a cwd
cwd = os.getcwd()

#Import root files
tree_ew = uproot.open(cwd+f'/data/ew_{file_id}.root:hitdumper/hitdumpertree;1')
tree_fb = uproot.open(cwd+f'/data/fb_{file_id}.root:hitdumper/hitdumpertree;1')
keys = tree_ew.keys() #they have the same keys, don't worry

#print('Opened trees')

#Get op,crt,muon info
run_info_keys = ['run','subrun','event']

opkeys = [key for key in keys if 'op' in key]
opkeys.extend(run_info_keys)
opkeys.remove('nophits')
crtkeys = [key for key in keys if 'crt' in key]
crtkeys.extend(run_info_keys)
ctrkkeys = [key for key in keys if 'ctrk' in key]
ctrkkeys.extend(run_info_keys)
muonkeys = [key for key in keys if 'muon' in key]
muonkeys.extend(run_info_keys)
#Manually specify
hitkeys = ['nhits','hit_cryostat','hit_tpc','hit_plane','hit_wire','hit_channel',
          'hit_peakT','hit_charge','hit_ph','hit_width','hit_full_integral']
hitkeys.extend(run_info_keys)

#print('Making arrays')

op_df_ew = tree_ew.arrays(opkeys,library='pd')
op_df_ew = op_df_ew.set_index(run_info_keys)
op_df_fb = tree_fb.arrays(opkeys,library='pd')
op_df_fb = op_df_fb.set_index(run_info_keys)

#print('Op dataframe')

# crt_df_ew = tree_ew.arrays(crtkeys,library='pd')
# crt_df_ew = crt_df_ew.set_index(run_info_keys)
# crt_df_fb = tree_fb.arrays(crtkeys,library='pd')
# crt_df_fb = crt_df_fb.set_index(run_info_keys)

# ctrk_df_ew = tree_ew.arrays(ctrkkeys,library='pd')
# ctrk_df_ew = ctrk_df_ew.set_index(run_info_keys)
# ctrk_df_fb = tree_fb.arrays(ctrkkeys,library='pd')
# ctrk_df_fb = ctrk_df_fb.set_index(run_info_keys)

muon_df_ew = tree_ew.arrays(muonkeys,library='pd')
muon_df_ew = muon_df_ew.set_index(run_info_keys)
muon_df_fb = tree_fb.arrays(muonkeys,library='pd')
muon_df_fb = muon_df_fb.set_index(run_info_keys)

# hit_df_ew = tree_ew.arrays(hitkeys,library='pd')
# hit_df_ew = hit_df_ew.set_index(run_info_keys)
# hit_df_fb = tree_fb.arrays(hitkeys,library='pd')
# hit_df_fb = hit_df_fb.set_index(run_info_keys)

#print('Made dataframes')

#Add tpc info to op info
op_df_ew.loc[:,'op_tpc'] = op_df_ew.loc[:,'ophit_opch'].values%2 #if the channel is odd, the tpc is 1, this works out nicely
op_df_fb.loc[:,'op_tpc'] = op_df_fb.loc[:,'ophit_opch'].values%2

#Employ early cuts

#Save dataframes
op_df_ew.to_pickle(f'data/op_df_ew__precut{file_id}.pkl')
op_df_fb.to_pickle(f'data/op_df_fb__precut{file_id}.pkl')
# crt_df_ew.to_pickle('data/crt_df_ew.pkl')
# crt_df_fb.to_pickle('data/crt_df_fb.pkl')
# ctrk_df_ew.to_pickle('data/ctrk_df_ew.pkl')
# ctrk_df_fb.to_pickle('data/ctrk_df_fb.pkl')
muon_df_ew.to_pickle(f'data/muon_df_ew__precut{file_id}.pkl')
muon_df_fb.to_pickle(f'data/muon_df_fb__precut{file_id}.pkl')
# hit_df_ew.to_pickle('data/hit_df_ew.pkl')
# hit_df_fb.to_pickle('data/hit_df_fb.pkl')




