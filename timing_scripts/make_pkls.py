import sys
import os
sys.path.append('/sbnd/app/users/brindenc/mypython') #My utils path
from bc_utils.pmtutils import pic as pmtpic
from bc_utils.pmtutils import plotters as pmtplotters
from bc_utils.utils import pic,plotters
import uproot

#You can also specify a cwd
cwd = os.getcwd()

#Import root files
tree_tpc0 = uproot.open(cwd+'/data/hitdumper_tb0.root:hitdumper/hitdumpertree;1')
tree_tpc1 = uproot.open(cwd+'/data/hitdumper_tb1.root:hitdumper/hitdumpertree;1')
keys = tree_tpc0.keys() #they have the same keys, don't worry


#Get op,crt,muon info
run_info_keys = ['run','subrun','event']

opkeys = [key for key in keys if 'op' in key]
opkeys.extend(run_info_keys)
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


op_df_tpc0 = tree_tpc0.arrays(opkeys,library='pd')
op_df_tpc0 = op_df_tpc0.set_index(run_info_keys)
op_df_tpc1 = tree_tpc1.arrays(opkeys,library='pd')
op_df_tpc1 = op_df_tpc1.set_index(run_info_keys)

crt_df_tpc0 = tree_tpc0.arrays(crtkeys,library='pd')
crt_df_tpc0 = crt_df_tpc0.set_index(run_info_keys)
crt_df_tpc1 = tree_tpc1.arrays(crtkeys,library='pd')
crt_df_tpc1 = crt_df_tpc1.set_index(run_info_keys)

ctrk_df_tpc0 = tree_tpc0.arrays(ctrkkeys,library='pd')
ctrk_df_tpc0 = ctrk_df_tpc0.set_index(run_info_keys)
ctrk_df_tpc1 = tree_tpc1.arrays(ctrkkeys,library='pd')
ctrk_df_tpc1 = ctrk_df_tpc1.set_index(run_info_keys)

muon_df_tpc0 = tree_tpc0.arrays(muonkeys,library='pd')
muon_df_tpc0 = muon_df_tpc0.set_index(run_info_keys)
muon_df_tpc1 = tree_tpc1.arrays(muonkeys,library='pd')
muon_df_tpc1 = muon_df_tpc1.set_index(run_info_keys)

hit_df_tpc0 = tree_tpc0.arrays(hitkeys,library='pd')
hit_df_tpc0 = hit_df_tpc0.set_index(run_info_keys)
hit_df_tpc1 = tree_tpc1.arrays(hitkeys,library='pd')
hit_df_tpc1 = hit_df_tpc1.set_index(run_info_keys)

#Add tpc info to op info
op_df_tpc0.loc[:,'op_tpc'] = op_df_tpc0.loc[:,'ophit_opch'].values%2 #if the channel is odd, the tpc is 1, this works out nicely
op_df_tpc1.loc[:,'op_tpc'] = op_df_tpc1.loc[:,'ophit_opch'].values%2



#Save dataframes
op_df_tpc0.to_pickle('data/op_df_tpc0.pkl')
op_df_tpc1.to_pickle('data/op_df_tpc1.pkl')
crt_df_tpc0.to_pickle('data/crt_df_tpc0.pkl')
crt_df_tpc1.to_pickle('data/crt_df_tpc1.pkl')
ctrk_df_tpc0.to_pickle('data/ctrk_df_tpc0.pkl')
ctrk_df_tpc1.to_pickle('data/ctrk_df_tpc1.pkl')
muon_df_tpc0.to_pickle('data/muon_df_tpc0.pkl')
muon_df_tpc1.to_pickle('data/muon_df_tpc1.pkl')
hit_df_tpc0.to_pickle('data/hit_df_tpc0.pkl')
hit_df_tpc1.to_pickle('data/hit_df_tpc1.pkl')




