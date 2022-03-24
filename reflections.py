#Get new PMTs dataframe
import pandas as pd
import uproot
import sys
import os
sys.path.append('/sbnd/app/users/brindenc/mypython')
from bc_utils.refutils import pic,plotters
cwd = os.getcwd()+'/'


pmts = pd.read_pickle(cwd+'PMT_info.pkl')
pmts = pmts.drop(['opdet_area','f','distance'],axis=1)
new_pmts = pic.get_ref_PMTs(pmts)
new_pmts.to_pickle('PMT_info_ref.pkl')
