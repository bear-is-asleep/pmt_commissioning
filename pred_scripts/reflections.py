#Get new PMTs dataframe
import pandas as pd
import uproot
import sys
import os
sys.path.append('/sbnd/app/users/brindenc/mypython')
from bc_utils.refutils import pic,plotters
import getopt
cwd = os.getcwd()+'/'

#Handle args passed
argList = sys.argv[1:]

#Options
options = 'n:'
longoptions = ''

if len(argList) > 0:
  args,values = getopt.getopt(argList,options,longoptions)
  for opt,arg in args:
    curArg = opt
    curVal = arg
    #print(curArg,curVal,args,values)
    if 'n' in curArg:
      n_ref = int(curVal)
      
else:
  n_ref = 1

pmts = pd.read_pickle(cwd+'PMT_info.pkl')
pmts = pmts.drop(['opdet_area','f','distance'],axis=1) #These are old, I'm going to drop them
new_pmts = pic.get_ref_PMTs(pmts,n_ref=n_ref)
new_pmts.to_pickle('PMT_info_ref.pkl')
