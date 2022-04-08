import os
import getopt
import sys

sbnd_version = 'v09_43_00'
scripts_dir = f'/sbnd/app/users/brindenc/analyze_sbnd/PMT/scripts/'

#Handle args passed
argList = sys.argv[1:]

#Options
options = 'f:'
longoptions = ''

if len(argList) > 0:
  args,values = getopt.getopt(argList,options,longoptions)
  for opt,arg in args:
    curArg = opt
    curVal = arg
    #print(curArg,curVal,args,values)
    if 'f' in curArg:
      file = curVal
      file_path = f'{os.getcwd()}/{file}'
      cwd = f'/sbnd/app/users/brindenc/analyze_sbnd/PMT/{sbnd_version}/analysis_pipeline'
      print(f'Uploaded file {file_path}')
      
else:
  file = 'hitdumper_ew_filtered_one.root'
  file_path = f'/sbnd/app/users/brindenc/analyze_sbnd/PMT/v09_43_00/analysis_pipeline/{file}'
  cwd = f'/sbnd/app/users/brindenc/analyze_sbnd/PMT/{sbnd_version}/analysis_pipeline'
  os.system(f'mkdir -p {cwd}')
  print(f'No args passed, defaulting to {sbnd_version} and file: {file_path}')

#Analyze in the right directory
os.system(f'mkdir -p {cwd}')
os.chdir(cwd)
os.system(f'cp -rf {scripts_dir}* .') #Copy all the scripts and directories into current directory


os.system(f'python3.9 PMT_analyze_fit.py -a {file_path}') #Run python
os.chdir('read_model')
os.system('python3.9 PMT_results.py')
os.system('python3.9 predict_working_PMTs.py')
print('*********************\n\n*********************')
print(f'Finished with analysis! Look for your stinking files in:\n{os.getcwd()}/')


