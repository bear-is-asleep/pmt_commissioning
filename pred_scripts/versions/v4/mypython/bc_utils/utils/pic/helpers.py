import numpy as np
from random import choice
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from datetime import date
import sys
import seaborn as sns
import matplotlib
from scipy import optimize
from sklearn.linear_model import LinearRegression
from numpy import isin, sqrt,exp,arctan,cos,sin,pi
from time import time

#Constants
hc = 1240 #eV nm
r_e = 2.818e-15 #m
alpha = 1/137 #alpha
m_e = 511e3 #eV/c^2
n_A = 1.3954/6.6335209e-23 #n_A/cm^3
m_u = 105.658 #MeV
efficiency = 0.25 #Efficiency of light collection
efficiency_TPB = 0.5 #Efficiency of TPB conversion
l_att = 66 #cm https://arxiv.org/pdf/1611.02481.pdf 

#Constants

#Centers of tpc
tpc0_xyz = np.array([-100,0,250])
tpc1_xyz = np.array([100,0,250])
r = r_pmt = 19/2 #cm https://www.hamamatsu.com/us/en/product/type/R5912/index.html
V = 4*5*2 #TPC volume m^3
V_cm = V*1e6
A = 4*5*2 + 2*4*2 + 2*5*2 #TPC area m^2
A_cm = A*1e4
L = 6*V_cm/A_cm #Characteristic length of detector cm

l_R = 55#Rayleigh scattering length cm +-5 https://pure.royalholloway.ac.uk/portal/files/29369028/2018GraceEPHD.pdf
l = 63 #Attenuation length cm +- 3 https://www.mpi-hd.mpg.de/gerda/public/2020/phd2020_BirgitZatschler.pdf_compressed.pdf
QE = 0.25 # cathode at 390 nm https://www.hamamatsu.com/us/en/product/type/R5912/index.html
Tw = 0.92 #Transmissivity of PMT https://www.azom.com/article.aspx?ArticleID=4765
Rw = 0 #Reflectivity of PMT

#Rs = 0.886 #Reflectivity of surface coating

E_ws = 0.5 #Assume same as PMT coating
e_ws = 0.5 #Shifting efficiency TPB
Ref_cov = 4*5/A #Area of inner surface covered by TPB coated reflective surfaces
dEdx = 2.2 #MeV/cm check later

#Gaussian fit
def gaussian(x,a, mean, stddev):
  #Ex: popt, pcov = optimize.curve_fit(gaussian, x, data)
  return a*np.exp(-((x - mean) / sqrt(2) / stddev)**2)

def fit_gaussian(data):
  #Fits gaussian to 1d data and returns parameters
  steps=len(data)
  x = np.linspace(min(data),max(data),steps)
  popt, pcov = optimize.curve_fit(gaussian, x, data)
  return x,popt,pcov
def fit_gauss(x,y):
  return optimize.curve_fit(gaussian,x,y,sigma=sqrt(y))

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


def truncate_df(df,keys_to_sum = {''},keys_to_median = {''},mod_chs=False):
  #If mod_chs, we take the ch%1000, so that reflected channels convert to the same channel 
  #Truncate df for average values across all events and runs, len(df) = number of PMTs
  cols = df.shape[1]
  col_options = np.ones(cols)
  if mod_chs:
    df.loc[:,'ophit_ch'] = df.loc[:,'ophit_ch'].values%1000 #Convert all channels into one

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

#Equations
def N_0 (length,dEdx):
  #Number of ionized photons given in SBND proposal paper
  ionization_rate = 2.4e4 #gamma/MeV
  return length*dEdx*ionization_rate 
def d_calc(x):
  #Dispersion of photons do to radiating spherically
  if np.isscalar(x): #Calculating dispersion for scalar
    if x < sqrt(r_pmt**2/2):
      #zz = 0
      return 1/2
    elif x >= sqrt(r_pmt**2/2):
      return r_pmt**2/(4*x**2)
      #zz = 0
    #return 1 #We're going to try this for now
  else: #Calculating dispersion for vector x input, ds are dispersion values
    ds = []
    for val in x:
      if val < sqrt(r_pmt**2/2):
        ds.append(1/2)
      elif val >= sqrt(r_pmt**2/2):
        ds.append(r_pmt**2/(4*val**2))
    return np.asarray(ds)
    #return np.ones(x.shape) #Try this for now

def d_calc_lazy(x,f=pi*r_pmt**2/A_cm):
  #Dispersion of photons, radiating in a box
  #Account for infinite reflections by multiplying by constant Rs
  #f is fractional coverage of PMT
  #R is reflectivity, this will be the area of the reflective surface divided by the total area
  #coating is type of PMT coating

  if np.isscalar(x): #Calculating dispersion for scalar
    return f #We're going to try this for now
  else: #Calculating dispersion for vector x input, ds are dispersion values
    return np.full(x.shape,f) #Try this for now

def I_calc(x,l):
  #Photons lost due to scattering
  if not np.isscalar(x):
    vals = []
    for val in x:
      vals.append(exp(-val/l))
    return np.asarray(vals)
  else:
    return exp(-x/l)
def N_calc(x,efficiency,length,dEdx,l,photon_interactions=True,dispersion=True):
  #Photons that reach PMT as a function of x (distance from PMT)
  if photon_interactions and dispersion:
    return N_0(length,dEdx)*d_calc(x)*I_calc(x,l)*efficiency
  elif photon_interactions and not dispersion:
    return N_0(length,dEdx)*I_calc(x,l)*efficiency
  elif not photon_interactions and dispersion:
    return N_0(length,dEdx)*d_calc(x)*efficiency
  else: 
    return N_0(length,dEdx)*efficiency
def N_calc_lazy(x,efficiency,length,dEdx,l):
  #Photons that reach PMT as a function of x (distance from PMT)
  #Use lazy d method
  return N_0(length,dEdx)*d_calc_lazy(x)*I_calc(x,l)*efficiency

def total_hits_PMT(df):
  #Returns dataframe which combines arrays to find total ophits
  arr = df.values #Get values into numpy array
  col_options = []
  for key in df.keys():
    if key == 'nophits':
      col_options.append(0)
    else:
      col_options.append(2)

  new_arr = sum_combine(arr,every_other=arr.shape[0],col_options=col_options)

  return pd.DataFrame(new_arr,columns=df.keys())

def distance(xyz1,xyz2):
  ds = np.zeros(3)
  for i in range(3):
    ds[i] = (xyz1[i]-xyz2[i])**2
  return sqrt(sum(ds))

def angles(xyz1,xyz2):
  ds = np.zeros(3)
  for i in range(3):
    ds[i] = abs(xyz2[i]-xyz1[i])
  theta_xy = arctan(ds[1]/ds[0])
  theta_xz = arctan(ds[2]/ds[0])
  return theta_xy,theta_xz  

def A_calc(a,b):
  return pi*a*b
def a_calc(d,r,theta_xz):
  return d*r*cos(theta_xz)/(d+r*sin(theta_xz))
def b_calc(d,r,theta_xy):
  return d*r*cos(theta_xy)/(d+r*sin(theta_xy))

def get_pmts(opdet):
  #Take opdet in and return pmt information
  opdet = opdet.drop(['nophits','ophit_peakT','ophit_width','ophit_area','ophit_amplitude','ophit_pe'],axis=1)
  opdet = opdet.drop_duplicates()
  pmts_coated = opdet[opdet['ophit_opdet_type'] == 0] #0:coated 1:uncoated
  pmts_uncoated = opdet[opdet['ophit_opdet_type'] == 1] #0:coated 1:uncoated
  pmts = pd.concat((pmts_coated,pmts_uncoated))

  temp_df = pd.DataFrame(index=pmts.index,columns=['opdet_tpc','opdet_area','f','distance'])

  for row,line in pmts.iterrows():
    x = line['ophit_opdet_x']
    y = line['ophit_opdet_y']
    z = line['ophit_opdet_z']
    xyz = np.array([x,y,z])
    if x < 0:
      tpc = 0
      xyz2 = tpc0_xyz
    elif x > 0:
      tpc = 1
      xyz2 = tpc1_xyz
    temp_df.loc[row]['opdet_tpc'] = tpc
    d = distance(xyz,xyz2) #Distance from center to PMTs
    phi_xy,phi_xz = angles(xyz,xyz2) #Angle from PMT
    theta_xy = abs(phi_xy)
    theta_xz = abs(phi_xz) 
    a = a_calc(d,r,theta_xz) #Compression in z-axis, right left
    b = b_calc(d,r,theta_xy) #Compression in y-axis, up down
    A = abs(A_calc(a,b)) #Area of resulting elipse
    temp_df.loc[row]['opdet_area'] = A
    temp_df.loc[row]['f'] = A/A_cm
    temp_df.loc[row]['distance'] = d
  print(pmts.shape)
  pmts = pd.concat((pmts,temp_df),axis=1)
  #pmts = pmts.loc[:,~pmts.columns.duplicated()] #Drop duplicate columns
  return pmts

def err_calc(true,pred):
  if isinstance(true,float) or isinstance(true,int) and isinstance(pred,float):
    if true < 1e-10:
      return pred #Handle zeros
    else:
      return (pred-true)/true
  else:
    err = np.zeros(true.shape)
    for i in range(len(true)): #Should be same size
      if true[i] < 1e-10:
        err[i] = pred[i] #Handle zeros
      else:
        err[i] = (pred[i]-true[i])/true[i]
    return err


def distances_to_PMT(hits,x1,y1,z1,x2,y2,z2,xp,yp,zp):
	#Distance of all track points to PMT given track
  hits = int(hits) #Convert to hits
  d = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2) #Total distance of track
  lengths = np.array([x2-x1,y2-y1,z2-z1]) #length of line segment each dimension
  step_sizes = lengths/hits
  d_to_PMT = [] #distance to PMT for each hit
  #Hit locations along line
  x_h = []
  y_h = []
  z_h = []
  for i in range(hits+1):
    x_h.append(x1+i*step_sizes[0])
    y_h.append(y1+i*step_sizes[1])
    z_h.append(z1+i*step_sizes[2])
    d_to_PMT.append(sqrt((x_h[i]-xp)**2+(y_h[i]-yp)**2+(z_h[i]-zp)**2)) #3D distance to PMT
  
  d_to_PMT = np.asarray(d_to_PMT)
  return d_to_PMT,x_h,y_h,z_h,np.linalg.norm(step_sizes) #Option to return hit loctations 

def print_stars():
  print('\n******************************\n')

def get_muon_tracks(pmt_df):
  #Return dataframe for plotting muon tracks in both TPCs
  columns=['muontrk_x1_0', 'muontrk_y1_0', 'muontrk_z1_0', 'muontrk_x2_0',
       'muontrk_y2_0', 'muontrk_z2_0','muontrk_x1_1', 'muontrk_y1_1', 'muontrk_z1_1', 'muontrk_x2_1',
       'muontrk_y2_1', 'muontrk_z2_1','event','subrun','run']
  mtrks_data = np.full((pmt_df.shape[0],len(columns)),-999)
  cnter = 0 #Use this to index
  for row,line in pmt_df.iterrows():
    i = cnter #event
    mtrks_data[i,12] = row[2] #Mark event
    mtrks_data[i,13] = row[1] #Mark subrun
    mtrks_data[i,14] = row[0] #Mark run

    if line['muontrk_tpc'] == 0: #Append for x1_0 - z2_0
        mtrks_data[i,0] = line['muontrk_x1']
        mtrks_data[i,1] = line['muontrk_y1']
        mtrks_data[i,2] = line['muontrk_z1']
        mtrks_data[i,3] = line['muontrk_x2']
        mtrks_data[i,4] = line['muontrk_y2']
        mtrks_data[i,5] = line['muontrk_z2']
    if line['muontrk_tpc'] == 1: #Append for x1_1 - z2_1
        mtrks_data[i,6] = line['muontrk_x1']
        mtrks_data[i,7] = line['muontrk_y1']
        mtrks_data[i,8] = line['muontrk_z1']
        mtrks_data[i,9] = line['muontrk_x2']
        mtrks_data[i,10] = line['muontrk_y2']
        mtrks_data[i,11] = line['muontrk_z2']
    cnter+=1


  mtrks = pd.DataFrame(mtrks_data,columns=columns).drop_duplicates()
  #Don't set indeces for this dataframe, it'll break the plotter code. Maybe I'll fix it tho
  return mtrks

def reflect_xz(y,flip=2,ref=1,L=400):
  #Returns new y, reflected across face in xz plane in detector
  #flip divisible by 2 is off bottom, otherwise off top
  #ref is number of reflections
  #L Length of y-coordinate
  y_prime=0

  if flip % 2 == 0:
    y = y+200 #Adjust such that xz plane starts at y=0
    y_prime = -y
    for n in range(2,ref+1):
      y_prime = -(y_prime+(n-1)*L)-(n-1)*L
  else: 
    y = y+200 #Adjust such that xz plane starts at y=0
    y_prime=y
    for n in range(1,ref+1):
      y_prime = -(y_prime-n*L)+n*L
  return y_prime-200

def reflect_yz(x,flip=2,ref=1):
  #Returns new y, reflected across face in xz plane in detector
  #flip divisible by 2 is off cpa, otherwise off apa (where PDS is)
  #ref is number of reflections
  L = 200 #Length of x-coordinate
  x_prime=0

  if flip % 2 == 0: #cpa
    cpa_ref = int((ref+1)/2)
    x_prime = x-L*ref #Since x is always the same
  else: #apa
    cpa_ref = int(ref/2)
    x_prime = x+L*ref #Since x is always the same
  return x_prime,cpa_ref,ref-cpa_ref #Also return number of cpa reflections, it has reflective coating

def reflect_xy(z,flip=2,ref=1,L=500):
  #Returns new z, reflected across face in xy plane in detector
  #flip divisible by 2 is off z=0, otherwise off z=500
  #ref is number of reflections
  #L Length of z-coordinate
  z_prime=0

  if flip % 2 == 0: #z=0
    z_prime = -z
    for n in range(2,ref+1):
      z_prime = -(z_prime+(n-1)*L)-(n-1)*L
  else: 
    z_prime=z
    for n in range(1,ref+1):
      z_prime = -(z_prime-n*L)+n*L
  return z_prime

def single_plane_reflection(x,y,z,flip_coord='xyz',x_refs=1,y_refs=1,z_refs=1,initialize=False,ch=0):
  #For single reflections only, this can happen 6 different ways:
  #Front-back,left-right,top-bottom
  #Default to returning all 6 types of reflections for set of coordinates, return cpa/apa ref. for each
  #Also default to 1 reflection along each dimension
  if initialize:
    if ch == 0:
      coord_ref = [[x,y,z,0,0,0]] #6 indexes for 3 coords. total # reflections, cpa ref., apa ref.
    else:
      coord_ref = [[x,y,z,0,0,0,ch]] #7 indexes for 3 coords. total # reflections, cpa ref., apa ref., ch #
  else: 
    coord_ref = []
  if 'xyz' in flip_coord: #Build other cases later if needed
    for x_ref in range(1,x_refs+1):
      for y_ref in range(1,y_refs+1):
        for z_ref in range(1,z_refs+1):
          x1,cpa_ref1,apa_ref1 = reflect_yz(x,ref=x_ref)
          x2,cpa_ref2,apa_ref2 = reflect_yz(x,flip=3,ref=x_ref)
          z1 = reflect_xy(z,ref=z_ref)
          z2 = reflect_xy(z,flip=3,ref=z_ref)
          y1 = reflect_xz(y,ref=y_ref)
          y2 = reflect_xz(y,flip=3,ref=y_ref)
          #Now there's 8+10+6=26 cases x,y,z,xy,xz,yz,xyz
          xs = np.array([x,x1,x2])
          ys = np.array([y,y1,y2])
          zs = np.array([z,z1,z2])
          for i,xal in enumerate(xs): 
            for j,yal in enumerate(ys):
              for k,zal in enumerate(zs): #Find all possible combinations of 3 xs, 3 ys, 3 zs
                if i == 0 and j == 0 and k == 0:
                  continue #skip events that don't reflect
                tot_ref = 0 #Get total number of reflections
                cpa_ref = 0 #cpa ref
                apa_ref = 0 #apa ref
                if i != 0:
                  tot_ref += x_ref
                  cpa_ref = 0
                  apa_ref = 0
                if j != 0:
                  tot_ref += y_ref
                if k != 0:
                  tot_ref += z_ref
                if i == 1:
                  cpa_ref = cpa_ref1
                  apa_ref = apa_ref1
                if i == 2:
                  cpa_ref = cpa_ref2
                  apa_ref = apa_ref2
                #6 indexes for 3 coords. total # reflections, cpa ref., apa ref.
                if ch == 0:
                  coord_ref.append([xal,yal,zal,tot_ref,cpa_ref,apa_ref])
                else: 
                  ch+=1000
                  coord_ref.append([xal,yal,zal,tot_ref,cpa_ref,apa_ref])
  if ch == 0:
    df = pd.DataFrame(coord_ref,columns=['x','y','z','tot_ref','cpa_ref','apa_ref']).drop_duplicates()
  else: 
    df = pd.DataFrame(coord_ref,columns=['ophit_opdet_x','ophit_opdet_y','ophit_opdet_z',
    'tot_ref','cpa_ref','apa_ref','ophit_opdet']).drop_duplicates()
  return df

def lin_fit_3d(pmt_hits_df,xkey='mhits',ykey='mean_distance',zkey='ophit_obs'):
  #Fit x,y to linear fit in 3 dimenstions
  pmt_hits_df_coated0 = pmt_hits_df[pmt_hits_df.loc[:,'ophit_opdet_type'] == 0]
  pmt_hits_df_coated1 = pmt_hits_df[pmt_hits_df.loc[:,'ophit_opdet_type'] == 1]
  dfs = [pmt_hits_df_coated0,pmt_hits_df_coated1]
  coefs = [] #xkey is first, then ykey
  bs = [] #z intercept of 3d plot 
  scores = [] #Model score on itself
  models = []
  for df in dfs:
    xx = df.loc[:,[xkey,ykey]].values
    y = df.loc[:,zkey].values
    model = LinearRegression()
    model.fit(xx,y)
    coefs.append(model.coef_)
    bs.append(model.intercept_)
    scores.append(model.score(xx,y))
    models.append(model)
  return coefs,bs,scores,models #Returns all of these u know..

def check_functionality(df,std=4):
  #Accept values up to std outside of 0 (no error)
  obs = df.loc[:,'ophit_obs']
  pred = df.loc[:,'ophit_pred']
  err = df.loc[:,'pred_err']
  calc_err = err_calc(obs,pred)
  std = np.std(err)
  temp = np.zeros(df.shape[0])
  for i,val in enumerate(calc_err):
    if abs(val) > 4*std:
      temp[i] = 0
    else: 
      temp[i] = 1

  df.loc[:,'functioning'] = temp
  return df

def fake_bad_pmt(df,chs,reduce=0.0):
  #Multiply given channels by reduction factor
  temp_df = df.copy() #Make a copy, not to lose info of original df
  bad_indeces = []
  for ch in chs:
    bad_indeces.append(temp_df[temp_df['ophit_ch'] == ch].index.values[0])
  temp_df.loc[bad_indeces,'ophit_obs'] *= reduce
  pred = temp_df.loc[:,'ophit_pred'].values
  obs = temp_df.loc[:,'ophit_obs'].values
  new_err = err_calc(obs,pred)
  temp_df.loc[:,'pred_err'] = new_err
  return temp_df

def make_analysis_df(pmts,muontracks,op_hits,columns,events,subruns,runs,cwd,
fit_results=False,move_dfs=True):
  #pmts: Locations and coating of pmts
  #muontracks: Locations and #hits for muon tracks
  #op_hits: # of hits and waveforms of pmts
  #columns: Columns of output dataframe
  #.
  #.
  #.
  #fit_results: Fit results to op_obs information
  #move_dfs: Move results into read_model folder (recommended for organization)

  #Fit parameters
  if fit_results:
    #print('Set dumby fit values to initialize')
    md0 = mh0 = b0 = md1 = mh1 = b1 = 0
  else:
    print('Loading fit parameters')
    fit_df_loaded = pd.read_pickle(cwd+'fit_params.pkl') #This file should exist after one fit
    fit_df0 = fit_df_loaded[fit_df_loaded.loc[:,'coating'] == 0]
    fit_df1 = fit_df_loaded[fit_df_loaded.loc[:,'coating'] == 1]

    #Coating = 0
    md0 = fit_df0.loc[:,'coef_mean_distance'].values[0] #coefficient for linear fit to distance from track
    mh0 = fit_df0.loc[:,'coef_mhits'].values[0] #coefficient for linear fit to track hits
    b0  = fit_df0.loc[:,'intercept'].values[0] #z intercept for fit

    #Coating = 1
    md1 = fit_df1.loc[:,'coef_mean_distance'].values[0] #coefficient for linear fit to distance from track
    mh1 = fit_df1.loc[:,'coef_mhits'].values[0] #coefficient for linear fit to track hits
    b1  = fit_df1.loc[:,'intercept'].values[0] #z intercept for fit

  pmt_info = np.zeros((pmts.shape[0],muontracks.shape[0],len(columns)))
  j=0 #Counter for 2nd dim. (muon tracks)

  op_hits_index = op_hits.index.drop_duplicates()
  op_hits_index = op_hits_index.values

  print('Starting analysis')
  print_stars()

  total_trks = 0 #Keep track of total number of cosmics processed
  for run in runs: 
    for subrun in subruns:
      for event in events:
        
        start = time()
        index = (run,subrun,event) #Index PMT_hits

        #print(index in op_hits_index,index in muontracks.index.values,index,op_hits_index[-1],type(index))

        #Check if event is in both muon and op information
        in_op = False
        in_muon = False
        for tup in op_hits_index:
          if tup == index:
            in_op = True
        for tup in muontracks.index.values:
          if tup == index:
            in_muon = True
        if in_op and in_muon:# and subrun < 101 and run < 2: #temporary subspace to test
          #print(f'Processing Run: {run} Subrun: {subrun} Event: {event}')
          op_event_info = op_hits.loc[pd.IndexSlice[(index)]]
          muon_temp = muontracks.loc[pd.IndexSlice[(index)]]
        else:
          #print(f'Skipped Run: {run} Subrun: {subrun} Event: {event}')
          continue

        #Divide hits to tracks based on length of tracks
        for i_muon,row_muon in muon_temp.iterrows():
          #Coordinates for muon tracks
          x1 = row_muon['muontrk_x1']
          y1 = row_muon['muontrk_y1']
          z1 = row_muon['muontrk_z1']
          x2 = row_muon['muontrk_x2']
          y2 = row_muon['muontrk_y2']
          z2 = row_muon['muontrk_z2']
          length = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
          muon_temp.loc[i_muon,'muontrk_length'] = length
        #Get fraction of hits belonging to each track
        muon_temp.loc[:,'frac_hits'] = muon_temp.loc[:,'muontrk_length'].values/muon_temp.loc[:,'muontrk_length'].values.sum()
        for i_muon,row_muon in muon_temp.iterrows():
          total_trks+=1 #iterate counter
          i=0 #Counter for second dim. 
          for i_op,row_op in pmts.iterrows():
            #Skip muon type of other (mtype = 5)
            mtype = row_muon['muontrk_type']

            #Coordinates for muon tracks
            x1 = row_muon['muontrk_x1']
            y1 = row_muon['muontrk_y1']
            z1 = row_muon['muontrk_z1']
            x2 = row_muon['muontrk_x2']
            y2 = row_muon['muontrk_y2']
            z2 = row_muon['muontrk_z2']

            mtrks = row_muon['nmuontrks']
            mhits = row_muon['nhits']*row_muon['frac_hits'] #Muon track hits, normalized to track length
            tpc = row_muon['muontrk_tpc'] #Muon tpc
            
            #Coordinates for PMT
            ch = row_op['ophit_opdet']
            xp = row_op['ophit_opdet_x']
            yp = row_op['ophit_opdet_y']
            zp = row_op['ophit_opdet_z']
            ch_hits = float(len(op_event_info[op_event_info['ophit_opdet'] == ch]))*row_muon['frac_hits'] #Account for double counting
            
            det_type = row_op['ophit_opdet_type']
            d_s,x_h,y_h,z_h,dx = distances_to_PMT(mhits,x1,y1,z1,x2,y2,z2,xp,yp,zp) #Distance array to PMT from each hit
            if det_type == 0: #Coated
              N_pred = md0*np.mean(d_s)+mh0*mhits+b0 #Fit prediction
            elif det_type == 1: #Uncoated
              N_pred = md1*np.mean(d_s)+mh1*mhits+b1 #Fit prediction
            

            #Run model - old methods
            #N_pred_arr = N_calc_lazy(d_s,efficiency,length,dEdx,l,det_type)
            #N_pred_arr = N_calc(d_s,efficiency,dx,dEdx,l)
            #N_pred = np.sum(N_pred_arr)
            err = err_calc(ch_hits,N_pred)
              

            #Fill in arrays
            pmt_info[i][j][1] = N_pred
            pmt_info[i][j][0] = ch
            pmt_info[i][j][2] = ch_hits
            pmt_info[i][j][3] = err
            pmt_info[i][j][4] = run
            pmt_info[i][j][5] = subrun
            pmt_info[i][j][6] = event
            pmt_info[i][j][7] = pmts.loc[ch,'ophit_opdet_type']
            pmt_info[i][j][8] = pmts.loc[ch,'distance']
            pmt_info[i][j][9] = pmts.loc[ch,'ophit_opdet_x']
            pmt_info[i][j][10] = pmts.loc[ch,'ophit_opdet_y']
            pmt_info[i][j][11] = pmts.loc[ch,'ophit_opdet_z']
            pmt_info[i][j][12] = pmts.loc[ch,'opdet_tpc']
            pmt_info[i][j][13] = row_muon['muontrk_tpc']
            pmt_info[i][j][14] = row_muon['muontrk_x1']
            pmt_info[i][j][15] = row_muon['muontrk_y1']
            pmt_info[i][j][16] = row_muon['muontrk_z1']
            pmt_info[i][j][17] = row_muon['muontrk_x2']
            pmt_info[i][j][18] = row_muon['muontrk_y2']
            pmt_info[i][j][19] = row_muon['muontrk_z2']
            pmt_info[i][j][20] = mhits
            pmt_info[i][j][21] = mtype
            pmt_info[i][j][22] = np.mean(d_s)
            pmt_info[i][j][23] = dx
            pmt_info[i][j][24] = row_muon['muontrk_length']
            
            i+=1 #Count up for PMTs
          j+=1 #Count up for muon tracks
        end = time()
        elapsed = end-start
        print(f'Run: {run} Subrun: {subrun} Event: {event} NTracks: {int(mtrks)} Time: {elapsed:.2f}s')
  print_stars()
  print(f'Ended analysis, Total tracks: {total_trks}')
  #Resulting DF
  pmt_info = pmt_info.reshape((pmts.shape[0]*muontracks.shape[0],len(columns)))
  pmt_hits_df = pd.DataFrame(pmt_info,columns=columns)
  pmt_hits_df = pmt_hits_df.set_index(['run','subrun','event'])
  pmt_hits_df_trunc = truncate_df(pmt_hits_df,keys_to_sum={'ophit_pred','ophit_obs'})

  #Remove extra zeros that may show up...
  pmt_hits_df = pmt_hits_df[pmt_hits_df.loc[:,'ophit_ch'] != 0]
  pmt_hits_df_trunc = pmt_hits_df_trunc[pmt_hits_df_trunc.loc[:,'ophit_ch'] != 0]

  #Make muon tracks DF
  mtrks = get_muon_tracks(pmt_hits_df)
    
  #Make fit params DF
  if fit_results:
    coefs,bs,scores,models = lin_fit_3d(pmt_hits_df)
    coef_mhits0 = coefs[0][0] #coefficient for mhits for coated pmts
    coef_mhits1 = coefs[1][0] #coefficient for mhits for uncoated pmts
    coef_trkdistance0 = coefs[0][0] #coefficient for mean track distance for coated pmts
    coef_trkdistance1 = coefs[1][0] #coefficient for mean track distance for uncoated pmts
    b0 = bs[0] #z-intercept for coated pmts
    b1 = bs[1] #z-intercept for uncoated pmts
    score0 = scores[0] #score for coated pmts
    score1 = scores[1] #score for uncoated pmts
    data = [[coef_mhits0,coef_trkdistance0,b0,score0,0],
            [coef_mhits1,coef_trkdistance1,b1,score1,1]] #coefs, intercepts, score, coating
    fit_df = pd.DataFrame(data,columns=['coef_mhits','coef_mean_distance','intercept','score','coating'])
    fit_df.to_pickle(cwd+'fit_params.pkl')
    print_stars()
    print('Finished fitting!')
    if move_dfs:
      os.system('mkdir -p read_model')
      os.system(f'cp {cwd}fit_params.pkl {cwd}read_model')
  else:
    #Save DFs
    mtrks.to_pickle(cwd+'muon_tracks.pkl')
    pmt_hits_df.to_pickle(cwd+'pmt_hits_df.pkl')
    pmt_hits_df_trunc.to_pickle(cwd+'pmt_hits_df_trunc.pkl')
    print_stars()
    print('Finished analyzing fit!')
    #Move DFs
    if move_dfs:
      os.system('mkdir -p read_model')
      os.system(f'mv {cwd}muon_tracks.pkl {cwd}pmt_hits_df* {cwd}read_model')

def make_analytical_df(pmts,muontracks,op_hits,events,subruns,runs,cwd,move_dfs=True):
  #pmts: Locations and coating of pmts
  #muontracks: Locations and #hits for muon tracks
  #op_hits: # of hits and waveforms of pmts
  #columns: Columns of output dataframe
  #.
  #.
  #.
  #move_dfs: Move results into read_model folder (recommended for organization)
  #This function uses analytical calculations

  #Columns for pmt info
  columns = ['ophit_ch','ophit_pred','ophit_obs','pred_err','run','subrun',
          'event','ophit_opdet_type','distance','ophit_opdet_x','ophit_opdet_y',
          'ophit_opdet_z','opdet_tpc','muontrk_tpc','muontrk_x1', 'muontrk_y1', 'muontrk_z1', 'muontrk_x2',
          'muontrk_y2', 'muontrk_z2','mhits','mtype','mean_distance','muontrk_dx','muontrk_length','ophit_pred_noatt',
          'pred_err_noatt','ophit_pred_nodisp','pred_err_nodisp'] 

  pmt_info = np.zeros((pmts.shape[0],muontracks.shape[0],len(columns)))
  j=0 #Counter for 2nd dim. (muon tracks)

  op_hits_index = op_hits.index.drop_duplicates()
  op_hits_index = op_hits_index.values

  print('Starting analysis')
  print_stars()

  total_trks = 0 #Keep track of total number of cosmics processed
  for run in runs: 
    for subrun in subruns:
      for event in events:
        
        start = time()
        index = (run,subrun,event) #Index PMT_hits

        #print(index in op_hits_index,index in muontracks.index.values,index,op_hits_index[-1],type(index))

        #Check if event is in both muon and op information
        in_op = False
        in_muon = False
        for tup in op_hits_index:
          if tup == index:
            in_op = True
        for tup in muontracks.index.values:
          if tup == index:
            in_muon = True
        if in_op and in_muon:# and subrun < 101 and run < 2: #temporary subspace to test
          #print(f'Processing Run: {run} Subrun: {subrun} Event: {event}')
          op_event_info = op_hits.loc[pd.IndexSlice[(index)]]
          muon_temp = muontracks.loc[pd.IndexSlice[(index)]]
        else:
          #print(f'Skipped Run: {run} Subrun: {subrun} Event: {event}')
          continue

        #Divide hits to tracks based on length of tracks
        for i_muon,row_muon in muon_temp.iterrows():
          #Coordinates for muon tracks
          x1 = row_muon['muontrk_x1']
          y1 = row_muon['muontrk_y1']
          z1 = row_muon['muontrk_z1']
          x2 = row_muon['muontrk_x2']
          y2 = row_muon['muontrk_y2']
          z2 = row_muon['muontrk_z2']
          length = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
          muon_temp.loc[i_muon,'muontrk_length'] = length
        #Get fraction of hits belonging to each track
        muon_temp.loc[:,'frac_hits'] = muon_temp.loc[:,'muontrk_length'].values/muon_temp.loc[:,'muontrk_length'].values.sum()
        for i_muon,row_muon in muon_temp.iterrows():
          total_trks+=1 #iterate counter
          i=0 #Counter for second dim. 
          for i_op,row_op in pmts.iterrows():
            #Skip muon type of other (mtype = 5)
            mtype = row_muon['muontrk_type']

            #Coordinates for muon tracks
            x1 = row_muon['muontrk_x1']
            y1 = row_muon['muontrk_y1']
            z1 = row_muon['muontrk_z1']
            x2 = row_muon['muontrk_x2']
            y2 = row_muon['muontrk_y2']
            z2 = row_muon['muontrk_z2']

            mtrks = row_muon['nmuontrks']
            mhits = row_muon['nhits']*row_muon['frac_hits'] #Muon track hits, normalized to track length
            tpc = row_muon['muontrk_tpc'] #Muon tpc
            
            #Coordinates for PMT
            ch = row_op['ophit_opdet']%1000
            xp = row_op['ophit_opdet_x']
            yp = row_op['ophit_opdet_y']
            zp = row_op['ophit_opdet_z']
            ch_hits = float(len(op_event_info[op_event_info['ophit_opdet'] == ch]))*row_muon['frac_hits'] #Account for double counting

            det_type = row_op['ophit_opdet_type']
            d_s,x_h,y_h,z_h,dx = distances_to_PMT(mhits,x1,y1,z1,x2,y2,z2,xp,yp,zp) #Distance array to PMT from each hit
            if det_type == 0: #Coated
              efficiency = e_ws*QE
            elif det_type == 1: #Uncoated
              efficiency = QE
            #Run model - old methods
            
            #Everything
            N_pred_arr = N_calc(d_s,efficiency,dx,dEdx,l)
            N_pred = np.sum(N_pred_arr)
            err = err_calc(ch_hits,N_pred)
              
            #No attenuation
            N_pred_arr_noatt = N_calc(d_s,efficiency,dx,dEdx,l,photon_interactions=False)
            N_pred_noatt = np.sum(N_pred_arr_noatt)
            err_noatt = err_calc(ch_hits,N_pred_noatt)

            #No dispersion
            N_pred_arr_nodisp = N_calc(d_s,efficiency,dx,dEdx,l,dispersion=False)
            N_pred_nodisp = np.sum(N_pred_arr_nodisp)
            err_nodisp = err_calc(ch_hits,N_pred_nodisp)

            #Fill in arrays
            pmt_info[i][j][1] = N_pred
            pmt_info[i][j][0] = ch
            pmt_info[i][j][2] = ch_hits
            pmt_info[i][j][3] = err
            pmt_info[i][j][4] = run
            pmt_info[i][j][5] = subrun
            pmt_info[i][j][6] = event
            pmt_info[i][j][7] = pmts.loc[ch,'ophit_opdet_type']
            pmt_info[i][j][8] = pmts.loc[ch,'distance']
            pmt_info[i][j][9] = pmts.loc[ch,'ophit_opdet_x']
            pmt_info[i][j][10] = pmts.loc[ch,'ophit_opdet_y']
            pmt_info[i][j][11] = pmts.loc[ch,'ophit_opdet_z']
            pmt_info[i][j][12] = pmts.loc[ch,'opdet_tpc']
            pmt_info[i][j][13] = row_muon['muontrk_tpc']
            pmt_info[i][j][14] = row_muon['muontrk_x1']
            pmt_info[i][j][15] = row_muon['muontrk_y1']
            pmt_info[i][j][16] = row_muon['muontrk_z1']
            pmt_info[i][j][17] = row_muon['muontrk_x2']
            pmt_info[i][j][18] = row_muon['muontrk_y2']
            pmt_info[i][j][19] = row_muon['muontrk_z2']
            pmt_info[i][j][20] = mhits
            pmt_info[i][j][21] = mtype
            pmt_info[i][j][22] = np.mean(d_s)
            pmt_info[i][j][23] = dx
            pmt_info[i][j][24] = row_muon['muontrk_length']
            pmt_info[i][j][25] = N_pred_noatt
            pmt_info[i][j][26] = err_noatt
            pmt_info[i][j][27] = N_pred_nodisp
            pmt_info[i][j][28] = err_nodisp

            i+=1 #Count up for PMTs
          j+=1 #Count up for muon tracks
        end = time()
        elapsed = end-start
        print(f'Run: {run} Subrun: {subrun} Event: {event} NTracks: {int(mtrks)} Time: {elapsed:.2f}s')
  print_stars()
  print(f'Ended analysis, Total tracks: {total_trks}')
  #Resulting DF
  pmt_info = pmt_info.reshape((pmts.shape[0]*muontracks.shape[0],len(columns)))
  pmt_hits_df = pd.DataFrame(pmt_info,columns=columns)
  pmt_hits_df = pmt_hits_df.set_index(['run','subrun','event'])
  pmt_hits_df_trunc = truncate_df(pmt_hits_df,keys_to_sum={'ophit_pred','ophit_obs'})

  #Remove extra zeros that may show up...
  pmt_hits_df = pmt_hits_df[pmt_hits_df.loc[:,'ophit_ch'] != 0]
  pmt_hits_df_trunc = pmt_hits_df_trunc[pmt_hits_df_trunc.loc[:,'ophit_ch'] != 0]

  #Make muon tracks DF
  mtrks = get_muon_tracks(pmt_hits_df)
  

  #Save DFs
  mtrks.to_pickle(cwd+'muon_tracks.pkl')
  pmt_hits_df.to_pickle(cwd+'pmt_hits_df.pkl')
  pmt_hits_df_trunc.to_pickle(cwd+'pmt_hits_df_trunc.pkl')
  print_stars()
  print('Finished analyzing fit!')
  #Move DFs
  if move_dfs:
    os.system('mkdir -p read_model')
    os.system(f'mv {cwd}muon_tracks.pkl {cwd}pmt_hits_df* {cwd}read_model')


def make_ref_df(pmts,muontracks,op_hits,events,subruns,runs,cwd,move_dfs=True):
  #pmts: Locations and coating of pmts
  #muontracks: Locations and #hits for muon tracks
  #op_hits: # of hits and waveforms of pmts
  #columns: Columns of output dataframe
  #.
  #.
  #.
  #move_dfs: Move results into read_model folder (recommended for organization)
  #This function uses analytical calculations with reflections from PMT_info_ref.pkl

  #Columns for pmt info
  columns = ['ophit_ch','ophit_pred','ophit_obs','pred_err','run','subrun',
          'event','ophit_opdet_type','ophit_opdet_x','ophit_opdet_y',
          'ophit_opdet_z','opdet_tpc','muontrk_tpc','muontrk_x1', 'muontrk_y1', 'muontrk_z1', 'muontrk_x2',
          'muontrk_y2', 'muontrk_z2','mhits','mtype','mean_distance','muontrk_dx','muontrk_length','ophit_pred_noatt',
          'pred_err_noatt','ophit_pred_nodisp','pred_err_nodisp','Reflectivity','ophit_pred_noref','pred_err_noref',
          'maxreflectivity','ophit_pred_consdisp','pred_err_consdisp'] 

  pmt_info = np.zeros((pmts.shape[0],muontracks.shape[0],len(columns)))
  j=0 #Counter for 2nd dim. (muon tracks)

  op_hits_index = op_hits.index.drop_duplicates()
  op_hits_index = op_hits_index.values

  print('Starting analysis')
  print_stars()

  total_trks = 0 #Keep track of total number of cosmics processed
  for run in runs: 
    for subrun in subruns:
      for event in events:
        
        start = time()
        index = (run,subrun,event) #Index PMT_hits

        #print(index in op_hits_index,index in muontracks.index.values,index,op_hits_index[-1],type(index))

        #Check if event is in both muon and op information
        in_op = False
        in_muon = False
        for tup in op_hits_index:
          if tup == index:
            in_op = True
        for tup in muontracks.index.values:
          if tup == index:
            in_muon = True
        if in_op and in_muon:# and subrun < 101 and run < 2: #temporary subspace to test
          #print(f'Processing Run: {run} Subrun: {subrun} Event: {event}')
          op_event_info = op_hits.loc[pd.IndexSlice[(index)]]
          muon_temp = muontracks.loc[pd.IndexSlice[(index)]]
        else:
          #print(f'Skipped Run: {run} Subrun: {subrun} Event: {event}')
          continue

        #Divide hits to tracks based on length of tracks
        for i_muon,row_muon in muon_temp.iterrows():
          #Coordinates for muon tracks
          x1 = row_muon['muontrk_x1']
          y1 = row_muon['muontrk_y1']
          z1 = row_muon['muontrk_z1']
          x2 = row_muon['muontrk_x2']
          y2 = row_muon['muontrk_y2']
          z2 = row_muon['muontrk_z2']
          length = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
          muon_temp.loc[i_muon,'muontrk_length'] = length
        #Get fraction of hits belonging to each track
        muon_temp.loc[:,'frac_hits'] = muon_temp.loc[:,'muontrk_length'].values/muon_temp.loc[:,'muontrk_length'].values.sum()
        for i_muon,row_muon in muon_temp.iterrows():
          total_trks+=1 #iterate counter
          i=0 #Counter for second dim.
          for i_op,row_op in pmts.iterrows():
            #Quick check for reflectivity
            reflectivity = row_op['Reflectivity']
            max_reflectivity = row_op['maxreflectivity'] 
            ch = row_op['ophit_opdet']
            if reflectivity == 0 and ch%1000 != row_op['ophit_opdet']:
              continue #skip events with no possible photons, but keep ones with location info
            mtype = row_muon['muontrk_type']

            #Coordinates for muon tracks
            x1 = row_muon['muontrk_x1']
            y1 = row_muon['muontrk_y1']
            z1 = row_muon['muontrk_z1']
            x2 = row_muon['muontrk_x2']
            y2 = row_muon['muontrk_y2']
            z2 = row_muon['muontrk_z2']

            mtrks = row_muon['nmuontrks']
            mhits = row_muon['nhits']*row_muon['frac_hits'] #Muon track hits, normalized to track length
            tpc = row_muon['muontrk_tpc'] #Muon tpc
            
            #Coordinates for PMT
            xp = row_op['ophit_opdet_x']
            yp = row_op['ophit_opdet_y']
            zp = row_op['ophit_opdet_z']
            tpcp = row_op['opdet_tpc'] #tpc for pmt
            ch_hits = 0
            if ch%1000 == row_op['ophit_opdet']: #Avoid double counting
              ch_hits = float(len(op_event_info[op_event_info['ophit_opdet'] == ch]))*row_muon['frac_hits'] #Account for double counting

            det_type = row_op['ophit_opdet_type']
            d_s,x_h,y_h,z_h,dx = distances_to_PMT(mhits,x1,y1,z1,x2,y2,z2,xp,yp,zp) #Distance array to PMT from each hit
            if det_type == 0: #Coated
              efficiency = e_ws*QE*reflectivity
            elif det_type == 1: #Uncoated
              efficiency = QE*reflectivity
            

            #Run model - old methods
            #N_calc(x,efficiency,length,dEdx,l,photon_interactions=True,dispersion=True)
            #N_pred_arr = N_calc_lazy(d_s,efficiency,length,dEdx,l,det_type)
            
            #Everything
            N_pred_arr = N_calc(d_s,efficiency,dx,dEdx,l)
            N_pred = np.sum(N_pred_arr)
            err = err_calc(ch_hits,N_pred)
              
            #No attenuation
            N_pred_arr_noatt = N_calc(d_s,efficiency,dx,dEdx,l,photon_interactions=False)
            N_pred_noatt = np.sum(N_pred_arr_noatt)
            err_noatt = err_calc(ch_hits,N_pred_noatt)

            #No dispersion
            N_pred_arr_nodisp = N_calc(d_s,efficiency,dx,dEdx,l,dispersion=False)
            N_pred_nodisp = np.sum(N_pred_arr_nodisp)
            err_nodisp = err_calc(ch_hits,N_pred_nodisp)

            #No reflections
            if det_type == 0 and reflectivity == max_reflectivity or det_type == 1 and reflectivity == max_reflectivity: #Only keep channel with max reflectivity (1 for no ref.)
              N_pred_arr_noref = N_calc(d_s,efficiency,dx,dEdx,l)
              N_pred_noref = np.sum(N_pred_arr_noref)
              err_noref = err_calc(ch_hits,N_pred_noref)
            
            #Constant dispersion (equal to A_pmt/A_det)
            N_pred_arr_consdisp = N_calc_lazy(d_s,efficiency,dx,dEdx,l)
            N_pred_consdisp = np.sum(N_pred_arr_consdisp)
            err_consdisp = err_calc(ch_hits,N_pred_consdisp)
            
            #Set reflected locations to 0, so they don't affect final positions
            if ch%1000 != row_op['ophit_opdet']:
              xp = yp = zp = 0


            #Fill in arrays
            pmt_info[i][j][1] = N_pred
            pmt_info[i][j][0] = ch
            pmt_info[i][j][2] = ch_hits
            pmt_info[i][j][3] = err
            pmt_info[i][j][4] = run
            pmt_info[i][j][5] = subrun
            pmt_info[i][j][6] = event
            pmt_info[i][j][7] = det_type
            pmt_info[i][j][8] = xp
            pmt_info[i][j][9] = yp
            pmt_info[i][j][10] = zp
            pmt_info[i][j][11] = tpcp
            pmt_info[i][j][12] = row_muon['muontrk_tpc']
            pmt_info[i][j][13] = row_muon['muontrk_x1']
            pmt_info[i][j][14] = row_muon['muontrk_y1']
            pmt_info[i][j][15] = row_muon['muontrk_z1']
            pmt_info[i][j][16] = row_muon['muontrk_x2']
            pmt_info[i][j][17] = row_muon['muontrk_y2']
            pmt_info[i][j][18] = row_muon['muontrk_z2']
            pmt_info[i][j][19] = mhits
            pmt_info[i][j][20] = mtype
            pmt_info[i][j][21] = np.mean(d_s)
            pmt_info[i][j][22] = dx
            pmt_info[i][j][23] = row_muon['muontrk_length']
            pmt_info[i][j][24] = N_pred_noatt
            pmt_info[i][j][25] = err_noatt
            pmt_info[i][j][26] = N_pred_nodisp
            pmt_info[i][j][27] = err_nodisp
            pmt_info[i][j][28] = reflectivity
            pmt_info[i][j][29] = N_pred_noref
            pmt_info[i][j][30] = err_noref
            pmt_info[i][j][31] = max_reflectivity
            pmt_info[i][j][32] = N_pred_consdisp
            pmt_info[i][j][33] = err_consdisp

            i+=1 #Count up for PMTs
          j+=1 #Count up for muon tracks
        end = time()
        elapsed = end-start
        print(f'Run: {run} Subrun: {subrun} Event: {event} NTracks: {int(mtrks)} Time: {elapsed:.2f}s')
  print_stars()
  print(f'Ended analysis, Total tracks: {total_trks}')
  #Resulting DF
  pmt_info = pmt_info.reshape((pmts.shape[0]*muontracks.shape[0],len(columns)))
  pmt_hits_df = pd.DataFrame(pmt_info,columns=columns)
  pmt_hits_df = pmt_hits_df.set_index(['run','subrun','event'])
  pmt_hits_df_trunc = truncate_df(pmt_hits_df,keys_to_sum={'ophit_pred','ophit_obs','ophit_opdet_x','ophit_opdet_y',
          'ophit_opdet_z'},mod_chs=True) #mod channels

  #Remove extra zeros that may show up...
  pmt_hits_df = pmt_hits_df[pmt_hits_df.loc[:,'ophit_ch'] != 0]
  pmt_hits_df_trunc = pmt_hits_df_trunc[pmt_hits_df_trunc.loc[:,'ophit_ch'] != 0]

  #Make muon tracks DF
  mtrks = get_muon_tracks(pmt_hits_df)
  

  #Save DFs
  mtrks.to_pickle(cwd+'muon_tracks.pkl')
  pmt_hits_df.to_pickle(cwd+'pmt_hits_df.pkl')
  pmt_hits_df_trunc.to_pickle(cwd+'pmt_hits_df_trunc.pkl')
  print_stars()
  print('Finished analyzing fit!')
  #Move DFs
  if move_dfs:
    os.system('mkdir -p read_model')
    os.system(f'mv {cwd}muon_tracks.pkl {cwd}pmt_hits_df* {cwd}read_model')
