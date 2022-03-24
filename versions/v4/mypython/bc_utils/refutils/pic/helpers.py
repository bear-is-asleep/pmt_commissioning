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
sys.path.append('/sbnd/app/users/brindenc/mypython')
import seaborn as sns
import matplotlib
from scipy import optimize
from sklearn.linear_model import LinearRegression
from numpy import isin, sqrt,exp,arctan,cos,sin,pi
from time import time
from bc_utils.utils import pic
from math import comb
import scipy

#Constants
R_cpa = 0 #VUV
R_apa = R_fc = 0.26 #VUV
R_cpa_vis = 0.92 #vis
R_apa_vis = R_fc_vis = 0.585 #vis


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
  L = 213.4 #Length of x-coordinate
  x_prime=0

  if flip % 2 == 0: #cpa
    if x > 0: #cpa is on other side for opposite tpc
      cpa_ref = int((ref+1)/2)
    else:
      cpa_ref = int(ref/2)
    x_prime = x-L*ref #Since x is always the same
  else: #apa
    if x > 0: #cpa is on other side for opposite tpc
      cpa_ref = int(ref/2)
    else:
      cpa_ref = int((ref+1)/2)
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
      coord_ref = [[x,y,z,0,0,0,0,0,0]] #6 indexes for 3 coords. total # reflections, cpa ref., apa ref.
    else:
      coord_ref = [[x,y,z,0,0,0,ch,0,0,0]] #7 indexes for 3 coords. total # reflections, cpa ref., apa ref., ch #
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
                x_ref_val = 0 #x ref for current element
                y_ref_val = 0 #y ref for current element
                z_ref_val = 0 #z ref for current element
                if i != 0:
                  x_ref_val += x_ref
                  tot_ref += x_ref
                  cpa_ref = 0
                  apa_ref = 0
                if j != 0:
                  y_ref_val += y_ref
                  tot_ref += y_ref
                if k != 0:
                  z_ref_val += z_ref
                  tot_ref += z_ref
                if i == 1:
                  cpa_ref = cpa_ref1
                  apa_ref = apa_ref1
                if i == 2:
                  cpa_ref = cpa_ref2
                  apa_ref = apa_ref2
                #6 indexes for 3 coords. total # reflections, cpa ref., apa ref.
                if ch == 0:
                  coord_ref.append([xal,yal,zal,tot_ref,cpa_ref,apa_ref,x_ref_val,y_ref_val,z_ref_val])
                else: 
                  ch+=1000
                  coord_ref.append([xal,yal,zal,tot_ref,cpa_ref,apa_ref,ch,x_ref_val,y_ref_val,z_ref_val])
  if ch == 0:
    df = pd.DataFrame(coord_ref,columns=['x','y','z','tot_ref','cpa_ref','apa_ref','x_ref','y_ref','z_ref']).drop_duplicates()
  else: 
    df = pd.DataFrame(coord_ref,columns=['ophit_opdet_x','ophit_opdet_y','ophit_opdet_z',
    'tot_ref','cpa_ref','apa_ref','ophit_opdet','x_ref','y_ref','z_ref']).drop_duplicates()
  return df

def flatten_list(_2d_list):
    flat_list = []
    # Iterate through the outer list
    for element in _2d_list:
        if type(element) is list:
            # If the element is of type list, iterate through the sublist
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list

def convert(list):
      
    # Converting integer list to string list
    s = [str(i) for i in list]
      
    # Join list items using join()
    res = int("".join(s))
      
    return(res)
def sortstrings_numerically(strings,sort=True,drop_ints=True):
  #Return tuple, matching string to int
  ints = []
  for s in strings:
    ints.append([int(x) for x in s if x.isdigit()])
  int_new = []
  for l in ints:
    int_new.append(convert(l))
  tups = [[]]
  for i in range(len(strings)):
    tups.append([strings[i],int_new[i]])
  #tups.pop(0)
  tups = tups[1:]
  if sort:
    tups = sorted(tups,key=lambda x: x[1])
  else:
    tups = tups
  if drop_ints:
    for j in tups:
      del j[1]
    tups = flatten_list(tups)
  return tups

def calc_prob_vis(cpa_ref):
  #Calculate the probability of a photon being vis given # of cpa_ref
  return 1 - scipy.stats.binom.pmf(1,cpa_ref,pic.e_ws)

  

def get_ref_coefs(singref,coating):
  #singref is a single PMT with all of its reflections (see make_single_reflection)
  Rs = np.zeros(singref.shape[0]) #calculates reflectivity based on the surface in question
  Ref_coefs = np.array([R_cpa,R_apa,R_fc]) #Reflectivity of cpa,apa,field cage VUV
  Ref_coefs_vis = np.array([R_cpa_vis,R_apa_vis,R_fc_vis]) #Reflectivity of cpa,apa,field cage VUV
  for row,line in singref.iterrows():
    #Get number of reflections in each plane
    x_ref = line['x_ref']
    y_ref = line['y_ref']
    z_ref = line['z_ref']
    cpa_ref = line['cpa_ref']
    apa_ref = line['apa_ref']
    
    refs = [x_ref,y_ref,z_ref,cpa_ref,apa_ref]

    #Map to Rs array
    R=1
    vis_prob = calc_prob_vis(cpa_ref) #Probability of being visible photon (using binomial theorem)
    vuv_prob = 1 - vis_prob
    for i,ref in enumerate(refs):
      if ref != 0: #There is a reflection
        if i == 1 or i == 2: #fc reflection
          R*=Ref_coefs[2]**ref*vuv_prob+Ref_coefs_vis[2]**ref*vis_prob
        elif i == 3: #cpa ref
          if coating == 0:
            R*=(Ref_coefs[0])**ref*vuv_prob+(Ref_coefs_vis[0])**ref*vis_prob
          else:
            R*=(Ref_coefs[0])**ref*vuv_prob+(Ref_coefs_vis[0])**ref*vis_prob
        elif i == 4: #apa ref
          R*=Ref_coefs[1]**ref*vuv_prob+Ref_coefs_vis[1]**ref*vis_prob
    if coating == 1 and x_ref == 0 and y_ref == 0 and z_ref == 0: #set coated PMTs initial R to 0
      Rs[row] = 0
    else:
      Rs[row] = R
  singref.loc[:,'Reflectivity'] = Rs
  return singref

def get_ref_PMTs(pmts,n_ref=1):
  #pmts is PMT_info.pkl dataframe

  new_pmts = pd.DataFrame(columns=['ophit_opdet_x','ophit_opdet_y','ophit_opdet_z','ophit_opdet_type',
  'opdet_tpc','tot_ref','cpa_ref','apa_ref','ophit_opdet','x_ref','y_ref','z_ref','Reflectivity'])

  for row,line in pmts.iterrows():
    ch = row
    x = line['ophit_opdet_x']
    y = line['ophit_opdet_y']
    z = line['ophit_opdet_z']
    coating = line['ophit_opdet_type']
    tpc = line['opdet_tpc']
    singref = single_plane_reflection(x,y,z,flip_coord='xyz',x_refs=n_ref,y_refs=n_ref,
    z_refs=n_ref,initialize=True,ch=ch)
    singref = get_ref_coefs(singref,coating)
    singref.loc[:,'ophit_opdet_type'] = coating
    singref.loc[:,'opdet_tpc'] = tpc
    singref.loc[:,'maxreflectivity'] = max(singref.loc[:,'Reflectivity'])
    frames = [new_pmts,singref]
    new_pmts = pd.concat(frames) #Merge frames
  return new_pmts