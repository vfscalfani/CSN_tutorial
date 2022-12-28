#!/usr/bin/env python
# coding: utf-8

# # Chemical Space Network Calculations (More Efficient Memory Use)
# 
# **V.F. Scalfani, V.D. Patel, and A.M. Fernandez** \
# v. December 20, 2022
# 
# This is a revised version of the CSN calculations that uses less memory.
# In the initial version we were unnecessarily duplicating mol objects in the
# subsets variable for the same SMILES, so now we create a much smaller
# separate dictionary of mol objects.

### 1. Import RDKit, Networkx, and other libraries

# for timing of the calculations
from time import time
start_time = time()

from pprint import pprint

# RDKit stuff
from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFMCS
from rdkit import DataStructs
from rdkit.Chem import rdmolops

# numpy
import numpy as np

# pandas
import pandas as pd

# networkx
import networkx as nx

# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt

# Print versions of libraries used
print('RDKit version: ',rdBase.rdkitVersion)
print('Numpy version:', np.__version__)
print('Pandas version:', pd.__version__)
print('Networkx version',nx.__version__)
print('MatplotLib version:', mpl.__version__)

### 2. Load ChEMBL Dataset

# pd.options.display.max_rows = 30
df = pd.read_csv('glucocorticoid_receptor_2034_2.csv', sep =';')

### 3. Data Preparation and Checks

# Note that we are using the chemical structures as-is (i.e., using ChEMBL's standardization workflow).
# See supporting manuscript for more details. 
# Create a new dataframe with only Chembl ID, Smiles, and Standard Value (Ki, for this example)
# We are ignoring the Standard relation here for Ki values (i.e., >, <, etc. are treated as =)

df1 = df[['Molecule ChEMBL ID','Smiles','Standard Value']].copy()

# drop any rows with NaN (missing Standard Values)
df1.dropna(inplace=True)
print("length of df1: ",len(df1))

# Check for presence of disconnected SMILES notation via string matching 
df2 = df1[~df1['Smiles'].str.contains("\.")]
print("length of df2: ",len(df2))

# We should double check for disconnected fragments in the event that
# the dot disconnect bond is used with ring-closures
# see: http://www.dalkescientific.com/writings/diary/archive/2004/12/12/library_generation_with_smiles.html

smis = df1['Smiles'].tolist()
num_frags = []
for smi in smis:
    mol = Chem.MolFromSmiles(smi)
    num_frags.append(len(Chem.GetMolFrags(mol))) # returns number of fragments

# now check that all molecules have only one fragment
print("all fragments equal to 1:")
print(all(frag == 1 for frag in num_frags))

# Group by ChEMBLID/Smiles rows, and then take the mean of the standard Ki value to account for duplicates.
# See: Zhang, B. et al. J Comput Aided Mol Des 2015, 29 (10), 937–950. 
# https://doi.org/10.1007/s10822-015-9872-1.
df3 = df2.groupby(['Molecule ChEMBL ID', 'Smiles'])['Standard Value'].mean().reset_index()

# rename the Standard Value Column to Ki
df3.rename(columns={'Standard Value': 'Ki'}, inplace=True)

# Double check that all ChEMBL IDs are now unique
chembl_ids = df3['Molecule ChEMBL ID'].tolist()

# create a set, if the set length matches the list, there are no duplicates
set_chembl_ids = set(chembl_ids)
print("length of set_chembl_ids is equal to length chembl_ids")
print(len(set_chembl_ids) == len(chembl_ids))

# Double check that all SMILES are unique (different) compounds
# To be on the safe side, we can parse the SMILES as RDKit mol objects
# then write out canonical smiles and check

smis = df3['Smiles'].tolist()
rdkit_can_smiles = []
for smi in smis:
    mol = Chem.MolFromSmiles(smi)       
    rdkit_can_smiles.append(Chem.MolToSmiles(mol, canonical=True)) # default is true
    
set_rdkit_can_smiles = set(rdkit_can_smiles)
print("length of set_rdkit_can_smiles is equal to length rdkit_can_smiles")
print(len(set_rdkit_can_smiles) == len(rdkit_can_smiles))

# We'll use the original SMILES as unique dictionary keys, so we should verify that the
# ChEMBL SMILES are unique strings too.
set_smis = set(smis)

print("length of set_smis is equal to length smis")
print(len(set_smis) == len(smis))

### 4. Compile Node data

# We will use pKi for coloring nodes.
from math import log10

# Ki values are in nM units
# 1. convert to M
df3.loc[:,"Ki_M"] = (df3.loc[:,"Ki"] * (10**-9))

# 2. then compute -log10[Ki]
def minuslog(x):
    return -log10(x)

df3.loc[:,"pKi"] = (df3.loc[:,"Ki_M"].apply(minuslog))

# drop Ki and Ki_M columns (no longer needed)
df3.drop(["Ki","Ki_M"],axis=1,inplace=True)

# get the max/min Ki values, which will be used for a color bar later
print("max and min pKi values:")
print(df3.pKi.max())
print(df3.pKi.min())

# set the dataframe index as Smiles (we already verified they are all unique from eachother)
df4 = df3.set_index('Smiles')

pd.set_option("expand_frame_repr", False)

print("Fist 20 data rows: ")
print(df4.head(20)) # view first 20

# save to a dictionary
node_data = df4.to_dict('index')

# SMILES are the keys
list(node_data.keys())[0]

# ChEMBL ID and pKi are the associated values
list(node_data.values())[0]

# print(node_data)

### 5. Compute and Compile Edge Data

# We first need to create subset pairs of the SMILES
smis = [] # using ChEMBL provided SMILES
for key,value in node_data.items():
    smis.append(key)

from itertools import combinations
smis_subsets = list(combinations(smis,2))
print("length of smis_subsets: ")
print(len(smis_subsets))

# View first 5
# print(smis_subsets[0:5])

# create a dictionary of smiles and corresponding mol object
molObjects = {}
for smi in smis:
    molObjects.update({smi: Chem.MolFromSmiles(smi)})

# View first 5
# print(list(molObjects.items())[0:5])

# create a dictionary, subsets
subsets = {}
for i, (smi1,smi2) in enumerate(smis_subsets):
    field = {}
    field["smi1"] = smi1
    subsets[i] = field
    
    field["smi2"] = smi2
    subsets[i] = field
len(subsets)

# Indexing examples
# get first key
list(subsets)[0]

# get first value
list(subsets.values())[0]

# get smi1
list(subsets.values())[0]['smi1']

list(subsets.keys())[0]

list(subsets.values())[0]

# ### Compute RDKit fingerprint-based Tanimoto Similarity (multiprocessing)

# get number of processors
import multiprocessing
print("multiprocessing cpu count:", multiprocessing.cpu_count())

# From the Python docs, this below is number of usable CPUs (works on Unix/Linux)
# https://docs.python.org/3/library/multiprocessing.html
# we subtracted 2 from total number, so that we can still easily use computer for other tasks
import os
num_cpus = len(os.sched_getaffinity(0)) - 2
print("num_cpus used: ", num_cpus)

# create a list of mol1, mol2, and their dictionary key as tuples
mol_tuples = []
for key, value in subsets.items():
    mol_tuples.append((molObjects[value["smi1"]],molObjects[value["smi2"]], key))

print("Starting computing Tanimoto Similarity with RDKit fingerprints...")

# compute and add Tanimoto Similarity using default RDKit fingerprints
def t_sim(mol1,mol2,key):
    fp1 = Chem.RDKFingerprint(mol1)
    fp2 = Chem.RDKFingerprint(mol2)
    tan_sim = round(DataStructs.TanimotoSimilarity(fp1,fp2), 3)
    return key, tan_sim

# run multiprocessing on the t_sim function
from multiprocessing import Pool
if __name__ == '__main__':
    with Pool(num_cpus) as p:
        star_map = p.starmap(t_sim, mol_tuples)
    for key, tan_sim in star_map:
        subsets[key].update({"tan_similarity": tan_sim})

pprint(list(subsets.values())[0:5])

# ### Compute MCS-based Tanimoto Coefficient (multiprocessing)

# get number of processors
import multiprocessing
print("multiprocessing cpu count:", multiprocessing.cpu_count())

# From the Python docs, this below is number of usable CPUs (works on Unix/Linux)
# https://docs.python.org/3/library/multiprocessing.html
# we subtracted 2 from total number, so that we can still easily use computer for other tasks
import os
num_cpus = len(os.sched_getaffinity(0)) - 2
print("num_cpus used: ", num_cpus)

# Add maximum common substructure (MCS)-based Tanimoto Coefficient
# See Vogt, M. et al. J Comput Aided Mol Des 2016, 30 (3), 191–208. 
# https://doi.org/10.1007/s10822-016-9906-3.

############
############

# Compute and add MCS-based similarity

# Here are benchmark times with the 10 second timeout in FindMCS and the 81,000 compound gluco pairs:

    # Intel Core i9-12900K using 22 of the 24 cores: ~25 minutes
    
#############
#############

print("Starting MCS calculations...")

def tc_mcs(mol1,mol2,key):
    # get maximum common substructure instance
    mcs = rdFMCS.FindMCS([mol1,mol2],timeout=10) # adding a 10 second timeout
    
    # get number of common bonds
    mcs_bonds = mcs.numBonds
    
    # get number of bonds for each
    # default is only heavy atom bonds
    mol1_bonds = mol1.GetNumBonds()
    mol2_bonds = mol2.GetNumBonds()
    
    # compute MCS-based Tanimoto
    tan_mcs = mcs_bonds / (mol1_bonds + mol2_bonds - mcs_bonds)
    return key, tan_mcs

# run multiprocessing on the tc_mcs function
from multiprocessing import Pool

if __name__ == '__main__':
    with Pool(num_cpus) as p:
        star_map = p.starmap(tc_mcs, mol_tuples)
    for key, tan_mcs in star_map:
        subsets[key].update({"tan_mcs": round(tan_mcs,3)})

pprint(list(subsets.values())[0:5])

# Keys are integers
list(subsets.keys())[0:5]

# Here is what the star_map variable looks like
star_map[0:5]

### 6. Save Data

# Save the data, so you don't have to re-compute the tanimoto and MCS similarity again.
# Save the subsets data as a pickle
import pickle
with open('subsets.pickle', 'wb') as outfile:
    pickle.dump(subsets, outfile, pickle.HIGHEST_PROTOCOL)

# Save the node data
with open('node_data.pickle', 'wb') as outfile:
    pickle.dump(node_data, outfile, pickle.HIGHEST_PROTOCOL)

print("pickle files saved!")
print(f"Final Time: {time()-start_time} seconds")

# single core processing code below (not used)

# compute and add Tanimoto Similarity using default RDKit fingerprints
# for key,value in subsets.items():
#     fp1 = Chem.RDKFingerprint(molObjects[value["smi1"]])
#     fp2 = Chem.RDKFingerprint(molObjects[value["smi2"]])
#     tan_sim = round(DataStructs.TanimotoSimilarity(fp1,fp2), 3)
#     subsets[key].update({"tan_similarity": tan_sim})

# Add maximum common substructure (MCS)-based Tanimoto Coefficient
# See Vogt, M. et al. J Comput Aided Mol Des 2016, 30 (3), 191–208. 
# https://doi.org/10.1007/s10822-016-9906-3.

############
############

# Compute and add MCS-based similarity
# This may take a very long time...

#############
#############

#def tc_mcs(mol1,mol2):
#    # get maximum common substructure instance
#    mcs = rdFMCS.FindMCS([mol1,mol2],timeout=10) # adding a 10 second timeout for now
    
#    # get number of common bonds
#    mcs_bonds = mcs.numBonds
    
#    # get number of bonds for each
#    # default is only heavy atom bonds
#    mol1_bonds = mol1.GetNumBonds()
#    mol2_bonds = mol2.GetNumBonds()
    
#    # compute MCS-based Tanimoto  
#    tan_mcs = mcs_bonds / (mol1_bonds + mol2_bonds - mcs_bonds)
#    return tan_mcs


# loop through subsets and compute tc_mcs

#for key,value in subsets.items():
#    tan_mcs_value = round(tc_mcs(molObjects[value["smi1"]], molObjects[value["smi2"]]), 3)
#    print(key) # to watch progress
#    subsets[key].update({"tan_mcs": tan_mcs_value})


