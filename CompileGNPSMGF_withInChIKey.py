#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 10:23:22 2018

@author: madeleineernst
"""
#############################################################################################################
#                                                                                                           #
#  compile .mgf file of unique InChIKeys and associated MS/MS spectra retrieved from GNPS Libraries         #
#                                                                                                           #
#############################################################################################################

from pyteomics import mgf, auxiliary
import pandas as pd
import numpy as np

# load GNPS library matches downloaded from https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=03fba62d93cb4cbfa3f72106d18f7d2c
lib = pd.read_csv("ProteoSAFe-MOLECULAR-LIBRARYSEARCH-03fba62d-view_all_annotations_DB/MOLECULAR-LIBRARYSEARCH-03fba62d-view_all_annotations_DB-main.tsv", sep = "\t")
# load SMILES and associated InChIKeys 
# SMILES were converted to InChIKeys using Marvin Beans MolConverter: https://chemaxon.com/marvin-archive/3.3.3/marvin/doc/user/molconvert.html
smiles = pd.read_csv("SMILES_GNPSLibraries.csv",sep=',', index_col = 0)
ikeys = pd.read_csv("InchiKeys_GNPSLibraries.txt",  sep='\t',header = None)

lib = lib.dropna(subset=['Smiles'])  
lib = lib[lib.Smiles != ' ']
lib['Smiles'] = lib['Smiles'].str.strip() # remove white spaces

ikeys = [j for i in ikeys.values.tolist() for j in i]
ikeys = [w.replace('InChIKey=', '') for w in ikeys]
smiles["inchikey"] = ikeys
smiles = smiles.rename(columns = {'SMILES':'Smiles'})
smiles = smiles.drop_duplicates(subset=['Smiles']) 

libcomb = pd.merge(lib, smiles,how="left",on="Smiles")
libcomb = libcomb.dropna(subset=['inchikey'])
libcomb = libcomb.drop_duplicates(subset='inchikey', keep='first', inplace=False) # remove duplicate InChIKeys

# load GNPS .mgf file downloaded from https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=6e22f85aeb0744208e872d1640f508d9
mgf_file = 'ProteoSAFe-METABOLOMICS-SNETS-6e22f85a-download_cluster_buckettable/METABOLOMICS-SNETS-6e22f85a-download_clustered_spectra-main_ChargeReplaced.mgf'
scans = libcomb.Scan.tolist()

counter=0
with mgf.read(mgf_file) as reader:
    for spectrum in reader:
        for idx, scan in enumerate(scans):
            if spectrum['params']['scans'] == str(scan):
                file_name = '{}.mgf'.format("GNPSLibraries_uniqueSMILES_withFeatureIDs")
                spectrum['params']['SMILES'] = libcomb.Smiles.tolist()[idx]
                spectrum['params']['InchiKey'] = libcomb.inchikey.tolist()[idx]
                counter+=1
                spectrum['params']['FEATURE_ID'] = counter
                mgf.write((spectrum,), file_name)