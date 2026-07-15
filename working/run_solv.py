#!/usr/bin/env python

# *********************************************************************
# MAIN PROGRAM TO INVERT OBSERVED SURFACE DISPLACEMENTS FOR SURFACE-
#  LOAD DISTRIBUTION
#
# Copyright (c) 2026
#
#    LoadSolv is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    LoadSolv is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with LoadSolvC.  If not, see <https://www.gnu.org/licenses/>.
#
# *********************************************************************

# IMPORT PRINT FUNCTION
from __future__ import print_function

# IMPORT MPI MODULE
from mpi4py import MPI

# MODIFY PYTHON PATH TO INCLUDE 'LoadDef' DIRECTORY
import sys
import os
sys.path.append(os.getcwd() + "/../")

# IMPORT PYTHON MODULES
import numpy as np
import scipy as sc
import datetime
import netCDF4 
from solv import load_inversion
from utility import read_loadDesignMatrix

# --------------- SPECIFY USER INPUTS --------------------- #

### Design Matrix file (can be created using "LoadDef" software or other softwares which predict the elastic deformation to the Earth to a surface load (e.g., PyLith, SPOTL))
dm_file = ("../input/DesignMatrixLoad/designmatrix_cf_PREM_cells_-56.0_13.0_275.0_328.0_0.5_commonMesh_regional_-56.0_13.0_275.0_328.0_0.01_0.01_oceanmask_gps_locations_SA_UNR.nc") 
### Data directory 
###  Datafile Format: Station, Latitude[+N], Longitude[+E], East-Displacement[mm], North-Displacement[mm], Up-Displacement[mm]
###      If only vertical-component data are available, please set "uonly" flag below to True. 
data_dir = ("../input/GNSS/ngl20_sa_new/")

### Datafile prefix
data_pre = ("201710")

### Does the datafile contain only vertical-component data? If so, set uonly = True. Default is uonly = False.
### Datafile Format: Station, Latitude[+N], Longitude[+E], Up-Displacement[mm]
### [NOTE: It is always assumed that the Jacobian is built based on three components; this parameter applies only to the data!]
uonly = True

### A Priori Data File
### Datafile format: Latitude[+N], Longitude[+E], Amplitude[m]
### [NOTE: It is assumed that the A Priori data is aligned with the load cells being solved for.]
apriori_file = ("../input/APriori/GRACE_edge_0517.txt")

### Is A Priori data to be used? If so, set apriori = True. Default is apriori = False.
apriori = False

### Weighted inversion? If so, set weighted = True. Default is weighted = False.
### [NOTE: You must provide GNSS uncertainty with input data file]
weighted = False

### Tikhonov regularization parameter and order (options for order: 'zeroth' or 'second')
### Beta is to be supplied if using zeroth- and second-order Tikhonov regularization. Default is beta = 0.0
alpha = 1.0
beta = 0.0
tikhonov = 'second'

### Reference Height (m) and Density (kg/m^3) of Load Used to Compute Design Matrix
### [NOTE: It is assumed that a design matrix is built by placing a uniform 1 meter freshwater load in each pixel of the model grid.]
### [NOTE: If you have used a different height of water or density these parameters must be altered to reflect this.]
ref_height = 1.
ref_density = 1000.

### Output file suffix
outfile = ("_TikReg-" + tikhonov + "-" + str(beta) + "-" + str(alpha) + ".txt")

# -------------------- SETUP MPI --------------------------- #

# Get the Main MPI Communicator That Controls Communication Between Processors
comm = MPI.COMM_WORLD
# Get My "Rank", i.e. the Processor Number Assigned to Me
rank = comm.Get_rank()
# Get the Total Number of Other Processors Used
size = comm.Get_size()

# ---------------------------------------------------------- #

# -------------------- BEGIN CODE -------------------------- #

# Ensure that the Output Directories Exist
if (rank == 0):
    if not (os.path.isdir("../output/")):
        os.makedirs("../output/")
    if not (os.path.isdir("../output/Inversion/")):
        os.makedirs("../output/Inversion/")
    if not (os.path.isdir("../output/Inversion/SurfaceLoad/")):
        os.makedirs("../output/Inversion/SurfaceLoad/")

# Determine Number of Available Data Files
data_files = []
print(data_dir)
if os.path.isdir(data_dir):
    for mfile in os.listdir(data_dir): # Filter by Data Directory
        if mfile.startswith(data_pre): # Filter by Datafile Prefix
            data_files.append(data_dir + mfile) # Append File to List
else:
    sys.exit('Error: Data directory not found.')

# Test for Data Files
if not data_files:
    sys.exit('Error: Data files not found.')

# Determine Number of Files Read In
if isinstance(data_files,float) == True: # only 1 file
    numel = 1 
else:
    numel = len(data_files)

# Determine the Chunk Sizes for the Inversion
total_files = len(data_files)
nominal_load = total_files // size # Floor Divide
# Final Chunk Might Be Different in Size Than the Nominal Load
if rank == size - 1:
    procN = total_files - rank * nominal_load
else:
    procN = nominal_load

# Perform the Inversion(s)
# Primary Rank
if (rank == 0):
    model_vector = load_inversion.main(dm_file,data_files,apriori_file,rank,procN,comm,reference_height=ref_height,reference_density=ref_density,alpha=alpha,beta=beta,tikhonov=tikhonov,outfile=outfile,uonly=uonly,apriori=apriori,weighted=weighted)
# Worker Ranks
else:
    load_inversion.main(dm_file,data_files,apriori_file,rank,procN,comm,reference_height=ref_height,reference_density=ref_density,alpha=alpha,beta=beta,tikhonov=tikhonov,outfile=outfile,uonly=uonly,apriori=apriori,weighted=weighted)

# --------------------- END CODE --------------------------- #


