LoadSolv

This software is used to invert GNSS observed surface displacements to infer changes in water storage within a chosen area.

To use this code you need the following ingredients:
        - Design Matrix (Jacobian): Dataset describing the predicted displacement response of each GNSS station within a study area (e.g., Western US) to a unit load in each pixel of
          your model domain (e.g., 0.25 degree pixel grid in the Western US)
        - GNSS Displacement Data: A file (or series of files) that contain the displacement of GNSS stations on a given epoch (e.g., January 2006). This file should include the following:
          Station Name, Latitude, Longitude, E Displacement [mm], N Displacement [mm], V Displacement [mm], E Sigma [mm], N Sigma [mm], V Sigma [mm]

This code base relies upon common Python packages such as: NumPy, SciPy, Pandas, datetime, Netcdf4, mpi4py, sys, os, matplotlib

Results presented in "Extreme Winter Precipitation Drives Recharge of Deep Mountain Groundwater" by Swarr et al. (2026) were produced using the following setting in the script run_solv.py:

# Design Matrix File
dm_file = ("input/DesignMatrixLoad/designmatrix_cf_PREM_wusa_grid_0.25_06-22_gps_wusa_AR.nc")

# Data Files
# Example: MEaSUREs Monthly GNSS displacements between Jan, 2006 and Sep, 2022
data_dir = ("input/GNSS/meas_rmgrace_five/")
# Example: MEaSUREs daily GNSS displacements between Oct 1, 2022 and June 28, 2024
data_dir = ("input/GNSS/meas_rmgldas_five/")

# Data Prefix
# This will find all data files with the prefix "20" and solve the inverse problem using each data file
data_pre = ("20")

# We are only using vertical displacements
uonly = True

# We are not solving a weighted inverse problem. However, we do read in the observational uncertainties to compute formal parameter uncertainties
weighted = False

# Regularization parameter values
alpha = 3.5
beta = 1.0

# Choice of regularization
# We want to include second order Tikhonov regularization so this flag ensure that the 2-D finite difference Laplacian is constructed for our use
tikhonov = 'second'

# Constants used for scaling our output estimates of water thickness returned after solving the inverse problem
# These do not need to be changed unless different values were used in creating the design matrix of the inverse problem
# These are the defaults used when constructing a design matrix in LoadDef
ref_height = 1.
ref_density = 1000.

To solve the inverse problem for all months between Jan, 2006 and Sep, 2022 using MEaSUREs displacement data you would use the following inputs in run_solv.py
data_dir = ("input/GNSS/meas_rmgrace_five/")
data_pre = ("20")

You would then use the following command in your terminal "python run_solv.py"
