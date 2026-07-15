# *********************************************************************
# FUNCTION TO INVERT SURFACE DISPLACEMENTS FOR SURFACE LOAD
# 
# Copyright (c) 2026
#
# This file is part of LoadSolv.
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
#    along with LoadSolv.  If not, see <https://www.gnu.org/licenses/>.
#
# *********************************************************************

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
from scipy.sparse import linalg
from utility import read_datafile
from utility import read_datafile_uonly
from utility import read_datafile_weighted
from utility import read_datafile_uonly_weighted
from utility import read_apriori_datafile

def main(datafile,fid,apriori_file,design_matrix,sta_ids,sta_comp_ids,sta_comp_lat,sta_comp_lon,load_cell_ids,load_cell_lat,load_cell_lon,alpha,beta,tikhonov,nine_point=False,uonly=False,apriori=False,weighted=False):

    # Read the Datafile
    if (weighted == True):
        if (uonly == True): 
            #   Format: Station, Latitude[+N], Longitude[+E], Up-Displacement[mm], Up-Sigma[mm]
            sta,slat,slon,udisp,usig = read_datafile_uonly_weighted.main(datafile)
        else:
            #   Format: Station, Latitude[+N], Longitude[+E], East-Displacement[mm], North-Displacement[mm], Up-Displacement[mm], East-Sigma[mm], North-Sigma[mm], Up-Sigma[mm]
            sta,slat,slon,edisp,ndisp,udisp,esig,nsig,usig = read_datafile_weighted.main(datafile)
    else:
        if (uonly == True):
            #   Format: Station, Latitude[+N], Longitude[+E], Up-Displacement[mm]
            # Using this utility file to read in uncertainties here as we will use them later
            sta,slat,slon,udisp,usig = read_datafile_uonly_weighted.main(datafile)
        else:
            #   Format: Station, Latitude[+N], Longitude[+E], East-Displacement[mm], North-Displacement[mm], Up-Displacement[mm]
            sta,slat,slon,edisp,ndisp,udisp = read_datafile.main(datafile)
    
    # Find Stations that Match Design Matrix File
    # Otherwise, Remove Rows from Design Matrix where No Data Exists for that Station
    # Otherwise, Remove Observations where No Design Matrix Information Exists for that Station
    # Also, ensure correct order for stations in Design Matrix and in Data Vector
    common_sta,data_idx,dm_idx = np.intersect1d(sta,sta_ids,return_indices=True)
    slat = slat[data_idx]
    slon = slon[data_idx]
    if (len(data_idx) != len(dm_idx)):
        sys.exit(':: Error: Mismatch in stations. [perform_inversion.py]')
    # Build G matrix, data vector, and weighting matrix differently depending on number of components
    if (weighted ==True):
        if (uonly == True): # one component
            # Initialize updated Jacobian
            G = np.empty((len(dm_idx),len(load_cell_ids)))
            # Initialize datavector
            d = np.empty((len(data_idx),))
            # Initialize weighting vector
            w = np.empty((len(data_idx),))
            for bb in range(0,len(dm_idx)):
                crow_dm = dm_idx[bb]
                crow_data = data_idx[bb]
                G[bb,:] = design_matrix[crow_dm+2,:] # up (note: it is assumed that the Jacobian was built for 3 components -- e,n,u; hence, the +2)
                d[bb] = udisp[crow_data] # up
                w[bb] = 1/usig[crow_data]**2 # up-sigma
            # Create Weighting Matrix and Weight the Design Matrix and Data Vector
            W = np.diag(w)
            A = G.copy()
            dat = d.copy()
            G = np.dot(W,G)
            d = np.dot(W,d)
        else: # three components (e,n,u)
            # Initialize updated Design Matrix
            G = np.empty((len(dm_idx)*3,len(load_cell_ids)))
            # Initialize datavector
            d = np.empty((len(data_idx)*3,))
            # Initialize weighting vector
            w = np.empty((len(data_idx)*3,))
            for bb in range(0,len(dm_idx)):
                crow_dm = dm_idx[bb]
                crow_data = data_idx[bb]
                G[(bb*3),:] = design_matrix[crow_dm,:] # east
                G[(bb*3)+1,:] = design_matrix[crow_dm+1,:] # north
                G[(bb*3)+2,:] = design_matrix[crow_dm+2,:] # up
                d[(bb*3)] = edisp[crow_data] - edisp[crow_data].mean() # east
                d[(bb*3)+1] = ndisp[crow_data] - ndisp[crow_data].mean() # north
                d[(bb*3)+2] = udisp[crow_data] # up
                w[(bb*3)] = 1/esig[crow_data]**2 # east-sigma
                w[(bb*3)+1] = 1/nsig[crow_data]**2 # north-sigma 
                w[(bb*3)+2] = 1/usig[crow_data]**2 # up-sigma
            # Create Weighting Matrix and Weight the Design Matrix and Data Vector
            W = np.diag(w)
            A = G.copy()
            dat = d.copy()
            G = np.dot(W,G)
            d = np.dot(W,d)
    else:
        if (uonly == True): # one component
            # Initialize updated Jacobian
            G = np.empty((len(dm_idx),len(load_cell_ids)))
            # Initialize datavector
            d = np.empty((len(data_idx),))
            for bb in range(0,len(dm_idx)):
                crow_dm = dm_idx[bb]
                crow_data = data_idx[bb]
                G[bb,:] = design_matrix[crow_dm+2,:] # up (note: it is assumed that the Jacobian was built for 3 components -- e,n,u; hence, the +2)
                d[bb] = udisp[crow_data] # up
                A = G.copy()
                dat = d.copy()
        else: # three components (e,n,u)
            # Initialize updated Design Matrix
            G = np.empty((len(dm_idx)*3,len(load_cell_ids)))
            # Initialize datavector
            d = np.empty((len(data_idx)*3,))
            for bb in range(0,len(dm_idx)):
                crow_dm = dm_idx[bb]
                crow_data = data_idx[bb]
                G[(bb*3),:] = design_matrix[crow_dm,:] # east
                G[(bb*3)+1,:] = design_matrix[crow_dm+1,:] # north
                G[(bb*3)+2,:] = design_matrix[crow_dm+2,:] # up
                d[(bb*3)] = edisp[crow_data] # east
                d[(bb*3)+1] = ndisp[crow_data] # north
                d[(bb*3)+2] = udisp[crow_data] # up
                A = G.copy()
                dat = d.copy()

    # Bias the Solution with Tikhonov Regularization
    if (tikhonov == 'zeroth'): # Zeroth-Order Tikhonov Regularization
        # Create identity matrix 
        identity = np.identity(len(load_cell_ids))
        # ** Option D: Solve by redefining the design matrix and data vector to include regularization (Aster, Borchers, & Thurber [2013], Eq. 4.5)
        C = np.concatenate((G,identity),axis=0)
        zero_vec = np.zeros(len(load_cell_ids))
        b = np.concatenate((d,zero_vec))
        # Solve Using Conjugate Gradient Method
        CTC = np.dot(C.T,C)
        CTb = np.dot(C.T,b)
        mvec,info = linalg.cg(CTC, CTb)

    elif (tikhonov == 'second' or tikhonov == 'zeroth_second'): # Second-Order Tikhonov Regularization (in two dimensions)             # Search for load cells that are surrounded on all four sides by other load cells (n,e,w,s)
        coordinates = np.empty((len(load_cell_lat),2))
        coordinates[:,0] = load_cell_lat
        coordinates[:,1] = load_cell_lon
        unq_lat = np.unique(load_cell_lat)
        unq_lon = np.unique(load_cell_lon)
        lat_sep = (max(unq_lat) - min(unq_lat)) / (len(unq_lat)-1)
        lon_sep = (max(unq_lon) - min(unq_lon)) / (len(unq_lon)-1)
        print(':: Note: Computing longitude separation value and neareset neighbors may fail if working across the Prime Meridian. Please use caution and verify results. [perform_inversion.py]')
        print(':: Note: For small or unusual cell separation distances, you may need to verify that cells are being found using numpy.where (tolerance of zero). [perform_inversion.py]')
        # Initialize L matrix (see Aster, Borchers, & Thurber (2013), Chapter 4, Exercise 4.3)
        L = []
        for cc in range(0,len(load_cell_ids)):
            ccell_id = load_cell_ids[cc]
            ccelllat = load_cell_lat[cc]
            ccelllon = load_cell_lon[cc]
            # Search for nearest neighbors (south, north, west, and east)
            nn_slat_minus = ccelllat - lat_sep  - 0.1*lat_sep #+ 0.1*lat_sep
            nn_slat_plus = ccelllat - lat_sep  + 0.1*lat_sep #- 0.1*lat_sep
            nn_nlat_minus = ccelllat + lat_sep - 0.1*lat_sep #+ 0.1*lat_sep
            nn_nlat_plus = ccelllat + lat_sep + 0.1*lat_sep #- 0.1*lat_sep
            nn_wlon_minus = ccelllon - lon_sep - 0.1*lon_sep #+ 0.1*lon_sep
            nn_wlon_plus = ccelllon - lon_sep + 0.1*lon_sep #- 0.1*lon_sep
            nn_elon_minus = ccelllon + lon_sep - 0.1*lat_sep #+ 0.1*lon_sep
            nn_elon_plus = ccelllon + lon_sep + 0.1*lat_sep #- 0.1*lon_sep
            nn_south = np.where((load_cell_lat >= nn_slat_minus) & (load_cell_lat <= nn_slat_plus) & (load_cell_lon == ccelllon)); nn_south = nn_south[0]
            nn_north = np.where((load_cell_lat >= nn_nlat_minus) & (load_cell_lat <= nn_nlat_plus) & (load_cell_lon == ccelllon)); nn_north = nn_north[0]
            nn_west = np.where((load_cell_lon >= nn_wlon_minus) & (load_cell_lon <= nn_wlon_plus) & (load_cell_lat == ccelllat)); nn_west = nn_west[0]
            nn_east = np.where((load_cell_lon >= nn_elon_minus) & (load_cell_lon <= nn_elon_plus) & (load_cell_lat == ccelllat)); nn_east = nn_east[0]
            ## Nine Point Regularization## 
            if (nine_point == True): # search for southwest, southeast, northwest, and northeast points (e.g., Rosser 1975 Comp. & Maths. with Appls.)
                # See https://www.sciencedirect.com/science/article/pii/0898122175900358
                nn_sw = np.where(np.isin(load_cell_lat, nn_south) & np.isin(load_cell_lon, nn_west))[0] 
                nn_se = np.where(np.isin(load_cell_lat, nn_south) & np.isin(load_cell_lon, nn_east))[0]
                nn_nw = np.where(np.isin(load_cell_lat, nn_north) & np.isin(load_cell_lon, nn_west))[0]
                nn_ne = np.where(np.isin(load_cell_lat, nn_north) & np.isin(load_cell_lon, nn_east))[0] 
                # Current cell has all eight neighbors; add to L matrix
                # See Aster, Borchers, & Thurber (2013), Exercise 4.3 (Section 4.9, Page 126) for an example of how to
                #   generate the L matrix for second-order Tikhonov regularization in two dimensions
                # See Rosser (1975) Comp. & Maths. with Appls., Vol 1, pp. 351-360 for 9-point stencil
                crow = np.zeros((len(load_cell_ids),))
                crow[cc] = -20
                if nn_south is None: # no south neighbor
                    crow[cc] += 4
                else:
                    crow[nn_south] = 4
                if nn_north is None: # no north neighbor
                    crow[cc] += 4
                else:
                    crow[nn_north] = 4
                if nn_west is None: # no west neighbor
                    crow[cc] += 4
                else:
                    crow[nn_west] = 4
                if nn_east is None: # no east neighbor
                    crow[cc] += 4
                else:
                    crow[nn_east] = 4
                if nn_sw is None: # no southwest neighbor
                    crow[cc] += 1
                else:
                    crow[nn_sw] = 1
                if nn_se is None: # no southeast neighbor
                    crow[cc] += 1
                else:
                    crow[nn_se] = 1
                if nn_nw is None: # no northwest neighbor
                    crow[cc] += 1
                else:
                    crow[nn_nw] = 1
                if nn_ne is None: # no northeast neighbor
                    crow[cc] += 1
                else:
                    crow[nn_ne] = 1
                L.append(crow)
            else:
                # Current cell has all four neighbors; add to L matrix
                # See Aster, Borchers, & Thurber (2013), Exercise 4.3 (Section 4.9, Page 126) for an example of how to 
                #   generate the L matrix for second-order Tikhonov regularization in two dimensions
                crow = np.zeros((len(load_cell_ids),))
                crow[cc] = -4
                if nn_south is None: # no south neighbor
                    crow[cc] += 1
                    print("no south neighbor")
                else: 
                    crow[nn_south] = 1
                if nn_north is None: # no north neighbor
                    crow[cc] += 1
                    #print("no north neighbor")
                else:
                    crow[nn_north] = 1
                if nn_west is None: # no west neighbor
                    crow[cc] += 1
                    #print("no west neighbor")
                else:
                    crow[nn_west] = 1
                if nn_east is None: # no east neighbor
                    crow[cc] += 1
                    #print("no east neighbor")
                else:
                    crow[nn_east] = 1
                L.append(crow)
                #print(crow[cc])
        L = np.asarray(L)

    # No Regularization Applied
    elif (tikhonov == 'none'):
        L = np.eye(load_cell_ids)
        C = G.copy()
        b = d.copy()

    else:
        sys.exit(':: Error: Invalid Tikhonov regularization code. [perform_inversion.py]')

    C = np.concatenate((G,(alpha*L)),axis=0)
    zero_vec_rows = L.shape[0]
    zero_vec = np.zeros(zero_vec_rows)
    b = np.concatenate((d,zero_vec))
    LTL = np.dot(L.T,L)
        
    # Option K: Applying Zeroth- and Second-Order Regularization. 
    # Requires a square matrix if we wish to use conjugate gradient
    identity = np.identity(len(load_cell_ids))
    O = np.concatenate((G,alpha*L,beta*identity),axis=0)
    # If we wished to use an a-priori model the first zero_vec would be replaced by our a-priori model
    # could set a-priori model to grace around edge and get rid of zeroth
    # order.
    e = np.concatenate((d,zero_vec,zero_vec))
    OTO = np.dot(O.T,O)
    OTe = np.dot(O.T,e)
    mvec,info = linalg.cg(OTO,OTe)
    print(min(mvec))
    print(max(mvec))
        
    # Residual Norm and Solution Seminorm
    res_norm = LA.norm(np.dot(G,mvec)-d)
    sol_norm = LA.norm(np.dot(L,mvec))
    print("Residual Norm:",res_norm)
    print("Solution Seminorm:",sol_norm)

    # Output Statistics
    # Model Covariance Matrix
    if (weighted == True):
        covm = np.linalg.inv(G.T @ G + alpha * (L.T @ L) + beta * np.eye(G.shape[1]))
    else:
        # Compute Hessian
        H_inv = np.linalg.inv(G.T @ G + alpha * (L.T @ L) + beta * np.eye(G.shape[1]))
        covm = H_inv @ G.T @ (usig[:, None]**2 * G) @ H_inv
    # Computing 95% confidence intervals
    confidence_intervals = 1.96 * np.sqrt(np.diag(covm))
 
    plt.scatter(load_cell_lon,load_cell_lat,c=mvec)
    plt.show()

    # Return Model Vector
    return mvec,confidence_intervals
