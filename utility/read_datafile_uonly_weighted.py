# *********************************************************************
# FUNCTION TO READ IN OBSERVATIONS of SURFACE DISPLACEMENT
#
# Copyright (c) 2014-2019: HILARY R. MARTENS, LUIS RIVERA, MARK SIMONS
#
# This file is part of LoadDef.
#
#    LoadDef is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    LoadDef is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with LoadDef.  If not, see <https://www.gnu.org/licenses/>.
#
# *********************************************************************

import numpy as np
import sys

def main(filename):

    # Read in station names
    sta = np.loadtxt(filename,skiprows=1,usecols=(0,),unpack=True,dtype='U')

    # Read in data
    lat,lon,udisp,usig = np.loadtxt(filename,skiprows=1,usecols=(1,2,5,8),unpack=True)

    # Return Parameters
    return sta,lat,lon,udisp,usig
