# *********************************************************************
# FUNCTION TO READ IN A PRIORI MODEL
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
import sys 

def main(filename):

    # Read in data
    lat,lon,height = np.loadtxt(filename,usecols=(0,1,2),unpack=True)

    # Return Parameters
    return lat,lon,height
