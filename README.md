# LoadSolv

A Python based software used to invert GNSS observations of surface displacement to infer changes in water storage.

## History
This software's original strucutre and scripts were developed at the University of Montana by Dr. Hilary R. Martens in 2021 for use with the elastic deformation modeling software LoadDef. Refinement and alterations of this version of the software were made by Dr. Matthew J Swarr over the course of his PhD at the University of Montana. The software has been archived in various open-source repositories associated with publications over the years, but we have it here to serve as a 'living' version that can be worked on as times goes on.

## General Structure
To use this code you need the following ingredients:
        - Design Matrix (Jacobian): Dataset describing the predicted horizontal and vertical displacement response of each GNSS station within a study area (e.g., Western US) to a unit load in each 
          pixel of your model domain (e.g., 0.25 degree pixel grid in the Western US). The design matrix can be computed using common semi-analytical models that simualted the solid Earth's elastic 
          deformation response (e.g., LoadDef, SPOTL). However, design matrices may also be computed using fully-numerical models of Earth deformation (e.g., PyLith, Citcom, Abaqus).
        - GNSS Displacement Data: A file (or series of files) that contain the displacement of GNSS stations on a given epoch (e.g., January 2006). This file should include the following:
          Station Name, Latitude, Longitude, E Displacement [mm], N Displacement [mm], V Displacement [mm], E Sigma [mm], N Sigma [mm], V Sigma [mm]

Outputs:
        - Gridded estimates of water storage and associated formal uncertainties for each epoch in which data files were provided.

## Requirements
This code base relies upon common Python packages such as: NumPy, SciPy, Pandas, datetime, Netcdf4, mpi4py, sys, os, matplotlib

## Publications
This software has been used to produce results presented in several peer-reviewed publications. The following list serves an extensive (but perhaps not up to date) list of publications:

- Swarr, M. J., Argus, D. F., Martens, H. R., Hoylman, Z. H., Young, Z. M., Borsa, A. A., et al. (2026). Extreme winter precipitation drives recharge of mountain groundwater in the Sierra Nevada and cascades of the western United States. Water Resources Research, 62, e2025WR040910. https://doi.org/10.1029/2025WR040910
- Martens, H. R., Lau, N., Swarr, M. J., Argus, D. F., Cao, Q., Young, Z. M., et al. (2024). GNSS geodesy quantifies water-storage gains and drought improvements in California spurred by atmospheric rivers. Geophysical Research Letters, 51, e2023GL107721. https://doi.org/10.1029/2023GL107721
- Swarr, M. J., Martens, H. R., & Fu, Y. (2024). Sensitivity of GNSS-derived estimates of terrestrial water storage to assumed Earth structure. Journal of Geophysical Research: Solid Earth, 129, e2023JB027938. https://doi.org/10.1029/2023JB027938
- White, A. M., Lajoie, L. J., Knappe, E., Martens, H. R., Swarr, M. J., Khatiwada, A., et al. (2023). High-density integrated GNSS and hydrologic monitoring network for short-scale hydrogeodesy in high mountain watersheds. Earth and Space Science, 10, e2022EA002678. https://doi.org/10.1029/2022EA002678
