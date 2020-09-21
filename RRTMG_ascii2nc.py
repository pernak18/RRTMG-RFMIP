#!/usr/bin/env python

from __future__ import print_function

import os, sys
import numpy as np

# RC GitLab repo
# git clone git@lex-gitlab.aer.com:RC/common_modules.git
sys.path.append('common')
import utils

# ESGF-compliant RFMIP results netCDF template so we can extract
# profile ('site') locations and weights
ESGF = '/Users/rpernak/Work/RC/RFMIP_RRTMGP/PI_O3_trop_PD_strat/' + \
  'rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc'

def findProfiles(inDir):
  """
  Extract the RRTMG flux files for given spectral domain

  inDir -- string, directory with RRTMG output files for a given
    Buehler experiment
  """

  import glob

  inFiles = sorted(glob.glob('{}/OUTPUT_RRTM_*'.format(inDir)))
  if len(inFiles) == 0: sys.exit('No profiles found, returning')

  return inFiles
# end findProfiles()

class rrtmg():
  def __init__(self, inFile, upwelling=False):
    """
    Extract the RRTMG broadband flux files for a given RFMIP profile
    and Buehler experiment
    """

    self.inFile = str(inFile)
    self.doUpwelling = bool(upwelling)

    # band wavenumber limits used in calculations
    # RRTMG defines these limits
    self.bandWN = np.array([10, 3250])
  # end constructor

  def readASCII(self):
    """
    Read a single RRTMG ASCII flux file and return the model output
    in a dictionary to be used in output netCDF creation (writeNC())

    Input

    Output
      None, but outDict is saved as the `profile` attribute for the
        class and is a dictionary with the following fields:

        level_pressure: pressure at layer boundaries (nLevel array)
        flux: either upwelling or downwelling [W/m2], depending on
          the user specification; nLevel elements
    """

    # for the Buehler exeriments, we're only modeling broadband
    # fluxes, so the load of the data is pretty simple -- skip the
    # header, ignore the footer, keep everything from 1 band
    level, pLev, upFlux, downFlux, netFlux, hr = \
      np.loadtxt(self.inFile, unpack=True, skiprows=3, max_rows=61)

    outDict = {}

    # mbar to Pa conversion, and staying with TOA to surface convention
    outDict['level_pressures'] = np.array(pLev) * 100

    outFlux = upFlux if self.doUpwelling else downFlux
    outDict['flux'] = np.array(outFlux)
    self.profile = dict(outDict)
  # end readASCII()
# end rrtmg()

def writeNC(experiments, fluxes, pressures, esgfNC=ESGF, \
  outNC='rrtmg_buehler.nc', upwelling=False):
  """
  Write broadband RRTMG fluxes for all experiments into a single
  netCDF file that is ESGF-compliant.

  Inputs
    experiments -- string list; names of experiments (present day and
      leave-one-out)
    fluxes -- float array, nExp x nProf x nLev RRTMG broadband fluxes
    pressures -- float array of profile pressures

  Keywords
    esgfNC -- string, path to netCDF that contains ESGF-compliant
      RFMIP profile information (fluxes, locations, weights)
    outNC -- string, name of output netCDF file to which model results
      are written
    upwelling -- boolean, upwelling instead of downwelling fluxes are
      provided
  """

  import netCDF4 as nc

  utils.file_check(ESGF)
  with nc.Dataset(ESGF, 'r') as ncObj:
    # get profile location and weights for each "site" (profile)
    weights = np.array(ncObj.variables['profile_weight'])
    lats = np.array(ncObj.variables['lat'])
    lons = np.array(ncObj.variables['lon'])
    times = np.array(ncObj.variables['time'])
  # endwith

  # output dimension information
  nExp, nProf, nLev = fluxes.shape
  dimNames = ['expt', 'site', 'level']
  dimVals = [nExp, nProf, nLev]

  # ESGF abbreviations for longwave up and down
  fluxStr = 'rlu' if upwelling else 'rld'

  # write variables and dimensions to netCDF, trying to conform to
  # ESGF convention (but probably failing)
  with nc.Dataset(outNC, 'w') as ncObj:
    ncObj.set_fill_on()
    ncObj.description = 'RRTMG fluxes for Buehler RFMIP experiments'
    ncObj.source = 'Atmospheric and Environmental Research (AER)'

    for name, val in zip(dimNames, dimVals):
      ncObj.createDimension(name, val)

    ncVar = ncObj.createVariable('lat', float, ('site'), \
      fill_value=np.nan)
    ncVar[:] = lats
    ncVar.units = 'degrees north [-90, 90]'
    ncVar.description = 'Geographic latitude of profile'
    ncVar.valid_range = (-90, 90)

    ncVar = ncObj.createVariable('lon', float, ('site'), \
      fill_value=np.nan)
    ncVar[:] = lons
    ncVar.units = 'degrees east [0, 360]'
    ncVar.description = 'Geographic longitude of profile'
    ncVar.valid_range = (0, 360)

    ncVar = ncObj.createVariable('plev', float, ('site', 'level'), \
      fill_value=np.nan)
    ncVar[:] = pressures
    ncVar.units = 'Pa'
    ncVar.description = 'Pressure at layer edge'
    ncVar.valid_range = (0, 1e6)

    ncVar = ncObj.createVariable('profile_weight', float, ('site'), \
      fill_value=np.nan)
    ncVar[:] = weights
    ncVar.units = ''
    ncVar.description = 'Profile weight to recover global mean'
    ncVar.valid_range = (0, 1)

    ncVar = ncObj.createVariable(fluxStr, float, \
      ('expt', 'site', 'level'), fill_value=-1000)
    ncVar[:] = fluxes
    ncVar.units = 'W m-2'
    ncVar.description = ''
    ncVar.valid_range = (-1000, 1e6)

    ncVar = ncObj.createVariable('time', float, ('site'), \
      fill_value=np.nan)
    ncVar[:] = times
    ncVar.units = 'Gregorian'
    ncVar.description = 'Days since 2014-01-01 00:00:00'
    ncVar.valid_range = (-1000, 1e6)
  # endwith

  print('Wrote {}'.format(outNC))
# end writeNC

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(\
    description='Convert ASCII RRTMG output to netCDF format ' + \
    'that is ESGF-compliant. Only done for LW Buehler experiments ' + \
    '(https://github.com/pernak18/RFMIP_Buehler_Experiments). ' + \
    'This only works for RRTMG broadband calculations.')
  parser.add_argument('--rrtmg_results_dir', '-r', type=str, \
    default='OUTPUT_RRTMG', \
    help='Directory with RRTMG LW results, with separate ' + \
    'directories for each experiment.')
  parser.add_argument('--upwelling', '-u', action='store_true', \
    help='Extract upwelling fluxes instead of the downwelling ones.')
  parser.add_argument('--output_nc', '-o', type=str, \
    default='rld_Efx_RRTMG_rad_irf_r1i1p1f1_gn.nc', \
    help='Name of output netCDF that contains all downwelling or ' + \
    'upwelling fluxes for all profiles and experiments.')
  args = parser.parse_args()

  resultsRRTMG = args.rrtmg_results_dir
  utils.file_check(resultsRRTMG)

  # list to store flux arrays for each experiment
  # will be converted to nExp x nProfile x nLevel array
  fluxRRTMG = []

  # sorts by time created..i think
  experiments = os.listdir(resultsRRTMG)
  expNames = []
  for exp in experiments:
    fullPath = '{}/{}'.format(resultsRRTMG, exp)
    if not os.path.isdir(fullPath): continue
    #print(exp)
    expNames.append(exp)
    ascii = findProfiles(fullPath)

    # profile output loop -- 1 object per profile
    flux, pLev = [], []
    for txt in ascii:
      rObj = rrtmg(txt, upwelling=args.upwelling)
      rObj.readASCII()
      flux.append(rObj.profile['flux'])
      pLev.append(rObj.profile['level_pressures'])
    # end loop over text files

    fluxRRTMG.append(flux)
  # end experiment loop

  # pressures do not change with experiment
  pressures = np.array(pLev)
  fluxRRTMG = np.array(fluxRRTMG)

  writeNC(expNames, fluxRRTMG, pressures,
    outNC=args.output_nc, upwelling=args.upwelling)
# end main()
