#!/usr/bin/env python

# for Python 3 compatibility
from __future__ import print_function

import os, sys, glob
import numpy as np

sys.path.append('/home/rpernak/python/GIT_python_modules')
import utils

TOPDIR = '/rd47/scratch/RRTMGP/RRTMGP_Git/examples/rfmip-clear-sky'

class timingGPTL():
  def __init__(self, inDir, rrtmgp=False):
    """
    Read in timing.*.*.* files that were generated in either the 
    RRTMG or RRTMGP drivers, then calculate the moments (over the 
    number of iterations) for a given block size

    We assume a naming convention of timing.TRIAL.BLOCK_SIZE.MODEL 
    underneath inDir, where TRIAL is the trial number, BLOCK_SIZE is 
    the block size used in the computation, and MODEL is either 
    RRTMG or RRTMGP
    """

    utils.file_check(inDir)
    self.modelStr = 'RRTMGP' if rrtmgp else 'RRTMG'
    self.inFiles = sorted(glob.glob('%s/timing.*.*.%s' % \
      (inDir, self.modelStr)))

    self.allBlockSizes = \
      np.array([os.path.basename(inFile).split('.')[2] \
        for inFile in self.inFiles]).astype(int)
    self.blockSizes = np.unique(self.allBlockSizes)
  # end constructor

  def extractTime(self):
    """
    Extract wall clock time for all model calculations from timing 
    file
    """

    # making some assumptions about GPTL output here -- i.e., that it
    # is uniform over all trials and block sizes
    iRec = [40, 41] if self.modelStr == 'RRTMGP' else [40]
    wallClock = []
    for inFile in self.inFiles:
      dat = np.array(open(inFile).read().splitlines())[iRec]
      if self.modelStr == 'RRTMGP':
        iWallClock = 4
        wallClock.append(float(dat[0].split()[iWallClock]) + \
          float(dat[1].split()[iWallClock-1]))
      else:
        iWallClock = 4
        wallClock.append(float(dat[0].split()[iWallClock]))
      # endif G or GP
    # end loop over files

    self.times = np.array(wallClock)
    return self
  # end extractTime()

  def calcMoments(self):
    """
    Calculate means and standard deviations of the wallclock times 
    for model runs over all trials for a given block size
    """

    timingMean, timingSD = [], []
    for bSize in self.blockSizes:
      iBlock = np.where(self.allBlockSizes == bSize)[0]
      if iBlock.size == 0:
        timingMean.append(np.nan)
        timingSD.append(np.nan)
      else:
        timingMean.append(self.times[iBlock].mean())
        timingSD.append(self.times[iBlock].std(ddof=1))
      # endif mean/SD calc
    # end block size loop
    self.means = np.array(timingMean)
    self.sigmas = np.array(timingSD)

    return self
  # end calcMoments()
# end timingGPTL()

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(\
    description='Read in GPTL standard output from RRTMG and ' + \
    'RRTMGP drivers and calculate statistics for each block size ' + \
    'used. Do this for RRTMGP and RRTMG and LW and SW ' + \
    '(separately, so that there are 4 sets of statistics). This ' + \
    'is a very ad-hoc script and should be used in concert with ' + \
    'RFMIP_wrapper.py.')
  parser.add_argument('-gl', '--rrtmg_lw_dir', type=str, \
    default='%s/LW' % TOPDIR, \
    help='Directory with RRTMG LW results.')
  parser.add_argument('-gs', '--rrtmg_sw_dir', type=str, \
    default='%s/SW' % TOPDIR, \
    help='Directory with RRTMG SW results.')
  parser.add_argument('-gpl', '--rrtmgp_lw_dir', type=str, \
    default='%s/rrtmgp_LW' % TOPDIR, \
    help='Directory with RRTMGP LW results.')
  parser.add_argument('-gps', '--rrtmgp_sw_dir', type=str, \
    default='%s/rrtmgp_SW' % TOPDIR, \
    help='Directory with RRTMGP SW results.')
  args = parser.parse_args()

  from pandas import DataFrame as DF

  allDir = [args.rrtmg_lw_dir, args.rrtmg_sw_dir, \
    args.rrtmgp_lw_dir, args.rrtmgp_sw_dir]
  gpBool = [False, False, True, True]
  swBool = [False, True, False, True]
  meanDict, sigmaDict = {}, {}

  for runDir, gpBoo, swBoo in zip(allDir, gpBool, swBool):
    runStr = '%s_%s' % \
      ('RRTMGP' if gpBoo else 'RRTMG', 'SW' if swBoo else 'LW')
    runObj = timingGPTL(runDir, rrtmgp=gpBoo)
    runObj.extractTime()
    runObj.calcMoments()
    meanDict[runStr] = runObj.means
    sigmaDict[runStr] = runObj.sigmas
  # end run loop

  # making the assumption that all runs had the same block size arrays
  meanDict['block_size'] = np.array(runObj.blockSizes)
  sigmaDict['block_size'] = np.array(runObj.blockSizes)
  meanCSV = 'timing_RRTMG_RRTMGP_mean.csv'
  sdCSV = 'timing_RRTMG_RRTMGP_sigma.csv'

  DF.from_dict(meanDict).to_csv(meanCSV, index=False)
  DF.from_dict(sigmaDict).to_csv(sdCSV, index=False)
# end main()

