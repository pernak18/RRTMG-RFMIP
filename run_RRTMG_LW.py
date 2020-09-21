#!/usr/bin/env python

import os, sys, glob
import subprocess as sub

sys.path.append('common')
import utils

class driverRRTMGLW:
  def __init__(self, inFile, exe):
    """
    `run_RRTMG_LW.py -h`
    """

    self.inFile = str(inFile)
    self.exe = str(exe)
    for path in paths: utils.file_check(path)
  # end constructor

  def runRRTMG(self):
    """
    Run RRTMG LW on RFMIP profiles
    """

    os.symlink(inFile, 'INPUT_RRTM')
    sub.call([self.exe])
  # end runRRTGM()

  def rename(self):
    """
    Rename OUTPUT_RRTM file from RRTMG to something more unique
    """

    oRRTMG = 'OUTPUT_RRTM'
    os.rename(oRRTMG, self.inFile.replace('INPUT_RRTM', oRRTMG))
  # end rename
# end driverRRTMGLW

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(\
    description='Run RRTMG_LW executable (that allows for ' + \
    'additional halocarbons, so probably v4.86 or >5.0) on ' + \
    'inputs generated with TAPE5_to_INPUT_RRTM.py')
  parser.add_argument('input_dir', type=str,
    help='Directory with INPUT_RRTM files for all RFMIP ' + \
    'profiles. This directory is likely for an entire experiment.')
  parser.add_argument('--exe_path', '-x', type=str, \
    default='rrtmg_lw', help='Path to RRTMG_LW executable.')
  args = parser.parse_args()

  inDir = args.input_dir; utils.file_check(inDir)
  inFiles = glob.glob('{}/INPUT_RRTM_*'.format(inDir))
  if len(inFiles) == 0:
    print('Found no RRTMG_LW inputs in {}'.format())
    sys.exit(1)
  # endif nFiles

  for inFile in inFiles:
    rObj = driverRRTMGLW(inFile, args.exe_path)
    rObj.getInputs()
    rObj.runRRTMG()
    rObj.rename()
  # end inFile loop
# endif main()
