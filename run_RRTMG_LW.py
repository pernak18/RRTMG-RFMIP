#!/usr/bin/env python

import os, sys, glob
import subprocess as sub

sys.path.append('common')
import utils

class driverRRTMGLW:
  def __init__(self, inFile, exe, outDir):
    """
    `run_RRTMG_LW.py -h`
    """

    self.inFile = str(inFile)
    self.exe = str(exe)
    paths = [self.inFile, self.exe]
    for path in paths: utils.file_check(path)

    self.outDir = str(outDir)
    if not os.path.exists(self.outDir): os.makedirs(self.outDir)
  # end constructor

  def runRRTMG(self):
    """
    Run RRTMG LW on RFMIP profiles
    """

    iRRTMG = 'INPUT_RRTM'
    if os.path.exists(iRRTMG): os.remove(iRRTMG)
    os.symlink(inFile, iRRTMG)
    sub.call([self.exe])
    os.remove(iRRTMG)
  # end runRRTGM()

  def rename(self):
    """
    Rename OUTPUT_RRTM file from RRTMG to something more unique
    """

    oRRTMG = 'OUTPUT_RRTM'
    base = os.path.basename(self.inFile).replace('INPUT_RRTM', oRRTMG)
    outFile = '{}/{}'.format(self.outDir, base)
    os.rename(oRRTMG, outFile)
    print('Finished {}'.format(outFile))
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
  parser.add_argument('--outdir', '-o', type=str, \
    default='OUTPUT_RRTMG/PD', \
    help='Directory into which OUTPUT_RRTM files are moved.')
  args = parser.parse_args()

  inDir = args.input_dir; utils.file_check(inDir)
  inFiles = sorted(glob.glob('{}/INPUT_RRTM_*'.format(inDir)))
  if len(inFiles) == 0:
    print('Found no RRTMG_LW inputs in {}'.format())
    sys.exit(1)
  # endif nFiles

  for inFile in inFiles:
    rObj = driverRRTMGLW(inFile, args.exe_path, args.outdir)
    rObj.runRRTMG()
    rObj.rename()
  # end inFile loop
# endif main()
