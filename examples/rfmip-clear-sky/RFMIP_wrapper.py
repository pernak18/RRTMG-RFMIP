#!/usr/bin/env python

# for Python 3 compatibility
from __future__ import print_function

import os, sys
import subprocess as sub

RRTMGP_ROOT = os.getenv('RRTMGP_ROOT')
if RRTMGP_ROOT is None:
  print('Please set RRTMGP_ROOT environment variable')
  sys.exit(1)
# endif RRTMGP_ROOT

sys.path.append('/home/rpernak/python/GIT_python_modules')
import utils

class driverRFMIP():
  def __init__(self, inVars, trial=1, blockSize=100):
    """
    Run a driver got a given block size and rename its GPTL timing 
    output. Block size by default is the number of RFMIP profiles
    """
    utils.file_check(inVars['exe'])

    self.driverExe = str(inVars['exe'])
    self.runDir = os.path.dirname(self.driverExe)
    self.iter = int(trial)
    self.blockSize = int(blockSize)
    self.topDir = os.getcwd()

    # making some assumptions here regarding driver name
    self.model = 'RRTMGP' if 'rrtmgp' in self.driverExe else 'RRTMG'
  # end constructor()

  def runDriver(self):
    os.chdir(self.runDir)
    sub.call([self.driverExe, str(self.blockSize)])

    # the driver should produce a "timing.block_size" file with GPTL
    # output that we need to rename so we don't overwrite in other 
    # iterations
    os.rename('timing.%d' % self.blockSize, \
      'timing.%d.%d.%s' % (self.iter, self.blockSize, self.model))

    # cd back into original directory
    os.chdir(self.topDir)
  # end runDriver()
# end driverRFMIP

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(\
    description='Run a driver (RRTMG or RRTMGP, e.g.) for a ' + \
    'given number of iterations and block sizes. It is assumed ' + \
    'the driver implements GPTL and writes timing information to ' + \
    'timing.block_size files.')
  parser.add_argument('-e', '--exe', type=str, \
    default='%s/examples/rfmip-clear-sky/SW/rrtmg_rfmip_sw' % \
    RRTMGP_ROOT, help='Full path to driver executable.')
  parser.add_argument('-b', '--block_size', type=int, nargs='+', \
    default=[4, 8, 100, 1800], \
    help='Size of memory blocks for computation.')
  parser.add_argument('-n', '--niter', type=int, default=10, \
    help='Number of iterations over which to run driver.')
  args = parser.parse_args()

  for iTrial in range(1, args.niter+1):
    for iBlock in args.block_size:
      print('Trial, Block Size = %d, %d' % (iTrial, iBlock) )
      oRFMIP = driverRFMIP(vars(args), trial=iTrial, blockSize=iBlock)
      oRFMIP.runDriver()
    # end block loop
  # end iter loop

# end main()

