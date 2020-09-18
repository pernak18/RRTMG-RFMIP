#!/usr/bin/env python

import os, sys
import numpy as np
import netCDF4 as nc

class convertT5:
  def __init__(self, inDict, iProf):
    """
    `TAPE5_to_INPUT_RRTM.py -h`
    """

    # chemical names for XS species, their indices in a given record
    # CF4 is another acceptable cross section for RRTMG, but initially
    # we are treating it as a line-by-line gas in LBLRTM
    self.xsLBL = ['CCL3F', 'CCL2F2', 'CHClF2', 'CHF3', 'CH2F2', \
      'CHF2CF3', 'CFH2CF3', 'CF3CH3', 'CCL4']
    self.alias = ['CFC11', 'CFC12', 'CFC22', 'HFC23', 'HFC32', \
      'HFC-125', 'HFC-134a', 'HFC-143a', 'CCL4']

    self.inT5 = inDict['tape5']
    if not os.path.exists(self.inT5):
      sys.exit('Cannot find {}'.format(self.inT5))

    self.ncFile = inDict['nc_file']
    self.iProf = iProf
    self.dropXS = inDict['drop']

    with nc.Dataset(self.ncFile, 'r') as ncObj: self.surfaceT = \
      np.array(ncObj.variables['surface_temperature'])[0, iProf]

    self.outFile = \
      '{:s}_profile{:03d}'.format(inDict['outfile_prefix'], iProf+1)
  # end constructor

  def omit(self):
    """
    Omit a single XS species from list of LBLRTM XS species (chemical)
    names and their associated aliases
    """

    try:
      # where does input species match list of RRTMG aliases?
      iXS = self.alias.index(self.dropXS)

      # in-place removal of RRTMG alias and LBL name
      self.xsLBL.pop(iXS)
      self.alias.pop(iXS)
      print('{} dropped from RRTMG alias list'.format(self.dropXS))
    except:
      print('Could not find {} in allowable names'.format(self.dropXS))
      print('  --drop should be in [{}]'.format(', '.join(self.alias)))
      raise
    # end exception
  # endif omit

  def inputRRTMG(self):
    """
    Read in TAPE5 specifications and write the necessary data to an
    input file for RRTMG
    """

    inDat = open(self.inT5).read().splitlines()

    self.outDat = []
    xsNames = []

    # indices in records 2.2.4 and 2.2.5 (see lblrtm_instructions.html)
    # of species to keep
    idx224, idx225 = [], []
    for iLine, line in enumerate(inDat):
      # do not need any of the RADSUM blocks, which are at the end of
      # the TAPE5
      if iLine >= 608 and iLine < 650: continue

      # keep LBL records 1.1, 1.2, 2.1 (number of molecules in 2.1
      # goes from )
      if iLine == 0: self.outDat.append(line)

      if iLine == 1:
        # some extra options for RRTMG
        line = ' HI=1 F4=1 CN=1 AE 0 EM=0 SC=0 FI=0 PL=0 TS=0 ' + \
          'AM=0 MG=1 LA=0    0 XS=1   00   00    3    0   00'
        self.outDat.append(line)

        # RRTMG record 1.4 is surface temperature and emissivity
        # this is not a record in the LBL TAPE5, so we have to add it
        # emissivity is just 1 since we're assuming blackbody
        line = '{:10.3E}'.format(self.surfaceT)
        self.outDat.append(line)
      # endif line 1

      if iLine == 3:
        self.outDat.append(line.replace('42','{:2d}'.format(7)))

      # profile information for line-by-line molecules
      # layer P, T, Z
      if (iLine < 424) and ((iLine-4) % 7 == 0):
        self.outDat.append(line)

      # "Big 7" molecules + broadening
      if (iLine < 424) and ((iLine-5) % 7 == 0):
        self.outDat.append(line)

      # number of XS species
      if iLine == 424: self.outDat.append(line)

      # species names: only keep the 10 in RRTMG v4.86
      if iLine == 425:
        split = line.split()
        for iXS, xs in enumerate(split):
          if xs in self.xsLBL:
            # keep species name and corresponding index for usage
            # with records 2.2.4 and 2.2.5
            xsNames.append(self.alias[self.xsLBL.index(xs)])
            if iLine == 425: idx224.append(iXS)
          # endif XS
        # end split loop
      # endif 425

      if iLine == 426:
        split = line.split()
        for iXS, xs in enumerate(split):
          if xs in self.xsLBL:
            # keep species name and corresponding index for usage
            # with records 2.2.4 and 2.2.5
            xsNames.append(self.alias[self.xsLBL.index(xs)])
            if iLine == 426: idx225.append(iXS)
          # endif XS
        # end split loop

        # write to output file the first 7 XS names
        # record 2.2.1
        record221 = ['{:10s}'.format(name) for name in xsNames[:7]]
        record221 = ''.join(record221)
        self.outDat.append(record221)

        # write the species names (RRTMG allows at most 10 XS
        # species, so record 2.2.1 spans 2 lines at most)
        record221 = ['{:10s}'.format(name) for name in xsNames[7:]]
        record221 = ''.join(record221)
        self.outDat.append(record221)
      # endif 426

      # profile information for cross section species
      if iLine == 427: self.outDat.append(line)

      # layer P, T, Z
      if (iLine > 427) and ((iLine-428) % 3 == 0):
        self.outDat.append(line)

      # species concentrations
      if (iLine > 427) and ((iLine-428) % 3 == 1):
        xsConc = []
        split = line.split()
        broad = split[7]

        for i224 in idx224:
          xsConc.append(split[i224])
          nXS = len(xsConc)
          if nXS == 7:
            xsConc.append(broad)
            strConc = ['{:15s}' for xs in xsConc]
            self.outDat.append(''.join(strConc))
            xsConc = []
          # endif nXS
        # end i224 loop
      # endif conc line 1

      if (iLine > 427) and ((iLine-428) % 3 == 2):
        split = line.split()
        for i225 in idx225:
          xsConc.append(split[i225])
          nXS = len(xsConc)
          if nXS == 7:
            xsConc.append(broad)
            strConc = ['{:15s}'.format(xs) for xs in xsConc]
            self.outDat.append(''.join(strConc))
            xsConc = []
          # endif nXS
        # end i225 loop

        strConc = ['{:15s}'.format(xs) for xs in xsConc]
        self.outDat.append(''.join(strConc))
      # endif conc line 2
    # end line loop
  # end inputRRTMG()

  def write(self):
    """
    Write the file generated with inputRRTMG()
    """

    # write RRTM input file
    with open(self.outFile, 'w') as outFP:
      for line in self.outDat: outFP.write('{:s}\n'.format(line))

    print('Wrote {}'.format(self.outFile))
  # end write()
# end convertT5

if __name__ == '__main__':
  import argparse, glob

  parser = argparse.ArgumentParser(\
    description='Distill an LBLRTM TAPE5 down to what is needed ' + \
    'for RRTMG input. The two inputs are very similar, but there ' + \
    'is less information needed for RRTMG. Very ad-hoc, but ' + \
    'flexibility might not be necessary.')
  parser.add_argument('tape5_dir', type=str, \
    help='Directory with LBLRTM TAPE5 files for ' + \
    'RFMIP profiles, generated by pyscripts/make_TAPE5.py in ' + \
    'https://github.com/pernak18/RFMIP_Buehler_Experiments')
  parser.add_argument('--nc_file', '-n', type=str, \
    default='input_profiles_1-2.nc', \
    help='RFMIP specifications netCDF.')
  parser.add_argument('--outfile_prefix', '-o', type=str, \
    default='INPUT_RRTM_PD', help='Prefix for output file, which ' + \
    'will have the profile number appended to it. All outupt ' + \
    'files are written to the working directory.')
  parser.add_argument('--drop', '-d', type=str, \
    help='RRTMG XS species name to omit. All species are kept ' + \
    'if this is not set. Case-sensitive. See ALIAS global variable.')
  args = parser.parse_args()

  dirT5 = args.tape5_dir
  if not os.path.exists(dirT5):
    print('Could not find {}, returning'.format(dirT5))
    sys.exit(1)
  # endif dirT5

  allT5 = sorted(glob.glob('{}/TAPE5*'.format(dirT5)))
  if len(allT5) == 0:
    print('Found no TAPE5s in {}'.format(dirT5))
    sys.exit(1)
  # endif nT5

  inArgs = vars(args)
  for iProf, inT5 in enumerate(allT5):
    # we're making some assumptions that the TAPE5s are named in some
    # kind of sequential order
    inArgs['tape5'] = inT5
    cObj = convertT5(inArgs, iProf=iProf)
    if cObj.dropXS is not None: cObj.omit()
    cObj.inputRRTMG()
    cObj.write()
  # end profile loop

# endif main()
