#!/usr/bin/env python
import sys,getopt
import ICPS2MRSClasses as hu
import DiscoClasses as di
import numpy as np
from niceout import *


DEBUG=0 #3 show pylab plots, 2 all, 1 some, 0 none
hu.setDEBUG(DEBUG)
hu.setNSIDE(16)
di.setDEBUG(DEBUG)
import ROOT

ROOT.gROOT.SetBatch()



def main(argv):
    opts, args = getopt.getopt(argv, "l:s:f:c:m:")
    for opt, arg in opts:
        if opt in "-l":
            luminosity=float(arg)
        if opt in "-s":
            series=int(arg)
        if opt in "-f":
            folder=str(arg)
        if opt in "-c":
            chunksize = int(arg)
        if opt in "-m":
            M = float(arg)
        #if opt in "-b":
            #basename=arg


    maps=[]
    hu.setDDir(folder)

    maps.append(hu.TwoMRSMap(hu.Make2MRSMap(0.00, 0.02, 16)))
    
    print PF(),"maps loaded."

    inj=hu.Injector(maps,sim=True, pthresh=np.power(10., -3.7), m=M)
    #DDisco=di.disco(inj, significance , power,basename,ddir)
    #DDisco.AnalyzeDiscoveryPot()
    hu.DoOneLuminosity(luminosity, inj, series, chunksize=chunksize)


#-----------------------------------------------------------------------------
if __name__ == "__main__":
    main(sys.argv[1:])

