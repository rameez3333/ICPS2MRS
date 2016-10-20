#!/usr/bin/env python
import sys,getopt
import ICPS2MRSClasses as hu
import DiscoClasses as di

from niceout import *


DEBUG=3 #3 show pylab plots, 2 all, 1 some, 0 none
hu.setDEBUG(DEBUG)
hu.setNSIDE(16)
di.setDEBUG(DEBUG)
import ROOT

ROOT.gROOT.SetBatch()



def main(argv):
    #opts, args = getopt.getopt(argv, "d:p:s:b:")
    #for opt, arg in opts:
        #if opt in "-d":
            #ddir=arg
        #if opt in "-p":
            #power=float(arg)
        #if opt in "-s":
            #significance=float(arg)
        #if opt in "-b":
            #basename=arg


    maps=[]
    hu.setDDir("Try")

    maps.append(hu.TwoMRSMap(hu.Make2MRSMap(0.00, 0.02, 16)))
    
    print PF(),"maps loaded."

    inj=hu.Injector(maps,sim=True)
    #DDisco=di.disco(inj, significance , power,basename,ddir)
    #DDisco.AnalyzeDiscoveryPot()
    hu.DoOneLuminosity(1e41, inj)


#-----------------------------------------------------------------------------
if __name__ == "__main__":
    main(sys.argv[1:])

