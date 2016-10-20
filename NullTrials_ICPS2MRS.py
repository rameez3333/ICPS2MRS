#!/usr/bin/env python
import sys,getopt
import ICPS2MRSClasses as hu
import DiscoClasses as di
from ICPS2MRSClasses import DoOneScramble
from niceout import *


DEBUG=0 #3 show pylab plots, 2 all, 1 some, 0 none
hu.setDEBUG(DEBUG)
hu.setNSIDE(16)
di.setDEBUG(DEBUG)
import ROOT

ROOT.gROOT.SetBatch()



def main(argv):
    opts, args = getopt.getopt(argv, "d:p:s:b:")
    for opt, arg in opts:
        if opt in "-d":
            ddir=arg
        if opt in "-p":
            power=float(arg)
        if opt in "-s":
            significance=float(arg)
        if opt in "-b":
            basename=arg


    maps=[]
    hu.setDDir(ddir)

    maps.append(hu.TwoMRSMap(hu.Make2MRSMap(0.02, 16)))
    
    print PF(),"maps loaded."

    inj=hu.Injector(maps,sim=True)
    inj.SetSeed()
    fout = open(basename+'output.txt', "w")
    for i in range(0, 1000):
        TS,nsfit,nserr,nsig,gr=DoOneScramble(None,(0,inj,0,30))
        fout.write(str(TS)+' '+str(nsfit)+' '+str(nserr)+' '+str(nsig)+'\n')
    fout.close()

#-----------------------------------------------------------------------------
if __name__ == "__main__":
    main(sys.argv[1:])

