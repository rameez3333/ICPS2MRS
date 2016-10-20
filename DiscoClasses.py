#!/usr/bin/env python
import os,copy
from scipy.stats import chisqprob,poisson

from niceout import *
from ICPS2MRSClasses import DoOneScramble

import ROOT
TH1D       =ROOT.TH1D
TGraph     =ROOT.TGraph
TCanvas    =ROOT.TCanvas

import numpy  as np
log10  =np.log10
fabs   =np.fabs
sqrt   =np.sqrt

from array import array as arr

def setDEBUG(debug):
    global DEBUG
    DEBUG=debug

class disco:
    def __init__(self,inj,signif,power,outbasename,ddir,chi2=1):
        self.inj=inj
        self.signif=signif
        self.power=power
        self.outbasename=outbasename
        self.DDIR=ddir
        self.chi2=chi2
        self.CollectedTrials={} #organized by ns injected, contains TS,ns fitted,prob from chi-squared
        self.CollectedTrialsFile=self.DDIR+"/keeptrials.dat"
        if not os.path.isfile(self.CollectedTrialsFile):
            f=open(self.CollectedTrialsFile,"w")
            f.write("ns_inj\tnsig\tTS\tns_fit\tprob\n")
            f.close()
        if not os.path.isfile(self.DDIR+"/rezults_"+self.outbasename+".dat"):
            f=open(self.DDIR+"/rezults_"+self.outbasename+".dat","w")
            f.write("grpah val.\tfit val.\t+-fit err.\t+-uncertainty\n")
            f.close()
        self.maxns=0
        for CR in ["IsotropicHS"]:
            if self.inj.Det[CR]:
                self.maxns+=len(self.inj.Det[CR].EvList)
        print PF(),"self.maxns",self.maxns


    def AnalyzeDiscoveryPot(self):
        self.ReUseStoredSim(self.CollectedTrialsFile)

        #FIRST PASS: Do trials until first detection
        self.FindStartPoint()
        last_loop=self.FindLastLoop()

        mew_mean,mew_mean_err=self.GetNextMean(iternum=last_loop)
        print PF(),"mew_mean", mew_mean," +- ",mew_mean_err
        self.DrawCollectedTrials(mew_mean,max_ns=self.maxns,cname=self.DDIR+"/c_"+self.outbasename+"_"+str(last_loop)+".root")

        for i in range(1000):
            last_loop+=1
            self.InjectWithPoissonMean(nInjections=20,mean=mew_mean)
            mew_mean,mew_mean_err=self.GetNextMean(step=1./(i+1),iternum=last_loop)
            print PF(),"mew_mean",mew_mean," +- ",mew_mean_err
            self.DrawCollectedTrials(mew_mean,max_ns=self.maxns,cname=self.DDIR+"/c_"+self.outbasename+"_"+str(i)+".root")

    def FindLastLoop(self):
        for iternum in range(1000):
            if not os.path.isfile(self.DDIR+"/fcnscan_"+self.outbasename+"_iter_"+str(iternum)+".root"):
                return iternum

    def ReUseStoredSim(self,simfile):
            f=open(simfile)
            f.readline()
            for line in f.readlines():
                ns_inj=int(line.split()[0])
                nsig=int(line.split()[1])
                TS=float(line.split()[2])
                ns_fit=float(line.split()[3])
                self.AddTrial(ns_inj,nsig,ns_fit,TS,fwrite=False)
            f.close()

    def FindStartPoint(self):
        #Do trials until first detection, return ns for first detection
        if len(self.CollectedTrials) !=0:
            nlow_detection=99999
            for i in self.CollectedTrials.keys():
                for k in range(len(self.CollectedTrials[i])):
                    if self.CollectedTrials[i][k][2]<self.signif:
                        print PF(),"DETECTION in stored sims!!!"
                        if i<nlow_detection:
                            nlow_detection=i
            if nlow_detection<99999:
                return nlow_detection

        detection=False
        print PF(),"DETECTION not in stored sims!!! , find injecting signal"
        n_detection=0
        while not detection:
            n_detection+=1
            inj_copy=copy.deepcopy(self.inj)
            TS,nsfit,nserr,nsig,gr=DoOneScramble(None,(n_detection,inj_copy,n_detection,len(inj_copy.Det["IsotropicHS"].EvList)))
            if DEBUG>1:
                print PF(),TS,nsfit,nserr,nsig
            prob=self.AddTrial(n_detection,nsig,nsfit,TS)
            if prob < self.signif:
                detection=True
                print PF(),"DETECTION!!!"
        return n_detection #returning the ns for the first discovery, can happen that jobs with lower possible ne are still running, doesn't matter, rough estimate is goo enough

    def AddTrial(self,ns_inj,nsig,ns_fit,TS,fwrite=True):
        if ns_inj not in self.CollectedTrials.keys():
            self.CollectedTrials[ns_inj]=[]
        prob=(chisqprob(TS, self.chi2))/2. #one-sided chi-sq prob
        self.CollectedTrials[ns_inj].append([TS,ns_fit,prob,nsig])
        if fwrite:
            f=open(self.CollectedTrialsFile,"a")
            f.write(str(ns_inj)+"\t"+str(nsig)+"\t"+str(TS)+"\t"+str(ns_fit)+"\t"+str(prob)+"\n")
            f.close()
        return prob

    def DrawCollectedTrials(self,mean,max_ns,cname):
        c_CT=TCanvas()
        c_CT.Divide(2,2)

        h_ns_inj=TH1D("h_ns_inj","h_ns_inj",max_ns+1,-0.5,max_ns+0.5)
        h_ns_inj.SetLineColor(2)
        h_ns_fit=TH1D("h_ns_fit","h_ns_fit",max_ns+1,-0.5,max_ns+0.5)

        h_nsW_inj=TH1D("h_nsW_inj","h_nsW_inj",max_ns+1,-0.5,max_ns+0.5)
        h_nsW_inj.SetLineColor(2)

        h_nsW_fit=TH1D("h_nsW_fit","h_nsW_fit",max_ns+1,-0.5,max_ns+0.5)

        h_TS=TH1D("TS","TS",400,0,100)
        h_TSW=TH1D("TSW","TSW",400,0,100)

        weights=self.GetWeightsForPoisson(mean)

        for i in self.CollectedTrials.keys():
            for k in range(len(self.CollectedTrials[i])):
                h_ns_inj.Fill(i)

                h_ns_fit.Fill(self.CollectedTrials[i][k][1])
                h_nsW_inj.Fill(i,weights[i])


                h_nsW_fit.Fill(self.CollectedTrials[i][k][1],weights[i])

                h_TS.Fill(self.CollectedTrials[i][k][0])
                h_TSW.Fill(self.CollectedTrials[i][k][0],weights[i])

        c_CT.cd(1)
        h_TS.Draw()
        c_CT.cd(2)
        h_TSW.Draw()
        c_CT.cd(3)
        h_ns_inj.Draw()
        h_ns_fit.Draw("same")
        c_CT.cd(4)
        h_nsW_inj.Draw()
        h_nsW_fit.Draw("same")
        c_CT.Update()
        c_CT.SaveAs(cname)

    def GetWeightsForPoisson(self,mean=1.):
        pp=poisson(mean)
        weights={}
        cutfactor=pp.cdf(self.maxns)
        for i in self.CollectedTrials.keys():
            try:
                weights[i]=pp.pmf(i)/(len(self.CollectedTrials[i])*cutfactor)
            except Exception,e:
                weights[i]=0.
        return weights

    def fcn(self,npar, gin, f, par, iflag):
        weights=self.GetWeightsForPoisson(par[0])

        o_values=[]
        o_weigts=[]
        for i in self.CollectedTrials.keys():
            for k in range(len(self.CollectedTrials[i])):
                o_values.append(self.CollectedTrials[i][k][0])
                o_weigts.append(weights[i])
        valweights = zip(o_values, o_weigts)
        valweights.sort()
        indexv=0
        ssum=0
        for i in range(len(valweights)):
            ssum+=valweights[i][1]
            if ssum>=(1-self.power):
                indexv=i
                break
        percentile=valweights[indexv][0]
        if chisqprob(percentile, 1)/2.==0.5:
            f[0]=1.
        else:
            f[0]=log10(fabs(self.signif-(chisqprob(percentile, 1))/2.))

    def uncertaintyEstimator(self,nSuccess, nTrials):
        if nTrials<=0:
            return 1.
        pEst = float(nSuccess)/nTrials
        return sqrt(4.*nTrials*pEst*(1.-pEst) + 1.)/(nTrials+1.)

    def GetNextMean(self,step=1.,iternum=0):
        fcnscan=TH1D("fcnscan","fcnscan",self.maxns+1,-0.5,self.maxns+0.5)
        npoints=10000
        fcnscan=TGraph(npoints)
        f=[-99.]
        minmean=99999.
        minprobval=99999.
        for i in range(npoints):
            meanval=(self.maxns*(i+1))/float(npoints)
            self.fcn(None, None, f, [meanval], None)
            if minprobval>f[0]:
                minprobval=f[0]
                minmean=meanval
            fcnscan.SetPoint(i+1,meanval,f[0])
        fcnscan.SaveAs(self.DDIR+"/fcnscan_"+self.outbasename+"_iter_"+str(iternum)+".root")
        fcnscan.Draw("apl")

        startv=minmean
        print "startv",startv

        gMinuit = ROOT.TMinuit(3)
        gMinuit.SetFCN(self.fcn)
        #gMinuit.SetPrintLevel(-1)
        arglist = arr('d', 2*[0.01])
        ierflg = ROOT.Long(0)

        arglist[0] = 0.5
        gMinuit.mnexcm("SET ERR", arglist ,1,ierflg)

        # Set initial parameter values for fit
        vstart = arr( 'd', (startv,) )
        # Set step size for fit
        step =   arr( 'd', (step,) )

        # Define the parameters for the fit
        gMinuit.mnparm(0, "mean",  vstart[0], step[0], 0. ,self.maxns,ierflg)

        arglist[0] = 6000 # Number of calls to FCN before giving up.
        arglist[1] = 0.3  # Tolerance
        gMinuit.mnexcm("MIGRAD", arglist ,1,ierflg)

        nsfit,nserr = ROOT.Double(0), ROOT.Double(0)

        gMinuit.GetParameter(0,nsfit,nserr )

        weights=self.GetWeightsForPoisson(nsfit)
        sumSqError=0

        pp=poisson(nsfit)
        cutfactor=pp.cdf(self.maxns)
        for i in range(self.maxns):
            n=0
            nSuccess=0
            try:
                w=pp.pmf(i)/(1.*cutfactor)
            except:
                w=0.
            if i in self.CollectedTrials.keys():
                n=len(self.CollectedTrials[i])
                w=weights[i]
                for k in range(n):
                    if self.CollectedTrials[i][k][2]>=self.signif:
                        nSuccess+=1
            error = self.uncertaintyEstimator(nSuccess, n)
            try:
                sumSqError+=w*w*error*error
            except:
                pass

        uncertainty = sqrt(sumSqError)
        ffout=open(self.DDIR+"/rezults_"+self.outbasename+".dat","a+")
        ffout.write(str(minmean)+"\t"+str(nsfit)+"\t+-"+str(nserr)+"\t+-"+str(uncertainty)+"\n")
        ffout.close()
        return nsfit,nserr

    def InjectWithPoissonMean(self,nInjections,mean):
        nstartedtot=0
        nsig=[]
        for i in range(nInjections):
            print PF(), "trial",i,"of",nInjections
            val=np.random.poisson(mean, 1)[0]
            inj_copy=copy.deepcopy(self.inj)
            TS,nsfit,nserr,nsig,gr=DoOneScramble(None,(i,inj_copy,val,len(inj_copy.Det["IsotropicHS"].EvList)))
            self.AddTrial(val,nsig,nsfit,TS)
