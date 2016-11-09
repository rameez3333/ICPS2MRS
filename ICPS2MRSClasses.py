#!/usr/bin/env python
from niceout import *
import sys,random,glob
import PAO_ev_list231_final as PAO_ev_list
import healpy as hp
import scipy.integrate as integrate
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy  as np
from astropy.cosmology import Planck15 as cosmo


comoving_distance = InterpolatedUnivariateSpline(np.linspace(0, 3, 1000), cosmo.comoving_distance(np.linspace(0, 3, 1000)).value)

radians=np.radians
degrees=np.degrees
pi     =np.pi
def exp10(v):
    return np.power(10,v)
log10  =np.log10
sin    =np.sin
cos    =np.cos
arccos =np.arccos
sqrt   =np.sqrt
power  =np.power
arcsin =np.arcsin

import matplotlib.pyplot as plt


import ROOT
TCanvas    =ROOT.TCanvas
TH1D       =ROOT.TH1D
TH2D       =ROOT.TH2D
TEllipse   =ROOT.TEllipse
TMarker    =ROOT.TMarker
TFile      =ROOT.TFile

from array import array as arr

def setDEBUG(debug):
    global DEBUG
    DEBUG=debug
    if DEBUG>0:
        global pl
        import pylab as pl
        
def setNSIDE(nside):
    global NSIDE
    NSIDE=nside
    if DEBUG:
        print PF(),"NSIDE set to",NSIDE

def setDDir(ddir):
    global DDIR
    DDIR=ddir
    if not os.path.exists(DDIR):
        try:
            os.makedirs(DDIR)
        except :
            print "Weird"
    if DEBUG:
        print PF(),"the output directory DDIR set to",DDIR

HUBBLEParam = cosmo.H0.value #kms^-1Mpc^-1
MPctocm = 3.085677581e+24
ergtoGeV = 624.151        

def LumToFlux(Lum, znow):
    d = comoving_distance(znow)*MPctocm
    flin = Lum/(4.*np.pi*d**2.)
    phi = flin*ergtoGeV
    return phi

def LumFluxToZ(Lum, Flux):
    flin  = Flux/ergtoGeV
    d = np.power(Lum/(flin*4.*np.pi), 0.5)
    return d/MPctocm*HUBBLEParam/299792.0
     
    
    

def N0Scaler(m, zlim):
    integ = integrate.quad(lambda x: (1+x)**m*x**2, 0, zlim)
    return 3.*integ[0]/zlim**3.
     

#######################################################################################
#############Modify the following according to my preferred conventions################
#######################################################################################


def IndexToDeclRa(index):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -degrees(theta-pi/2.),degrees(pi*2.-phi)

def DeclRaToIndex(decl,RA):
    return hp.pixelfunc.ang2pix(NSIDE,radians(-decl+90.),radians(360.-RA))

def getPointD(decl,RA,m, interp=True):
    if interp:
        return hp.get_interp_val(m, radians(-decl+90.),radians(360.-RA))
    else:
        return m[DeclRaToIndex(decl,RA)]


#######################################################################################
#Verify the RA, DEC to healpy convention here is the same as others.
def scattomap(dec,ra, nside=16):
    hmap = np.zeros(hp.nside2npix(nside))
    hmap = hmap + np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(360. - ra)), minlength=hp.nside2npix(nside))
    return hmap


def Make2MRSMap(redshiftcutlow=0.00, redshiftcuthigh=0.173, nside=16):
    twomrsarr = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(1,2,24)).transpose()
    twomrsarr = np.vstack((twomrsarr[0], twomrsarr[1], twomrsarr[2]/299792.458))
    twomrsarr = twomrsarr.transpose()[twomrsarr[2]<redshiftcuthigh].transpose()
    twomrsarr = twomrsarr.transpose()[twomrsarr[2]>redshiftcutlow].transpose()
    return scattomap(twomrsarr[1], twomrsarr[0], nside)

#====================================================================================================================================================================
#====================================================================================================================================================================

def GetUnblindedResults(inj):
    print PF(), "calculating unblinnded results: Using real data!!!"
    inj.LoadRealDetSet("IsotropicHS")
    llh=LLH(inj,maxfilen=None)
    nsfit,nserr,TS,pvalue,gr=llh.Minimize()
    return nsfit,nserr,TS,pvalue,gr

#====================================================================================================================================================================
#====================================================================================================================================================================
#====================================================================================================================================================================
def DoOneScramble(results, args, post_trial=False,inj_bkg_pos_input=None,inj_sig_pos_input=None):
    n_tr,inj,nsig,nHS=args

    prestring="SCRAMBLE (nsig="+str(nsig)+") "+str(n_tr)+": "
    print prestring

    #inj.SetSeed()

    #inj.ChoseSigPosiotions()
    nsigHS=nsig

    #print "nsig="+str(nsig)+" split to "+str(nsigPAO)+" PAO and "+str(nsigTA)+" TA"
    inj.CRSet['IsotropicHS']=[]
    inj.Det["IsotropicHS"]=Detector(Name="IceCubePSSkyMap", EvList=IsoGenerator((inj.NU+inj.NL+inj.SU+inj.SL)/2, inj.NU+inj.NL, inj.SU+inj.SL))
    if post_trial:
        maxfilen,inj_bkg_pos=inj.CreateDetSet("IsotropicHS",n_sig=nsig,maxfilen=0       ,inj_bkg_pos_input=inj_bkg_pos_input,inj_sig_pos_input=inj_sig_pos_input)
    else:
        maxfilen,inj_bkg_posPAO=inj.CreateDetSet("IsotropicHS")

    llh=LLH(inj,maxfilen)
    nsfit,nserr,TS,pvalue,gr=llh.Minimize()

    if results==None and not post_trial:
        return TS,nsfit,nserr,nsig,gr
    elif results==None:
        return TS,nsfit,nserr,nsig,gr


def DoOneLuminosity(luminosity, inj, series= 1, chunksize=1):
    chunkscaler = False
    if not os.path.exists(DDIR+"/"+str(luminosity)):
        try:
            os.makedirs(DDIR+"/"+str(luminosity))
        except:
            print "Weird2"
    fout = open(DDIR+'/'+str(luminosity)+'/Trial'+str(series)+'.txt', "w")
    fout.write(str(inj.NL)+"|"+ str(inj.NO)+"|"+str(inj.NU)+"|"+str(inj.SL)+"|"+str(inj.SO)+"|"+str(inj.SU)+"\n")
    i=0
    optscale =0
    
    absdetectcount=0
    
    cnt = 0
    while True:
        cond, line = inj.AddOneSource(luminosity, chunksize=chunksize, counter=cnt)
        cnt = cnt +1
        if cond :
            fout.write(str(i-1)+"|"+line)
            inj.SetSeed()
            maxfilen,inj_bkg_pos = inj.CreateDetSet("IsotropicHS")
            llh=LLH(inj,maxfilen)
            nsfit,nserr,TS,pvalue,gr=llh.Minimize()
            print i, nsfit, nserr, TS, pvalue, gr
            print inj.Nsrc, inj.Ndensity, inj.TotDiffuseFlux, len(inj.Det["IsotropicHS"].EvList), inj.PTDisc
            fout.write(str(i)+"|"+str(inj.Nsrc)+"|"+str(inj.Ndensity)+"|"+str(inj.TotDiffuseFlux)+"|"+str(len(inj.Det["IsotropicHS"].EvList))+"|"+str(nsfit)+"|"+str(nserr)+"|"+str(TS[0])+"|"+str(pvalue)+"|"+str(inj.PTDisc)+"|"+str(inj.Det["IsotropicHS"].EvList[-1][-1])+"\n")
            if (inj.TotDiffuseFlux > 3.65e-3) or (len(inj.Det["IsotropicHS"].EvList)-(inj.NL+inj.SL+inj.NU+inj.SU)/2.) > 2*(inj.NO+inj.SO-(inj.NL+inj.SL+inj.NU+inj.SU)/2.):
                break
            if i%10 == 0:
                fout.write(str(i)+"|"+str(inj.Nsrc)+"|"+str(inj.Ndensity)+"|"+str(inj.TotDiffuseFlux)+"|"+str(len(inj.Det["IsotropicHS"].EvList))+"|"+str(np.nan)+"|"+str(np.nan)+"|"+str(np.nan)+"|"+str(np.nan)+"|"+str(inj.PTDisc)+"|"+str(inj.Det["IsotropicHS"].EvList[-1][-1])+"\n")
            i=i+1
            if inj.PTDisc and not chunkscaler:
                chunksize = chunksize*10
                chunkscaler = True
            if inj.PTDisc:
                absdetectcount=absdetectcount+1
    abscountnow = absdetectcount
    print "Switching at count:", cnt
    chunksize = chunksize*10.
    for j in range(0,30):
        cond, line = inj.AddOneSource(luminosity, chunksize=chunksize, counter=cnt+j)
        if cond:
            if (len(inj.Det["IsotropicHS"].EvList)-(inj.NL+inj.SL+inj.NU+inj.SU)/2.) > 10*(inj.NO+inj.SO-(inj.NL+inj.SL+inj.NU+inj.SU)/2.) and (inj.TotDiffuseFlux/3.65e-4 < 100. ) and not(inj.TotDiffuseFlux/3.65e-4 > 10. ) and (chunksize < 50000000):
                chunksize = chunksize*10.
                chunkscaler = True
            print inj.Nsrc, inj.Ndensity, inj.TotDiffuseFlux, len(inj.Det["IsotropicHS"].EvList), inj.PTDisc
            fout.write(str(i)+"|"+str(inj.Nsrc)+"|"+str(inj.Ndensity)+"|"+str(inj.TotDiffuseFlux)+"|"+str(len(inj.Det["IsotropicHS"].EvList))+"|"+str(np.nan)+"|"+str(np.nan)+"|"+str(np.nan)+"|"+str(np.nan)+"|"+str(inj.PTDisc)+"|"+str(inj.Det["IsotropicHS"].EvList[-1][-1])+"\n")
            if (inj.TotDiffuseFlux > 3.65e-3) and (len(inj.Det["IsotropicHS"].EvList)-(inj.NL+inj.SL+inj.NU+inj.SU)/2.) > 3*(inj.NO+inj.SO-(inj.NL+inj.SL+inj.NU+inj.SU)/2.):
                break
            i=i+1
            if inj.PTDisc:
                absdetectcount=absdetectcount+1

            
    print inj.Nsrc, inj.Ndensity, inj.TotDiffuseFlux, len(inj.Det["IsotropicHS"].EvList), inj.PTDisc
    fout.write(str(i)+"|"+str(inj.Nsrc)+"|"+str(inj.Ndensity)+"|"+str(inj.TotDiffuseFlux)+"|"+str(len(inj.Det["IsotropicHS"].EvList))+"|"+str(np.nan)+"|"+str(np.nan)+"|"+str(np.nan)+"|"+str(np.nan)+"|"+str(inj.PTDisc)+"|"+str(inj.Det["IsotropicHS"].EvList[-1][-1])+"\n")
    fout.close()
    arr = np.asarray(inj.Det["IsotropicHS"].EvList).transpose()[0]
    np.savetxt(DDIR+'/'+str(luminosity)+'/HSdist'+str(series)+'.txt', arr.transpose(), delimiter="|")
    
    
    


#====================================================================================================================================================================
class Detector:
    def __init__(self, Name="IceCubePSSkyMap", EvList=None):
        #self.b=b # detector latitude in degrees
        #self.tm=tm #detector max zenith angle in degrees
        #self.emin=emin
        #self.emax=emax
        #self.gamma=gamma
        self.Name=Name
        self.EvList=EvList
        self.ExposureNormMap=self.FillDetExposureNormMap()

        if DEBUG>2:
            hp.mollview(self.ExposureNormMap, title="Exposure Map "+self.Name)
            pl.savefig(DDIR+"/Exposures/Exposure_map_"+self.Name+"_moll.png")
            pl.show()
            pl.close('all')

            if not os.path.exists(DDIR+"/Exposures"):
                try:
                    os.makedirs(DDIR+"/Exposures")
                except:
                    print "Weird4"
            fig = plt.figure()
            plt.xlabel('Decl $[\circ]$)')
            plt.title("Exposure Function "+self.Name)
            plt.ylim([0,2])
            ax = fig.gca()
            ax.set_xticks(np.arange(-90.,90.,10.))
            ax.set_yticks(np.arange(0.,2.0,0.2))
            plt.grid()
            xdata=np.linspace(-90.,90.,num=1000)
            ydata=map(lambda v: self.DetExposure(v) ,xdata)
            ydataMapInterp=map(lambda v: getPointD(v,180.,self.ExposureNormMap)/self.ExpNorm,xdata)
            ydataMap180=map(lambda v: getPointD(v,180.,self.ExposureNormMap,False)/self.ExpNorm,xdata)
            ydataMap180_1=map(lambda v: getPointD(v,180.1,self.ExposureNormMap,False)/self.ExpNorm,xdata)
            plt.plot(xdata,ydata,"r-",xdata,ydataMapInterp,"ro",xdata,ydataMap180_1,"b--",xdata,ydataMap180,"b--")
            fig.savefig(DDIR+"/Exposures/Exposure_Function_"+self.Name+".png" )

            fig2 = plt.figure()
            plt.ylim([0,1.2])
            plt.xlabel('Decl $[\circ]$)')
            plt.title("Exposure Function "+self.Name+", bin size in RA taken into account")
            ax = fig.gca()
            ax.set_xticks(np.arange(-90.,90.,10.))
            ax.set_yticks(np.arange(0.,2.0,0.2))
            plt.grid()
            ydata=map(lambda v: self.DetExposure(v)*cos(radians(v)) ,xdata)
            plt.plot(xdata,ydata,"r-")
            fig2.savefig(DDIR+"/Exposures/Exposure_Function_decl_bin_size_"+self.Name+".png" )

    def GetRandomPosition(self, size=1):
        ra = np.random.uniform(0, 360, size=size)
        dec = np.rad2deg(np.arcsin(np.random.uniform(0, 0.9961946981, size=size))*np.random.choice([1.,-1.], size=size))
        if size==1:
            raD = ra[0]
            dec = dec[0]
        
        return dec,raD

    def GetRandomGamma(self, size=1):
        gamma = np.random.uniform(-1., -4.0, size=size)
        if size==1:
            gamma = gamma[0]
        return gamma

    def DetExposure(self,decl):
        if ((decl > 85.0) or (decl < -85.0)):
            return 1.0
        else:
            return 1.0

    def FillDetExposureNormMap(self):
        m = np.arange(hp.nside2npix(NSIDE))
        m=m*0.
        m=np.asarray(map(lambda x: self.DetExposure(IndexToDeclRa(x)[0]) , range(len(m))))
        if DEBUG>2:
            self.ExpNorm=1./sum(m)
        return m/sum(m)

#====================================================================================================================================================================
#====================================================================================================================================================================
#====================================================================================================================================================================

class PS_hs:
    def __init__(self, decl,RA,Gamma,sigma,HESEMapsSum,ExposureNormMap,sim):
        self.decl=decl
        self.RA=RA
        self.sigma=sigma
        self.ExposureRez=self.GetExposureRez(ExposureNormMap)
        if sim:
            self.StackedRez=self.GetStackedLLH(HESEMapsSum)
        else:
            self.StackedRez=self.GetStackedLLH(HESEMapsSum,ExposureNormMap)

    def GetExposureRez(self,ExposureNormMap):
        return getPointD(self.decl,self.RA,ExposureNormMap)

    def GetStackedLLH(self,HESEMapsSum,ExposureNormMap=None):
        sys.stdout = open(os.devnull, "w")
        sigmaR=radians(self.sigma)
        msmap=hp.smoothing(HESEMapsSum, sigma=sigmaR,verbose=False,lmax=64)
        sys.stdout = sys.__stdout__
        if ExposureNormMap!=None:
            msmap=msmap*ExposureNormMap
            msmap=msmap/sum(msmap)
        return getPointD(self.decl,self.RA,msmap)

#====================================================================================================================================================================

def IsoGenerator(size, N=1., S=1.):
    arr = []
    for i in range(0,size):
        ra = np.random.uniform(0, 360, size=1)
        dec = np.rad2deg(np.arcsin(np.random.uniform(0, 1, size=1))*np.random.choice([1.,-1.], size=1, p=[(float(N)/float(N+S)), (float(S)/float(N+S))]))
        en = np.random.uniform(2, 4, size=1)
        arr.append([dec, ra, en, np.nan])
    return arr

#====================================================================================================================================================================
#====================================================================================================================================================================
#====================================================================================================================================================================
class Injector:
    def __init__(self,TwoMRSMaps=None,sim=True, pthresh=3.e-3, zslices = 16, smoothing=0.0, m=3, zmax=3.):
        self.HESEMaps=TwoMRSMaps
        self.sim=sim
        self.pthresh = pthresh
        self.obs = ICPSObservations()
        self.perf = ICPSPerformance('DiscoSens.txt')
        self.NU, self.NL, self.NO, self.SU, self.SL, self.SO = self.obs.GetHSCounts(pthresh)
        self.PTDisc=False
        self.Det={}
        if sim:
            print 'Total Median Hotspots above threshold p value ',pthresh,':', (self.NU+self.NL+self.SU+self.SL)/2,'N:S', (self.NU+self.NL)/2, (self.SU+self.SL)/2
            print 'Total Observed Hotspots', self.NO + self.SO
            self.Det["IsotropicHS"]=Detector(Name="IceCubePSSkyMap", EvList=IsoGenerator((self.NU+self.NL+self.SU+self.SL)/2, self.NU+self.NL, self.SU+self.SL))
        self.c_evWeights={}
        self.c_evWeights["IsotropicHS"]=TCanvas()
        self.h_injPattern={}
        self.h_injPattern["IsotropicHS"]=None

        self.TwoMRSMapsSum={}
        self.TwoMRSMapsSum["IsotropicHS"]=self.HESEMapsSumf("IsotropicHS")
        self.CRSet={}
        self.CRSet["IsotropicHS"]=[]
        self.sigWsum={}
        self.TomoMaps={}
        self.Zslices = zslices
        self.Zarr = np.linspace(0., 0.15, zslices)
        
        for i in range(len(self.Zarr)-1):
            print PF(), 'Loading map in range ', self.Zarr[i], ' to ', self.Zarr[i+1]
            self.TomoMaps[self.Zarr[i]] = TwoMRSMap(Make2MRSMap(self.Zarr[i],  self.Zarr[i+1], 16), "SourceDist", smoothing)
        self.Nsrc = 0
        self.Ndensity = 0.
        self.M = m
        self.TotDiffuseFlux=0.
        self.Fluxes=[]
        self.Zmax=zmax
        self.scaler = N0Scaler(self.M, self.Zmax)
        self.DCMRmax=cosmo.comoving_distance(self.Zmax).value
        if self.M:
            self.Evoprobx = np.linspace(0, self.Zmax, 1500)
            self.Evoproby = np.power(1.+self.Evoprobx, self.M)*np.power(self.Evoprobx, 2.)
            self.backpolate = InterpolatedUnivariateSpline(self.Evoproby, self.Evoprobx)
        

    def GenerateDecRaReal(self, z):
        if z<self.Zarr[-1]:
            znearest = self.Zarr[np.abs(self.Zarr - z).argmin()]
            znearest = float('%.1g' % znearest)
            samplemap = self.TomoMaps[znearest].EqNormMap
        else:
            samplemap = np.ones(hp.nside2npix(16))/float(hp.nside2npix(16))
        
        index = np.random.choice(np.arange(hp.nside2npix(16)), size=1, p = samplemap)
        
        return IndexToDeclRa(index[0])
        
        
    def AddOneSource(self, lum, CR="IsotropicHS", chunksize=1, counter=0):
        self.PTDisc = False
        if self.M:
            rs = np.random.uniform(0, np.power((1.+self.Zmax), self.M)*np.power(self.Zmax,2), size=chunksize)
            z = self.backpolate(rs)
        else:
            rs = np.random.uniform(0, np.power(self.Zmax, 2), size=chunksize)
            z = np.power(rs, 0.5)
        
        lbefore = str(self.Nsrc)+"|"+str(self.Ndensity)+"|"+str(self.TotDiffuseFlux)+"|"+str(len(self.Det["IsotropicHS"].EvList))+"|"+str(np.nan)+"|"+str(np.nan)+"|"+str(np.nan)+"|"+str(np.nan)+"|"+str(self.PTDisc)+"|"+str(self.Det["IsotropicHS"].EvList[-1][-1])+"\n"
        
        flux = LumToFlux(lum, z)
        
        minsens =self.perf.minsens
        mindisc =self.perf.mindisc
        
        self.TotDiffuseFlux = self.TotDiffuseFlux+flux.sum()
        

        
        z = z[(flux>(minsens+mindisc)/2.)]
        flux = flux[(flux>(minsens+mindisc)/2.)]
        if len(z):
            print PF(),counter,':' ,len(z), "Sources survive"
        self.Nsrc=self.Nsrc+chunksize
        self.Ndensity = float(self.Nsrc)*3./(4.*np.pi*np.power(self.DCMRmax, 3.)*self.scaler)
        
        #print PF(), len(z), 'sources survive'
        
        if len(z)>0:
            print 'Current Density : ', self.Ndensity
        
        
        
        count = 0
        for znow, fl in zip(z, flux):
            dec, ra = self.GenerateDecRaReal(znow)
            sens = self.perf.Sensit(dec)
            disc = self.perf.Disco(dec)
            ptdisc = self.perf.PTSensit(dec)
            #print 'Diag:', dec, ra, fl ,sens, disc, ptdisc
            if fl > (sens+disc)/2.:
                self.Det[CR].EvList.append([dec, ra, 2., fl])
                count = count+1
            if fl > ptdisc:
                self.PTDisc = True
        return count, lbefore                  

        
        

    def CreateDetSet(self,CR,decl_sig=-9999.,RA_sig=-9999.,maxfilen=0,inj_bkg_pos_input=None,inj_sig_pos_input=None):
        print PF(),"CreateDetSet"
        #h_injBKG=TH2D("h_"+CR+"injBKG","h_"+CR+"injBKG",300,0.,360,300,-90,90)
        #h_injSIG=TH2D("h_"+CR+"injSIG","h_"+CR+"injSIG",300,0.,360,300,-90,90)
        #h_injEne=TH1D("h_"+CR+"injGamma","h_"+CR+"injGamma",40,-4.0,-1.0)

        inj_bkg_pos=[]
        #inj_sig_pos=[]
        for i in range(len(self.CRSet[CR]), len(self.Det[CR].EvList)):
            if inj_bkg_pos_input==None:
                pos=self.Det[CR].EvList[i]
                ene=pos[2]
                inj_bkg_pos.append(pos)
            else:
                pos=inj_bkg_pos_input[i][:2]
                ene=inj_bkg_pos_input[i][2]
            #h_injBKG.Fill(pos[1],pos[0])
            #h_injEne.Fill(ene)
            #print PF(),"creating "+CR+" (n_sig="+str(n_sig)+")    ",i,"/",len(self.Det[CR].EvList),ene,pos
            #print "Test:", pos[0]
            self.CRSet[CR].append(PS_hs(pos[0],pos[1],ene,0.,self.TwoMRSMapsSum[CR],self.Det[CR].ExposureNormMap,True))
        if inj_bkg_pos_input!=None:
            inj_bkg_pos=inj_bkg_pos_input

        #for i in range(n_sig):
            #if inj_sig_pos_input==None:
                #ene=self.Det[CR].GetRandomGamma()
                #decl_sig,RA_sig=self.GetRandomHESE(ene,CR)
                #inj_sig_pos.append([decl_sig,RA_sig,ene])
            #else:
                #decl_sig=inj_bkg_pos_input[i][0]
                #RA_sig=inj_bkg_pos_input[i][1]
                #ene=inj_bkg_pos_input[i][2]
            #h_injSIG.Fill(RA_sig,decl_sig)
            #h_injEne.Fill(ene)
            #print PF(),"creating "+CR+" (n_sig="+str(n_sig)+") sig",i,"/",n_sig,ene,decl_sig,RA_sig
            #self.CRSet[CR].append(PS_hs(decl_sig,RA_sig,ene,0.,self.TwoMRSMapsSum[CR],self.Det[CR].ExposureNormMap,True))
        #if inj_sig_pos_input!=None:
            #inj_sig_pos=inj_sig_pos_input

        #if maxfilen==0:
            #for ff in glob.glob(DDIR+"/injection_*_"+CR+".root"):
                #if int(ff.split("_")[-2])>maxfilen:
                    #maxfilen=int(ff.split("_")[-2])
            #maxfilen+=1

        #f=TFile(DDIR+"/injection_"+str(maxfilen)+"_"+CR+".root","RECREATE")
        #h_injBKG.Write()
        #h_injSIG.Write()
        #h_injEne.Write()
        #f.Close()
        print PF(),"CreateDetSet done"
        return maxfilen,inj_bkg_pos

    #def GetRandomHESE(self,ene,CR):
        #rn=np.random.uniform(0, self.sigWsum[CR])
        #for i in range(len(self.HESEMaps)):
            #rn-=self.HESEMaps[i].SigWeight[CR]
            #decl=-999
            #if rn<=0:
                #s=0.0
                #while True:
                    #decl=random.gauss(self.HESEMaps[i].SigDecl, s)
                    #if decl<-90.:
                        #decl=-180-decl
                    #if decl>90.:
                        #decl=180-decl
                    #if self.Det[CR].DetExposure(decl)>0.:
                        #break

                #s_scaled=s/cos(np.radians(decl))
                #ra=random.gauss(self.HESEMaps[i].SigRA, s_scaled)
                #ra=ra%360.
                #if DEBUG:
                    #print PF(),"selected HESE position",decl,ra
                #return decl,ra

    def SetSeed(self):
        np.random.seed()
        random.seed()

    #def ChoseSigPosiotions(self):
        #self.sigWsum["IsotropicHS"]=0.
        #for i in range(len(self.HESEMaps)):
            #self.HESEMaps[i].SetSignalPosiotion(self.Det)
            #for CR in ["IsotropicHS"]:
                #self.sigWsum[CR]+=self.HESEMaps[i].SigWeight[CR]


    def HESEMapsSumf(self,CR):
        if not os.path.exists(DDIR+"/SourcesPlots"):
            try:
                os.makedirs(DDIR+"/SourcesPlots")
            except:
                print "Weird3"
        m = np.arange(hp.nside2npix(NSIDE))
        m=m*0.

        for i in range(len(self.HESEMaps)):
            m=m+self.HESEMaps[i].EqNormMap
        if self.sim:
            m=m*self.Det[CR].ExposureNormMap
        m=m/sum(m)

        self.h_injPattern[CR]=TH2D("h_injPattern_"+CR,"h_injPattern_"+CR,300,0.,360,300,-90,90)
        for x in range(1,self.h_injPattern[CR].GetNbinsX()):
            for y in range(1,self.h_injPattern[CR].GetNbinsY()):
                self.h_injPattern[CR].SetBinContent(x,self.h_injPattern[CR].GetNbinsY()-y,hp.get_interp_val(m, (self.h_injPattern[CR].GetYaxis().GetBinCenter(y)+90.)*pi/180.,(360.-self.h_injPattern[CR].GetXaxis().GetBinCenter(x))*pi/180.))
        #self.h_injPattern[CR].SaveAs(DDIR+"/SourcesPlots/h_injPattern_"+CR+".root")
        self.c_evWeights[CR].cd()
        self.h_injPattern[CR].Draw("col")
        self.c_evWeights[CR].Update()
        if DEBUG>0:
            hp.mollview(m,title="",cbar=True,margins=(0.1,0.05,0.05,0.1))
            hp.graticule(color='white')
            hp.projtext(179, 30., '30.       ', lonlat=True, coord='E',color="Black",horizontalalignment='right')
            hp.projtext(179, 60., '60.       ', lonlat=True, coord='E',color="Black",horizontalalignment='right')
            hp.projtext(179, -30., '-30.       ', lonlat=True, coord='E',color="Black",horizontalalignment='right')
            hp.projtext(179., -60., '-60.       ', lonlat=True, coord='E',color="Black",horizontalalignment='right')
            hp.projtext(179., 0., '360.       ', lonlat=True, coord='E',color="Black",horizontalalignment='right')
            rvals=[30.*v for v in range(12)]
            rvals[rvals.index(180.0)]=180.01
            for i in rvals:
                if i!=179.99:
                    hp.projtext(360.-i, 1., '{:.0f}'.format(i), lonlat=True, coord='E',color="White")
                else:
                    hp.projtext(360.-i, 1., "180", lonlat=True, coord='E',color="White")
            pl.savefig(DDIR+"/SourcesPlots/h_injPattern_"+CR+"_moll.png")
            #pl.show()
            pl.close('all')
        return m


    def LoadRealDetSet(self,CR):
        self.c_evWeights[CR].cd()
        elps=[]
        mrks=[]

        allw=np.zeros(len(self.Det[CR].EvList))
        for cd in range(len(self.Det[CR].EvList)):
            RA_sig=self.Det[CR].EvList[cd][0]
            decl_sig=self.Det[CR].EvList[cd][1]
            ene=self.Det[CR].EvList[cd][2]*1E18
            if DEBUG>2:
                print PF(),"decl_sig",decl_sig
                print PF(),"RA_sig",RA_sig
                print PF(),"ene",ene
            self.CRSet[CR].append(PS_hs(decl_sig,RA_sig,ene,0.,self.TwoMRSMapsSum[CR],self.Det[CR].ExposureNormMap,False))
            elps.append(TEllipse(RA_sig,decl_sig,self.CRSet[CR][-1].StackedRez*3E6))
            elps[-1].SetFillColor(1)
            elps[-1].SetFillStyle(3001)
            elps[-1].Draw("same")
            mrks.append(TMarker(RA_sig,decl_sig,29))
            mrks[-1].SetMarkerSize(2.3)
            mrks[-1].SetMarkerColor(2)
            mrks[-1].Draw("same")
            allw[cd]=self.CRSet[CR][-1].StackedRez

        h_w=TH1D("h_w_"+CR,"h_w_"+CR,100,0,max(allw))
        h_w.FillN(len(allw),allw,np.ones(len(self.Det[CR].EvList)))
        #self.c_evWeights[CR].SaveAs(DDIR+"/c_evWeights_"+CR+".root")
        #h_w.SaveAs(DDIR+"/h_w_"+CR+".root")

#====================================================================================================================================================================
#====================================================================================================================================================================


class TwoMRSMap:
    def __init__(self, rawmap, mtype="SourceDist", smoothing = 0.0):
        if mtype=="SourceDist":
            self.EqNormMap=self.MapToEqNorm(rawmap, smoothing)
        else:
            print PF("err"), "mtype==\""+mtype+"\" unknown"
            exit(1)
        self.SigWeight={}

    def MapToEqNorm(self,eqmap, smoothing):
        if NSIDE!=16:
            if DEBUG:
                print PF("warn"),"rebinning EqNormMap to ",NSIDE
            eqmap=hp.pixelfunc.ud_grade(eqmap,NSIDE)
        
        if smoothing:
            eqmap = hp.smoothing(eqmap, sigma=np.deg2rad(smoothing),verbose=False,lmax=64)
        #print 'wtf', eqmap, np.sum(eqmap) 
        return eqmap/np.sum(eqmap)

    #def SetSignalPosiotion(self,Det):
        #s=np.random.uniform(0, 1)
        #cumulmap=np.cumsum(self.EqNormMap)
        #idx = (np.abs(cumulmap - s)).argmin()
        #self.SigDecl,self.SigRA=IndexToDeclRa(idx)

        #self.SigWeight["IsotropicHS"]=Det["IsotropicHS"].DetExposure(self.SigDecl)

#====================================================================================================================================================================


class LLH():
    def __init__(self,inj,maxfilen):
        self.hss=inj.CRSet["IsotropicHS"]
        self.SetMinimizer()
        self.ntot=len(self.hss)
    def LLHFunc(self,ns,histo=None):
        retval=0

        for i in range(self.ntot):
            try:
                retval+=np.log(ns*self.hss[i].StackedRez+(self.ntot-ns)*self.hss[i].ExposureRez)
            except Exception,e:
                print e
                print "StackedRez",self.hss[i].StackedRez,"ExposureRez",self.hss[i].ExposureRez,(self.ntot-ns),ns
                print "retval",retval
        return retval
        
    def LLHFuncRatio(self,ns):
        return -self.LLHFunc(ns)+self.LLHFunc(0)

    def SetMinimizer(self):
        def fcn(npar, gin, f, par, iflag):
            f[0]=self.LLHFuncRatio(par[0])

        self.gMinuit = ROOT.TMinuit(3)
        self.gMinuit.SetFCN(fcn)
        self.gMinuit.SetPrintLevel(-1)

    def Minimize(self,plotgr=False):
        arglist = arr('d', 2*[0.01])
        ierflg = ROOT.Long(0)

        arglist[0] = 0.5
        self.gMinuit.mnexcm("SET ERR", arglist ,1,ierflg)

        # Set initial parameter values for fit
        vstart = arr( 'd', (1.0,) )
        # Set step size for fit
        step =   arr( 'd', (0.01,) )

        # Define the parameters for the fit
        self.gMinuit.mnparm(0, "ns",  vstart[0], step[0], 0. ,self.ntot,ierflg)

        arglist[0] = 6000 # Number of calls to FCN before giving up.
        arglist[1] = 0.3  # Tolerance
        self.gMinuit.mnexcm("MIGRAD", arglist ,1,ierflg)

        nsfit,nserr = ROOT.Double(0), ROOT.Double(0)

        self.gMinuit.GetParameter(0,nsfit,nserr )

        TS=-2.*self.LLHFuncRatio(nsfit)
        nsfit=float(nsfit)
        nserr=float(nserr)
        pvalue=ROOT.TMath.Prob(TS,1)/2.
        print PF(),"nsfit=",nsfit,"nserr=",nserr,"ts=",TS, "pvalue=",pvalue

        if plotgr:
            gr=ROOT.TGraph((self.ntot)*2)
            for j in range((self.ntot)*2):
                gr.SetPoint(j,j/2.,self.LLHFuncRatio(j/2.))
            gr.Draw("apl")
        else:
            gr=None
        return nsfit,nserr,TS,pvalue,gr
    
class ICPSPerformance():
    def __init__(self, fname):
        self.sindec = np.genfromtxt(fname, delimiter=",", usecols=(0))
        self.fomflux =  np.genfromtxt(fname, delimiter=",", usecols=(1))
        self.sensdiscratio = self.fomflux[42]/self.fomflux[25]  
        self.sensdiscratiopt = self.fomflux[43]/self.fomflux[25]
        self.minsens = np.min(self.fomflux[0:42])*self.sensdiscratio*1.e3
        self.mindisc = self.minsens/self.sensdiscratio
    
    def Disco(self, decl):
        decl = np.deg2rad(decl)
        return np.interp([decl], np.arcsin(self.sindec)[0:42], self.fomflux[0:42])[0]*1.e3
    
    def Sensit(self, decl):
        return self.Disco(decl)*self.sensdiscratio
    
    def PTSensit(self, decl):
        return self.Disco(decl)*self.sensdiscratiopt
    
class ICPSObservations():
    def __init__(self):
        print PF(), "Loading IceCube Observations"
        self.NorthUcount = np.genfromtxt('HSNorthU.txt', usecols=(1), delimiter=",")
        self.NorthUthresh = np.genfromtxt('HSNorthU.txt', usecols=(0), delimiter=",")
        
        self.SouthUcount = np.genfromtxt('HSSouthU.txt', usecols=(1), delimiter=",")
        self.SouthUthresh = np.genfromtxt('HSSouthU.txt', usecols=(0), delimiter=",")

        self.NorthLcount = np.genfromtxt('HSNorthL.txt', usecols=(1), delimiter=",")
        self.NorthLthresh = np.genfromtxt('HSNorthL.txt', usecols=(0), delimiter=",")
        
        self.SouthLcount = np.genfromtxt('HSSouthL.txt', usecols=(1), delimiter=",")
        self.SouthLthresh = np.genfromtxt('HSSouthL.txt', usecols=(0), delimiter=",")

        self.NorthOcount = np.genfromtxt('HSNorthO.txt', usecols=(1), delimiter=",")
        self.NorthOthresh = np.genfromtxt('HSNorthO.txt', usecols=(0), delimiter=",")
        
        self.SouthOcount = np.genfromtxt('HSSouthO.txt', usecols=(1), delimiter=",")
        self.SouthOthresh = np.genfromtxt('HSSouthO.txt', usecols=(0), delimiter=",")

    def GetHSCounts(self, pval):
        xint = -1*np.log10(pval)
        
        NorthU = int(np.interp([xint], self.NorthUthresh, self.NorthUcount)[0])
        NorthL = int(np.interp([xint], self.NorthLthresh, self.NorthLcount)[0])
        NorthO = int(np.interp([xint], self.NorthOthresh, self.NorthOcount)[0])
        
        SouthU = int(np.interp([xint], self.SouthUthresh, self.SouthUcount)[0])
        SouthL = int(np.interp([xint], self.SouthLthresh, self.SouthLcount)[0])
        SouthO = int(np.interp([xint], self.SouthOthresh, self.SouthOcount)[0])      
        
        return NorthU, NorthL, NorthO, SouthU, SouthL, SouthO
    
