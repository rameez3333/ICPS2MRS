from glob import glob
import os, sys
import numpy as np
import time
from scipy.stats import linregress
import matplotlib.pyplot as plt


def converternow(input):                                                                                 
    if input=='True':
        return 1.0
    else:
        return 0.0

def converterpval(input):
    if float(input):
        return float(input)
    else:
        return np.nan

M = sys.argv[1]

pslimlistdict={}
diffmeaslowlistdict={}
diffmeashighlistdict={}
our5sigmadict={}
our3sigmadict={}
oursensdict={}
countsens=[]
countsensdict={}
pslim=[]
our5sigma=[]
our3sigma=[]
oursens=[]

diffmeaslow=[]
diffmeashigh=[]
flistdict = {}

ourtestlums=[]
ourtestdens=[]

allums=[]
alldens = []
countlums=[]

foldname = 'Whew'+str(M)

base = np.linspace(39.0, 45, 25)

for lum in base:
    pslimlistdict[lum] = []
    diffmeaslowlistdict[lum] = []
    diffmeashighlistdict[lum]=[]
    our3sigmadict[lum]=[]
    oursensdict[lum]=[]
    countsensdict[lum]=[]
    flistdict[lum] = sorted(glob(foldname+"/"+str(np.power(10., lum))+"/Trial*.txt"))
    
for lum in base:
    for f in flistdict[lum]:
        print "Handling Now:", f
        try:
            arr = np.genfromtxt(f, delimiter="|", skip_header=1, converters={8:converterpval, 9:converternow, 10:lambda x: float(x.replace('\n', ''))})
            densit = arr.transpose()[2][np.where(arr.transpose()[9]==1.0)[0][0]-1]
            pslimlistdict[lum].append(densit)
        except:
            print "Probably unfinished job"
        try:
            difflow = arr.transpose()[2][np.where(arr.transpose()[3]>9.2e-7)[0][0]]
            diffmeaslowlistdict[lum].append(difflow)
        except:
            print 'wtf1'
        try:
            diffhigh = arr.transpose()[2][np.where(arr.transpose()[3]<3.7e-4)[0][-1]]
            diffmeashighlistdict[lum].append(diffhigh)
        except:
            print 'wtf2'

        try:
            ours = arr.transpose()[2]
            for den in ours:
                allums.append(lum)
                alldens.append(den)
        except:
            print 'wtf3'     
        
        
        try:
            ours = arr.transpose()[2][np.where(arr.transpose()[8]<0.003)]
            for den in ours:
                ourtestlums.append(lum)
                ourtestdens.append(den)
            
        except:
            print 'wtf4' 
        
        try:
            #print np.where(arr.transpose()[8]<0.01)
            ourlow = arr.transpose()[2][np.where(arr.transpose()[8]<0.003)[0][0]]
            print ourlow
            our3sigmadict[lum].append(ourlow)
        except:
            print 'wtf5'
        try:
            oursensval = arr.transpose()[2][np.where(arr.transpose()[8]<0.48)[0][0]]
            oursensdict[lum].append(oursensval)
        except:
            print 'wtf6'
        try:
            countsensval = arr.transpose()[2][np.where(arr.transpose()[4]>=48)[0][0]]
            countsensdict[lum].append(countsensval)
        except:
            print 'wtf7'
            
        #print densit, difflow, diffhigh
        
        
        
            
        
        
disclums=[]
senslums=[]
difflowlums=[]
diffhighlums=[]
pslums=[]
        
for lum in base:
    print lum, np.asarray(oursensdict[lum])
    try:
        pslim.append(np.percentile(np.asarray(pslimlistdict[lum]), 50.))
        pslums.append(lum)
    except:
        print "PSlim not valid"
    try:
        diffmeaslow.append(np.percentile(np.asarray(diffmeaslowlistdict[lum]), 10.))
        difflowlums.append(lum)
    except:
        print "Hmm"
    try:
        diffmeashigh.append(np.percentile(np.asarray(diffmeashighlistdict[lum]), 90.))
        diffhighlums.append(lum)
    except:
        print "Hmm"
    try:
        countsens.append(np.percentile(np.asarray(countsensdict[lum]), 90.))
        countlums.append(lum)
    except:
        print "Fuck"
print oursensdict    
for lum in base:
    try:
        oursens.append(np.percentile(np.asarray(oursensdict[lum]), 10.))
        senslums.append(lum)
    except:
        print 'Sens not valid at lum', lum
        break

print our3sigmadict    
for lum in base:
    try:
        our3sigma.append(np.percentile(our3sigmadict[lum], 50.))
        disclums.append(lum)
    except:
        print 'Disc not valid at lum', lum
        break
#print pslim


#print pslums, pslim
#print senslums, oursens

print "Lets See", countlums, countsens, countsensdict

print linregress(np.asarray(pslums)[5:], np.log10(np.asarray(pslim))[5:])

plt.plot(np.asarray(pslums), np.asarray(pslim), color='black', linewidth=2, label='90%C.L UL from PS non detection')
plt.plot(np.asarray(countlums), np.asarray(countsens), color='brown', linewidth=2, label='90%C.L UL from PS Warmspot count')
plt.plot(np.asarray(diffhighlums), np.asarray(diffmeaslow), color='red', linewidth=2)
plt.plot(np.asarray(difflowlums), np.asarray(diffmeashigh), color='blue', linewidth=2)
plt.plot(np.asarray(senslums), np.asarray(oursens), color='green', linewidth=2, label='Sensitivity')
plt.plot()
#plt.plot(disclums[0:12], np.asarray(our3sigma)[0:12], color='blue', linewidth=2, label= '3 sigma discovery potential')


#plt.scatter(allums, alldens, color='black')
#plt.scatter(ourtestlums, ourtestdens, color='red')
plt.ylabel('Source Density at z=0 (MPc^-3)')
plt.xlabel('log10(Neutrino Luminosity, [erg/s])')
plt.title('Evolution, (1+z)^'+str(M))
plt.yscale('log')
plt.legend(loc='best', fontsize=18.)

plt.show()






