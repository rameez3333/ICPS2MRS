from glob import glob
import os, sys
import numpy as np
import time
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

foldname = 'Try'+str(M)

for lum in np.linspace(39, 45, 25):
    pslimlistdict[lum] = []
    diffmeaslowlistdict[lum] = []
    diffmeashighlistdict[lum]=[]
    our3sigmadict[lum]=[]
    oursensdict[lum]=[]
    flistdict[lum] = sorted(glob(foldname+"/"+str(np.power(10., lum))+"/Trial*.txt"))
    
for lum in np.linspace(39, 45, 25):
    for f in flistdict[lum]:
        print "Handling Now:", f
        arr = np.genfromtxt(f, delimiter="|", skip_header=1, converters={8:converterpval, 9:converternow, 10:lambda x: float(x.replace('\n', ''))})
        densit = arr.transpose()[2][np.where(arr.transpose()[9]==1.0)[0][0]-1]
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
            oursensval = arr.transpose()[2][np.where(arr.transpose()[8]<0.49999980976949693)[0][0]]
            oursensdict[lum].append(oursensval)
        except:
            print 'wtf6'
        
        #print densit, difflow, diffhigh
        
        
        
        if densit:
            pslimlistdict[lum].append(densit)
        
        
disclums=[]
senslums=[]

        
for lum in np.linspace(39, 45, 25)[0:len(np.linspace(39, 45, 25))-3]:
    print lum, np.asarray(oursensdict[lum])
    pslim.append(np.percentile(np.asarray(pslimlistdict[lum]), 90.))
    diffmeaslow.append(np.percentile(np.asarray(diffmeaslowlistdict[lum]), 10.))
    diffmeashigh.append(np.percentile(np.asarray(diffmeashighlistdict[lum]), 90.))


print oursensdict    
for lum in np.linspace(39, 45, 25)[0:len(np.linspace(39, 45, 25))-3]:
    try:
        oursens.append(np.percentile(np.asarray(oursensdict[lum]), 10.))
        senslums.append(lum)
    except:
        print 'Sens not valid at lum', lum
        break

print our3sigmadict    
for lum in np.linspace(39, 45, 25)[0:len(np.linspace(39, 45, 25))-3]:
    try:
        our3sigma.append(np.percentile(our3sigmadict[lum], 50.))
        disclums.append(lum)
    except:
        print 'Disc not valid at lum', lum
        break
#print pslim

plt.plot(np.linspace(39, 45, 25)[0:len(np.linspace(39, 45, 25))-3], np.asarray(pslim), color='black', linewidth=2, label='90\%C.L UL from PS non detection')
#plt.plot(np.linspace(39, 45, 25)[0:len(np.linspace(39, 45, 25))-3], np.asarray(diffmeaslow), color='red', linewidth=2)
#plt.plot(np.linspace(39, 45, 25)[0:len(np.linspace(39, 45, 25))-3], np.asarray(diffmeashigh), color='green', linewidth=2)
plt.plot(senslums, np.asarray(oursens), color='green', linewidth=2, label='Sensitivity')
plt.plot(disclums[0:12], np.asarray(our3sigma)[0:12], color='blue', linewidth=2, label= '3 sigma discovery potential')


#plt.scatter(allums, alldens, color='black')
#plt.scatter(ourtestlums, ourtestdens, color='red')
plt.ylabel('Source Density at z=0 (MPc^-3)')
plt.xlabel('log10(Neutrino Luminosity, [erg/s])')
plt.title('Evolution, (1+z)^0')
plt.yscale('log')
plt.legend(loc='best', fontsize=18.)

plt.show()






