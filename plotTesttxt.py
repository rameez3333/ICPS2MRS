#!/usr/bin/env python
from pylab import *
from ROOTtoPython import HistToHist
import ROOT
import glob
import sys
import numpy as np
from scipy import special

rootfile = sys.argv[1]
#rootfile2 = sys.argv[2]

ax = subplot(111)
ax.cla()

file = open(rootfile)

h = ROOT.TH1D("hTestStatistic", ";2 ln #lambda;trials", 10000, -1,20);
ls=[]
count =0

for line in file:
    #if len(line.split()) == 2:
	print line.split()[0]
        h.Fill(float(line.replace("[ ","").replace("]","").replace("[","").split()[0]))
        print float(line.replace("[","").replace("]","").replace("[","").split()[0])
	ls.append(float(line.replace("[","").replace("]","").replace("[","").split()[0]))
        count = count+1

print "Total ", count, "trials"

print 'Median TS', np.median(np.asarray(ls))
        

#f2 = ROOT.TFile(rootfile2)

#h = f.Get("hTestStatistic")

x=(2.67482,2.67482)

observed = (1,90)





#h2 =  f2.Get("hTestStatistic")

#h.Rebin()
h.Rebin()
h.Rebin()
h.Rebin()
h.Rebin()
h.Rebin()
h.Rebin()
h.Rebin()
h.Rebin()

#h2.Rebin()
#h2.Rebin()
#h2.Rebin()
#h2.Rebin()
#h2.Rebin()
#h2.Rebin()

ax.set_yscale("log")
#HistToHist(h2, ax, color="green", label = "TS Distribution (ns = 25)")
HistToHist(h, ax, color="Green", histtype = 'stepfilled', alpha = 0.5, label = "Null TS Distribution (ns = 0)")
#ax.plot(x, observed, color="black", lw=5, label="Observed $2\log\lambda$ (Cygnus X3, Best of 8) = 2.67482")

x = array([h.GetBinCenter(i) for i in range(0,h.GetNbinsX()+1)])
print x

x = array(arange(0., 20.,0.01))


k = 3.0


nTrials = h.GetSum()
xRange = h.GetXaxis().GetXmax() - h.GetXaxis().GetXmin()
binsPerUnit = h.GetNbinsX() / xRange
norm = (nTrials/2.)/binsPerUnit


chi2 = norm * pow(x, k/2 - 1)*exp(-x*0.5)/(pow(2,k/2)*special.gamma(k/2))

ax.plot(x, chi2, color="red", lw=2, label="chi square ndof = 3")

k = 1.0

chi2 = norm * pow(x, k/2 - 1)*exp(-x*0.5)/(pow(2,k/2)*special.gamma(k/2))

ax.plot(x, chi2, color="black", lw=2, label="ndof = 1")


k = 2.0

chi2 = norm * pow(x, k/2 - 1)*exp(-x*0.5)/(pow(2,k/2)*special.gamma(k/2))

ax.plot(x, chi2, color="blue", lw=2, label="ndof = 2")


#ax.text(10, 30, 'Post Trial p value = 80.5%', style='normal', bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

ax.set_ylabel("Number of trials", fontsize=15)
ax.set_xlabel("$2\log\lambda$", fontsize=15)
ax.set_ylim(1e0, 1e4)
print h.GetBinLowEdge(2)
print h.GetBinCenter(2)
ax.set_xlim(0 - h.GetBinWidth(1)/2., 20)



outfile = rootfile.replace(".txt", ".png")
legend()
grid()
savefig(outfile)
show()
    
    
