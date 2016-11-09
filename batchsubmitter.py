from glob import glob
import os, sys
import numpy as np
import time
nside = str(32)

flines = open('submittemplate.sh').readlines()


M = sys.argv[1]


def lumtochunksize(lum, Mloc):
    if Mloc==3:
        if lum >44:
            return 1
        else:
            return int(np.power(32, 44. - lum))
    elif Mloc==0:
        if lum >44:
            return 1
        else:
            return int(np.power(4, 44. - lum))

jname="ICPS2MRS"
foldname = 'Whew'+str(M)
for lum in np.linspace(45, 39, 25):
    chunksize = lumtochunksize(lum, int(M))
    print 'Chunksize : ', chunksize
    for ser in range(0,100):
        #elif lum >=38.75:
            #chunksize = 2370000
        #elif lum >=38.5:
            #chunksize = 5620000
        #elif lum >=38.25
            #chunksize = 13620000
        #chunksize = 1
        jobname=jname+str(lum)+'lum'+str(ser)+'ser'+str(M)+'job'
        print 'Lum:', lum
        print 'Chunksize : ', chunksize
        #raw_input("Press Enter to Continue")
        fout = open(jobname+'.slurm', "w")
        jobline = 'python Scan_ICPS2MRS.py -l '+str(np.power(10., lum)) + ' -s ' +str(ser) + ' -f '+foldname+ ' -c ' + str(chunksize) + ' -m '+str(M)
        for line in flines:
            fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
        fout.close()
        os.system('chmod +x ' + jobname+'.slurm')
        os.system('sbatch -p icecube '+ jobname+'.slurm')
        time.sleep(0.1)
        
