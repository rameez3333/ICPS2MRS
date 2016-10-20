from glob import glob
import os, sys
import numpy as np
import time
nside = str(32)

flines = open('submittemplate.sh').readlines()


M = sys.argv[1]

jname="ICPS2MRS"
foldname = 'TryNlim'+str(M)
for lum in np.linspace(39, 45, 25):
    for ser in range(0,100):
        #if lum < 40.:
            #chunksize = 100
        #elif lum < 41:
            #chunksize = 10
        #else:
            #chunksize = 1
        chunksize = 1
        jobname=jname+str(lum)+'lum'+str(ser)+'ser'+str(M)+'job'
        fout = open(jobname+'.slurm', "w")
        jobline = 'python Scan_ICPS2MRS.py -l '+str(np.power(10., lum)) + ' -s ' +str(ser) + ' -f '+foldname+ ' -c ' + str(chunksize) + ' -m '+str(M)
        for line in flines:
            fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
        fout.close()
        os.system('chmod +x ' + jobname+'.slurm')
        os.system('sbatch -p icecube '+ jobname+'.slurm')
        time.sleep(0.1)
        #raw_input("Press Enter to Continue")
