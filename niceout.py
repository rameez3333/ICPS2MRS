#!/usr/bin/env python
import inspect,os

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

def PF(state=0):
    if state not in [0,"err","warn"]:
        print bcolors.FAIL+"state "+state+" unknown."+bcolors.ENDC
        exit(1)
    callerframerecord = inspect.stack()[1]    # 0 represents this line, 1 represents line at caller
    frame = callerframerecord[0]
    info = inspect.getframeinfo(frame)
    cole=bcolors.ENDC
    if state==0:
        col=bcolors.HEADER
    if state=="err":
        col=bcolors.FAIL
    if state=="warn":
        col=bcolors.WARNING
    retstr=col+os.path.basename(info.filename).ljust(20)+":"+str(info.lineno)+" in "+info.function+cole
    retstr=retstr.ljust(50)
    return retstr



