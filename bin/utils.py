import os
import optparse 
import os.path
import optparse
import subprocess
import sys
import glob
import utils
import commands
import shutil
import errno

def formSamples(proclist):
    samp=[]
    procs= proclist.split(",")
    for proc in procs:
        isAllProcesses = proc=="All"
        if proc=="ST" or isAllProcesses:
            from samplesST import samples as stemp
            samp.extend(stemp)
        if proc=="ST_sd" or isAllProcesses:
            from samplesST_sd import samples as stemp
            samp.extend(stemp)
        if proc=="ST_p_sd" or isAllProcesses:
            from samplesST_p_sd import samples as stemp
            samp.extend(stemp)
        if proc=="ST_sd_reshaped" or isAllProcesses:
            from samplesST_sd_reshaped import samples as stemp
            samp.extend(stemp)
        if proc=="TT" or isAllProcesses:
            from samplesTT import samples as stemp
            samp.extend(stemp)
        if proc=="SystsTop" or isAllProcesses:
            from samplesSystsTop import samples as stemp
            samp.extend(stemp)
        if proc=="VJ" or isAllProcesses:
            from samplesVJ import samples as stemp
            samp.extend(stemp)
        if proc=="VV" or isAllProcesses:
            from samplesVV import samples as stemp
            samp.extend(stemp)
        if proc=="QCDMu" or isAllProcesses:
            from samplesQCDMu import samples as stemp
            samp.extend(stemp)
        if proc=="QCDEle" or isAllProcesses:
            from samplesQCDEle import samples as stemp
            samp.extend(stemp)
        if proc=="SingleMuon" or isAllProcesses:
            from samplesSingleMuon import samples as stemp
            samp.extend(stemp)
        if proc.startswith("_"):
            samp.append(proc[1:])
    return samp

def remove_folder(path):
    # check if folder exists
    if os.path.exists(path):
         # remove if exists
         print '*** Removing the folder : ',path
         shutil.rmtree(path)
         
def copyfiles(src, dest):
    os.system("cp "+src+"/*.txt "+dest)

def copytree(src, dest):
    try:
        shutil.copytree(src, dest)
    except OSError as e:
        # If the error was caused because the src wasn't a directory
        if e.errno == errno.ENOTDIR:
            shutil.copytree(src, dest)
        else:
            print('Directory not copied. Error: %s' % e)

def mymkdir(path):
    if commands.getstatusoutput("ls "+path)[0]!=0:
        os.system("mkdir " +path )
