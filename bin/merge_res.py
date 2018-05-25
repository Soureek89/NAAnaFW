import os
import optparse 
import os.path
import optparse
import subprocess
import sys
import glob
import commands
from utils import *
#import colorama
#from colorama import Fore, Back, Style
#text = "The quick brown fox jumps over the lazy dog"
#print(Fore.RED + text)

#More complex working example: this will run only the ST_T_tch and the V+Jets samples and split them into batches of 10 files, taking them from the remote folder on Orso's public:
#Usage: python merge_res.py -l ./res -P ST,TT,VJ,VV --rm True

from os.path import join,exists
print 'Python version', sys.version_info
print '--------------------'
if sys.version_info < (2, 7):
    raise "Must use python 2.7 or greater. Have you forgotten to do cmsenv?"

workdir = 'work'
fileListDir = join(workdir,'files')

usage = 'usage: python merge_res.py -l ./res -P ST,TT,VJ,VV --rm True'
parser = optparse.OptionParser(usage)
#Input files:
parser.add_option('-P', '--process',        dest='process',  type='string',     default = 'all', help="samples to add, according to the convention in 'script_rename.py'. Options are: 'All','ST','VJ','VV','QCDMu','QCDEle', or '_'+specific process, e.g. _ST_T_tch to add only the _ST_T_tch. Accepts also multiple processes, separated by comma, e.g. -P ST,_TT,VJ will select the V+jets samples and single top sample sets, as well as the one sample named TT.")

#Details of the analysis step:
parser.add_option('-n', '--dryrun',   dest='dryrun',  action='store_true', default=False)
#Splitting options:
parser.add_option('-l', '--localpath',        dest='localpath',  type='string',     default = 'res/', )
parser.add_option('-c', '--channel',        dest='channel',  type='string',     default = 'muon', )
parser.add_option('-o', '--outpath',        dest='outpath',  type='string',     default = '', )
parser.add_option('--merge',   dest='domerge',  action='store_true', default=True)
parser.add_option('--retryFailed',   dest='doretry',  action='store_true', default=True)
parser.add_option('--rm',   dest='doremove',  action='store_true', default=False)
parser.add_option('--rmf',   dest='doforceremove', action='store_true', default=False)
parser.add_option('-t','--doTrees', dest='doTrees', action= 'store_true', default=False)

(opt, args) = parser.parse_args()



#define samples, one folder for each mass value
samples = []


samplesfound = formSamples(opt.process)

path = opt.localpath
if not path.endswith("/"): path = path+"/"

pathout = path
if not opt.outpath=='':
    pathout = opt.outpath
if not pathout.endswith("/"): pathout = pathout+"/"
#channelstorun=["muon","electron"]
channelstorun=[]
#pogchamp = []
for f in (opt.channel).split(","): channelstorun.append(f) 
if opt.doTrees: 
    for index, item in enumerate(samplesfound): 
        s = samplesfound[index]
        samplesfound[index]  = 'trees_'+s
        
print samplesfound
#print " pogchamp ", pogchamp
#channelstorun.append(opt.channel)
#channelstorun=["muonantiiso"]
syststorun = ["noSyst"]

toretry = opt.doretry

def domerge(samples, channels, toRetry={} ):
    samplesToRetry = {}
    systs=[]
    if(len(toRetry)>0): print "retrying ", toRetry
    for s in samples:
        systToRetry=[]
        if(len(toRetry)>0):
            keys = toRetry.keys()
            if not s in keys:
                print "s is ",s, " not to be retried"
                continue

        for c in channels:
            systs=[]
            cmdlsbase = "ls "+path+s+"_part*"+c+"*.root"
            
            linesall = commands.getstatusoutput(cmdlsbase)[1].splitlines()
            for l in linesall:
                systappend=(l.split("_")[-1])[:-5]
                if systappend==c: systappend= "noSyst"
                if systappend in systs: continue
                systs.append(systappend)
                print "appended syst",systs[-1]
            for sys in systs:
                if(len(toRetry)>0):
                    if not sys in toRetry[s]:
                        print "sys is ",sys, " not to be retried now"
                        continue
           
                print "check syst",sys
                cmdls = "ls "+path+s+"_part*"+c+"_"+sys+".root"
                if "noSys" in sys:
                    cmdls = "ls "+path+s+"_part*"+c+".root"

                if (commands.getstatusoutput(cmdls)[0]==0):
                    linessys = commands.getstatusoutput(cmdls)[1].splitlines()
                    print "looking for files with the command: ", cmdls
                    print linessys[0]
                    halves= linessys[0].split("_part")
                    secondhalf=halves[1]
                    secondhalfnumber=len(secondhalf.split("_")[0])
                    halves[1]=secondhalf[secondhalfnumber:]
                    print halves
                    mergeoutput=False
                    if(len(halves)==2 and opt.domerge):
                        cmdmerge="hadd -f "+halves[0].replace(path,pathout)+halves[1]+" "+path+s+"_part*"+c+"_"+sys+".root"
                        if "noSys" in sys:
                            cmdmerge="hadd -f "+halves[0].replace(path,pathout)+halves[1]+" "+path+s+"_part*"+c+".root"
                        print "merging with the command: ",cmdmerge
                        mergeoutput=True
                        statusoutput = (commands.getstatusoutput(cmdmerge))
                        mergeoutput=statusoutput[0] ==0
                        if (mergeoutput): print "merge successful!"
                        else: 
                            print "error! status is: ", statusoutput[0], "\n output is: \n",statusoutput[1]
                            finalline =(statusoutput[1].split("\n"))[-1]
                            print "final line is ", finalline
                            if "due to error in " in finalline and ".root" in finalline.split()[-1]:
                                toremove = finalline.split()[-1]
                                print "removing the problematic file:" , toremove, " trying again after other jobs finish"
                                os.system("rm "+ toremove)
                                systToRetry.append(sys)
                    if (mergeoutput and opt.doremove and not toretry):
                        cmdrm = "rm "+path+s+"_part*"+c+"_"+sys+".root"
                        if "noSys" in sys:
                            cmdrm = "rm "+path+s+"_part*"+c+".root"
                        print " now removing with the command: ", cmdrm
                        os.system(cmdrm)
	
                    if opt.doforceremove :
                        cmdrm = "rm "+path+s+"_part*"+c+"_"+sys+".root"
                        if "noSys" in sys:
                            cmdrm = "rm "+path+s+"_part*"+c+".root"
                        print " forcing removal with with the command: ", cmdrm
                        os.system(cmdrm)
	                    
	                
                    else: 
                        if len(halves)!=2: print "careful! didn't find part 0 for channel ",c," sys ",sys
        if len(systToRetry)>0:
            samplesToRetry[s]=systToRetry
    if(len(samplesToRetry)>0):
        print "retrying samples ",samplesToRetry
        domerge(samples, channels, samplesToRetry)
        
domerge(samplesfound, channelstorun)
