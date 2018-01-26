import os
import commands
import os.path
import optparse
import subprocess
import sys
import glob
import commands
import optparse
from subprocess import call


os.system("rm -r listfiles.txt")
os.system("rm core.*")


usage = ''
parser = optparse.OptionParser(usage)
parser.add_option('-P', '--process',dest='process',  type='string',     default = None, help="samples to add, according to the convention in 'script_rename.py'. Options are: 'All','ST','VJ','VV','QCDMu','QCDEle', or '_'+specific process, e.g. _ST_T_tch to add only the _ST_T_tch. Accepts also multiple processes, separated by comma, e.g. -P ST,_TT,VJ will select the V+jets samples and single top sample sets, as well as the one sample named TT.")

#Details of the analysis step:                                                                                                                                                                             
parser.add_option('-n', '--dryrun',   dest='dryrun',  action='store_true', default=False)
parser.add_option('-m', '--mode',   dest='mode', type='string',     default = 'resubmit', )

parser.add_option('-c', '--channel',        dest='channel',  type='string',     default = 'muon', )
parser.add_option('-d', '--dir',        dest='directory',  type='string',     default = '/afs/cern.ch/work/o/oiorio/public/xWajid/', )
parser.add_option('-S', '--size',        dest='size',  type=int,     default = 1000, )
parser.add_option('-I', '--inclusivelaunch',        dest='inc',  action='store_true', default=False)
#parser.add_option('-p', '--process',        dest='process',  type='string',     default = None )


(opt, args) = parser.parse_args() 

proc = opt.process
channel = opt.channel

initialDir = "./"
initialDir = "/afs/cern.ch/work/o/oiorio/public/xWajid/"
initialDir = opt.directory
#os.system("ls -lS "+initialDir+"res/*muonantiiso.root > listfilesantiiso.txt")
#os.system("ls -lS "+initialDir+"res/*.root > listfilesall.txt")
#os.system("ls -lS "+initialDir+"res/*electron.root > listfilesele.txt")
#allfilesele =open("listfilesele.txt","r")

if proc != None:
    samples =[]
    from utils import *
    samples = formSamples(opt.process)


os.system("ls -lS "+initialDir+"trees/*"+channel+".root > listfiles.txt")
allfiles =open("listfiles.txt","r")
    
allines= allfiles.readlines()
resubmitlines=[]

for l in allines:
    passesSample=True
    doPass=False
    if proc!=None:
        for s in samples:
            if "trees_"+s in l:
                doPass=True
        if not doPass: passesSample=False
    if not passesSample:continue
#    cmdlsbase= "yes yes | bpeek "+ l.split("\n")[0]
    size=int(l.split()[4]) 
    name=(l.split()[8]).split("trees/trees_")[1]
    
    #    print name
    print "size", size," name trees_"+name
    chan = "_"+channel
    if opt.mode=="resubmit":
        if(size<2000):
            print name
            print "smallsize! command is "
            namefolder= name[:-5]
            if opt.inc:
                namefolder=namefolder.replace(channel,"muonANDelectron")
            runfile = "work/qexe/"+namefolder+"_CMVAT_noSys/run.sh"
            runfolder = "work/qexe/"+namefolder+"_CMVAT_noSys"
            print "bsub -q 8nh "+ runfile
            finalpartsize=len(chan)+5
            cmdrun="./run.py -w work -q 8nh -e bsub "+namefolder+"_CMVAT_noSys"+" SingleTopAnalysis "+name[:-finalpartsize]+" work/files/"+name[:-finalpartsize]+".txt "+ chan[1:]+" CMVAT noSys noSync DATA trees noMVA ."
            print cmdrun
            if opt.dryrun:
                print 'Dry Run (command will not be executed)'
                continue
            else: os.system(cmdrun)
    if opt.mode=="tmpcp":
        if(size<2000):
            print name
            print "smallsize! command is "
            namefolder= name[:-5]
            finalpartsize=len(chan)+5
            cmdrun = "python new_singletop.py -m tmpcp -f work/splitfiles/ -P "+"_"+name[:-finalpartsize]+ " -t $WORKSPACE/public/xWajid/trees -o $WORKSPACE/public/xWajid/ -c "+channel
            print  cmdrun
            os.system(cmdrun)
#    cmdlsbase= "ls "+ l.split("\n")[0]
#    print cmdlsbase
#    linesout = commands.getstatusoutput(cmdlsbase)[1].splitlines()
##    linesout = ['hello','Error in <TNetXNGFile::Open>: [ERROR] Operation expired']
##    print linesout
#    nerr = 0
#    haserror=False
#    errorline = False
#    for lerr in linesout:
#        if "<< output from stderr >>" in lerr:
#            errorline=True
#        if errorline:
#            nerr = nerr+1
#            print "error line # ",nerr," is: ",lerr
#    if (nerr<4):
#        os.system("brequeue "+l.split("\n")[0])
         
#for r in resubmitlines:
#    print "Job with error is",r,"Resubmitting ..."
    #os.system("brequeue "+r.split("\n")[0])

