import os
import optparse 
import os.path
import optparse
import subprocess
import sys
import glob
from utils import *
#import colorama
#from colorama import Fore, Back, Style
#text = "The quick brown fox jumps over the lazy dog"
#print(Fore.RED + text)
#Usage: python new_singletop.py -c muon -s noSys --t3batch 
#Usage: python new_singletop.py -c muon -s noSys -m t3se 
#Usage: python new_singletop.py -c muon -s noSys -m local -l test/ -P _test
#Usage: python new_singletop.py -c muon -s noSys -m local -l test/data/ -P _test -d DATA
#Usage: python new_singletop.py -c muonantiiso -s noSys -m local -P _test -d DATA
#Usage: python new_singletop.py -c muonantiiso -s noSys -m local -P _test -l mc
#Usage: python new_singletop.py --t3batch -f files/final/ -P SingleMuon -d DATA -S 40 -c munonantiiso
#Usage: python new_singletop.py -c muonantiiso -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 40 -t trees
#python new_singletop.py -m tmpcp --tmpdir "/tmp/wajid" -f files/pisa/renamed/ -P _ST_T_tch_small -S 5 -n 1


#More complex working example: this will run only the ST_T_tch and the V+Jets samples and split them into batches of 10 files, taking them from the remote folder on Orso's public:
#Usage: python new_singletop.py -c muon -s noSys -m local -S 10 -P _ST_T_tch,VJ --t3batch

from os.path import join,exists
print 'Python version', sys.version_info
print '--------------------'
if sys.version_info < (2, 7):
    raise "Must use python 2.7 or greater. Have you forgotten to do cmsenv?"

workdir = 'work'
fileListDir = join(workdir,'files')

#define samples paths
#pathlocal = "/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_16/src/Analysis/NAAnaFW/test/crab_projects/crab_st_top/results/ST/" 
pathlocal = "/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin/ST/"
filepath='/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin/mc/'


usage = ''
parser = optparse.OptionParser(usage)
#Output files/directory
parser.add_option('-o', '--outputpath',        dest='outputpath',  type='string',     default = '.',  help='output for root files')

#Input files:
parser.add_option('-f', '--filepath',        dest='filepath',  type='string',     default = './files/final/',      help='files with the trees location, default to the afs area where the latest version is always located')
parser.add_option('-l', '--localpath',        dest='localpath',  type='string',     default = '/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin/ST/')
#Target process
parser.add_option('-P', '--process',        dest='process',  type='string',     default = 'all', help="samples to add, according to the convention in 'script_rename.py'. Options are: 'All','ST','VJ','VV','QCDMu','QCDEle', or '_'+specific process, e.g. _ST_T_tch to add only the _ST_T_tch. Accepts also multiple processes, separated by comma, e.g. -P ST,_TT,VJ will select the V+jets samples and single top sample sets, as well as the one sample named TT.")

#Details of the analysis step:
parser.add_option('-c', '--channel',  dest='channel', type='string',     default = 'muon', help='Channel to analyze: singleH or singleZ')
parser.add_option('-C', '--cat',      dest='cat',     type='string',     default = 'CMVAT',    help='Category to analyze: cat0 or cat1 or cat2')
parser.add_option('-s', '--sys',      dest='sys',     type='string',     default = 'noSys',   help='Systematics: jesUp, jesDown, jerUp, jerDown')
parser.add_option('',   '--sync',     dest='sync',    type='string',     default = 'noSync',  help='Synchro exercise')
parser.add_option('-d', '--isData',   dest='isData',  type='string',     default = 'MC',      help='is Data or MC?')
parser.add_option('-t', '--treesDir',        dest='treesDir',  type='string',  default = "noTrees",        help='trees directory, if blank no tree is generated')
parser.add_option('--mva',        dest='mva',  type='string',  default = "noMVA",        help='which mva to train')
parser.add_option('--tmpdir',     dest='tmpdir',  type='string',  default = "/tmp/oiorio/",        help='tmp directory for local tmp cp mode')
#Running mode details:
parser.add_option('-g', '--gdb',      dest='gdb',     action='store_true', default=False)
parser.add_option('-n', '--dryrun',   dest='dryrun',  action='store_true', default=False)
parser.add_option('-m', '--mode',     dest='mode',    default='t3se', choices=['local','t3se','tmpcp'])
parser.add_option('--t3batch',        dest='t3batch', action='store_true', default=False)
#Splitting options:
parser.add_option('-S', '--split',        dest='split',  type=int,     default = 0, help="Splitting each channel by batches of n jobs running simultaneously, where n= opt.splitMultiplicity. If s==0, no splitting is done(default). NOTA BENE! Doesn't work on local files atm.")


isData="MC"
(opt, args) = parser.parse_args()


if opt.sys not in ["noSys", "jesUp", "jesDown", "jerUp", "jerDown", "metUnclUp", "metUnclDown"]:
    parser.error('Please choose an allowed value for sys: "noSys", "jesUp", "jesDown", "jerUp", "jerDown","metUnclUp", "metUnclDown"')



#define samples, one folder for each mass value
samples = []

samples = formSamples(opt.process)

print "samples are:",samples

# Create working area if it doesn't exist
if not exists(fileListDir):
    os.makedirs(fileListDir)

# Create test/ directory if it doesn't exist
if not exists('test'):
    os.makedirs('test')

filePath= opt.filepath
pathlocal = opt.localpath

    
if opt.split!=0:#Modify the samples and txt list!
    samplesOld=samples
    samples=[]
    os.system("mkdir "+workdir+"/splitfiles")
    filePath= workdir+"/splitfiles"
    for s in samplesOld:
        sT2Path = join(opt.filepath,s+'.txt')
        print sT2Path 
        f = open(sT2Path,'r')
        listing = f.read()                
        files = listing.split()
        splitfiles=[]
        last =0.0
        while last < (len(files)):
            splitfiles.append(files[int(last):int(last + opt.split)])
#            print files[int(last)]
#            print files[int(last+opt.split)]
            last +=opt.split
        print splitfiles
        nsplit = len(splitfiles)
        print nsplit
        os.system("rm "+filePath+"/"+s+"_part*")
        for n in xrange(0,nsplit):
            partname=s+"_part"+str(n)
            samples.append(partname)
            fnew = open(filePath+"/"+partname+".txt",'w')
            for fl in splitfiles[n]:
                fnew.write(fl+"\n")
            fnew.close()
            
    for s in samples:
        print s 
        print join(filePath,s+'.txt')

for s in samples:
#    if opt.dry: continue
    #f (s.startswith("JetHT") or s.startswith("SingleMu") or s.startswith("SingleEl") or  s.startswith("MET") or opt.isData): isData="DATA"
    isData=opt.isData
    tmpPath =""
    if opt.mode == "tmpcp":
        tmpPath= opt.tmpdir+"/"+s
        os.system("mkdir "+tmpPath)
        sT2Path = join(filePath,s+'.txt')
        print sT2Path 
        f = open(sT2Path,'r')
        listing = f.read()
        fileslist = listing.split()
        issuccessful=False
        for f in fileslist[:-1]:
            localcpcommand="xrdcp "+f+ " "+ tmpPath + " > xrdout.log &"
 
#            localcpcommand="gfal-ls  srm://stormfe1.pi.infn.it//cms/"+(fileslist[-1]).replace("root://cms-xrd-global.cern.ch//","")
            #localcpcommand="gfal-copy  --force srm://stormfe1.pi.infn.it//cms/"+(f).replace("root://cms-xrd-global.cern.ch//","") + " " + tmpPath + "> gfalout.log &"
            #        localcpcommand="gfal-ls  srm://stormfe1.pi.infn.it//cms/"+(fileslist[-1]).replace("root://cms-xrd-global.cern.ch//","")
            print "copying with command ",localcpcommand
            #            print "final command is ", localcpcommand
            os.system(localcpcommand)
        localcpcommand="xrdcp "+fileslist[-1]+ " "+ tmpPath +" > xrdout.log "
#        localcpcommand="gfal-copy  --force srm://stormfe1.pi.infn.it//cms/"+(fileslist[-1]).replace("root://cms-xrd-global.cern.ch//","") + " " + tmpPath +" > gfalout.log "
#        localcpcommand="gfal-ls  srm://stormfe1.pi.infn.it//cms/"+(fileslist[-1]).replace("root://cms-xrd-global.cern.ch//","")
#        localcpcommand="gfal-ls  srm://stormfe1.pi.infn.it//cms/"+(fileslist[-1]).replace("root://cms-xrd-global.cern.ch//","")+ " "+ tmpPath 
        print "final command is ", localcpcommand
        os.system(localcpcommand)
        sPath = join(tmpPath,'/*.root')
        files = glob.glob(sPath)
        print 'Info: Sample',s,'Files found',len(files)
        issuccessufl= (commands.getstatusoutput(localcpcommand)[0]==0)
        
    if opt.mode == 'local':
        print 'Info: Running in local mode ...'
        sPath = join(pathlocal,'*.root')
        print 'Info: Looking for the *.root files at: ',sPath



        # listing = subprocess.check_output(lLs.split()+[sPath])
        files = glob.glob(sPath)
        print 'Info: Sample',s,'Files found',len(files)

    elif opt.mode == 't3se':
        print 'Info: Running on t3se in interactive mode ...'
        sPath = join(pathlocal,'*.root')
        sT2Path = join(filePath,s+'.txt')
        print sT2Path 
        f = open(sT2Path,'r')
        listing = f.read()
        files = listing.split()
        print 'Info: Sample',s,'Files found',len(files) 
        f.close() 

    # Save it to a semi-temp file
    sampleFileList = join(fileListDir,s+'.txt')
    
    print 'Info :',sampleFileList
    with open(sampleFileList,'w') as sl:
        sl.write('\n'.join(files))
    
    if opt.outputpath:
      #opath='/afs/cern.ch/work/w/wajid/public/xWajid'
      outDirs = [opt.outputpath + '/' + 'res', opt.outputpath + '/' + 'trees']
    else:
      outDirs = ['res','trees']
    
    #outDirs = ['res','trees']

    if (opt.treesDir!="./" and opt.treesDir!='noTrees'): outDirs.append(opt.treesDir)
    for d in outDirs:
        if exists(d): continue
        os.makedirs(d)
      
#    treescheck = (not (commands.getstatusoutput('ls '+ opt.treesDir)[0]==0))
#    print "treescheck, command is ",("ls "+opt.treesDir)," result ",treescheck
    
        
    cmd = 'SingleTopAnalysis '+ s + ' ' + sampleFileList  + ' ' + opt.channel + ' ' + opt.cat + ' ' + opt.sys + ' ' + opt.sync + ' ' + isData + ' ' + opt.treesDir + ' ' + opt.mva + ' ' + opt.outputpath
    channels= opt.channel    
    if "," in channels:
        cmd =''
        for c in channels.split(","):
            cmd = cmd + 'SingleTopAnalysis '+ s + ' ' + sampleFileList  + ' ' + c + ' ' + opt.cat + ' ' + opt.sys + ' ' + opt.sync + ' ' + isData + ' ' + opt.treesDir + ' ' + opt.mva + ' ' + opt.outputpath + " & \n "
    print "cmd is: ",cmd
       
#    if opt.mode=="tmpcp":
#        import time
#        time.sleep(30)
#        alldone=False
#        nretries=0
#        f = open(sT2Path,'r')
#        listing = f.read()
#        fileslistsize = len(listing.split())
#        while not alldone:
#            nretries = nretries+1
#            sizesVsFiles = commands.getstatusoutput("ls "+tmpPath+" -s")[1].splitlines()
#            if nretries >60: 
#                print "retrying for more than 60 times, the copy probably failed! Aborting..."
#                exit
#            if len(sizesVsFiles[1:])!=fileslistsize:
#                print " n files:",fileslistsize," files copying ", len(sizesVsFiles[1:]), "waiting 10 seconds before retrying"
#                time.sleep(10)
#                continue
#            print " n files:",fileslistsize," files copying ", len(sizesVsFiles[1:]), "waiting 10 seconds before retrying"
#            fileDone=True
#            for f in sizesVsFiles[1:]:
#                print " file ", f.split()[1], " size ", f.split[0]
#                if not( f.split()[0]>100):
#                    time.sleep(20)
#                fileDone=fileDone and ( f.split()[0]>100)
#            
    if opt.gdb:
        cmd = 'gdb --args '+cmd
    
    elif opt.t3batch:
        #usage: qexe.py [-h] [-w WORKDIR] [-q QUEUE] [-n] [-e {qsub,bsub}]
        print 'Info: Running on t3se in batch mode ...'
        sPath = join(pathlocal,'*.root')
#        que = '-q 2nw'
        que = '-q 8nh'
#        que = '-q 1nd'
        dry = '-n'
        batch = '-e bsub'
        jid = '%s_%s_%s_%s' % (s,opt.channel.replace(",","AND"),opt.cat,opt.sys)
        #cmd = './run.py -w ' + workdir + ' ' + que + ' ' + dry + ' ' + batch + ' ' + jid + ' '+cmd
        cmd = './run.py -w ' + workdir + ' ' + que + ' ' + batch + ' ' + jid + ' " '+cmd+ ' " '
        print "now cmd is ", cmd
    print 'Info:',cmd
 
    if opt.dryrun:
        print 'Dry Run (command will not be executed)'
        continue
    if opt.mode == "tmpcp":
        print "copy mode: no running will be done!"
        continue
    print '--------------------'
    subprocess.call(cmd,shell=True)

