import os
import optparse 
import os.path
import optparse
import subprocess
import sys
import glob
import commands
#from utils import *
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

(opt, args) = parser.parse_args()

#define samples, one folder for each mass value
namesAndCfgs = {}  
discriminants=[]
#discriminants.append("ST_vs_TT")
#discriminants.append("ST_vs_VJ")
#discriminants.append("STsd_vs_VJ")
#discriminants.append("STsd_vs_TT")
#discriminants.append("STsd_vs_ST")
discriminants.append("STsd_vs_VJ")
#discriminants.append("STsd_vs_ST_test")
for d in discriminants:
    namesAndCfgs[d]="cfg"+d+".txt"
    os.system("root -l -b -q 'trainTMVA.C(\""+d+"\",\""+namesAndCfgs[d]+"\")' >& train"+d+".log &")
    #os.system("root -l 'trainTMVA.C(\""+d+"\",\""+namesAndCfgs[d]+"\")' >& train"+d+".log &")
