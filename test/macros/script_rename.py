import os, sys
import shutil
import errno
os.seteuid(os.geteuid())
import optparse 
from utils import *
#usage: python script_rename.py -s ~oiorio/public/xWajid/files/trees/nov09/ -d files/renamed 

usage = 'usage: python script_rename.py -s ~oiorio/public/xWajid/files/trees/nov09/ -d files/renamed '
parser = optparse.OptionParser(usage)
parser.add_option('-f', '--src',        dest='src',  type='string',     default = './', help="directory to rename files")
parser.add_option('-o', '--dest',       dest='dest',  type='string',    default = './renamed/', help="directory to rename files")
parser.add_option('-l', '--lepton',       dest='lep',  type='string',    default = 'SingleMuon', help="data trigger file")

(opt, args) = parser.parse_args()
#src='./'
#dest='./renamed'

src = opt.src
dest= opt.dest
mymkdir(dest)
#remove_folder(dest)

print '*** Copying the list of files from : ',src, 'to : ', dest 
#copytree(src,dest)
copyfiles(src,dest)
os.chdir(dest)
print "Current directory is: %s" %os.getcwd()

# listing directories
list=os.listdir(os.getcwd())

for i in list: 
  if 'listST_t-channel_top_4f' in i:
    print '*** Renaming file', i 
    os.renames(i,"ST_T_tch.txt")
    
  if 'listST_t-channel_antitop' in i:
    print '*** Renaming file', i 
    os.renames(i,"ST_Tbar_tch.txt")

  if 'listST_s-channel_4f' in i:
    print '*** Renaming file', i 
    os.renames(i,"ST_T_sch.txt")
  
  if 'listST_tW_top_5f' in i:
    print '*** Renaming file', i 
    os.renames(i,"ST_T_tW.txt")
  
  if 'listST_tW_antitop_5f' in i:
    print '*** Renaming file', i 
    os.renames(i,"ST_Tbar_tW.txt")

  if 'listTT' in i:
    print '*** Renaming file', i 
    os.renames(i,"TT.txt")
  
  if 'listWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt' in i or "listWToLNu" in i:
    print '*** Renaming file', i 
    firstpart=i.split("_")[0]+i.split("_")[1]
    os.renames(i, firstpart.replace("list","")+".txt")

#  if 'listWJetsTo' in i:
#    print '*** Renaming file', i 
#    os.renames(i,"WJetsToL.txt")
  
  if 'listDYJetsToLL' in i:
    print '*** Renaming file', i 
    os.renames(i,"DYJetsToLL.txt")

  if 'listWWTo1L1Nu2Q' in i:
    print '*** Renaming file', i 
    os.renames(i,"WWTo1L1Nu2Q.txt")
  
  if 'listWWTo2L2Nu' in i:
    print '*** Renaming file', i 
    os.renames(i,"WWTo2L2Nu.txt")

  if 'listWZTo1L1Nu2Q' in i:
    print '*** Renaming file', i 
    os.renames(i,"WZTo1L1Nu2Q.txt")
   
  if 'listWZTo2L2Q' in i:
    print '*** Renaming file', i 
    os.renames(i,"WZTo2L2Q.txt")
    
  if 'listZZTo2L2Q' in i:
    print '*** Renaming file', i 
    os.renames(i,"ZZTo2L2Q.txt")

  if 'listQCD_Pt-20toInf_MuEnrichedPt15' in i:
    print '*** Renaming file', i 
    os.renames(i,"QCDMuPt20toInf.txt")

  if 'listRun2016' in i:
    print '*** Renaming file', i 
    os.renames(i,i.replace("listRun2016",opt.lep+"_Run2016").replace("_B2GAnaFW_80X_V3p1","").replace("_June",""))
#    os.renames(i,i.replace("listRun2016","SingleMuon_Run2016").replace("_B2GAnaFW_80X_V3p1","").replace("_June",""))
#    if '_B2GAnaFW_80X_V2p1' in i:
#      os.renames(i,i.replace("_B2GAnaFW_80X_V2p1",""))


print '*** Renaming Done ***'
