import os, commands
import subprocess


crabStatus = "crab status "
crabResubmit = "crab resubmit "
cmdls = "ls "
import optparse 

usage = 'usage: %prog -d dirName'
parser = optparse.OptionParser(usage)
parser.add_option("-d","--dir",dest="dir",type="string",help="Crab project folder")

(opt, args) = parser.parse_args()


if opt.dir is None:
    parser.error('Please define the crab project folder with the option -d')

path = str(opt.dir + "/")

cmd = cmdls + path
status,ls_la = commands.getstatusoutput( cmd )


if(not os.path.isdir("Logs")): os.system('mkdir Logs')

d= "/"
cmd = crabStatus + path 
if(path.endswith("/")): 
    print path[:-1]
    path = path[:-1]
writer =open("Logs/"+path.replace("/","")+".log", 'w') 
process = subprocess.call(cmd, shell = True, stdout=writer)
writer.close()
ifile = open("Logs/"+path.replace("/","")+".log", 'r')
lines = ifile.readlines()
isFailed = False
for l in lines:
    if "failed" in l:
        isFailed = True
        break
    if ("Jobs status"in l) and ("100.0" in l): print " 100% jobs for "+d+" project are Finished"
    if ("Publication status"in l) and ("100.0" in l): print "---> 100% jobs for "+d+" project are Published"
if (isFailed):
    cmdResub = crabResubmit + path 
    print "Resubmitting ", path
    writerResub =open("Logs/"+path.replace("/","")+"_resub.log", 'w')
    process = subprocess.call(cmdResub, shell = True, stdout=writerResub)
    writerResub.close()

