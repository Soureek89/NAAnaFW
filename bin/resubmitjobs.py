import os
import commands
import os.path
import optparse
import subprocess
import sys
import glob
import commands
import optparse

os.system("rm -r listjobs.txt")
os.system("bjobs -rp -o jobid | grep -v ^JOBID | cut -d ':' -f2 > listjobs.txt")
alljobs =open("listjobs.txt","r")
allines= alljobs.readlines()
resubmitlines=[]
for l in allines:
    cmdlsbase= "bpeek "+ l.split("\n")[0]
    print cmdlsbase
    linesout = commands.getstatusoutput(cmdlsbase)[1].splitlines()
#    print linesout
    haserror=False
    for lerr in linesout:
        #print lerr
        #if "output from stderr" in lerr:
        if haserror: continue
        if ( ("[FATAL]" in lerr) and ("Redirect" in lerr) ) or ( ("[ERROR]" in lerr) and ("Invalid" in lerr) ):
            print "job",l," has fatal redirect error: ",lerr
            resubmitlines.append(l)
            haserror=True
for r in resubmitlines:
    print "job with error is ",r,", resubmitting..."
    os.system("brequeue "+r.split("\n")[0])
#echo "Getting the list of jobId's from CERN Batch Cluster ..."
#
##bjobs -r -o jobid | grep -v ^JOBID | cut -d ':' -f2 > listjobs.txt
#
#for i in `cat listjobs.txt`;
#    do
#      bkill $i;
#    done

