### Decorrelate Q2 syst among the processes
import sys
import os, commands
import shutil
#from samples.toPlot import samples
from Analysis.NAAnaFW.MakePlot.samples.toPlot import samples as samplesDictionary
import optparse


usage = 'usage: %prog -d histosDir'
parser = optparse.OptionParser(usage)
parser.add_option('-d', '--dir', dest='path', type='string', default = './Plot', help='Histos folder')
(opt, args) = parser.parse_args()

samplesnames = []
samplesnames.append("ST_tch")
samplesnames.append("ST_tch_sd")
samplesnames.append("ST_tch_p_sd")
samplesnames.append("ST_sch")
samplesnames.append("ST_tW")
samplesnames.append("WJets")
samplesnames.append("VV")
samplesnames.append("TT")
samplesnames.append("TT_sd")
samplesnames.append("DYJets")
samplesnames.append("DDQCD")

samples = []
for sn in samplesnames:
    samples.append(samplesDictionary[sn])

channels = ["muon", "electron"]
systs = ["q2","psq2", "hdamp"]
vs = ["Up", "Down"]

newsysts = {"q2":["q2TT","q2ST_tch","q2ST_tch_sd", "q2WJets", "q2DYJets", "q2VV"], "psq2":["psq2TT","psq2ST_tch"], "hdamp":["hdampST_tch","hdampTT"]}
systsamples  =  {"q2":["TT","ST_tch","ST_tch_sd", "WJets", "DYJets", "VV"], "psq2":["TT","ST_tch"], "hdamp":["ST_tch","TT"]}
allsysts = ["", "btag", "mistag", "jes", "jer", "pu", "pdf_total"]
newlepSyst = {"lep":["lepMu","lepE"]}
for sy in systs:
    for s in newsysts[sy]:
        allsysts.append(s)
for s in newlepSyst["lep"]:
    allsysts.append(s)

def createSyst(versus, syst, newsyst, systsamp, verbose = True):
    for ch in channels:
        for s in samples:
            if ("QCDMu" in s.label)  and "electron" in ch : continue 
            pathch = opt.path+"/"+ch+"/"

            if s.label == "Data":continue
            src = s.label + "_" + ch + ".root"
            dst =  [s.label + "_" + ch + "_" + nSyst + versus + ".root" for nSyst in newsyst]
#            print dst

            if(verbose):
                for d in dst:
                    print "Copying " + pathch + src + " to "+ pathch + d
    
            [shutil.copyfile(pathch + src, pathch + d) for d in dst]

        for sysa in systsamp:
            s = samplesDictionary[sysa]
            srcSys = s.label + "_" + ch + "_" + syst + versus + ".root"  
            dstSys = s.label + "_" + ch + "_" + syst + s.label + versus + ".root" 
            
            if verbose: print "Copying " + pathch + srcSys + " to "+ pathch + dstSys
            shutil.copyfile(pathch + srcSys, pathch + dstSys)

def createLepSyst(versus, newlepsyst, verbose = True):
    for ch in channels:
        for s in samples:
            if ("QCDMu" in s.label)  and "electron" in ch : continue 
            pathch = opt.path+"/"+ch+"/"

            if s.label == "Data":continue
            src = s.label + "_" + ch + ".root"
            dst =  [s.label + "_" + ch + "_" + nlepSyst + versus + ".root" for nlepSyst in newlepsyst]
#            print dst

            if(verbose):
                for d in dst:
                    print "Copying " + pathch + src + " to "+ pathch + d
    
            [shutil.copyfile(pathch + src, pathch + d) for d in dst]

            if ch == "muon":
                srcSys = s.label + "_" + ch + "_lep" + versus + ".root"  
                dstSys = s.label + "_" + ch + "_lepMu" + versus + ".root" 
            else:
                srcSys = s.label + "_" + ch + "_lep" + versus + ".root"  
                dstSys = s.label + "_" + ch + "_lepE" + versus + ".root" 
            if verbose: print "Copying " + pathch + srcSys + " to "+ pathch + dstSys
            shutil.copyfile(pathch + srcSys, pathch + dstSys)

def invertname(samp_label, syst_, versus):
    for ch in channels:
        pathch = opt.path+"/"+ch+"/"
        srcSys = samp_label + "_" + syst_ + versus + "_" + ch + ".root"  
        dstSys = samp_label + "_" + ch + "_" + syst_ + versus + ".root" 
        shutil.copyfile(pathch + srcSys, pathch + dstSys)
        print "Renaming " + pathch + srcSys + " in " + pathch + dstSys

def rmdircont(namedir):
    os.system("rm " + opt.path + "/" + namedir + "/*.root")
    print "Deleting all previous files from " + opt.path + "/" + namedir + "/"

def copyall():
    for ch in channels:
        pathch = opt.path+"/"+ch+"/"
        pathout = opt.path+"/fit/"
        os.system("cp " + pathch + "*.root " + pathout)  
        print "Copying all root files from " + pathch + " to " + pathout
 
def merge(syst_, versus):
    for ch in channels:
        pathin = opt.path+"/"+ch+"/"
        pathout = opt.path+"/fit/"
        if syst_=="":
            srcST = "ST_tch_sd_" + ch + ".root"
            srcSTp = "ST_tch_p_sd_" + ch + ".root"
            srcTT = "TT_sd_" + ch + ".root"
        else:
            srcST = "ST_tch_sd_" + ch + "_" + syst_ + versus + ".root"
            srcSTp = "ST_tch_p_sd_" + ch + "_" + syst_ + versus + ".root"
            srcTT = "TT_sd_" + ch + "_" + syst_ + versus + ".root"
        os.system("hadd -f " + pathout + srcST + " " + pathin + srcST + " " + pathin + srcTT + " " + pathin + srcSTp)

def cureNaN(ch, sample, syst_, versus):
    pathin = opt.path+"/"+ch+"/"
    pathout = opt.path+"/fit/"
    os.system("cp " + pathin + sample + "_" + ch + ".root " + pathout + sample + "_" + ch + "_" + syst_ + versus + ".root")
    print "Copying "  + pathin + sample + "_" + ch + ".root" + " to " + pathout + sample + "_" + ch + "_"  + syst_ + versus + ".root"

samp = ["ST_tch", "TT"]
syst_name = ["hdamp","psq2"]
for v in vs:
    for sa in samp:
        for sy in syst_name:
            invertname(sa, sy, v)

for v in vs:
    for s in systs:
        createSyst(v, s, newsysts[s], systsamples[s])

for v in vs:
    createLepSyst(v, newlepSyst["lep"])

rmdircont("fit")
copyall()

for sy in allsysts:
    for v in vs:
        merge(sy, v)

for v in vs:
    for ch in channels:   
        cureNaN(ch, "ST_tW", "pdf_total", v)
        cureNaN(ch, "ST_sch", "pdf_total", v)
    cureNaN("electron", "DDQCD", "pdf_total", v)


 
