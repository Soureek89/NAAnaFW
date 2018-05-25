import os, commands
import optparse
from ROOT import * 
from variabile import variabile 
from Analysis.NAAnaFW.MakePlot.samples.VJets import *
from Analysis.NAAnaFW.MakePlot.samples.SingleTop import *
from Analysis.NAAnaFW.MakePlot.samples.TT import *
from Analysis.NAAnaFW.MakePlot.samples.VV import *
from Analysis.NAAnaFW.MakePlot.samples.toPlot import samples as samplesDictionary

#definizioni funzioni 
def checksyst(s, lepch, syst, histoname):
     print s.label+","+lepch+","+str(syst)+","+histoname
     infilenom = TFile.Open(pathin+"/"+s.label+"_"+str(lepch)+".root")
     infileDown = TFile.Open(pathin+"/"+s.label+"_"+str(lepch)+"_"+str(syst)+"Down.root")
     infileUp = TFile.Open(pathin+"/"+s.label+"_"+str(lepch)+"_"+str(syst)+"Up.root")
     
     infilenom.cd()
     h_nom = (TH1F)(infilenom.Get(histoname))
     int_nom = h_nom.Integral()

     infileUp.cd()
     h_Up = (TH1F)(infileUp.Get(histoname)) 
     int_Up = h_Up.Integral()
     infileDown.cd()
     h_Down = (TH1F)(infileDown.Get(histoname))
     int_Down = h_Down.Integral()

     return  str(int_Down)+";  "+str(int_nom)+";  "+str(int_Up)+"\n"

usage = 'python checksyst.py'
parser = optparse.OptionParser(usage)
parser.add_option('--pathin', dest='pathin', type='string', default='Plot/fit', help='')
parser.add_option('--pathout', dest='pathout', type='string', default='checksyst', help='')
parser.add_option('-L', '--lep', dest='lep', type='string', default='muon', help='Default checks integrals of muon systematics')
parser.add_option('-S', '--syst', dest='syst', type='string', default='all', help='Default checks integrals of all systematics')
parser.add_option('-P', '--process', dest='process', type='string', default = 'all', help="samples to add, according to the convention in 'script_rename.py'. Options are: 'All','ST','VJ','VV','QCDMu','QCDEle', or '_'+specific process, e.g. _ST_T_tch to add only the _ST_T_tch. Accepts also multiple processes, separated by comma, e.g. -P ST,_TT,VJ will select the V+jets samples and single top sample sets, as well as the one sample named TT.")
(opt, args) = parser.parse_args()

pathin = opt.pathin
if opt.pathin != "Plot/fit":
     pathin == opt.pathin

pathout = opt.pathout
if opt.pathout != "checksyst":
     pathout == opt.pathout

leptons=[]
if opt.lep != "muon":
     for lep in (opt.lep).split(","):
          leptons.append(lep)
else:
     leptons.append("muon")
     leptons.append("electron")
for lep in leptons:
     print lep

samplesnames = []
samples = []
if opt.process!="all":
    for sn in (opt.process).split(","): 
        samplesnames.append(sn)
else:
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
for sn in samplesnames:
    samples.append(samplesDictionary[sn])

components = []
for s in samples:
    components.extend(s.components)

systematics = []
if opt.syst!="all" and opt.syst!="noSyst":
    systematics.append((opt.syst).split(","))
elif opt.syst!="all" and opt.syst=="noSyst":
    systematics.append("") #di default per syst="" alla variabile si applica il peso w*w_nominal
else:
#    systematics.append("") #di default per syst="" alla variabile si applica il peso w*w_nominal
    systematics.append("btag")  
    systematics.append("mistag")
    systematics.append("pu")
#    '''
    systematics.append("lepMu")
    systematics.append("lepE")
    systematics.append("jes")
#    systematics.append("jer")
    systematics.append("q2ST_tch")
    systematics.append("q2ST_tch_sd")
    systematics.append("q2TT")
    systematics.append("q2DYJets")
    systematics.append("q2VV")
    systematics.append("q2WJets")
    systematics.append("psq2ST_tch")
    systematics.append("psq2TT")
    systematics.append("hdampST_tch")
    systematics.append("hdampTT")
    systematics.append("pdf_total")
#    '''

histonamesMu = []
histonamesMu.append("h1D_2j1t_cr_3838")
histonamesMu.append("h1D_2j1t_sr_3838")
histonamesMu.append("h1D_3j1t_cr_3838")
histonamesMu.append("h1D_3j1t_sr_3838")
histonamesMu.append("h_3j2t_BDT_ST_vs_TT_3j2t")

histonamesE = []
histonamesE.append("h1D_2j1t_cr_3838")
histonamesE.append("h1D_2j1t_sr_3838")
histonamesE.append("h1D_3j1t_cr_3838")
histonamesE.append("h1D_3j1t_sr_3838")
histonamesE.append("h_3j2t_BDT_ST_vs_TT_3j2t_E")

lumi=35.89

print "Deleting old text files"
os.system("rm checksyst/*")

for s in samples:
     for lep in leptons:
          filename = s.label +"_"+ lep + ".txt"
          f_out = open(pathout+"/"+filename, "a+")
          for syst in systematics:
               if (lep == "muon"):
                    for histoname in histonamesMu:
                         out = checksyst(s, lep, syst, histoname)
                         f_out.write(str(histoname)+"; "+str(syst)+"; "+out)
               if (lep == "electron"):
                    for histoname in histonamesE:
                         out = checksyst(s, lep, syst, histoname)
                         f_out.write(str(histoname)+"; "+str(syst)+"; "+out)
          f_out.close
