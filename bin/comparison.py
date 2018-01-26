import os, commands
import optparse
from ROOT import *

usage = 'python checksyst.py'
parser = optparse.OptionParser(usage)
parser.add_option('--pathin', dest='pathin', type='string', default='res', help='')
parser.add_option('--pathout', dest='pathout', type='string', default='hdamp', help='')
parser.add_option('-L', '--lep', dest='lep', type='string', default='muon', help='Default checks integrals of muon systematics')
(opt, args) = parser.parse_args()

filenames = []
'''
filenames.append("hdampup")
filenames.append("hdampdown")
filenames.append("nominal")
'''
#filenames.append("ST_tch_hdampUp_muon")
#filenames.append("ST_tch_hdampDown_muon")
filenames.append("ST_tch_psq2Up_muon")
filenames.append("ST_tch_psq2Down_muon")
filenames.append("ST_tch_muon")
#filenames.append("TT_hdampUp")
#filenames.append("TT_hdampDown")
#filenames.append("TT_psq2Up")
#filenames.append("TT_psq2Down")
#filenames.append("TT")

histos = []
'''
histos.append("h_nJets")
#histos.append("h_nJets_cr")
#histos.append("h_nJets_sr")
histos.append("h_nJets_0b")
#histos.append("h_nJets_0b_cr")
#histos.append("h_nJets_0b_sr")
histos.append("h_nJets_1b")
#histos.append("h_nJets_1b_cr")
#histos.append("h_nJets_1b_sr")
histos.append("h_nJets_2b")
#histos.append("h_nJets_2b_cr")
#histos.append("h_nJets_2b_sr")
histos.append("h_nbJets")
#histos.append("h_nbJets_cr")
#histos.append("h_nbJets_sr")
histos.append("h_2j1t_topMass")
histos.append("h_2j1t_mtwcut_topMass")
histos.append("h_2j1t_mtwcut_sr_topMass")
histos.append("h_2j1t_etajprime")
'''
histos.append("h1D_2j1t_cr_3838")
histos.append("h1D_2j1t_sr_3838")
histos.append("h1D_3j1t_cr_3838")
histos.append("h1D_3j1t_sr_3838")
histos.append("h_3j2t_BDT_ST_vs_TT_3j2t")
'''
histos.append("h_2j1t_BDT_ST_vs_TT_mtweta_mtw_G_50_AND_etajprime_L_2p5")
histos.append("h_2j1t_BDT_ST_vs_VJ_mtweta_mtw_G_50_AND_etajprime_L_2p5")
histos.append("h_2j1t_BDT_STsd_vs_ST_sr_mtw_G_50_AND_etajprime_G_2p5")
histos.append("h_2j1t_BDT_STsd_vs_All_sr_mtw_G_50_AND_etajprime_G_2p5")
#histos.append("h_2j1t_BDT_ST_vs_TT_mtweta_mtw_G_50_AND_etajprime_G_2p4")
#histos.append("h_3j1t_BDT_ST_vs_TT_mtweta_mtw_G_50_AND_etajprime_L_2p5")
#histos.append("h_3j1t_BDT_ST_vs_TT_mtweta_mtw_G_50_AND_etajprime_G_2p4")
histos.append("h_2j1t_topMass_mtw_G_50_AND_etajprime_L_2p5")
histos.append("h_2j1t_mtw_mtw_G_50_AND_etajprime_L_2p5")
histos.append("h_2j1t_mlb_mtw_G_50_AND_etajprime_L_2p5")
histos.append("h_2j1t_nextrajets_mtw_G_50_AND_etajprime_L_2p5")
'''
pathin = opt.pathin
if opt.pathin != "res":
     pathin == opt.pathin

pathout = opt.pathout
if opt.pathout != "hdamp":
     pathout == opt.pathout
infile = []

gStyle.SetPalette(1)
gStyle.SetOptStat(0)
gROOT.SetBatch()        # don't pop up canvases

scale = False
rebin = False

for fil in filenames:
    infile.append(TFile.Open(pathin+"/"+fil+".root"))
for hist in histos:
     h1 = TH1F()
     h2 = TH1F()
     h3 = TH1F()
     c1 = TCanvas(str(hist),"c1",50,50,700,600)
     leg = TLegend(0.7,0.7,0.9,0.9)
     for inf in infile:
          inf.cd()
          tmp = (TH1F)(inf.Get(hist))
          if scale:
               norm = (TH1F)(inf.Get("h_cutFlow"))
               print "norm",norm.GetBinContent(0)
          filename = (str)(inf.GetName())
#          print "filename" , filename
          if("Up" in filename):
               h1 = tmp.Clone("nominal") 
               if scale:
                    h1.Scale(1./norm.GetBinContent(0))
               if rebin:
                    h1.Rebin(4)
               h1.SetTitle(str(hist))
#               h1.SetName(str(hist))
          elif("Down" in filename):
               h2 = tmp.Clone("nominal") 
               if scale:
                    h2.Scale(1./norm.GetBinContent(0))
               if rebin:
                    h2.Rebin(4)
               h2.SetTitle("down")
          else:
               h3 = tmp.Clone("nominal")
#               print "i0ntegral ", h3.Integral()
               if scale:
                    h3.Scale(1./norm.GetBinContent(0))
#               h3.SetName("nominal")
               if rebin:
                    h3.Rebin(4)
               h3.SetTitle("nominal")
          tmp.Reset("ICES")
     c1.cd()
     maximum = max(h1.GetMaximum(),h2.GetMaximum(),h3.GetMaximum())
     minimum = min(h1.GetMinimum(),h2.GetMinimum(),h3.GetMinimum())
     h1.SetMaximum(maximum*1.5)
     h1.SetMinimum(0.0)
     h1.Draw("E")
     print "int h1 up",h1.Integral()
     h1.SetLineColor(kRed)
     h2.Draw("ESAME")
     h2.SetLineColor(kGreen)
     print "int h2 down",h2.Integral()
     h3.Draw("HIST SAME")
     h3.SetLineColor(kBlue)
     print "int h3 nominal",h3.Integral()
     leg.AddEntry(h1, "up", "l")
     leg.AddEntry(h3, "nominal", "l")
     leg.AddEntry(h2, "down", "l")
#     gPad.BuildLegend()
#     c1.SetTitle(str(hist))
     leg.Draw("SAME")
     c1.Update()
     c1.Print(pathout+"/"+hist+".pdf")
     c1.Print(pathout+"/"+hist+".png")
     c1.Print(pathout+"/"+hist+".root")
