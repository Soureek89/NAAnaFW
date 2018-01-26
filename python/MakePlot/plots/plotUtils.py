#!/bin/env python
########################################
###  Author: Annapaola de Cosa             
###  Code to produce pre e post fit plots  
###  from combine output (mlfit.root)
########################################
import sys
import os
import ROOT
from plots.services import Histo, Stack, Legend
from samples.toPlot import samples
import tdrstyle, CMS_lumi

ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)
from ROOT import SetOwnership

    
def raiseExep(s):
    raise RuntimeError('Failed to retrieve %s' %s )

def getHisto(samp):
    h = ROOT.gDirectory.Get(samp)
    if not isinstance(h, ROOT.TH1):
        raise RuntimeError('Failed to retrieve %s' %samp )
    return h
                 
def getSamples(folder, fitPhase, region, samps):
    ifilename= "/afs/cern.ch/work/a/adeiorio/CMSSW_7_4_7/src/ttDM/stat/"+folder+"/mlfit.root"
    ifile = ROOT.TFile.Open(ifilename)
    #print "-----> Opening input file ", ifilename
    ROOT.gDirectory.cd("shapes_"+fitPhase+"/"+region)
    #print "--- Looking into directory: " + "shapes_"+fitPhase+"/"+region
    histos = { s: getHisto(s) for s in samps}
    [h.Scale(h.GetBinWidth(1)) for h in histos.itervalues()]
    #print "Histo type after getting: ", histos["TT"]
    return histos

def getData(channel, var):
    ifilename= "/afs/cern.ch/work/a/adeiorio/CMSSW_8_0_26_patch1/src/Analysis/NAAnaFW/bin/Plot/"+channel+"/Data_"+channel+".root"
    print " file name data ", ifilename
    ifile = ROOT.TFile.Open(ifilename)
    h = ROOT.gDirectory.Get(var)
    leg = Legend()
    if not isinstance(h, ROOT.TH1):
        raise RuntimeError('Failed to retrieve Data')
    h_data = Histo.fromTH1(h)
    h_data.SetStyle(ROOT.kBlack, samples["Data"].style, samples["Data"].fill)
    leg.AddEntry(h_data, samples["Data"].leglabel, "lp")
    return h_data,leg
    
def setSampStyle(histos):
    histos_new = { s: Histo.fromTH1(h) for (s,h) in histos.iteritems()}
    [ h.SetStyle(samples[s].color, samples[s].style, samples[s].fill) if s in samples.keys() else  raiseExep(s) for (s,h) in histos_new.iteritems() ]
    return histos_new

    
def doStack(histos, leg):
    stack = Stack("met","met")
    h_sig = ROOT.TH1F()
    if "ST_tch_sd" in samples:
        h_sig = histos["ST_tch_sd"]
        h_sig.Scale(100)
        leg.AddEntry(h_sig, "t_sd(*100)", "l")
    for s,h in histos.iteritems():
        if not (samples[s].leglabel=="t, t-ch_sd"):
            stack.Add(h)
    for s,h in reversed(sorted(histos.iteritems())):
        print samples[s].leglabel
        if not (samples[s].leglabel=="t, t-ch_sd"):
            leg.AddEntry(h, samples[s].leglabel, "f")
    return stack,leg,h_sig

def getStack(folder, fitPhase, region, samps, leg):
    histos =  setSampStyle(getSamples(folder, fitPhase, region, samps))
    print "+ TT integral prefit: ", histos["TT"].Integral()
    stack,leg,h_sig = doStack(histos, leg)
    return stack,leg,h_sig


def plot(folder, fitPhase, region, samps, lumi = "35.89"):

    ### Setting canvas and pad style
    tdrstyle.setTDRStyle();
    H=600
    W=700
    #T = 0.08*H
    #B = 0.12*H
    L = 0.12*W
    R = 0.08*W
    c1 = ROOT.TCanvas("c1","c1",50,50,W,H)
    #c1 = ROOT.TCanvas()
    c1.SetFillColor(0)
    c1.SetBorderMode(0)
    c1.SetFrameFillStyle(0)
    c1.SetFrameBorderMode(0)
    c1.SetLeftMargin(0.12 )
    c1.SetRightMargin(0.9)
    c1.SetTopMargin( 1)
    c1.SetBottomMargin(-1)
    c1.SetTickx(1)
    c1.SetTicky(1)
    c1.cd()
    unit = " fb^{-1}"
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Preliminary"
    lumi_sqrtS = str(lumi)+ unit +" (13 TeV)"
    iPeriod = 0
    iPos = 11
    # writing the lumi information and the CMS "logo"

    pad1= ROOT.TPad("pad", "pad", 0, 0.31 , 1, 1)
    SetOwnership(pad1, 0)
    pad1.SetTopMargin(0.1)
    pad1.SetBottomMargin(0)
    pad1.SetLeftMargin(0.12)
    pad1.SetRightMargin(0.05)
    pad1.SetBorderMode(0)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.Draw()
    pad1.cd()

    ### Creating and drawing TStack and Legend
    print "region " , region[0]
    tmp = ROOT.TH1F()
    h_data = Histo.fromTH1(tmp)
    if(region[1].startswith("muon")):
        h_data,leg=getData("muon",region[0])
    if(region[1].startswith("electron")):
        h_data,leg=getData("electron",region[0])
    xmin=h_data._h.GetXaxis().GetXmin()
    xmax=h_data._h.GetXaxis().GetXmax()

    stack,leg,h_sig = getStack(folder, fitPhase, region[1], samps, leg)
    stack.SetMaximum(stack.GetMaximum()*1.55)
    stack.DrawStack()
    h_data.Draw("eSAMEpx0")
    leg.Draw("same")
    if(h_sig.Integral()>0):
        h_sig._h.SetOption("hist same")
        h_sig._h.SetLineStyle(9)
        h_sig._h.SetLineColor(ROOT.kBlue)
        h_sig._h.SetLineWidth(3)
#        hsig.SetFillColor(0)
        h_sig._h.SetMarkerSize(0.)
        h_sig._h.SetMarkerColor(ROOT.kBlue)
        h_sig.Draw("same L")
#        h_sig.Draw("same")
    
    CMS_lumi.CMS_lumi(pad1, lumi_sqrtS, iPos)
    
    c1.cd()
    h_ratio = ROOT.TH1F()
    ratio = ROOT.TH1F()
    h_ratio = stack.GetLast()
    pad2= ROOT.TPad("pad2", "pad2", 0, 0.01 , 1, 0.30)
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.45)
    pad2.SetLeftMargin(0.12)
    pad2.SetRightMargin(0.05)
    ROOT.gStyle.SetHatchesSpacing(2)
    ROOT.gStyle.SetHatchesLineWidth(2)
    c1.cd()
    pad2.Draw()
    pad2.cd()
    ratio = h_data._h.Clone("ratio")
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetMaximum(1.5)
    ratio.SetMinimum(0.5)
    ratio.Sumw2()
    ratio.SetStats(0)
    ratio.Divide(h_ratio)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(0.9)
    ratio.Draw("epx0e0")
    ratio.SetTitle("")

    f1 = ROOT.TLine(xmin, 1., xmax,1.)
    f1.SetLineColor(ROOT.kBlack)
    f1.SetLineStyle(ROOT.kDashed)
    f1.Draw("same")
    ratio.GetYaxis().SetTitle("Data / Bkg")
    ratio.GetYaxis().SetNdivisions(503)
    ratio.GetXaxis().SetLabelFont(42)
    ratio.GetYaxis().SetLabelFont(42)
    ratio.GetXaxis().SetTitleFont(42)
    ratio.GetYaxis().SetTitleFont(42)
    ratio.GetXaxis().SetTitleOffset(1.1)
    ratio.GetYaxis().SetTitleOffset(0.35)
    ratio.GetXaxis().SetLabelSize(0.15)
    ratio.GetYaxis().SetLabelSize(0.15)
    ratio.GetXaxis().SetTitleSize(0.16)
    ratio.GetYaxis().SetTitleSize(0.16)
    ratio.GetYaxis().SetRangeUser(0.5,1.5)
    ratio.GetXaxis().SetTitle("")
    ratio.GetXaxis().SetLabelOffset(0.04)
    ratio.GetYaxis().SetLabelOffset(0.01)
    c1.cd()
    c1.Update()
    ROOT.TGaxis.SetMaxDigits(3)
    c1.RedrawAxis()
    path = "./post_fit_"+folder+"/"
    if not os.path.isdir(path):
        os.makedirs(path)
    pdfname = fitPhase+"_" + region[1]+".pdf"   
    pngname = fitPhase+"_" + region[1]+".png"   
    c1.Print(path+pdfname)
    c1.Print(path+pngname)
 
