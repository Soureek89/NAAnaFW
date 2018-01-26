#!/bin/env python


########################################
###  Author: Annapaola de Cosa             
###  Code to produce pre e post fit plots  
###  from combine output (mlfit.root)
########################################

#Cooment and uncomment the different folder to produce the plots for muons, electrons, and/or their combination

from plots.plotUtils import *

bkgs_m2j1tcr = ["TT", "DDQCD", "ST_tch", "ST_tch_sd", "WJets","VV","DYJets","ST_tW","ST_sch"]
bkgs_m2j1tsr = ["TT", "DDQCD", "ST_tch", "ST_tch_sd", "WJets","VV","DYJets","ST_tW","ST_sch"]
bkgs_m3j1tcr = ["TT", "DDQCD", "ST_tch", "ST_tch_sd", "WJets","VV","DYJets","ST_tW","ST_sch"]
bkgs_m3j1tsr = ["TT", "DDQCD", "ST_tch", "ST_tch_sd", "WJets","VV","DYJets","ST_tW","ST_sch"]
bkgs_m3j2t   = ["TT", "DDQCD", "ST_tch", "ST_tch_sd", "WJets","VV","DYJets","ST_tW","ST_sch"]

regions_mu = {
    "m2j1tcr":("h1D_2j1t_cr_3838","muon_2j1t_central"),
    "m2j1tsr":("h1D_2j1t_sr_3838","muon_2j1t_forward"),
    "m3j1tcr":("h1D_3j1t_cr_3838","muon_3j1t_central"),
    "m3j1tsr":("h1D_3j1t_sr_3838","muon_3j1t_forward"),
    "m3j2t":("h_3j2t_BDT_ST_vs_TT_3j2t", "muon_3j2t"),
    }

regions_ele = {
    "e2j1tcr":("h1D_2j1t_cr_3838","electron_2j1t_central"),
    "e2j1tsr":("h1D_2j1t_sr_3838","electron_2j1t_forward"),
    "e3j1tcr":("h1D_3j1t_cr_3838","electron_3j1t_central"),
    "e3j1tsr":("h1D_3j1t_sr_3838","electron_3j1t_forward"),
    "e3j2t":("h_3j2t_BDT_ST_vs_TT_3j2t_E", "electron_3j2t"),
}

fitPhase = "fit_s"
folder = "muon"
if folder=="muon":    
    plot(folder, fitPhase, regions_mu["m2j1tcr"], bkgs_m2j1tcr)
    plot(folder, fitPhase, regions_mu["m2j1tsr"], bkgs_m2j1tsr)
    plot(folder, fitPhase, regions_mu["m3j1tcr"], bkgs_m3j1tcr)
    plot(folder, fitPhase, regions_mu["m3j1tsr"], bkgs_m3j1tsr)
    plot(folder, fitPhase, regions_mu["m3j2t"], bkgs_m3j2t)
#folder = "electron"
if folder=="electron":
    plot(folder, fitPhase, regions_ele["e2j1tcr"], bkgs_m2j1tcr)
    plot(folder, fitPhase, regions_ele["e2j1tsr"], bkgs_m2j1tsr)
    plot(folder, fitPhase, regions_ele["e3j1tcr"], bkgs_m3j1tcr)
    plot(folder, fitPhase, regions_ele["e3j1tsr"], bkgs_m3j1tsr)
    plot(folder, fitPhase, regions_ele["e3j2t"], bkgs_m3j2t)

fitPhase = "fit_s"
#folder = "mu_ele"
if folder=="mu_ele":    
    plot(folder, fitPhase, regions_mu["m2j1tcr"], bkgs_m2j1tcr)
    plot(folder, fitPhase, regions_mu["m2j1tsr"], bkgs_m2j1tsr)
    plot(folder, fitPhase, regions_mu["m3j1tcr"], bkgs_m3j1tcr)
    plot(folder, fitPhase, regions_mu["m3j1tsr"], bkgs_m3j1tsr)
    plot(folder, fitPhase, regions_mu["m3j2t"], bkgs_m3j2t)
    plot(folder, fitPhase, regions_ele["e2j1tcr"], bkgs_m2j1tcr)
    plot(folder, fitPhase, regions_ele["e2j1tsr"], bkgs_m2j1tsr)
    plot(folder, fitPhase, regions_ele["e3j1tcr"], bkgs_m3j1tcr)
    plot(folder, fitPhase, regions_ele["e3j1tsr"], bkgs_m3j1tsr)
    plot(folder, fitPhase, regions_ele["e3j2t"], bkgs_m3j2t)
