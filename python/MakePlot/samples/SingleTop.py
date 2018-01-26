######################################
#
# Annapaola de Cosa, January 2015
#
######################################


from utils import *

branch_frac_sd = 0.002

ST_T_tch = sample()
ST_T_tch.files = outlist (d,"ST_T_tch")
ST_T_tch.skimEff = 1
ST_T_tch.color = ROOT.kRed
#ST_T_tch.sigma =  70.29828*4.65068
#ST_T_tch.sigma =  216.99 #* 0.324
ST_T_tch.sigma =  136.02*(1.-branch_frac_sd) #* 0.324
ST_T_tch.jpref = jetLabel 
ST_T_tch.jp = jetLabel
ST_T_tch.style = 1
ST_T_tch.fill = 1001
ST_T_tch.leglabel = "ST"
ST_T_tch.label = "ST_T_tch"

ST_Tbar_tch = sample()
ST_Tbar_tch.files = outlist (d,"ST_Tbar_tch")
ST_Tbar_tch.skimEff = 1.
ST_Tbar_tch.color = ROOT.kRed
ST_Tbar_tch.sigma = 80.95*(1.-branch_frac_sd)
ST_Tbar_tch.jpref = jetLabel 
ST_Tbar_tch.jp = jetLabel
ST_Tbar_tch.style = 1
ST_Tbar_tch.fill = 1001
ST_Tbar_tch.leglabel = "ST"
ST_Tbar_tch.label = "ST_Tbar_tch"

ST_T_tch_sd = sample()
ST_T_tch_sd.files = outlist (d,"ST_T_tch_sd")
ST_T_tch_sd.skimEff = 1
ST_T_tch_sd.color = ROOT.kRed
ST_T_tch_sd.sigma =  136.02*branch_frac_sd 
ST_T_tch_sd.jpref = jetLabel 
ST_T_tch_sd.jp = jetLabel
ST_T_tch_sd.style = 1
ST_T_tch_sd.fill = 1001
ST_T_tch_sd.leglabel = "ST"
ST_T_tch_sd.label = "ST_T_tch_sd"

ST_Tbar_tch_sd = sample()
ST_Tbar_tch_sd.files = outlist (d,"ST_Tbar_tch_sd")
ST_Tbar_tch_sd.skimEff = 1.
ST_Tbar_tch_sd.color = ROOT.kRed
ST_Tbar_tch_sd.sigma = 80.95*branch_frac_sd 
ST_Tbar_tch_sd.jpref = jetLabel 
ST_Tbar_tch_sd.jp = jetLabel
ST_Tbar_tch_sd.style = 1
ST_Tbar_tch_sd.fill = 1001
ST_Tbar_tch_sd.leglabel = "ST"
ST_Tbar_tch_sd.label = "ST_Tbar_tch_sd"

ST_T_tch_p_sd = sample()
ST_T_tch_p_sd.files = outlist (d,"ST_T_tch_p_sd")
ST_T_tch_p_sd.skimEff = 1
ST_T_tch_p_sd.color = ROOT.kRed
ST_T_tch_p_sd.sigma =  0.62*(1-branch_frac_sd)
ST_T_tch_p_sd.jpref = jetLabel 
ST_T_tch_p_sd.jp = jetLabel
ST_T_tch_p_sd.style = 1
ST_T_tch_p_sd.fill = 1001
ST_T_tch_p_sd.leglabel = "ST"
ST_T_tch_p_sd.label = "ST_T_tch_p_sd"

ST_Tbar_tch_p_sd = sample()
ST_Tbar_tch_p_sd.files = outlist (d,"ST_Tbar_tch_p_sd")
ST_Tbar_tch_p_sd.skimEff = 1.
ST_Tbar_tch_p_sd.color = ROOT.kRed
ST_Tbar_tch_p_sd.sigma = 0.36*(1-branch_frac_sd)
ST_Tbar_tch_p_sd.jpref = jetLabel 
ST_Tbar_tch_p_sd.jp = jetLabel
ST_Tbar_tch_p_sd.style = 1
ST_Tbar_tch_p_sd.fill = 1001
ST_Tbar_tch_p_sd.leglabel = "ST"
ST_Tbar_tch_p_sd.label = "ST_Tbar_tch_p_sd"

ST_T_sch = sample()
ST_T_sch.files = outlist (d,"ST_T_sch")
ST_T_sch.skimEff = 1.0
ST_T_sch.color = 95 
ST_T_sch.sigma = 10.32
ST_T_sch.jpref = jetLabel 
ST_T_sch.jp = jetLabel
ST_T_sch.style = 1
ST_T_sch.fill = 1001
ST_T_sch.leglabel = "ST"
ST_T_sch.label = "ST_T_sch"

ST_Tbar_sch = sample()
ST_Tbar_sch.files = outlist (d,"ST_Tbar_sch")
ST_Tbar_sch.skimEff = 1.0
ST_Tbar_sch.color = 95 
ST_Tbar_sch.sigma = 3.97
ST_Tbar_sch.jpref = jetLabel 
ST_Tbar_sch.jp = jetLabel
ST_Tbar_sch.style = 1
ST_Tbar_sch.fill = 1001
ST_Tbar_sch.leglabel = "ST"
ST_Tbar_sch.label = "ST_Tbar_sch"

ST_T_tW = sample()
ST_T_tW.files = outlist (d,"ST_T_tW")
ST_T_tW.skimEff = 1.0
ST_T_tW.sigma = 35.6
ST_T_tW.jpref = jetLabel 
ST_T_tW.jp = jetLabel
ST_T_tW.color = ROOT.kOrange+2
ST_T_tW.style = 1
ST_T_tW.fill = 1001
ST_T_tW.leglabel = "ST"
ST_T_tW.label = "ST_T_tW"

ST_Tbar_tW = sample()
ST_Tbar_tW.files = outlist (d,"ST_Tbar_tW")
ST_Tbar_tW.skimEff = 1.
ST_Tbar_tW.sigma = 35.6
ST_Tbar_tW.jpref = jetLabel 
ST_Tbar_tW.jp = jetLabel
ST_Tbar_tW.color = ROOT.kOrange+2
ST_Tbar_tW.style = 1
ST_Tbar_tW.fill = 1001
ST_Tbar_tW.leglabel = "ST"
ST_Tbar_tW.label = "ST_Tbar_tW"

ST_T_tch_psq2Down = sample()
ST_T_tch_psq2Down.files = outlist (d,"ST_T_tch_psq2Down")
ST_T_tch_psq2Down.skimEff = 1.
ST_T_tch_psq2Down.sigma = 136.02*(1.-branch_frac_sd)
ST_T_tch_psq2Down.jpref = jetLabel
ST_T_tch_psq2Down.jp = jetLabel
ST_T_tch_psq2Down.color = ROOT.kRed
ST_T_tch_psq2Down.style = 1
ST_T_tch_psq2Down.fill = 1001
ST_T_tch_psq2Down.leglabel = ""
ST_T_tch_psq2Down.label = "ST_T_tch_psq2Down"

ST_Tbar_tch_psq2Down = sample()
ST_Tbar_tch_psq2Down.files = outlist (d,"ST_Tbar_tch_psq2Down")
ST_Tbar_tch_psq2Down.skimEff = 1.
ST_Tbar_tch_psq2Down.sigma = 136.02*(1.-branch_frac_sd)
ST_Tbar_tch_psq2Down.jpref = jetLabel
ST_Tbar_tch_psq2Down.jp = jetLabel
ST_Tbar_tch_psq2Down.color = ROOT.kRed
ST_Tbar_tch_psq2Down.style = 1
ST_Tbar_tch_psq2Down.fill = 1001
ST_Tbar_tch_psq2Down.leglabel = ""
ST_Tbar_tch_psq2Down.label = "ST_Tbar_tch_psq2Down"

ST_T_tch_psq2Up = sample()
ST_T_tch_psq2Up.files = outlist (d,"ST_T_tch_psq2Up")
ST_T_tch_psq2Up.skimEff = 1.
ST_T_tch_psq2Up.sigma = 136.02*(1.-branch_frac_sd)
ST_T_tch_psq2Up.jpref = jetLabel
ST_T_tch_psq2Up.jp = jetLabel
ST_T_tch_psq2Up.color = ROOT.kRed
ST_T_tch_psq2Up.style = 1
ST_T_tch_psq2Up.fill = 1001
ST_T_tch_psq2Up.leglabel = ""
ST_T_tch_psq2Up.label = "ST_T_tch_psq2Up"

ST_Tbar_tch_psq2Up = sample()
ST_Tbar_tch_psq2Up.files = outlist (d,"ST_Tbar_tch_psq2Up")
ST_Tbar_tch_psq2Up.skimEff = 1.
ST_Tbar_tch_psq2Up.sigma = 136.02*(1.-branch_frac_sd)
ST_Tbar_tch_psq2Up.jpref = jetLabel
ST_Tbar_tch_psq2Up.jp = jetLabel
ST_Tbar_tch_psq2Up.color = ROOT.kRed
ST_Tbar_tch_psq2Up.style = 1
ST_Tbar_tch_psq2Up.fill = 1001
ST_Tbar_tch_psq2Up.leglabel = ""
ST_Tbar_tch_psq2Up.label = "ST_Tbar_tch_psq2Up"

ST_T_tch_hdampDown = sample()
ST_T_tch_hdampDown.files = outlist (d,"ST_T_tch_hdampDown")
ST_T_tch_hdampDown.skimEff = 1.
ST_T_tch_hdampDown.sigma = 136.02*(1.-branch_frac_sd)
ST_T_tch_hdampDown.jpref = jetLabel
ST_T_tch_hdampDown.jp = jetLabel
ST_T_tch_hdampDown.color = ROOT.kRed
ST_T_tch_hdampDown.style = 1
ST_T_tch_hdampDown.fill = 1001
ST_T_tch_hdampDown.leglabel = ""
ST_T_tch_hdampDown.label = "ST_T_tch_hdampDown"

ST_Tbar_tch_hdampDown = sample()
ST_Tbar_tch_hdampDown.files = outlist (d,"ST_Tbar_tch_hdampDown")
ST_Tbar_tch_hdampDown.skimEff = 1.
ST_Tbar_tch_hdampDown.sigma = 136.02*(1.-branch_frac_sd)
ST_Tbar_tch_hdampDown.jpref = jetLabel
ST_Tbar_tch_hdampDown.jp = jetLabel
ST_Tbar_tch_hdampDown.color = ROOT.kRed
ST_Tbar_tch_hdampDown.style = 1
ST_Tbar_tch_hdampDown.fill = 1001
ST_Tbar_tch_hdampDown.leglabel = ""
ST_Tbar_tch_hdampDown.label = "ST_Tbar_tch_hdampDown"

ST_T_tch_hdampUp = sample()
ST_T_tch_hdampUp.files = outlist (d,"ST_T_tch_hdampUp")
ST_T_tch_hdampUp.skimEff = 1.
ST_T_tch_hdampUp.sigma = 136.02*(1.-branch_frac_sd)
ST_T_tch_hdampUp.jpref = jetLabel
ST_T_tch_hdampUp.jp = jetLabel
ST_T_tch_hdampUp.color = ROOT.kRed
ST_T_tch_hdampUp.style = 1
ST_T_tch_hdampUp.fill = 1001
ST_T_tch_hdampUp.leglabel = ""
ST_T_tch_hdampUp.label = "ST_T_tch_hdampUp"

ST_Tbar_tch_hdampUp = sample()
ST_Tbar_tch_hdampUp.files = outlist (d,"ST_Tbar_tch_hdampUp")
ST_Tbar_tch_hdampUp.skimEff = 1.
ST_Tbar_tch_hdampUp.sigma = 136.02*(1.-branch_frac_sd)
ST_Tbar_tch_hdampUp.jpref = jetLabel
ST_Tbar_tch_hdampUp.jp = jetLabel
ST_Tbar_tch_hdampUp.color = ROOT.kRed
ST_Tbar_tch_hdampUp.style = 1
ST_Tbar_tch_hdampUp.fill = 1001
ST_Tbar_tch_hdampUp.leglabel = ""
ST_Tbar_tch_hdampUp.label = "ST_Tbar_tch_hdampUp"


#Grouping up the components:
ST_tch = sample()
ST_tch.color = ROOT.kRed
ST_tch.style = 1
ST_tch.fill = 1001
ST_tch.leglabel = "t, t-ch"
ST_tch.label = "ST_tch"
ST_tch.components = [ST_T_tch,ST_Tbar_tch]

ST_tch_sd = sample()
ST_tch_sd.color = ROOT.kRed
ST_tch_sd.style = 1
ST_tch_sd.fill = 1001
ST_tch_sd.leglabel = "t, t-ch_sd"
ST_tch_sd.label = "ST_tch_sd"
ST_tch_sd.components = [ST_T_tch_sd,ST_Tbar_tch_sd]

ST_tch_p_sd = sample()
ST_tch_p_sd.color = ROOT.kRed
ST_tch_p_sd.style = 1
ST_tch_p_sd.fill = 1001
ST_tch_p_sd.leglabel = "t, t-ch_p_sd"
ST_tch_p_sd.label = "ST_tch_p_sd"
ST_tch_p_sd.components = [ST_T_tch_p_sd,ST_Tbar_tch_p_sd]

ST_sch = sample()
ST_sch.color = 95 
ST_sch.style = 1
ST_sch.fill = 1001
ST_sch.leglabel = "t, s-ch"
ST_sch.label = "ST_sch"
ST_sch.components = [ST_T_sch]

ST_tW = sample()
ST_tW.color = ROOT.kOrange+2
ST_tW.style = 1
ST_tW.fill = 1001
ST_tW.leglabel = "tW"
ST_tW.label = "ST_tW"
ST_tW.components = [ST_T_tW,ST_Tbar_tW]

ST_tch_psq2Down = sample()
ST_tch_psq2Down.color = ROOT.kRed
ST_tch_psq2Down.style = 1
ST_tch_psq2Down.fill = 1001
ST_tch_psq2Down.leglabel = "t, t-ch"
ST_tch_psq2Down.label = "ST_tch_psq2Down"
ST_tch_psq2Down.components = [ST_T_tch_psq2Down,ST_Tbar_tch_psq2Down]

ST_tch_psq2Up = sample()
ST_tch_psq2Up.color = ROOT.kRed
ST_tch_psq2Up.style = 1
ST_tch_psq2Up.fill = 1001
ST_tch_psq2Up.leglabel = "t, t-ch"
ST_tch_psq2Up.label = "ST_tch_psq2Up"
ST_tch_psq2Up.components = [ST_T_tch_psq2Up,ST_Tbar_tch_psq2Up]

ST_tch_hdampDown = sample()
ST_tch_hdampDown.color = ROOT.kRed
ST_tch_hdampDown.style = 1
ST_tch_hdampDown.fill = 1001
ST_tch_hdampDown.leglabel = "t, t-ch"
ST_tch_hdampDown.label = "ST_tch_hdampDown"
ST_tch_hdampDown.components = [ST_T_tch_hdampDown,ST_Tbar_tch_hdampDown]

ST_tch_hdampUp = sample()
ST_tch_hdampUp.color = ROOT.kRed
ST_tch_hdampUp.style = 1
ST_tch_hdampUp.fill = 1001
ST_tch_hdampUp.leglabel = "t, t-ch"
ST_tch_hdampUp.label = "ST_tch_hdampUp"
ST_tch_hdampUp.components = [ST_T_tch_hdampUp,ST_Tbar_tch_hdampUp]

