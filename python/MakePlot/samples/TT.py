from utils import *

branch_frac_sd = 0.002

TTbar = sample()
TTbar.files = outlist (d,"TT")
#TTbar.files = TTbar.files[:1]
TTbar.jpref = jetLabel 
TTbar.jp = jetLabel 
TTbar.skimEff = 1.
TTbar.sigma = 831.76*(1.-branch_frac_sd)*(1.-branch_frac_sd)
TTbar.color = ROOT.kOrange
TTbar.style = 1
TTbar.fill = 1001
TTbar.leglabel = "t#bar{t}"
TTbar.label = "TT"

TTbar_psq2Down = sample()
TTbar_psq2Down.files = outlist (d,"TT_psq2Down")
#TTbar_psq2Down.files = TTbar_psq2Down.files[:1]
TTbar_psq2Down.jpref = jetLabel 
TTbar_psq2Down.jp = jetLabel 
TTbar_psq2Down.skimEff = 1.
TTbar_psq2Down.sigma = 831.76*(1.-branch_frac_sd)*(1.-branch_frac_sd)
TTbar_psq2Down.color = ROOT.kOrange
TTbar_psq2Down.style = 1
TTbar_psq2Down.fill = 1001
TTbar_psq2Down.leglabel = "t#bar{t}"
TTbar_psq2Down.label = "TT_psq2Down"

TTbar_psq2Up = sample()
TTbar_psq2Up.files = outlist (d,"TT_psq2Up")
#TTbar_psq2Up.files = TTbar_psq2Up.files[:1]
TTbar_psq2Up.jpref = jetLabel 
TTbar_psq2Up.jp = jetLabel 
TTbar_psq2Up.skimEff = 1.
TTbar_psq2Up.sigma = 831.76*(1.-branch_frac_sd)*(1.-branch_frac_sd)
TTbar_psq2Up.color = ROOT.kOrange
TTbar_psq2Up.style = 1
TTbar_psq2Up.fill = 1001
TTbar_psq2Up.leglabel = "t#bar{t}"
TTbar_psq2Up.label = "TT_psq2Up"

TTbar_hdampDown = sample()
TTbar_hdampDown.files = outlist (d,"TT_hdampDown")
#TTbar_hdampDown.files = TTbar_hdampDown.files[:1]
TTbar_hdampDown.jpref = jetLabel 
TTbar_hdampDown.jp = jetLabel 
TTbar_hdampDown.skimEff = 1.
TTbar_hdampDown.sigma = 831.76*(1.-branch_frac_sd)*(1.-branch_frac_sd)
TTbar_hdampDown.color = ROOT.kOrange
TTbar_hdampDown.style = 1
TTbar_hdampDown.fill = 1001
TTbar_hdampDown.leglabel = "t#bar{t}"
TTbar_hdampDown.label = "TT_hdampDown"

TTbar_hdampUp = sample()
TTbar_hdampUp.files = outlist (d,"TT_hdampUp")
#TTbar_hdampUp.files = TTbar_hdampUp.files[:1]
TTbar_hdampUp.jpref = jetLabel 
TTbar_hdampUp.jp = jetLabel 
TTbar_hdampUp.skimEff = 1.
TTbar_hdampUp.sigma = 831.76*(1.-branch_frac_sd)*(1.-branch_frac_sd)
TTbar_hdampUp.color = ROOT.kOrange
TTbar_hdampUp.style = 1
TTbar_hdampUp.fill = 1001
TTbar_hdampUp.leglabel = "t#bar{t}"
TTbar_hdampUp.label = "TT_hdampUp"

TTbar_sd = sample()
TTbar_sd.files = outlist (d,"TT_sd")
#TTbar.files = TTbar.files[:1]
TTbar_sd.jpref = jetLabel 
TTbar_sd.jp = jetLabel 
TTbar_sd.skimEff = 1.
TTbar_sd.sigma = 831.76*2*branch_frac_sd*(1.-branch_frac_sd)
TTbar_sd.color = ROOT.kOrange
TTbar_sd.style = 1
TTbar_sd.fill = 1001
TTbar_sd.leglabel = "t#bar{t}_{sd}"
TTbar_sd.label = "TT_sd"

#grouping up the components
TT = sample()
TT.leglabel = "t#bar{t}"
TT.color = ROOT.kOrange
TT.label = "TT"
TT.skimEff = 1.
TT.fill = 1001
TT.style = 1
TT.components = [TTbar]

TT_psq2Down = sample()
TT_psq2Down.leglabel = "t#bar{t}"
TT_psq2Down.color = ROOT.kOrange
TT_psq2Down.label = "TT_psq2Down"
TT_psq2Down.skimEff = 1.
TT_psq2Down.fill = 1001
TT_psq2Down.style = 1
TT_psq2Down.components = [TTbar_psq2Down]

TT_psq2Up = sample()
TT_psq2Up.leglabel = "t#bar{t}"
TT_psq2Up.color = ROOT.kOrange
TT_psq2Up.label = "TT_psq2Up"
TT_psq2Up.skimEff = 1.
TT_psq2Up.fill = 1001
TT_psq2Up.style = 1
TT_psq2Up.components = [TTbar_psq2Up]

TT_hdampDown = sample()
TT_hdampDown.leglabel = "t#bar{t}"
TT_hdampDown.color = ROOT.kOrange
TT_hdampDown.label = "TT_hdampDown"
TT_hdampDown.skimEff = 1.
TT_hdampDown.fill = 1001
TT_hdampDown.style = 1
TT_hdampDown.components = [TTbar_hdampDown]

TT_hdampUp = sample()
TT_hdampUp.leglabel = "t#bar{t}"
TT_hdampUp.color = ROOT.kOrange
TT_hdampUp.label = "TT_hdampUp"
TT_hdampUp.skimEff = 1.
TT_hdampUp.fill = 1001
TT_hdampUp.style = 1
TT_hdampUp.components = [TTbar_hdampUp]

TT_sd = sample()
TT_sd.leglabel = "t#bar{t}_{sd}"
TT_sd.color = ROOT.kOrange
TT_sd.label = "TT_sd"
TT_sd.skimEff = 1.
TT_sd.fill = 1001
TT_sd.style = 1
TT_sd.components = [TTbar_sd]
