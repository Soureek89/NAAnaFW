
######################################
#
# Annapaola de Cosa, January 2015
#
######################################


from utils import *


WToLNu0J = sample()
WToLNu0J.files = outlist (d,"WToLNu0J")
WToLNu0J.skimEff = 1.
WToLNu0J.sigma =  49670.
WToLNu0J.jpref = jetLabel 
WToLNu0J.jp = jetLabel
WToLNu0J.color = ROOT.kGreen+2
WToLNu0J.style = 1
WToLNu0J.fill = 1001
WToLNu0J.leglabel = "WJets"
WToLNu0J.label = "WToLNu0J"

WToLNu1J = sample()
WToLNu1J.files = outlist (d,"WToLNu1J")
WToLNu1J.skimEff = 1.
WToLNu1J.sigma =  8264.
WToLNu1J.jpref = jetLabel 
WToLNu1J.jp = jetLabel
WToLNu1J.color = ROOT.kGreen+2
WToLNu1J.style = 1
WToLNu1J.fill = 1001
WToLNu1J.leglabel = "WJets"
WToLNu1J.label = "WToLNu1J"

WToLNu2J = sample()
WToLNu2J.files = outlist (d,"WToLNu2J")
WToLNu2J.skimEff = 1.
WToLNu2J.sigma =  2544.
WToLNu2J.jpref = jetLabel 
WToLNu2J.jp = jetLabel
WToLNu2J.color = ROOT.kGreen+2
WToLNu2J.style = 1
WToLNu2J.fill = 1001
WToLNu2J.leglabel = "WJets"
WToLNu2J.label = "WToLNu2J"

DYJetsToLL = sample()
DYJetsToLL.files = outlist (d,"DYJetsToLL")
DYJetsToLL.skimEff = 1.
#DYJetsToLL.sigma = 23.1363
DYJetsToLL.sigma =  1921.8*3.
DYJetsToLL.jpref = jetLabel 
DYJetsToLL.jp = jetLabel
DYJetsToLL.color = ROOT.kGreen
DYJetsToLL.style = 1
DYJetsToLL.fill = 1001
DYJetsToLL.leglabel = "DYJets"
DYJetsToLL.label = "DYJetsToLL"

WJets = sample()
#WJets.color = 6
WJets.color = ROOT.kGreen+2 
WJets.style = 1
WJets.fill = 1001
WJets.leglabel = "W + Jets"
WJets.label = "WJets"
WJets.components = [WToLNu0J,WToLNu1J,WToLNu2J]
#WJets.components = [WJets_HT100_200, WJets_HT200_400, WJets_HT400_600, WJets_HT600_800, WJets_HT800_1200, WJets_HT1200_2500, WJets_HT2500_Inf]
# WJets.components = [WJets_HT100_200, WJets_HT200_400,  WJets_HT800_1200, WJets_HT1200_2500 ]

DYJets = sample()
#DYJets.color = 6
DYJets.color = ROOT.kViolet-1
DYJets.style = 1
DYJets.fill = 1001
DYJets.leglabel = "Z + Jets"
DYJets.label = "DYJets"
DYJets.components = [DYJetsToLL]




