from utils import *

QCD_100_200 = sample()
QCD_100_200.files = outlist (d,"QCD_HT100_200")
QCD_100_200.jpref = jetLabel 
QCD_100_200.jp = jetLabel 
QCD_100_200.sf = 1
QCD_100_200.skimEff = 1.
QCD_100_200.sigma = 27990000
QCD_100_200.color = ROOT.kYellow -9
QCD_100_200.style = 1
QCD_100_200.fill = 1001
QCD_100_200.leglabel = "QCD"
QCD_100_200.label = "QCD_HT100_200"

QCD_200_300 = sample()
QCD_200_300.files = outlist (d,"QCD_HT200_300")
QCD_200_300.jpref = jetLabel 
QCD_200_300.jp = jetLabel 
QCD_200_300.sf = 1
QCD_200_300.skimEff = 1.
QCD_200_300.sigma = 1712000
QCD_200_300.color = ROOT.kYellow -9
QCD_200_300.style = 1
QCD_200_300.fill = 1001
QCD_200_300.leglabel = "QCD"
QCD_200_300.label = "QCD_HT200_300"


QCD_300_500 = sample()
QCD_300_500.files = outlist (d,"QCD_HT300_500")
QCD_300_500.jpref = jetLabel 
QCD_300_500.jp = jetLabel 
QCD_300_500.sigma = 347700
QCD_300_500.skimEff = 1
QCD_300_500.sf = 1
QCD_300_500.color = 9
QCD_300_500.style = 1
QCD_300_500.fill = 1001
QCD_300_500.leglabel = "QCD"
QCD_300_500.label = "QCD_HT300_500"


QCD_500_700 = sample()
QCD_500_700.files = outlist (d,"QCD_HT500_700")
QCD_500_700.jpref = jetLabel 
QCD_500_700.jp = jetLabel 
QCD_500_700.sf = 1
QCD_500_700.skimEff = 1.
QCD_500_700.sigma = 32100
QCD_500_700.color = ROOT.kYellow -9
QCD_500_700.style = 1
QCD_500_700.fill = 1001
QCD_500_700.leglabel = "QCD"
QCD_500_700.label = "QCD_HT500_700"


QCD_700_1000 = sample()
QCD_700_1000.files = outlist (d,"QCD_HT700_1000")
QCD_700_1000.jpref = jetLabel 
QCD_700_1000.jp = jetLabel 
QCD_700_1000.sf = 1
QCD_700_1000.skimEff = 1.
QCD_700_1000.sigma = 6831
QCD_700_1000.color = ROOT.kYellow -9
QCD_700_1000.style = 1
QCD_700_1000.fill = 1001
QCD_700_1000.leglabel = "QCD"
QCD_700_1000.label = "QCD_HT700_1000"


QCD_1000_1500 = sample()
QCD_1000_1500.files = outlist (d,"QCD_HT1000_1500")
QCD_1000_1500.jpref = jetLabel 
QCD_1000_1500.jp = jetLabel 
QCD_1000_1500.sf = 1
QCD_1000_1500.skimEff = 1.
QCD_1000_1500.sigma = 1207
QCD_1000_1500.color = ROOT.kYellow -9
QCD_1000_1500.style = 1
QCD_1000_1500.fill = 1001
QCD_1000_1500.leglabel = "QCD"
QCD_1000_1500.label = "QCD_HT1000_1500"

QCD_1500_2000 = sample()
QCD_1500_2000.files = outlist (d,"QCD_HT1500_2000")
QCD_1500_2000.jpref = jetLabel 
QCD_1500_2000.jp = jetLabel 
QCD_1500_2000.sf = 1
QCD_1500_2000.skimEff = 1.
QCD_1500_2000.sigma = 119.9
QCD_1500_2000.color = ROOT.kYellow -9
QCD_1500_2000.style = 1
QCD_1500_2000.fill = 1001
QCD_1500_2000.leglabel = "QCD"
QCD_1500_2000.label = "QCD_HT1500_2000"

QCD_2000_Inf = sample()
QCD_2000_Inf.files = outlist (d,"QCD_HT2000_Inf")
QCD_2000_Inf.jpref = jetLabel 
QCD_2000_Inf.jp = jetLabel 
QCD_2000_Inf.sf = 1
QCD_2000_Inf.skimEff = 1.
QCD_2000_Inf.sigma = 25.24
QCD_2000_Inf.color = ROOT.kYellow -9
QCD_2000_Inf.style = 1
QCD_2000_Inf.fill = 1001
QCD_2000_Inf.leglabel = "QCD"
QCD_2000_Inf.label = "QCD_HT2000_Inf"


QCD = sample()
#QCD.color = 4
#QCD.color = ROOT.kBlue+2
#QCD.color = 7
QCD.color = ROOT.kOrange +7
QCD.style = 1
QCD.fill = 1001
QCD.leglabel = "QCD"
QCD.label = "QCD"
QCD.components = [ QCD_100_200, QCD_200_300,  QCD_300_500, QCD_500_700, QCD_700_1000, QCD_1000_1500,  QCD_1500_2000,  QCD_2000_Inf]
#QCD.components = [ QCD_100_200, QCD_200_300,  QCD_300_500, QCD_500_700, QCD_700_1000, QCD_1000_1500]


QCDMu_ = sample()
QCDMu_.files = outlist (d,"QCDMuPt20toInf")
QCDMu_.jpref = jetLabel 
QCDMu_.jp = jetLabel 
QCDMu_.sf = 1
QCDMu_.skimEff = 1.
QCDMu_.sigma = 720648000*0.00042*0.4
QCDMu_.color = ROOT.kGray
QCDMu_.style = 1
QCDMu_.fill = 1001
QCDMu_.leglabel = "QCD"
QCDMu_.label = "QCDMuPt20toInf"

QCDMu = sample()
QCDMu.files = outlist (d,"QCDMuPt20toInf")
QCDMu.jpref = jetLabel 
QCDMu.jp = jetLabel 
QCDMu.sf = 1
QCDMu.skimEff = 1.
QCDMu.sigma = 720648000*0.00042*0.4
QCDMu.color = ROOT.kGray
QCDMu.style = 1
QCDMu.fill = 1001
QCDMu.leglabel = "QCD"
QCDMu.label = "QCDMuPt20toInf"
QCDMu.components = [QCDMu_]

QCD = sample()
#QCD.color = 4
#QCD.color = ROOT.kBlue+2
#QCD.color = 7
QCD.color = ROOT.kOrange +7
QCD.style = 1
QCD.fill = 1001
QCD.leglabel = "QCD"
QCD.label = "QCD"
QCD.components = [ QCD_100_200, QCD_200_300,  QCD_300_500, QCD_500_700, QCD_700_1000, QCD_1000_1500,  QCD_1500_2000,  QCD_2000_Inf]
#QCD.components = [ QCD_100_200, QCD_200_300,  QCD_300_500, QCD_500_700, QCD_700_1000, QCD_1000_1500]
