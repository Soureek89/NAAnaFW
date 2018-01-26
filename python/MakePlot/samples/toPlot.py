######################################
#
# Annapaola de Cosa, January 2015
#
######################################

import ROOT
import collections
import os, commands
from SingleTop import *
from VJets import *
from TT import *
from QCD import *
from DDQCD import *
from VV import *
from Data import *

samples = collections.OrderedDict()

samples["Data"] = Data
samples["ST_tch"] = ST_tch
samples["ST_tch_psq2Down"] = ST_tch_psq2Down
samples["ST_tch_psq2Up"] = ST_tch_psq2Up
samples["ST_tch_hdampDown"] = ST_tch_hdampDown
samples["ST_tch_hdampUp"] = ST_tch_hdampUp
samples["ST_tch_sd"] = ST_tch_sd
samples["ST_tch_p_sd"] = ST_tch_p_sd
samples["ST_sch"] = ST_sch
samples["ST_tW"] = ST_tW
samples["TT"] = TT
samples["TT_sd"] = TT_sd
samples["TT_psq2Down"] = TT_psq2Down
samples["TT_psq2Up"] = TT_psq2Up
samples["TT_hdampDown"] = TT_hdampDown
samples["TT_hdampUp"] = TT_hdampUp
samples["QCDMu"] = QCDMu
samples["DDQCD"] = DDQCD
samples["WJets"] = WJets
samples["DYJets"] = DYJets
samples["VV"] = VV
#samples["SingleTop_tchannel"] = ST_tch
#samples["SingleTop_schannel"] = ST_sch
#samples["SingleTop_tW"] = ST_tW







