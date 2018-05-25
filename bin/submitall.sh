#!/bin/sh
#Running over muon and muonantiiso channel for MC
python new_singletop.py -c muon -t trees -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMu -S 6 
python new_singletop.py -c muonantiiso -t trees -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMu -S 6

#Running if for Data and DDriven QCD
python new_singletop.py -c muon -t trees -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 6
python new_singletop.py -c muonantiiso -t trees -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 6

