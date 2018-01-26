#nohup cmsRun topplusdmTrees_skim_cfg.py isData=False doSynch=True doPreselection=False changeJECs=False outputLabel="treesTTDM_qcdnopres.root" sample="file:/afs/cern.ch/work/o/oiorio/public/xWajid/synch/mc_qcd_jecv4/B2GSynchMCQCD.root" > & qcdnopres.log &
#nohup cmsRun topplusdmTrees_skim_cfg.py isData=False doSynch=False doPreselection=True changeJECs=False outputLabel="treesTTDM_qcdpres.root" sample="file:/afs/cern.ch/work/o/oiorio/public/xWajid/synch/mc_qcd_jecv4/B2GSynchMCQCD.root" > & qcdpres.log &

nohup cmsRun topplusdmTrees_skim_cfg.py isData=False doSynch=False doPreselection=True changeJECs=True outputLabel="treesTTDM_ele.root" sample="file:/afs/cern.ch/work/o/oiorio/public/xWajid/synch/mc_jecv4/B2GSynchMC.root" maxEvts=100    > & stepmc.log &
#nohup cmsRun topplusdmTrees_skim_cfg.py isData=False doSynch=False doPreselection=False changeJECs=False outputLabel="treesTTDM_ele_nopres.root" sample="file:/afs/cern.ch/work/o/oiorio/public/xWajid/synch/mc_jecv4/B2GSynchMC.root" maxEvts=1000 > & stepmcnopres.log &

#nohup cmsRun topplusdmTrees_skim_cfg.py doSynch=True isData=False doSynch=True doPreselection=False changeJECs=False outputLabel="../bin/synch/mc_jecv4/treesTTDM.root" sample="file:/afs/cern.ch/work/o/oiorio/public/xWajid/synch/mc_jecv4/B2GSynchMC.root" > & stepmc.log &
#nohup cmsRun topplusdmTrees_skim_cfg.py isData=True changeJECs=False outputLabel="../bin/synch/runb/treesTTDM.root" sample="file:/afs/cern.ch/work/o/oiorio/public/xWajid/synch/runb_new/B2GSynchRunB.root" > & steprunb.log &
#nohup cmsRun topplusdmTrees_skim_cfg.py isData=True channel="DATA2016H" changeJECs=False outputLabel="../bin/synch/runh/treesTTDM.root" sample="file:/afs/cern.ch/work/o/oiorio/public/xWajid/synch/runh_new/B2GSynchRunH.root" > & steprunh.log &

#nohup cmsRun topplusdmTrees_skim_cfg.py isData=True doPreselection=False channel="DATA2016H" changeJECs=False outputLabel="../bin/treesTTDMRunHSynch.root" sample="file:/afs/cern.ch/work/o/oiorio/public/xWajid/synch/runh_new/B2GSynchRunH.root" > & steprunhdosynch.log &
#nohup cmsRun topplusdmTrees_skim_cfg.py isData=True doPreselection=True channel="DATA2016H" changeJECs=False outputLabel="../bin/treesTTDMRunHNoSynch.root" sample="file:/afs/cern.ch/work/o/oiorio/public/xWajid/synch/runh_new/B2GSynchRunH.root" > & steprunhdonotsynch.log &
