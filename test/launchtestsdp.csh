#nohup cmsRun topplusdmTrees_skim_cfg.py isData=False doSynch=True doPreselection=False changeJECs=True channel="tch_psd" outputLabel="treesTTDM_psd.root" sample="file:/afs/cern.ch/work/o/oiorio/public/xWajid/synch/mc_jecv4/B2GSynchMC.root" > & a.log &

#nohup cmsRun topplusdmTrees_skim_cfg.py isData=False doSynch=True doPreselection=False changeJECs=True outputLabel="treesTTDM_norm.root" sample="file:/afs/cern.ch/work/o/oiorio/public/xWajid/synch/mc_jecv4/B2GSynchMC.root" > & b.log &

nohup cmsRun topplusdmTrees_skim_cfg.py isData=False doSynch=True doPreselection=False changeJECs=True channel="tch_sd" outputLabel="treesTTDM_sd.root" sample="file:/afs/cern.ch/work/o/oiorio/public/xWajid/synch/mc_jecv4/B2GSynchMC.root" maxEvts=-1> & sdd.log &



