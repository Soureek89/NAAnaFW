from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.workArea = 'crab_13July'
config.General.requestName = '80xV2_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'
config.section_('JobType')
config.JobType.psetName = 'topplusdmTrees_skim_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['Summer16_23Sep2016BCDV4_DATA_UncertaintySources_AK4PFchs.txt', 'Summer16_23Sep2016V4_MC_L2L3Residual_AK4PFchs.txt', 'Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt', 'Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt', 'Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt', 'Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt', 'Summer16_23Sep2016BCDV4_DATA_L1RC_AK4PFchs.txt', 'Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt', 'Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt', 'Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt', 'Summer16_23Sep2016EFV4_DATA_UncertaintySources_AK4PFchs.txt', 'Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt', 'Summer16_23Sep2016EFV4_DATA_L1RC_AK4PFchs.txt', 'Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt', 'Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt', 'Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt', 'Summer16_23Sep2016GV4_DATA_UncertaintySources_AK4PFchs.txt', 'Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt', 'Summer16_23Sep2016GV4_DATA_L1RC_AK4PFchs.txt', 'Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt', 'Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt', 'Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt', 'Summer16_23Sep2016HV4_DATA_UncertaintySources_AK4PFchs.txt', 'Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt', 'Summer16_23Sep2016HV4_DATA_L1RC_AK4PFchs.txt', 'Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt', 'Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt', 'Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt', 'Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt', 'cMVAv2_Moriond17_B_H.csv', 'CSVv2_Moriond17_B_H.csv', 'btagging_cmva.root']
config.JobType.pyCfgParams = ['isData=False', 'changeJECs=True']
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
config.Data.inputDataset = '/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/oiorio-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V3p2_June-abe82725a22c385c5abbaab3e6182693/USER'
config.Data.outputDatasetTag = '80xV2_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'
config.Data.publication = False
config.Data.unitsPerJob = 5
config.Data.ignoreLocality = True
config.Data.splitting = 'FileBased'
config.Data.inputDBS = 'phys03'
config.Data.outLFNDirBase = '/store/user/oiorio/SingleTop/2017/July/13July/'
config.section_('Site')
config.Site.storageSite = 'T2_IT_Pisa'
