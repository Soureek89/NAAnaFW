import FWCore.ParameterSet.Config as cms
import copy
#from ttDM.treeDumper.topplusdmTrees_Data_V4_cfg import isData

leptonssize = cms.untracked.int32(10)
#jetssize = cms.untracked.int32(20)
jetssize = cms.untracked.int32(20)
genpartsize = cms.untracked.int32(20)

pholabel = cms.string("photons")
mulabel = cms.string("muons")
elelabel = cms.string("electrons")
jetlabel = cms.string("jetsAK4CHS")
jetak8label = cms.string("jetsAK8CHS")
subjetak8label = cms.string("subjetsAK8CHS")
metlabel = cms.string("metFull")
#metNoHFlabel = cms.string("metNoHF")
eventlabel = cms.string("eventShapePFVars")
eventJetlabel = cms.string("eventShapePFJetVars")
centralitylabel = cms.string("centrality")
genpartlabel = cms.string("genpart")


#Systematics:
systsToSave = ["noSyst","jes__up","jes__down","jer__up","jer__down"]
#systsToSave = ["jes__up","jes__down","jer__up","jer__down","unclusteredMet__up","unclusteredMet__down"]
#systsToSave = ["noSyst","jes__up","jes__down"]
#systsToSave = ["noSyst","jer__up","jer__down"]
#systsToSave = ["jer__up","jer__down","unclusteredMet__up","unclusteredMet__down"]
systsToSave = ["noSyst","unclusteredMet__up","unclusteredMet__down"]
systsToSave = ["noSyst","jes__up","jes__down"]
systsToSave = ["noSyst"]

#systsToSave = ["noSyst","unclusteredMet__up","unclusteredMet__down"]
#systsToSave = []
#systsToSave = ["jer__up"]

metFilters = ["Flag_CSCTightHaloFilter","Flag_goodVertices", "Flag_eeBadScFilter"]

catMu = ["Tight","TightAntiIso","Loose"]
catEl = ["Tight","TightAntiIso","Veto"]
catJet = ["Tight"]
#catMet = ["CorrT1"]
catMet = []

#                    muonsTight_Pt -> Iso <0.15 
#ttDM__noSyst->Draw("muonsTightIso04_LE_0p06_Pt") -> Iso < 0.06

scanMu = ["Iso04_0p06_LE","Iso04_0p15_LE","Iso04_0p06_GE","Iso04_0p15_GE"]
scanEl = []
scanJet = ["CorrPt_20"]

sysMu = [""]
sysEl = [""]
sysJet = ["JESUp","JESDown"]

#scanJet = ["CorrPt_30"]
#scanJet = ["Pt_30","Pt_40"]

#Triggers 50 ns
#SingleElTriggers = ["HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v1"]
#SingleMuTriggers = ["HLT_IsoMu20_v2","HLT_IsoMu20_eta2p1_v2",  "HLT_IsoTkMu20_eta2p1_v2", "HLT_IsoMu20_eta2p1_IterTrk02_v1"]

#hadronTriggers = ["HLT_PFMET170_NoiseCleaned_v2", "HLT_PFMET120_NoiseCleaned_BTagCSV0p72_v2", "HLT_PFMET120_PFMHT120_IDTight_v1", "HLT_PFMET90_PFMHT90_IDTight_v1", "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v2", "HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v2", "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_v2", "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v2", "HLT_AK8PFJet360_TrimMass30_v2", "HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v2", "HLT_PFHT450_SixJet40_PFBTagCSV0p72_v2", "HLT_PFHT800_v1", "HLT_PFHT750_4Jet_v1"]

#hadronTriggersMC = ["HLT_PFMET170_NoiseCleaned_v1","HLT_PFMET120_NoiseCleaned_BTagCSV07_v1", "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v1", "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v1", "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p41_v1","HLT_AK8PFJet360_TrimMass30_v1","HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1","HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v1","HLT_PFHT450_SixJet40_PFBTagCSV0p72_v2", "HLT_PFHT800_v1", "HLT_PFHT750_4Jet_v1"]

#"HLT_PFMET120_NoiseCleaned_BTagCSV07_v1", "", "HLT_PFMET120_PFMHT120_IDTight_v1", "HLT_PFMET90_PFMHT90_IDTight_v1", "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v2", "HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v2", "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_v2", "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v2", "HLT_AK8PFJet360_TrimMass30_v2", "HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v2", "HLT_PFHT450_SixJet40_PFBTagCSV0p72_v2", "HLT_PFHT800_v1"]

addCentrality=False

cutOnTriggers = False
doPreselectionCuts = False

#What to use for jets/other variables
saveBase = cms.untracked.bool(False)
saveNoCategory = False
j= "jetsAK4CHS"
jpref= "jetAK4CHS"

doAK8=False
doPhoton=False

j8= "jetsAK4CHS"
j8pref= "jetAK4CHS"

sj = "subjetsAK8CHS"
sjpref = "subjetAK8CHS"


#sj = "subjetsCmsTopTag"
#sjpref = "subjetsCmsTopTag"
subjetak8label = cms.string(sj)

jecVersion = cms.string("Spring16_25nsV10")

#Initializing the analyzer
DMTreesDumper = cms.EDAnalyzer(
    'DMAnalysisTreeMaker',
    lhes = cms.InputTag('source'),
    genprod = cms.InputTag('generator'),
    muLabel = mulabel,
    eleLabel = elelabel,
    jetsLabel = jetlabel,
    photonLabel = pholabel,
    boostedTopsLabel = jetak8label,
    boostedTopsSubjetsLabel = subjetak8label,
    metLabel = metlabel,
    genParticles =cms.InputTag( "filteredPrunedGenParticles"),
    genJets =cms.InputTag("genJetsAK4"),
#    metNoHFLabel = metNoHFlabel,
    eventLabel = eventlabel,
    JECVersion = jecVersion,
    physicsObjects = cms.VPSet(
        ),
    
    doPreselection = cms.untracked.bool(doPreselectionCuts), 

    partID = cms.InputTag( "genPart", "genPartID"),
    partStatus = cms.InputTag( "genPart", "genPartStatus"),
    partMom0ID = cms.InputTag( "genPart", "genPartMom0ID"),
    partPt = cms.InputTag( "genPart", "genPartPt"),
    partEta = cms.InputTag( "genPart", "genPartEta"),
    partPhi = cms.InputTag( "genPart", "genPartPhi"),
    partE = cms.InputTag( "genPart", "genPartE"),

    ak8jetSubjetIndex0 = cms.InputTag( "jetsAK8CHS", "jetAK8CHStopSubjetIndex0"),
    ak8jetSubjetIndex1 = cms.InputTag( "jetsAK8CHS", "jetAK8CHStopSubjetIndex1"),
    ak8jetSubjetIndex2 = cms.InputTag( "jetsAK8CHS", "jetAK8CHStopSubjetIndex2"),
    ak8jetSubjetIndex3 = cms.InputTag( "jetsAK8CHS", "jetAK8CHStopSubjetIndex3"),
    #trigger part:
    useTriggers = cms.untracked.bool(True),
    cutOnTriggers = cms.untracked.bool(cutOnTriggers),
    triggerBits = cms.InputTag("TriggerUserData","triggerBitTree"),
    triggerNames = cms.InputTag("TriggerUserData","triggerNameTree"),
    triggerPrescales = cms.InputTag("TriggerUserData","triggerPrescaleTree"),
    #met filters
    useMETFilters=cms.untracked.bool(False),
    metBits = cms.InputTag("METUserData","triggerBitTree"),
    metNames = cms.InputTag("METUserData","triggerNameTree"),
    #lumi,run,number
    lumiBlock = cms.InputTag("eventInfo","evtInfoLumiBlock"),
    runNumber = cms.InputTag("eventInfo","evtInfoRunNumber"),
    eventNumber = cms.InputTag("eventInfo","evtInfoEventNumber"),
    #HBHE
    #HBHEFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun1"),
    HBHEFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Loose"),
    HBHEIsoFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),
    #vertex
    jetKeysAK4CHS = cms.InputTag("jetKeysAK4CHS", ""),
    muonKeys = cms.InputTag("muonKeys", ""),

    vertexZ =  cms.InputTag("vertexInfo","z"),
    vertexChi2 =  cms.InputTag("vertexInfo","chi"),
    vertexNdof =  cms.InputTag("vertexInfo","ndof"),
    vertexRho =  cms.InputTag("vertexInfo","rho"),

    #resolved top part:
    doResolvedTopHad=cms.untracked.bool(False),
    doResolvedTopSemiLep=cms.untracked.bool(False),
    #cuts for the jet scan
    jetScanCuts=cms.vdouble(30), #Note: the order is important, as the jet collection with the first cut is used for the definition of mt2w.
    
    #Systematics trees to produce. Include:
    #jes__up,jes__down,jer__up,jer__down,unclusteredMet__up,unclusteredMet__down    
    systematics = cms.vstring(systsToSave), #cms.vstring("jes__up","jes__down","jer__up","jer__down","unclusteredMet__up","unclusteredMet__down"),

    channelInfo = cms.PSet(
        channel = cms.string("ttDM"),#Name of the channel, to use in the trees
        crossSection = cms.double(1),#Cross section in pb
        originalEvents = cms.double(1),#Number of events in the MC
        hadronicTriggers = cms.vstring(""),
        SingleElTriggers = cms.vstring(""),
        SingleMuTriggers = cms.vstring(""),
        PhotonTriggers = cms.vstring(""),
        metFilters = cms.vstring(metFilters),
        useLHE = cms.untracked.bool(False),#Whether one uses the weights from the LHE in order to get scale uncertainties
        useLHEWeights = cms.untracked.bool(False),#Whether one uses the weights from the LHE in order to get scale uncertaintiesxb
        addLHAPDFWeights = cms.untracked.bool(False), #Whether to add the PDF for uncertainty evaluation (time consuming)
        maxWeights = cms.untracked.int32(111),
        doTopDecayReshaping = cms.untracked.bool(False),
        doTopResweighting = cms.untracked.bool(False),
        )
    )

#DMTreesDumper.physicsObjects.append( 
#    cms.PSet(
#        label = metNoHFlabel,
#        prefix = cms.string("metNoHF"),
#        maxInstances =  cms.untracked.int32(1),
#        saveBaseVariables = cms.untracked.bool(True),
#        categories = cms.vstring(),
#        variablesD = cms.VInputTag(),
#        variablesF = cms.VInputTag(
#            cms.InputTag("metNoHF","metNoHFPt"),
#            cms.InputTag("metNoHF","metNoHFPhi"),
#            cms.InputTag("metNoHF","metNoHFPx"),
#            cms.InputTag("metNoHF","metNoHFPy"),
#            ),
#        variablesI = cms.VInputTag(
#            ),
#        singleD = cms.VInputTag(),
#        singleI = cms.VInputTag(),
#        singleF = cms.VInputTag(),
#        toSave = cms.vstring(),
#        )
#    )


#Now taking the other input objects:
DMTreesDumper.physicsObjects.append(  
    cms.PSet(
        label = mulabel,
        prefix = cms.string("mu"),
        maxInstances = leptonssize,
        saveBaseVariables = saveBase,
        saveNoCat = cms.untracked.bool(saveNoCategory),
        #categories = cms.vstring(),
        variablesD = cms.VInputTag(),
        variablesF = cms.VInputTag(
            cms.InputTag("muons","muE"),
            cms.InputTag("muons","muPt"),
            cms.InputTag("muons","muEta"),
            cms.InputTag("muons","muPhi"),
            cms.InputTag("muons","muCharge"),
            cms.InputTag("muons","muIsLooseMuon"),
            cms.InputTag("muons","muIsSoftMuon"),
            cms.InputTag("muons","muIsTightMuon"),
#            cms.InputTag("muons","muD0"),
#            cms.InputTag("muons","muD0err"),
#            cms.InputTag("muons","muDz"),
#            cms.InputTag("muons","muDzerr"),
            cms.InputTag("muons","muGenMuonCharge"),
            cms.InputTag("muons","muGenMuonEta"),
            cms.InputTag("muons","muGenMuonPt"),
            cms.InputTag("muons","muGenMuonE"),
            cms.InputTag("muons","muGenMuonPhi"),
#            cms.InputTag("muons","muGenMuonY"),
            cms.InputTag("muons","muGlbTrkNormChi2"),
            #cms.InputTag("muons","muHLTmuonDeltaR"),
            #cms.InputTag("muons","muHLTmuonE"),
            #cms.InputTag("muons","muHLTmuonEta"),
            #cms.InputTag("muons","muHLTmuonPt"),
            #cms.InputTag("muons","muHLTmuonPhi"),
            cms.InputTag("muons","muInTrkNormChi2"),
            cms.InputTag("muons","muIsGlobalMuon"),
            cms.InputTag("muons","muIsPFMuon"),
            cms.InputTag("muons","muIsTrackerMuon"),
            cms.InputTag("muons","muIso04"),
            cms.InputTag("muons","muNumberMatchedStations"),
            cms.InputTag("muons","muNumberOfPixelLayers"),
            cms.InputTag("muons","muNumberOfValidTrackerHits"),
            cms.InputTag("muons","muNumberTrackerLayers"),
            cms.InputTag("muons","muNumberValidMuonHits"),
            cms.InputTag("muons","muNumberValidPixelHits"),
            cms.InputTag("muons","muSumChargedHadronPt"),
            cms.InputTag("muons","muSumNeutralHadronPt"),
            cms.InputTag("muons","muSumPUPt"),
            cms.InputTag("muons","muSumPhotonPt"),
#            cms.InputTag("muons","muY"),
            
            ),
        variablesI = cms.VInputTag(),
        singleD = cms.VInputTag(),
        singleI = cms.VInputTag(),
        singleF = cms.VInputTag(),
        scanCuts = cms.vstring(scanMu),
        categories = cms.vstring(catMu),
        systCats = cms.vstring(sysMu),
        #        categories = cms.vstring("Tight","Loose"),
        toSave = cms.vstring("muE","muPt","muEta","muPhi","muIso04","muCharge","muIsTightMuon","muIsLooseMuon","allExtra"),
        )
    )
if doPhoton:
    DMTreesDumper.physicsObjects.append( 
        cms.PSet(
            label = pholabel,
            prefix = cms.string("pho"),
            maxInstances =  leptonssize,
            saveBaseVariables = cms.untracked.bool(True),
            categories = cms.vstring(),
            scanCuts = cms.vstring(),
            systCats = cms.vstring(),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(
                cms.InputTag("photons","phoEta"),
                cms.InputTag("photons","phoPt"),
                cms.InputTag("photons","phoHoverE"),
                cms.InputTag("photons","phoSigmaIEtaIEta"),
                cms.InputTag("photons","phoNeutralHadronIso"),
                cms.InputTag("photons","phoNeutralHadronIsoEAcorrected"),
                cms.InputTag("photons","phoChargedHadronIso"),
                cms.InputTag("photons","phoChargedHadronIsoEAcorrected"),
                cms.InputTag("photons","phoPhotonIso"),
                cms.InputTag("photons","phoPhotonIsoEAcorrected"),
                cms.InputTag("photons","phoHasPixelSeed"),
                cms.InputTag("photons","phoPassLooseID"),
                cms.InputTag("photons","phoPassMediumID"),
                cms.InputTag("photons","phoPassTightID"),
                ),
            variablesI = cms.VInputTag(
                ),
            singleD = cms.VInputTag(),
            singleI = cms.VInputTag(),
            singleF = cms.VInputTag(),
            toSave = cms.vstring("phoPt","phoHoverE","phoSigmaIEtaIEta","phoNeutralHadronIso","phoNeutralHadronIsoEAcorrected","phoChargedHadronIso","phoChargedHadronIsoEAcorrected","phoPhotonIso","phoPhotonIsoEAcorrected","phoHasPixelSeed","phoPassLooseID","phoPassMediumID","phoPassTightID","allExtra"),
            )
    )
DMTreesDumper.physicsObjects.append( 
    cms.PSet(
        label = elelabel,
        prefix = cms.string("el"),
        maxInstances = leptonssize,
        saveBaseVariables = saveBase,
        saveNoCat = cms.untracked.bool(saveNoCategory),
        variablesD = cms.VInputTag(),
        variablesF = cms.VInputTag(
            cms.InputTag("electrons","elE"),
            cms.InputTag("electrons","elPt"),
#            cms.InputTag("electrons","elMass"),
            cms.InputTag("electrons","elEta"),
            cms.InputTag("electrons","elPhi"),
            cms.InputTag("electrons","elCharge"),
            cms.InputTag("electrons","elEta"),
            cms.InputTag("electrons","elSCEta"),
            cms.InputTag("electrons","elHoE"),
            cms.InputTag("electrons","elIso03"),
#            cms.InputTag("electrons","elY"),
            cms.InputTag("electrons","eldEtaIn"),
            cms.InputTag("electrons","eldPhiIn"),
            cms.InputTag("electrons","elmissHits"),
            cms.InputTag("electrons","elfull5x5siee"),
            cms.InputTag("electrons","elooEmooP"),
            cms.InputTag("electrons","elhasMatchedConVeto"),
            cms.InputTag("electrons","elvidVeto"),
            cms.InputTag("electrons","elvidTight"),
            cms.InputTag("electrons","elvidTightnoiso"),
            cms.InputTag("electrons","elvidMedium"),
            cms.InputTag("electrons","elDz"),
            cms.InputTag("electrons","elDxy"),
        ),
        variablesI = cms.VInputTag(
            ),
        singleD = cms.VInputTag(),
        singleI = cms.VInputTag(),
        singleF = cms.VInputTag(),
        scanCuts = cms.vstring(scanEl),
        categories = cms.vstring(catEl),
        systCats = cms.vstring(sysEl),
#        categories = cms.vstring("Tight","Veto"),
        toSave = cms.vstring("elE","elPt","elEta","elPhi","elIso03","elisTight","elvidTightnoiso","elCharge","elisMedium","elisLoose","elisVeto","elSCEta","allExtra","elDz","elDxy"),
        )
    ) 


#if (False):
#    DMTreesDumper.physicsObjects.append( 
#        cms.PSet(
#            label = eventlabel,
#            prefix = cms.string(""),
#            maxInstances =  cms.untracked.int32(1),
#            saveBaseVariables = cms.untracked.bool(True),
#            scanCuts = cms.vstring(),
#            categories = cms.vstring(),
#            variablesF = cms.VInputTag(),
#            variablesD = cms.VInputTag(),
#            variablesI = cms.VInputTag(),
#            singleI = cms.VInputTag(),
#            singleF = cms.VInputTag(),
#            singleD = cms.VInputTag(
#                cms.InputTag("eventShapePFVars","isotropy"),
#                cms.InputTag("eventShapePFVars","C"),
#                cms.InputTag("eventShapePFVars","D"),
#                cms.InputTag("eventShapePFVars","aplanarity"),
#                cms.InputTag("eventShapePFVars","circularity"),
#                cms.InputTag("eventShapePFVars","sphericity"),
#                cms.InputTag("eventShapePFVars","thrust"),
#                cms.InputTag("eventShapePFVars","thrust"),
#                ),
#            toSave = cms.vstring(),
#            )
#        )
#
#DMTreesDumper.physicsObjects.append( 
#    cms.PSet(
#        label = eventJetlabel,
#        prefix = cms.string(""),
#        maxInstances =  cms.untracked.int32(1),
#        saveBaseVariables = cms.untracked.bool(True),
#        categories = cms.vstring(),
#        scanCuts = cms.vstring(),
#        variablesF = cms.VInputTag(),
#        variablesD = cms.VInputTag(),
#        variablesI = cms.VInputTag(),
#        singleI = cms.VInputTag(),
#        singleF = cms.VInputTag(),
#        singleD = cms.VInputTag(
#            cms.InputTag("eventShapePFJetVars","isotropy"),
#            cms.InputTag("eventShapePFJetVars","C"),
#            cms.InputTag("eventShapePFJetVars","D"),
#            cms.InputTag("eventShapePFJetVars","aplanarity"),
#            cms.InputTag("eventShapePFJetVars","circularity"),
#            cms.InputTag("eventShapePFJetVars","sphericity"),
#            cms.InputTag("eventShapePFJetVars","thrust"),
#            cms.InputTag("eventShapePFJetVars","thrust"),
#            ),
#        toSave = cms.vstring(),
#        )
#    )

if addCentrality:
    DMTreesDumper.physicsObjects.append( 
        cms.PSet(
            label = centralitylabel,
            prefix = cms.string(""),
            maxInstances =  cms.untracked.int32(1),
            saveBaseVariables = cms.untracked.bool(True),
            categories = cms.vstring(),
            scanCuts = cms.vstring(),
            variablesF = cms.VInputTag(),
            variablesD = cms.VInputTag(),
            variablesI = cms.VInputTag(),
            singleI = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleD = cms.VInputTag(
                cms.InputTag("centrality","centrality"),
                ),
            toSave = cms.vstring(),
            )
        )


DMTreesDumper.physicsObjects.append( 
    cms.PSet(
        label = jetlabel,
        prefix = cms.string(jpref),
        maxInstances = jetssize,
        saveBaseVariables = saveBase,
        saveNoCat = cms.untracked.bool(saveNoCategory),
#        saveNoCat = cms.untracked.bool(True),
        variablesD = cms.VInputTag(),
        variablesF = cms.VInputTag(
            cms.InputTag(j,jpref+"E"),
            cms.InputTag(j,jpref+"Pt"),
#            cms.InputTag(j,jpref+"Mass"),
            cms.InputTag(j,jpref+"Eta"),
            cms.InputTag(j,jpref+"Phi"),
#            cms.InputTag(j,jpref+"QGL"),
            cms.InputTag(j,jpref+"PartonFlavour"),
            cms.InputTag(j,jpref+"Phi"),
            cms.InputTag(j,jpref+"CSVv2"),
            cms.InputTag(j,jpref+"CMVAv2"),
            #cms.InputTag(j,jpref+"CSVV1"),
            cms.InputTag(j,jpref+"Charge"),
#            cms.InputTag(j,jpref+"ChargeMuEnergy"),
#            cms.InputTag(j,jpref+"ChargedHadronMultiplicity"),
#            cms.InputTag(j,jpref+"ElectronEnergy"),
            cms.InputTag(j,jpref+"GenJetCharge"),
            cms.InputTag(j,jpref+"GenJetE"),
            cms.InputTag(j,jpref+"GenJetEta"),
            cms.InputTag(j,jpref+"GenJetPhi"),
            cms.InputTag(j,jpref+"GenJetPt"),
#            cms.InputTag(j,jpref+"GenJetY"),
            cms.InputTag(j,jpref+"GenPartonCharge"),
            cms.InputTag(j,jpref+"GenPartonE"),
            cms.InputTag(j,jpref+"GenPartonEta"),
            cms.InputTag(j,jpref+"GenPartonPhi"),
            cms.InputTag(j,jpref+"GenPartonPt"),
#            cms.InputTag(j,jpref+"GenPartonY"),
#            cms.InputTag(j,jpref+"HFEMEnergy"),
#            cms.InputTag(j,jpref+"HFEMMultiplicity"),
#            cms.InputTag(j,jpref+"HFHadronEnergy"),
#            cms.InputTag(j,jpref+"HFHadronMultiplicity"),
            cms.InputTag(j,jpref+"HadronFlavour"),
 #           cms.InputTag(j,jpref+"SmearedE"),
 #           cms.InputTag(j,jpref+"SmearedPt"),
 #           cms.InputTag(j,jpref+"SmearedPEta"),
 #           cms.InputTag(j,jpref+"SmearedPhi"),
#            cms.InputTag(j,jpref+"Y"),
#            cms.InputTag(j,jpref+"muonEnergyFrac"),

            cms.InputTag(j,jpref+"neutralHadronEnergyFrac"),
            cms.InputTag(j,jpref+"neutralEmEnergyFrac"),


            cms.InputTag(j,jpref+"chargedEmEnergyFrac"),
            cms.InputTag(j,jpref+"chargedHadronEnergyFrac"),

#            cms.InputTag(j,jpref+"NumConstituents"),
#            cms.InputTag(j,jpref+"electronMultiplicity"),
#            cms.InputTag(j,jpref+"muonMultiplicity"),
 #           cms.InputTag(j,jpref+"neutralHadronMultiplicity"),
            cms.InputTag(j,jpref+"neutralMultiplicity"),
#            cms.InputTag(j,jpref+"photonMultiplicity"),
#            cms.InputTag(j,jpref+"numberOfDaughters"),
#            cms.InputTag(j,jpref+"chargedHadronEnergy"),
            cms.InputTag(j,jpref+"chargedMultiplicity"),
#            cms.InputTag(j,jpref+"chargedEmEnergy"),
            cms.InputTag(j,jpref+"neutralEmEnergy"),
#            cms.InputTag(j,jpref+"neutralHadronEnergy"),
            cms.InputTag(j,jpref+"jecFactor0"),
            cms.InputTag(j,jpref+"jetArea"),
#            cms.InputTag(j,jpref+"pileupJetIdRMS"),#
#            cms.InputTag(j,jpref+"pileupJetIdbeta"),
#            cms.InputTag(j,jpref+"pileupJetIdbetaClassic"),
#            cms.InputTag(j,jpref+"pileupJetIdbetaStar"),
#            cms.InputTag(j,jpref+"pileupJetIdbetaStarClassic"),
#            cms.InputTag(j,jpref+"pileupJetIddR2Mean"),
#            cms.InputTag(j,jpref+"pileupJetIddRMean"),
#            cms.InputTag(j,jpref+"pileupJetIddZ"),
#            cms.InputTag(j,jpref+"pileupJetIdfrac01"),
#            cms.InputTag(j,jpref+"pileupJetIdfrac02"),
#            cms.InputTag(j,jpref+"pileupJetIdfrac03"),
#            cms.InputTag(j,jpref+"pileupJetIdfrac04"),
#            cms.InputTag(j,jpref+"pileupJetIdjetR"),
#            cms.InputTag(j,jpref+"pileupJetIdjetRchg"),
#            cms.InputTag(j,jpref+"pileupJetIdmajW"),
#            cms.InputTag(j,jpref+"pileupJetIdminW"),
#            cms.InputTag(j,jpref+"pileupJetIdnCharged"),
#            cms.InputTag(j,jpref+"pileupJetIdnNeutrals"),
#            cms.InputTag(j,jpref+"pileupJetIdnParticles"),
#            cms.InputTag(j,jpref+"pileupJetIdptD"),
#            cms.InputTag(j,jpref+"pileupJetIdpull")
            ),
        variablesI = cms.VInputTag(),
        singleD = cms.VInputTag(),
        singleI = cms.VInputTag(),
        singleF = cms.VInputTag(),
        #toSave = cms.vstring(jpref+"Eta",jpref+"Phi","allExtra"),
        scanCuts = cms.vstring(scanJet),
        categories = cms.vstring(catJet),
        systCats = cms.vstring(sysJet),
#        categories = cms.vstring("Tight"),
#        toSave = cms.vstring(jpref+"E",jpref+"Pt",jpref+"Eta",jpref+"Phi",jpref+"GenJetPt",jpref+"GenJetEta",jpref+"CSVv2",jpref+"PartonFlavour",jpref+"QGL", jpref+"jecFactor0",jpref+"pileupJetIdRMS", jpref+"pileupJetIdbeta", jpref+"pileupJetIdbetaClassic", jpref+"pileupJetIdbetaStar", jpref+"pileupJetIdbetaStarClassic", jpref+"pileupJetIddR2Mean", jpref+"pileupJetIddRMean", jpref+"pileupJetIddZ", jpref+"pileupJetIdfrac01", jpref+"pileupJetIdfrac02", jpref+"pileupJetIdfrac03", jpref+"pileupJetIdfrac04", jpref+"pileupJetIdjetR", jpref+"pileupJetIdjetRchg", jpref+"pileupJetIdmajW", jpref+"pileupJetIdminW", jpref+"pileupJetIdnCharged", jpref+"pileupJetIdnNeutrals", jpref+"pileupJetIdnParticles", jpref+"pileupJetIdptD", jpref+"pileupJetIdpull","allExtra"),
        toSave = cms.vstring(jpref+"E",jpref+"Pt",jpref+"Eta",jpref+"Phi",jpref+"GenJetPt",jpref+"GenJetEta",jpref+"CMVAv2",jpref+"CSVv2",jpref+"HadronFlavour",jpref+"PartonFlavour",jpref+"QGL", jpref+"jecFactor0",jpref+"pileupJetIdRMS", jpref+"pileupJetIdbeta", jpref+"pileupJetIdbetaClassic", jpref+"pileupJetIdbetaStar", jpref+"pileupJetIdbetaStarClassic","allExtra"),
        ),
    )

if(doAK8):
    DMTreesDumper.physicsObjects.append( 
        cms.PSet(
            label = jetak8label,
            prefix = cms.string(j8pref),
            maxInstances = cms.untracked.int32(10),
            scanCuts = cms.vstring(),
            categories = cms.vstring(),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(
                cms.InputTag(j8,j8pref+"E"),
                cms.InputTag(j8,j8pref+"Pt"),
#                cms.InputTag(j8,j8pref+"Mass"),
                cms.InputTag(j8,j8pref+"Eta"),
                cms.InputTag(j8,j8pref+"Phi"),
                cms.InputTag(j8,j8pref+"Mass"),
                cms.InputTag(j8,j8pref+"minmass"),
                cms.InputTag(j8,j8pref+"nSubJets"),
                cms.InputTag(j8,j8pref+"topSubjetIndex0"),
                cms.InputTag(j8,j8pref+"topSubjetIndex1"),
                cms.InputTag(j8,j8pref+"topSubjetIndex2"),
                cms.InputTag(j8,j8pref+"topSubjetIndex3"),
                cms.InputTag(j8,j8pref+"prunedMass"),
                cms.InputTag(j8,j8pref+"tau1"),
                cms.InputTag(j8,j8pref+"tau2"),
                cms.InputTag(j8,j8pref+"tau3"),
                cms.InputTag(j8,j8pref+"topMass"),
                cms.InputTag(j8,j8pref+"trimmedMass"),
                cms.InputTag(j8,j8pref+"wMass"),
                ),
            variablesI = cms.VInputTag(),
            singleD = cms.VInputTag(),
            singleI = cms.VInputTag(),
            singleF = cms.VInputTag(),
            toSave = cms.vstring(j8pref+"E",j8pref+"Pt",j8pref+"Eta",j8pref+"Phi",j8pref+"Mass",j8pref+"minmass", j8pref+"nSubJets", j8pref+"prunedMass", j8pref+"tau1", j8pref+"tau2", j8pref+"tau3", j8pref+"topMass", j8pref+"trimmedMass", j8pref+"wMass",  "allExtra"),
            )
        )

    
    DMTreesDumper.physicsObjects.append( 
        cms.PSet(
            label = subjetak8label,
            prefix = cms.string(sjpref),
            maxInstances = cms.untracked.int32(10),
            scanCuts = cms.vstring(),
            categories = cms.vstring(),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(
                cms.InputTag(sj,sjpref+"E"),
                cms.InputTag(sj,sjpref+"Pt"),
#                cms.InputTag(sj,sjpref+"Mass"),
                cms.InputTag(sj,sjpref+"Eta"),
                cms.InputTag(sj,sjpref+"Phi"),
                #            cms.InputTag(sj,sjpref+"subjetCSV"),
                cms.InputTag(sj,sjpref+"CSVv2"),
                ),                         
            variablesI = cms.VInputTag(),
            singleD = cms.VInputTag(),
            singleI = cms.VInputTag(),
            singleF = cms.VInputTag(),
            toSave = cms.vstring(sjpref+"E",sjpref+"Pt",sjpref+"Eta",sjpref+"Phi", "allExtra"),
            )
        )
    

DMTreesDumper.physicsObjects.append(
        cms.PSet(
        label = metlabel,
            prefix = cms.string("metFull"),
            maxInstances =  cms.untracked.int32(1),
            saveBaseVariables = cms.untracked.bool(True),
            categories = cms.vstring(),
            scanCuts = cms.vstring(),
#            systCats = cms.vstring(),
            systCats = cms.vstring(),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(
                cms.InputTag("metFull","metFulluncorPt"),
                cms.InputTag("metFull","metFulluncorPhi"),
                cms.InputTag("metFull","metFullPt"),
                cms.InputTag("metFull","metFullPhi"),
                cms.InputTag("metFull","metFullPx"),
                cms.InputTag("metFull","metFullPy"),
                ),
            variablesI = cms.VInputTag(),
            singleD = cms.VInputTag(),
            singleI = cms.VInputTag(),
            singleF = cms.VInputTag(),
            toSave = cms.vstring(),
            )
)
#if(not isData):
#DMTreesDumper.physicsObjects.append(  
#    cms.PSet(
#        label = genpartlabel,
#        prefix = cms.string("genPart"),
#        maxInstances = genpartsize,
#        saveBaseVariables = saveBase,
#        scanCuts = cms.vstring(),
#        categories = cms.vstring(),
#        variablesD = cms.VInputTag(),
#        variablesF = cms.VInputTag(
#            cms.InputTag("genPart","genPartCharge"),
#            cms.InputTag("genPart","genPartE"),
#            cms.InputTag("genPart","genPartEta"),
#            cms.InputTag("genPart","genPartID"),
#            cms.InputTag("genPart","genPartMass"),
#            cms.InputTag("genPart","genPartMom0ID"),
#            cms.InputTag("genPart","genPartPhi"),
#            cms.InputTag("genPart","genPartPt"),
#            cms.InputTag("genPart","genPartStatus"),
##            cms.InputTag("genPart","genPartY"),
#            ),
#        variablesI = cms.VInputTag(),
#        singleD = cms.VInputTag(),
#        singleI = cms.VInputTag(),
#        singleF = cms.VInputTag(),
#        toSave = cms.vstring("genPartCharge", "genPartE", "genPartEta", "genPartID", "genPartMass", "genPartMom0ID", "genPartPhi", "genPartPt", "genPartStatus", "genPartY", "genPartSize"),
#        )
#    )
    
    

    
