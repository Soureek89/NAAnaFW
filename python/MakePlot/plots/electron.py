met_range_had = (200,480, True)
met_rebin_had = 2 
common_rebin = 5 
rebin_sr = 10
from common import met_range

settings = {
    #'cutFlow'                  : (None                           , 20, None, (-0.5,7.5) ),

    #    'h_nPV' : ('number of primary vertices' , )
    'h_nMu': ('Number of tight Muons', 13, 1, (-0.5,12.5) ),
    'h_nJets': ('Number of tight jets', 13, 1, (-0.5,12.5) ),

#    'h_MuPt': ''
#    'h_2j1t_bjetpt' 'Leading jet b jet Pt distribution'

#    h_MuEta;1Leading Muon Eta distribution
#    h_MuPhi;1Leading Muon Phi distribution
#    h_MuE;1Leading Muon E distribution
#    h_MuCharge;1Leading Muon Charge distribution
#    h_2j0t_jetpt40;12j0t Leading jet Pt distribution
#    h_3j1t_bjetpt;13j1t Leading jet b jet Pt distribution
#    h_3j1t_MuPt;13j1t Leading Muon Pt distribution
#    h_3j1t_MuEta;13j1t Leading Muon Eta distribution
#    h_3j1t_MuPhi;13j1t Leading Muon Phi distribution
#    h_3j1t_MuE;13j1t Leading Muon E distribution
#    h_3j1t_MuCharge;13j1t Leading Muon Charge distribution
#    h_3j2t_bjetpt;13j2t Leading jet b jet Pt distribution
#    h_3j2t_2ndbjetpt;13j2t 2nd Lead. jet b jet Pt distribution
#    h_3j2t_MuPt;13j2t Leading Muon Pt distribution
#    h_3j2t_MuEta;13j2t Leading Muon Eta distribution
#    h_3j2t_MuPhi;13j2t Leading Muon Phi distribution
#    h_3j2t_MuE;13j2t Leading Muon E distribution
#    h_3j2t_MuCharge;13j2t Leading Muon Charge distribution
#    h_2j1t_bjetpt;12j1t Leading jet b jet Pt distribution
#    h_2j1t_MuPt;12j1t Leading Muon Pt distribution
#    h_2j1t_MuEta;12j1t Leading Muon Eta distribution
#    h_2j1t_MuPhi;12j1t Leading Muon Phi distribution
#    h_2j1t_MuE;12j1t Leading Muon E distribution
#    h_2j1t_MuCharge;12j1t Leading Muon Charge distribution
#    h_nEl;1Number of tight Electrons
#    h_ElPt;1Leading Elec. Pt distribution
#    h_ElEta;1Leading Elec. Eta distribution
#    h_ElPhi;1Leading Elec. Phi distribution
#    h_ElE;1Leading Elec. E distribution
#    h_ElCharge;1Leading Elec. Charge distribution
#    h_nJets;1Number of tight jets
#    h_nbJets;1Number of tight b-jets
#    h_jet1Pt;1Leading jet Pt distribution
#    h_jet2Pt;1Trailing Jet Pt distribution
#    h_jet3Pt;1Third Jet Pt distribution
#    h_bjet1Pt;1Leading b-jet Pt distribution
#    h_bjet2Pt;1Trailing b-jet Pt distribution
#    h_bjet3Pt;1Third b-jet Pt distribution
#    h_bjetsPt;1B-Jets Pt distribution

#
#    #CR3nw
#    #'metFinal_CR3nw'         : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),# To fix cpp
#    'muonmetFinal_CR3nw'         : ('E^{miss}_{T} [GeV]'                       , 20, met_rebin_had, met_range_had ),# To fix cpp
#    'electronmetFinal_CR3nw'         : ('E^{miss}_{T} [GeV]'                   , 20, met_rebin_had, met_range_had ),# To fix cpp
#
#    'h_lep1Pt_CR3nw'          : ('p_{T} leading lep'           , 20, 30, (30,310) ),
#    'h_mu1Pt_CR3nw'          : ('p_{T} leading muon'           , 20, 30, (30,190) ),
#    'h_el1Pt_CR3nw'          : ('p_{T} leading electron'           , 20, 30, (30,190) ),
#
#    'dphi_CR3nw'          : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'dphi_6j_CR3nw'          : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'mt2w_CR3nw'        : ('M^{W}_{T2} [GeV]'                         , 20, 10, (50,450) ),
#    'mt_CR3nw'          : ('M_{T} [GeV]'                        , 20, 20, (160,500) ),
#    'muondphi_CR3nw'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_6j_CR3nw'      : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'muonmt2w_CR3nw'         : ('M^{W}_{T2} [GeV]'                         , 20, 10, (50,450) ),
#    'muonmt_CR3nw'           : ('M_{T} [GeV]'                        , 20, 20, (160,500) ),
#    'electrondphi_CR3nw'     : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'electrondphi_6j_CR3nw'  : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'electronmt2w_CR3nw'     : ('M^{W}_{T2} [GeV]'                         , 20, 10, (50,450) ),
#    'electronmt_CR3nw'       : ('M_{T} [GeV]'                        , 20, 20, (160,500) ),
#
#    'h_nbJets_CR3nw'        : ('Number of tight b-jets'       , 20, 1, (2,6)),
#    'h_nJets_CR3nw'         : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_muonnJets_CR3nw'      : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_muonnbJets_CR3nw'     : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_electronnJets_CR3nw'  : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_electronnbJets_CR3nw' : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#
    "h_2j1t_mtw": (     "2j1t mtw ",250,common_rebin*2,(0,500)),
    "h_2j1t_topPt": (   "2j1t top pt",250,common_rebin,(0,500)),
    "h_3j1t_mtw": (     "3j1t mtw ",100,common_rebin*2,(0,500)),
    "h_3j2t_mtw": (        "3j2t mtw ",250,common_rebin*2,(0,500)),

#    'h_jet1Pt_CR3nw'        : ('Leading jet Pt distribution' , 20, 20, (30,390) ),
#    'h_jet2Pt_CR3nw'        : ('Sub-leading jet Pt  distribution' , 20, 20, (30,270) ),
#    'h_jet3Pt_CR3nw'        : ('Third Jet Pt distribution'    , 20, 20, (30,190) ),
#    
#    'h_muonjet1Pt_CR3nw'     : ('Leading jet Pt distribution' , 20, 10, (30,490) ),
#    'h_muonjet2Pt_CR3nw'     : ('Sub-leading jet Pt distribution'  , 20, 10, (30,490)),
#    'h_muonjet3Pt_CR3nw'     : ('Third Jet Pt distribution'    , 20, 10, (30,230) ),
#
#    'h_electronjet1Pt_CR3nw' : ('Leading jet Pt distribution' , 20, 10, (30,490) ),
#    'h_electronjet2Pt_CR3nw' : ('Sub-leading jet Pt distribution'  , 20, 10, (30,490)),
#    'h_electronjet3Pt_CR3nw' : ('Third Jet Pt distribution'    , 20, 10, (30,230) ),
#
#    #CR4
#    'metFinal_outtop'          : ('E^{miss}_{T} [GeV]'                 , 20, met_rebin_had, met_range_had ),#
#
#    'dphi_CRqcd_had'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'dphi_6j_CRqcd_had'         : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    
#    'h_nbJets_CRqcd_had'       : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_nJets_CRqcd_had'        : ('Number of tight jets'         , 20, 1, (4,12)),
#
#    'h_jet1Pt_CRqcd_had'       : ('Leading jet Pt distribution' , 20, 5, (30,490) ),
#    'h_jet2Pt_CRqcd_had'       : ('Sub-leading jet Pt distribution'  , 20, 20, (30,490)),
#    'h_jet3Pt_CRqcd_had'       : ('Third Jet Pt distribution'    , 20, 20, (30,230) ),
#
#    #CR5
#    'metFinal_CR5'             : ('E^{miss}_{T} [GeV]'                          , 200, met_rebin_had, met_range_had ),
#
#    'dphi_CR5'     : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'dphi_6j_CR5'  : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#
#    'h_nJets_CR5'  : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_nbJets_CR5' : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_jet1Pt_CR5' : ('Leading jet Pt distribution' , 20, 5, (30,490) ),
#    'h_jet2Pt_CR5' : ('Sub-leading jet Pt distribution'  , 20, 5, (30,490)),
#    'h_jet3Pt_CR5' : ('Third Jet Pt distribution'    , 20, 5, (30,230) ),
#
#    #CR6
#    #'metFinal_CR6'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'muonmetFinal_CR6'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'electronmetFinal_CR6'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#
#    'dphi_CR6'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'dphi_6j_CR6'         : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_CR6'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_6j_CR6'      : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'lectrondphi_CR6'      : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'lectrondphi_6j_CR6'   : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#
#    'h_nJets_CR6'        : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_nbJets_CR6'       : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_muonnJets_CR6'      : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_muonnbJets_CR6'     : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_electronnJets_CR6'  : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_electronnbJets_CR6' : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#
#    'h_jet1Pt_CR6'       : ('Leading jet Pt distribution' , 20, 20, (30,490) ),
#    'h_jet2Pt_CR6'       : ('Sub-leading jet Pt distribution'  , 20, 20, (30,490)),
#    'h_jet3Pt_CR6'       : ('Third Jet Pt distribution'    , 20, 20, (30,230) ),
#    'h_muonjet1Pt_CR6'     : ('Leading jet Pt distribution' , 20, 5, (30,490) ),
#    'h_muonjet2Pt_CR6'     : ('Sub-leading jet Pt distribution'  , 20, 5, (30,490)),
#    'h_muonjet3Pt_CR6'     : ('Third Jet Pt distribution'    , 20, 5, (30,230) ),
#    'h_electronjet1Pt_CR6' : ('Leading jet Pt distribution' , 20, 5, (30,490) ),
#    'h_electronjet2Pt_CR6' : ('Sub-leading jet Pt distribution'  , 20, 5, (30,490)),
#    'h_electronjet3Pt_CR6' : ('Third Jet Pt distribution'    , 20, 5, (30,230) ),
#
#    #CR6nw
#    'metFinal_CR6nw'             : ('E^{miss}_{T} [GeV]'                          , 1000, met_rebin_had, met_range_had ),
#    'muonmetFinal_CR6nw'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'electronmetFinal_CR6nw'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#
#    'dphi_CR6nw'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'dphi_6j_CR6nw'         : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_CR6nw'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_6j_CR6nw'      : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'lectrondphi_CR6nw'      : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'lectrondphi_6j_CR6nw'   : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#
#    'h_nJets_CR6nw'        : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_nbJets_CR6nw'       : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_muonnJets_CR6nw'      : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_muonnbJets_CR6nw'     : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_electronnJets_CR6nw'  : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_electronnbJets_CR6nw' : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#
#    'h_jet1Pt_CR6nw'       : ('Leading jet Pt distribution' , 20, 20, (30,490) ),
#    'h_jet2Pt_CR6nw'       : ('Sub-leading jet Pt distribution'  , 20, 20, (30,490)),
#    'h_jet3Pt_CR6nw'       : ('Third Jet Pt distribution'    , 20, 20, (30,230) ),
#    'h_muonjet1Pt_CR6nw'     : ('Leading jet Pt distribution' , 20, 5, (30,490) ),
#    'h_muonjet2Pt_CR6nw'     : ('Sub-leading jet Pt distribution'  , 20, 5, (30,490)),
#    'h_muonjet3Pt_CR6nw'     : ('Third Jet Pt distribution'    , 20, 5, (30,230) ),
#    'h_electronjet1Pt_CR6nw' : ('Leading jet Pt distribution' , 20, 5, (30,490) ),
#    'h_electronjet2Pt_CR6nw' : ('Sub-leading jet Pt distribution'  , 20, 5, (30,490)),
#    'h_electronjet3Pt_CR6nw' : ('Third Jet Pt distribution'    , 20, 5, (30,230) ),
#
#    #CR7
#    #'metFinal_CR7'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'muonmetFinal_CR7'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'electronmetFinal_CR7'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    #'mixmetFinal_CR7'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#
#    'dphi_CR7'             : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'dphi_6j_CR7'          : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_CR7'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_6j_CR7'      : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'electrondphi_CR7'     : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'electrondphi_6j_CR7'  : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    #'mixdphi_CR7'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 1, None ),
#    #'mixdphi_6j_CR7'         : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 1, None ),
#    
#    'h_mu1Pt_CR7'     : ('Leading muon Pt distribution' , 20, 40, (30,230) ),
#    'h_mu2Pt_CR7'     : ('Sub-leading muon Pt distribution' , 20, 40, (30,230) ),
#    'h_el1Pt_CR7'     : ('Leading electron Pt distribution' , 20, 40, (30,230) ),
#    'h_el2Pt_CR7'     : ('Sub-leading electron Pt distribution' , 20, 40, (30,230) ),
#    
#    'h_nJets_CR7'          : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_nbJets_CR7'         : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_muonnJets_CR7'      : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_muonnbJets_CR7'     : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_electronnJets_CR7'  : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_electronnbJets_CR7' : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_mixnJets_CR7'        : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_mixnbJets_CR7'       : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    
#    'h_jet1Pt_CR7'         : ('Leading jet Pt distribution' , 20, 50, (30,490) ),
#    'h_jet2Pt_CR7'         : ('Sub-leading jet Pt distribution'  , 20, 50, (30,490)),
#    'h_jet3Pt_CR7'         : ('Third Jet Pt distribution'    , 20, 50, (30,230) ),
#    'h_muonjet1Pt_CR7'     : ('Leading jet Pt distribution' , 20, 50, (30,490) ),
#    'h_muonjet2Pt_CR7'     : ('Sub-leading jet Pt distribution'  , 20, 50, (30,490)),
#    'h_muonjet3Pt_CR7'     : ('Third Jet Pt distribution'    , 20, 50, (30,230) ),
#    'h_electronjet1Pt_CR7' : ('Leading jet Pt distribution' , 20, 50, (30,490) ),
#    'h_electronjet2Pt_CR7' : ('Sub-leading jet Pt distribution'  , 20, 50, (30,490)),
#    'h_electronjet3Pt_CR7' : ('Third Jet Pt distribution'    , 20, 50, (30,230) ),
#    'h_mixjet1Pt_CR7'       : ('Leading jet Pt distribution' , 20, 20, (30,490) ),
#    'h_mixjet2Pt_CR7'       : ('Sub-leading jet Pt distribution'  , 20, 20, (30,490)),
#    'h_mixjet3Pt_CR7'       : ('Third Jet Pt distribution'    , 20, 20, (30,230) ),
#
#    #CR7nw
#    'metFinal_CR7nw'             : ('E^{miss}_{T} [GeV]'                          , 1000, met_rebin_had, met_range_had ),
#    'muonmetFinal_CR7nw'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'electronmetFinal_CR7nw'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'mixmetFinal_CR7nw'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#
#    'dphi_CR7nw'             : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'dphi_6j_CR7nw'          : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_CR7nw'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_6j_CR7nw'      : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'electrondphi_CR7nw'     : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'electrondphi_6j_CR7nw'  : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    #'mixdphi_CR7nw'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 1, None ),
#    #'mixdphi_6j_CR7nw'         : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 1, None ),
#    
#    'h_mu1Pt_CR7nw'     : ('Leading muon Pt distribution' , 20, 40, (30,230) ),
#    'h_mu2Pt_CR7nw'     : ('Sub-leading muon Pt distribution' , 20, 40, (30,230) ),
#    'h_el1Pt_CR7nw'     : ('Leading electron Pt distribution' , 20, 40, (30,230) ),
#    'h_el2Pt_CR7nw'     : ('Sub-leading electron Pt distribution' , 20, 40, (30,230) ),
#    
#    'h_nJets_CR7nw'          : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_nbJets_CR7nw'         : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_muonnJets_CR7nw'      : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_muonnbJets_CR7nw'     : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_electronnJets_CR7nw'  : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_electronnbJets_CR7nw' : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_mixnJets_CR7nw'        : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_mixnbJets_CR7nw'       : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    
#    'h_jet1Pt_CR7nw'         : ('Leading jet Pt distribution' , 20, 50, (30,490) ),
#    'h_jet2Pt_CR7nw'         : ('Sub-leading jet Pt distribution'  , 20, 50, (30,490)),
#    'h_jet3Pt_CR7nw'         : ('Third Jet Pt distribution'    , 20, 50, (30,230) ),
#    'h_muonjet1Pt_CR7nw'     : ('Leading jet Pt distribution' , 20, 50, (30,490) ),
#    'h_muonjet2Pt_CR7nw'     : ('Sub-leading jet Pt distribution'  , 20, 50, (30,490)),
#    'h_muonjet3Pt_CR7nw'     : ('Third Jet Pt distribution'    , 20, 50, (30,230) ),
#    'h_electronjet1Pt_CR7nw' : ('Leading jet Pt distribution' , 20, 50, (30,490) ),
#    'h_electronjet2Pt_CR7nw' : ('Sub-leading jet Pt distribution'  , 20, 50, (30,490)),
#    'h_electronjet3Pt_CR7nw' : ('Third Jet Pt distribution'    , 20, 50, (30,230) ),
#    'h_mixjet1Pt_CR7nw'       : ('Leading jet Pt distribution' , 20, 20, (30,490) ),
#    'h_mixjet2Pt_CR7nw'       : ('Sub-leading jet Pt distribution'  , 20, 20, (30,490)),
#    'h_mixjet3Pt_CR7nw'       : ('Third Jet Pt distribution'    , 20, 20, (30,230) ),
#
}


store = [
    "h_nMu", #SR shape
    "h_nJets", #SR shape
    "h_2j1t_mtw",
    "h_2j1t_topPt",
    "h_3j1t_mtw",
    "h_3j2t_mtw",

    ]
