/**				\
 * DMAnalysisTreeMaker			\
 * \
 * Produces analysis trees from edm-ntuples adding extra variables for resolved and unresolved tops\
 * For custom systematics scenarios\
 * \
 * \\Author A. Orso M. Iorio\
 * \
 * \
 *\\version  $Id:\
 * \
 * \
*/ 

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h" 
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <Math/VectorUtil.h>
#include "./MT2Utility.h"
#include "./mt2w_bisect.h"
#include "./mt2bl_bisect.h"
#include "./Mt2Com_bisect.h"
#include "./DMTopVariables.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <vector>
#include <algorithm>
#include <TLorentzVector.h>
#include <TMVA/Reader.h>
#include <string>
#include <iostream>

//using namespace reco;
using namespace edm;
using namespace std;

namespace LHAPDF
{
void initPDFSet(int nset, const std::string &filename, int member = 0);
int numberPDF(int nset);
void usePDFMember(int nset, int member);
double xfx(int nset, double x, double Q, int fl);
double getXmin(int nset, int member);
double getXmax(int nset, int member);
double getQ2min(int nset, int member);
double getQ2max(int nset, int member);
void extrapolate(bool extrapolate = true);
}

class  DMAnalysisTreeMaker : public edm::EDAnalyzer 
{
public:
  explicit DMAnalysisTreeMaker( const edm::ParameterSet & );   

private:

  virtual void analyze(const edm::Event &, const edm::EventSetup & );
  virtual void beginRun(const edm::Run  &, const edm::EventSetup & );
  virtual void endJob();
  vector<string> additionalVariables(string);
  string makeName(string label,string pref,string var);
  string makeBranchNameCat(string label, string cat, string pref, string var);
  string makeBranchName(string label, string pref, string var);
  void initializePdf(string centralpdf,string variationpdf);
  void initTreeWeightHistory(bool useLHEW);
  void getEventPdf();
  int eventFlavour(bool getFlavour, int nb, int nc,int nudsg);
  bool flavourFilter(string ch, int nb, int nc,int nudsg); 

  void initCategoriesSize(string label);
  void setCategorySize(string label, string category, size_t size);
  void fillCategory(string label, string category, int pos_nocat, int pos_cat);

  double getWPtWeight(double ptW);
  double getZPtWeight(double ptZ);

  double getWEWKPtWeight(double ptW);
  double getZEWKPtWeight(double ptZ);
  double getAPtWeight(double ptA);
  double getTopPtWeight(double ptT,double ptTbar, bool extrap = false);
  bool getEventTriggers();
  bool getMETFilters();
  void getEventLHEWeights();

  bool isEWKID(int id);
  //
  // Set up MVA reader
  //
  // spectator variables, not used for MVA evaluation
  int isSig, b_mis, w_mis, wb_mis;
  float mtop;

  double jetUncertainty(double pt, double eta, string syst);
  double jetUncertainty8(double pt, double eta, string syst);
  double smear(double pt, double genpt, double eta, string syst);
  double getEffectiveArea(string particle, double eta);
  double resolSF(double eta, string syst);
  double getScaleFactor(double pt, double eta, double partonFlavour, string syst);
  double pileUpSF(string syst);
  double nInitEvents;
  float nTightJets;

  bool isInVector(std::vector<std::string> v, std::string s);
  bool isMCWeightName(std::string s);
  std::vector<edm::ParameterSet > physObjects;
  std::vector<edm::InputTag > variablesFloat, variablesInt, singleFloat,  singleInt;
  
  std::vector<edm::InputTag > variablesDouble, singleDouble;
  
  edm::EDGetTokenT< LHEEventProduct > t_lhes_;
  edm::EDGetTokenT <std::vector<std::vector<int>>> t_jetKeys_, t_muKeys_;
  edm::EDGetTokenT< GenEventInfoProduct > t_genprod_;

  edm::EDGetTokenT< std::vector<string> > t_triggerNames_;
  edm::EDGetTokenT< std::vector<float> > t_triggerBits_;
  edm::EDGetTokenT< std::vector<int> > t_triggerPrescales_;
  edm::EDGetTokenT<std::vector<string> > t_triggerNamesR_;
  
  edm::EDGetTokenT< unsigned int > t_lumiBlock_;
  edm::EDGetTokenT< unsigned int > t_runNumber_;
  edm::EDGetTokenT< ULong64_t > t_eventNumber_;

  edm::EDGetTokenT< bool > t_BadPFMuonFilter_;
  edm::EDGetTokenT< bool > t_BadChargedCandidateFilter_;
  
  edm::EDGetTokenT< std::vector<string> > t_metNames_;
  edm::EDGetTokenT< std::vector<float> > t_metBits_;

  edm::EDGetTokenT< std::vector<float> > t_pvZ_,t_pvChi2_,t_pvRho_;
  edm::EDGetTokenT< std::vector<int> > t_pvNdof_;

  edm::EDGetTokenT< double > t_Rho_;
  edm::EDGetTokenT<int> t_ntrpu_;

  edm::EDGetTokenT< std::vector<float> > jetAK8vSubjetIndex0;
  edm::EDGetTokenT< std::vector<float> > jetAK8vSubjetIndex1;

  edm::EDGetTokenT< std::vector<float> > genPartID ;
  edm::EDGetTokenT< std::vector<float> > genPartStatus;
  edm::EDGetTokenT< std::vector<float> > genPartMom0ID; 
  edm::EDGetTokenT< std::vector<float> > genPartPt;
  edm::EDGetTokenT< std::vector<float> > genPartPhi;
  edm::EDGetTokenT< std::vector<float> > genPartEta;
  edm::EDGetTokenT< std::vector<float> > genPartE;

  edm::LumiReWeighting LumiWeights_, LumiWeightsUp_, LumiWeightsDown_;
  
  edm::EDGetTokenT<reco::GenParticleCollection> t_genParticleCollection_;

  TH1D * nInitEventsHisto;
  TTree * treesBase;
  map<string, TTree * > trees;
  std::vector<string> names;
  std::vector<string> systematics;
  map< string , float[100] > vfloats_values;
  map< string , int[100] > vints_values;
  map< string , vector<string> > obj_to_floats,obj_to_ints, obj_to_doubles;
  map< string , string > obs_to_obj;
  map< string , string > obj_to_pref;
  map< string , std::vector<string> > obj_cats;

  map< string , double[100] > vdoubles_values;
  map< string , double[100] > vdouble_values;
  map< string , double > double_values;


  map< string , float > float_values;
  map< string , int > int_values;
  //  map< string , ULong64_t > int_values;
  map< string , int > sizes;

  map< string , bool > got_label; 
  map< string , int > max_instances; 
  map< int, int > subj_jet_map;

  map<string, edm::Handle<std::vector<float> > > h_floats;
  map<string, edm::Handle<std::vector<int> > > h_ints;
  map<string, edm::Handle<float> > h_float;
  map<string, edm::Handle<int> >h_int;
  
  map<string, edm::Handle<std::vector<double> > > h_doubles;
  map<string, edm::Handle<double> > h_double;
  
  map<string, edm::EDGetTokenT< std::vector<float> >  > t_floats;
  map<string, edm::EDGetTokenT< std::vector<int> > > t_ints;
  map<string, edm::EDGetTokenT<float>  > t_float;
  map<string, edm::EDGetTokenT<int> >t_int;
  map<string, edm::EDGetTokenT<std::vector<double> > > t_doubles;
  map<string, edm::EDGetTokenT<double> > t_double;
  

  string gen_label,  mu_label, ele_label, jets_label, boosted_tops_label, boosted_tops_subjets_label, met_label, photon_label;//metNoHF_label

  bool getPartonW, getPartonTop, doWReweighting, doTopReweighting;
  bool getParticleWZ, getWZFlavour;
  
  //Do resolved top measurement:
  bool doResolvedTopHad,doResolvedTopSemiLep;
  int max_leading_jets_for_top;
  int max_bjets_for_top;
  int max_genparticles;
  int n0;

  bool recalculateEA;

  //MC info:
  edm::ParameterSet channelInfo;
  std::string channel;
  double crossSection, originalEvents;
  bool useLHEWeights, useLHE, useTriggers,cutOnTriggers, useMETFilters, addPV;
  bool addLHAPDFWeights;
  string centralPdfSet,variationPdfSet;
  std::vector<string> SingleElTriggers, SingleMuTriggers, PhotonTriggers, hadronicTriggers,metFilters, HadronPFHT900Triggers, HadronPFHT800Triggers, HadronPFJet450Triggers;
  int maxPdf, maxWeights;
  edm::Handle<LHEEventProduct > lhes;
  edm::Handle<GenEventInfoProduct> genprod;

  edm::Handle<std::vector<std::vector<int>>> jetKeys;
  edm::Handle<std::vector<std::vector<int>>> muKeys;
  
  //Trigger info
  edm::Handle<std::vector<float> > triggerBits;
  edm::Handle<std::vector<string> > triggerNames;
  edm::Handle<std::vector<int> > triggerPrescales;
  edm::Handle<std::vector<string> > triggerNamesR ;

  edm::Handle<bool> BadPFMuonFilter;
  edm::Handle<bool> BadChargedCandidateFilter;
  
  edm::Handle<unsigned int> lumiBlock;
  edm::Handle<unsigned int> runNumber;
  edm::Handle<ULong64_t> eventNumber;

  edm::Handle<std::vector<float> > metBits;
  edm::Handle<std::vector<string> > metNames;
  edm::InputTag metNames_;
  
  edm::Handle<std::vector<float> > pvZ,pvChi2,pvRho;
  edm::Handle<std::vector<int> > pvNdof;

  edm::Handle<std::vector<float> > ak8jetvSubjetIndex0;
  edm::Handle<std::vector<float> > ak8jetvSubjetIndex1;
  
  float nPV;
  edm::Handle<int> ntrpu;

  //JEC info
  bool changeJECs;
  bool isData, applyRes;
  bool isV2;
  string EraLabel;
  edm::Handle<double> rho;
  double Rho;
  edm::Handle<reco::GenParticleCollection> genParticles;   
  
  //edm::Handle<double> Rho;
  std::vector<double> jetScanCuts;
  std::vector<JetCorrectorParameters> jecPars, jecParsL1_vect;
  JetCorrectorParameters *jecParsL1, *jecParsL1RC, *jecParsL2, *jecParsL3, *jecParsL2L3Residuals;
  JetCorrectionUncertainty *jecUnc;
  FactorizedJetCorrector *jecCorr, *jecCorr_L1;

  std::vector<JetCorrectorParameters> jecPars8,  jecParsL1_vect8, jecParsNoL1_vect8;
  JetCorrectorParameters *jecParsL18, *jecParsL1RC8, *jecParsL28, *jecParsL38, *jecParsL2L3Residuals8;
  JetCorrectionUncertainty *jecUnc8;
  FactorizedJetCorrector *jecCorr8,  *jecCorr_L18, *jecCorr_NoL18 ;

  bool isFirstEvent;
  //Do preselection
  bool doPreselection;
  
  class BTagWeight
  {
  private:
    int minTags;
    int maxTags;
  public:
    struct JetInfo
    {
      JetInfo(float mceff, float datasf) : eff(mceff), sf(datasf) {}
      float eff;
      float sf;
    };
    BTagWeight():
      minTags(0), maxTags(0)
    {
      ;
    }
    BTagWeight(int jmin, int jmax) :
      minTags(jmin) , maxTags(jmax) {}
    bool filter(int t);
    float weight(vector<JetInfo> jets, int tags);
    float weightWithVeto(vector<JetInfo> jetsTags, int tags, vector<JetInfo> jetsVetoes, int vetoes);
  };
  vector<BTagWeight::JetInfo> jsfscsvt, 
    jsfscsvt_b_tag_up, 
    jsfscsvt_b_tag_down, 
    jsfscsvt_mistag_up, 
    jsfscsvt_mistag_down;

  vector<BTagWeight::JetInfo> jsfscsvm, 
    jsfscsvm_b_tag_up, 
    jsfscsvm_b_tag_down, 
    jsfscsvm_mistag_up, 
    jsfscsvm_mistag_down;
  
  vector<BTagWeight::JetInfo> jsfscsvl, 
    jsfscsvl_b_tag_up, 
    jsfscsvl_b_tag_down, 
    jsfscsvl_mistag_up, 
    jsfscsvl_mistag_down;
  
  vector<BTagWeight::JetInfo> jsfscsvt_subj, 
    jsfscsvt_subj_b_tag_up, 
    jsfscsvt_subj_b_tag_down, 
    jsfscsvt_subj_mistag_up, 
    jsfscsvt_subj_mistag_down;

  vector<BTagWeight::JetInfo> jsfscsvm_subj, 
    jsfscsvm_subj_b_tag_up, 
    jsfscsvm_subj_b_tag_down, 
    jsfscsvm_subj_mistag_up, 
    jsfscsvm_subj_mistag_down;
  
  vector<BTagWeight::JetInfo> jsfscsvl_subj, 
    jsfscsvl_subj_b_tag_up, 
    jsfscsvl_subj_b_tag_down, 
    jsfscsvl_subj_mistag_up, 
    jsfscsvl_subj_mistag_down;

  //tight AK4 b-tagging weights
  
  BTagWeight b_csvt_0_tags= BTagWeight(0,0),
    b_csvt_1_tag= BTagWeight(1,10),
    b_csvt_2_tag= BTagWeight(2,10);
  double b_weight_csvt_0_tags,
    b_weight_csvt_1_tag,
    b_weight_csvt_2_tag;
  double b_weight_csvt_0_tags_mistag_up,
    b_weight_csvt_1_tag_mistag_up,
    b_weight_csvt_2_tag_mistag_up;
  double b_weight_csvt_0_tags_mistag_down,
    b_weight_csvt_1_tag_mistag_down,
    b_weight_csvt_2_tag_mistag_down;
  double b_weight_csvt_0_tags_b_tag_down,
    b_weight_csvt_1_tag_b_tag_down,
    b_weight_csvt_2_tag_b_tag_down;
  double b_weight_csvt_0_tags_b_tag_up,
    b_weight_csvt_1_tag_b_tag_up,
    b_weight_csvt_2_tag_b_tag_up;

  //medium AK4 b-tagging weights
    
  BTagWeight b_csvm_0_tags= BTagWeight(0,0),
    b_csvm_1_tag= BTagWeight(1,10),
    b_csvm_2_tag= BTagWeight(2,10);
  double b_weight_csvm_0_tags,
    b_weight_csvm_1_tag,
    b_weight_csvm_2_tag;
  double b_weight_csvm_0_tags_mistag_up,
    b_weight_csvm_1_tag_mistag_up,
    b_weight_csvm_2_tag_mistag_up;
  double b_weight_csvm_0_tags_mistag_down,
    b_weight_csvm_1_tag_mistag_down,
    b_weight_csvm_2_tag_mistag_down;
  double b_weight_csvm_0_tags_b_tag_down,
    b_weight_csvm_1_tag_b_tag_down,
    b_weight_csvm_2_tag_b_tag_down;
  double b_weight_csvm_0_tags_b_tag_up,
    b_weight_csvm_1_tag_b_tag_up,
    b_weight_csvm_2_tag_b_tag_up;
  //medium subjets b-tagging weights
    
  BTagWeight b_subj_csvm_0_tags= BTagWeight(0,0),
    b_subj_csvm_1_tag= BTagWeight(1,10),
    b_subj_csvm_0_1_tags= BTagWeight(0,10),
    b_subj_csvm_2_tags= BTagWeight(2,10);

  double b_weight_subj_csvm_0_tags,
    b_weight_subj_csvm_1_tag,
    b_weight_subj_csvm_0_1_tags,
    b_weight_subj_csvm_2_tags;
  double b_weight_subj_csvm_0_tags_mistag_up,
    b_weight_subj_csvm_1_tag_mistag_up,
    b_weight_subj_csvm_0_1_tags_mistag_up,
    b_weight_subj_csvm_2_tags_mistag_up;
  double b_weight_subj_csvm_0_tags_mistag_down,
    b_weight_subj_csvm_1_tag_mistag_down,
    b_weight_subj_csvm_0_1_tags_mistag_down,
    b_weight_subj_csvm_2_tags_mistag_down;
  double b_weight_subj_csvm_0_tags_b_tag_down,
    b_weight_subj_csvm_1_tag_b_tag_down,
    b_weight_subj_csvm_0_1_tags_b_tag_down,
    b_weight_subj_csvm_2_tags_b_tag_down;
  double b_weight_subj_csvm_0_tags_b_tag_up,
    b_weight_subj_csvm_1_tag_b_tag_up,
    b_weight_subj_csvm_0_1_tags_b_tag_up,
    b_weight_subj_csvm_2_tags_b_tag_up;

  //loose AK4 b-tagging weights
  
  BTagWeight b_csvl_0_tags= BTagWeight(0,0),
    b_csvl_1_tag= BTagWeight(1,10),
    b_csvl_2_tag= BTagWeight(2,10);
  
  double b_weight_csvl_0_tags_mistag_up,
    b_weight_csvl_1_tag_mistag_up,
    b_weight_csvl_2_tag_mistag_up;
  double b_weight_csvl_0_tags_mistag_down,
    b_weight_csvl_1_tag_mistag_down,
    b_weight_csvl_2_tag_mistag_down;
  double b_weight_csvl_0_tags_b_tag_down,
    b_weight_csvl_1_tag_b_tag_down,
    b_weight_csvl_2_tag_b_tag_down;
  double b_weight_csvl_0_tags_b_tag_up,
    b_weight_csvl_1_tag_b_tag_up,
    b_weight_csvl_2_tag_b_tag_up;
  double b_weight_csvl_0_tags,
    b_weight_csvl_1_tag,
    b_weight_csvl_2_tag;

  //loose subjets b-tagging weights

  BTagWeight b_subj_csvl_0_tags= BTagWeight(0,0),
    b_subj_csvl_1_tag= BTagWeight(1,10),
    b_subj_csvl_0_1_tags= BTagWeight(0,10),
    b_subj_csvl_2_tags= BTagWeight(2,10);
  
  double b_weight_subj_csvl_0_tags_mistag_up,
    b_weight_subj_csvl_1_tag_mistag_up,
    b_weight_subj_csvl_0_1_tags_mistag_up,
    b_weight_subj_csvl_2_tags_mistag_up;
  double b_weight_subj_csvl_0_tags_mistag_down,
    b_weight_subj_csvl_1_tag_mistag_down,
    b_weight_subj_csvl_0_1_tags_mistag_down,
    b_weight_subj_csvl_2_tags_mistag_down;
  double b_weight_subj_csvl_0_tags_b_tag_down,
    b_weight_subj_csvl_1_tag_b_tag_down,
    b_weight_subj_csvl_0_1_tags_b_tag_down,
    b_weight_subj_csvl_2_tags_b_tag_down;
  double b_weight_subj_csvl_0_tags_b_tag_up,
    b_weight_subj_csvl_1_tag_b_tag_up,
    b_weight_subj_csvl_0_1_tags_b_tag_up,
    b_weight_subj_csvl_2_tags_b_tag_up;
  double b_weight_subj_csvl_0_tags,
    b_weight_subj_csvl_1_tag,
    b_weight_subj_csvl_0_1_tags,
    b_weight_subj_csvl_2_tags;
  
  double MCTagEfficiency(string algo, int flavor, double pt, double eta);
  double MCTagEfficiencySubjet(string algo, int flavor, double pt, double eta); 
  double TagScaleFactor(string algo, int flavor, string syst,double pt);
  double TagScaleFactorSubjet(string algo, int flavor, string syst,double pt);
  
  //
  bool doBTagSF;
  bool doPU;
  
  string season;
  //  season = dataPUFile_;
  
  string distr;
  
  
};


DMAnalysisTreeMaker::DMAnalysisTreeMaker(const edm::ParameterSet& iConfig){
  
  mu_label = iConfig.getParameter<std::string >("muLabel");
  ele_label = iConfig.getParameter<std::string >("eleLabel");
  jets_label = iConfig.getParameter<std::string >("jetsLabel");
  photon_label = iConfig.getParameter<std::string >("photonLabel");
  boosted_tops_label = iConfig.getParameter<std::string >("boostedTopsLabel");
  boosted_tops_subjets_label = iConfig.getParameter<std::string >("boostedTopsSubjetsLabel");
  met_label = iConfig.getParameter<std::string >("metLabel");
  gen_label = iConfig.getParameter<std::string >("genLabel");
  physObjects = iConfig.template getParameter<std::vector<edm::ParameterSet> >("physicsObjects");
  
  channelInfo = iConfig.getParameter<edm::ParameterSet >("channelInfo"); // The physics of the channel, e.g. the cross section, #original events, etc.
  channel = channelInfo.getParameter<std::string>("channel");
  crossSection = channelInfo.getParameter<double>("crossSection");
  originalEvents = channelInfo.getParameter<double>("originalEvents");

  doPreselection = iConfig.getUntrackedParameter<bool>("doPreselection",true);
  doPU = iConfig.getUntrackedParameter<bool>("doPU",true);

  useLHEWeights = channelInfo.getUntrackedParameter<bool>("useLHEWeights",false);
  useLHE = channelInfo.getUntrackedParameter<bool>("useLHE",false);
  addLHAPDFWeights = channelInfo.getUntrackedParameter<bool>("addLHAPDFWeights",false);

  getPartonW = channelInfo.getUntrackedParameter<bool>("getPartonW",false);
  getParticleWZ = channelInfo.getUntrackedParameter<bool>("getParticleWZ",false);
  getPartonTop = channelInfo.getUntrackedParameter<bool>("getPartonTop",false);
  doWReweighting = channelInfo.getUntrackedParameter<bool>("doWReweighting",false);

  getWZFlavour = channelInfo.getUntrackedParameter<bool>("getWZFlavour",false);

  doResolvedTopSemiLep = iConfig.getUntrackedParameter<bool>("doResolvedTopSemiLep",false);
  doResolvedTopHad = iConfig.getUntrackedParameter<bool>("doResolvedTopHad",false);

  edm::InputTag genprod_ = iConfig.getParameter<edm::InputTag>( "genprod" );
  t_genprod_ = consumes<GenEventInfoProduct>( genprod_ );
  
  useTriggers = iConfig.getUntrackedParameter<bool>("useTriggers",true);
  cutOnTriggers = iConfig.getUntrackedParameter<bool>("cutOnTriggers",true);

  edm::InputTag ak8jetvSubjetIndex0_ = iConfig.getParameter<edm::InputTag>("ak8jetvSubjetIndex0");
  jetAK8vSubjetIndex0  = consumes< std::vector<float> >( ak8jetvSubjetIndex0_);
  edm::InputTag ak8jetvSubjetIndex1_ = iConfig.getParameter<edm::InputTag>("ak8jetvSubjetIndex1");
  jetAK8vSubjetIndex1 = consumes< std::vector<float> >( ak8jetvSubjetIndex1_);

  edm::InputTag lumiBlock_ = iConfig.getParameter<edm::InputTag>("lumiBlock");
  t_lumiBlock_ = consumes< unsigned int >( lumiBlock_ );
  edm::InputTag runNumber_ = iConfig.getParameter<edm::InputTag>("runNumber");
  t_runNumber_ = consumes< unsigned int >( runNumber_ );
  edm::InputTag eventNumber_ = iConfig.getParameter<edm::InputTag>("eventNumber");
  t_eventNumber_ = consumes< ULong64_t >( eventNumber_ );
  
  if(useTriggers){
     edm::InputTag triggerBits_ = iConfig.getParameter<edm::InputTag>("triggerBits");
    t_triggerBits_ = consumes< std::vector<float> >( triggerBits_ );
    edm::InputTag triggerPrescales_ = iConfig.getParameter<edm::InputTag>("triggerPrescales");
    t_triggerPrescales_ = consumes< std::vector<int> >( triggerPrescales_ );

    edm::InputTag triggerNamesR_ = iConfig.getParameter<edm::InputTag>("triggerNames");
    t_triggerNamesR_ = mayConsume< std::vector<string>, edm::InRun>(edm::InputTag("TriggerUserData","triggerNameTree"));

    SingleElTriggers= channelInfo.getParameter<std::vector<string> >("SingleElTriggers");
    SingleMuTriggers= channelInfo.getParameter<std::vector<string> >("SingleMuTriggers");
    PhotonTriggers= channelInfo.getParameter<std::vector<string> >("PhotonTriggers");
    hadronicTriggers= channelInfo.getParameter<std::vector<string> >("hadronicTriggers");
    HadronPFHT900Triggers= channelInfo.getParameter<std::vector<string> >("HadronPFHT900Triggers");
    HadronPFHT800Triggers= channelInfo.getParameter<std::vector<string> >("HadronPFHT800Triggers");
    HadronPFJet450Triggers= channelInfo.getParameter<std::vector<string> >("HadronPFJet450Triggers");
  }

  useMETFilters = iConfig.getUntrackedParameter<bool>("useMETFilters",true);

  if(useMETFilters){
    metFilters = channelInfo.getParameter<std::vector<string> >("metFilters");
    edm::InputTag BadChargedCandidateFilter_ = iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter");
    t_BadChargedCandidateFilter_ = consumes< bool >( BadChargedCandidateFilter_ );
    edm::InputTag BadPFMuonFilter_ = iConfig.getParameter<edm::InputTag>("BadPFMuonFilter");
    t_BadPFMuonFilter_ = consumes< bool >( BadPFMuonFilter_ );
    edm::InputTag metBits_ = iConfig.getParameter<edm::InputTag>("metBits");
    t_metBits_ = consumes< std::vector<float> >( metBits_ );
    metNames_ = iConfig.getParameter<edm::InputTag>("metNames");
    t_metNames_ = consumes< std::vector<string>, edm::InRun >( metNames_ );
  }
  
  addPV = iConfig.getUntrackedParameter<bool>("addPV",true);
  changeJECs = iConfig.getUntrackedParameter<bool>("changeJECs",false);
  recalculateEA = iConfig.getUntrackedParameter<bool>("recalculateEA",true);
  
  EraLabel = iConfig.getUntrackedParameter<std::string>("EraLabel");

  isData = iConfig.getUntrackedParameter<bool>("isData",false);
  applyRes = iConfig.getUntrackedParameter<bool>("applyRes",false);
  isV2= iConfig.getUntrackedParameter<bool>("isV2",false);
  
  t_Rho_ = consumes<double>( edm::InputTag( "fixedGridRhoFastjetAll" ) ) ;
  t_genParticleCollection_ = consumes<reco::GenParticleCollection>(edm::InputTag("filteredPrunedGenParticles"));

  if(addPV || changeJECs){
    edm::InputTag pvZ_ = iConfig.getParameter<edm::InputTag >("vertexZ");
    t_pvZ_ = consumes< std::vector<float> >( pvZ_ );
    edm::InputTag pvChi2_ = iConfig.getParameter<edm::InputTag >("vertexChi2");
    t_pvChi2_ = consumes< std::vector<float> >( pvChi2_ );
    edm::InputTag pvRho_ = iConfig.getParameter<edm::InputTag >("vertexRho");
    t_pvRho_ = consumes< std::vector<float> >( pvRho_ );
    edm::InputTag pvNdof_ = iConfig.getParameter<edm::InputTag >("vertexNdof");
    t_pvNdof_ = consumes< std::vector< int > >( pvNdof_ );
  }
  
  if (doPU)
    t_ntrpu_ = consumes< int >( edm::InputTag( "eventUserData","puNtrueInt" ) );

  maxWeights = 9;
  if(useLHEWeights){
    maxWeights = channelInfo.getUntrackedParameter<int>("maxWeights",9);//Usually we do have 9 weights for the scales, might vary depending on the lhe
  }
  if(addLHAPDFWeights){
    centralPdfSet = channelInfo.getUntrackedParameter<string>("pdfSet","NNPDF");
    variationPdfSet = channelInfo.getUntrackedParameter<string>("pdfSet","NNPDF");
    initializePdf(centralPdfSet,variationPdfSet);

  }
  if(doResolvedTopHad){
    max_leading_jets_for_top  = iConfig.getUntrackedParameter<int>("maxLeadingJetsForTop",8);//Take the 8 leading jets for the top permutations
  }
  if(doResolvedTopSemiLep){
    max_bjets_for_top  = iConfig.getUntrackedParameter<int>("maxBJetsForTop",2);//Take the 8 leading jets for the top permutations
  }
  if(!isData){
    max_genparticles  = iConfig.getUntrackedParameter<int>("maxGenPart",30);
   }
  systematics = iConfig.getParameter<std::vector<std::string> >("systematics");

  jetScanCuts = iConfig.getParameter<std::vector<double> >("jetScanCuts");

  std::vector<edm::ParameterSet >::const_iterator itPsets = physObjects.begin();

  bool addNominal=false;
  for (size_t s = 0; s<systematics.size();++s){
    if(systematics.at(s).find("noSyst")!=std::string::npos) {
      addNominal=true;
      break;
    }
  }
  if(systematics.size()==0){
    addNominal=true;
    systematics.push_back("noSyst");
  }//In case there's no syst specified, do the nominal scenario
  //addNominal=true;
  Service<TFileService> fs;
  TFileDirectory DMTrees;

  if(addNominal){
    DMTrees = fs->mkdir( "systematics_trees" );
  }

  trees["noSyst"] =  new TTree((channel+"__noSyst").c_str(),(channel+"__noSyst").c_str());

  nInitEventsHisto = new TH1D("initialEvents","initalEvents",10,0,10);

  edm::InputTag lhes_ = iConfig.getParameter<edm::InputTag>( "lhes" );
  t_lhes_ = consumes< LHEEventProduct >( lhes_ );
  
  edm::InputTag jetKeys_ = iConfig.getParameter<edm::InputTag >("jetKeysAK4CHS");
  edm::InputTag muKeys_ = iConfig.getParameter<edm::InputTag >("muonKeys");
  t_jetKeys_ = consumes<std::vector<std::vector<int>>> (jetKeys_);
  t_muKeys_ = consumes<std::vector<std::vector<int>>> (muKeys_);
  
  for (;itPsets!=physObjects.end();++itPsets){ 
    int maxI = itPsets->getUntrackedParameter< int >("maxInstances",10);
    variablesFloat = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesF"); 
    variablesInt = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesI");
    singleFloat = itPsets->template getParameter<std::vector<edm::InputTag> >("singleF"); 
    singleDouble = itPsets->template getParameter<std::vector<edm::InputTag> >("singleD"); 
    singleInt = itPsets->template getParameter<std::vector<edm::InputTag> >("singleI"); 
    string namelabel = itPsets->getParameter< string >("label");
    string nameprefix = itPsets->getParameter< string >("prefix");
    bool saveBaseVariables = itPsets->getUntrackedParameter<bool>("saveBaseVariables",true);
    bool saveNoCat = itPsets->getUntrackedParameter<bool>("saveNoCat",true);

    std::vector<std::string > categories = itPsets->getParameter<std::vector<std::string> >("categories");
    std::vector<std::string > toSave= itPsets->getParameter<std::vector<std::string> >("toSave");
    
    std::vector<edm::InputTag >::const_iterator itF = variablesFloat.begin();
    std::vector<edm::InputTag >::const_iterator itI = variablesInt.begin();
    std::vector<edm::InputTag >::const_iterator itsF = singleFloat.begin();
    std::vector<edm::InputTag >::const_iterator itsD = singleDouble.begin();
    std::vector<edm::InputTag >::const_iterator itsI = singleInt.begin();

    for(size_t sc = 0; sc< categories.size() ;++sc){
      string category = categories.at(sc);
      obj_cats[namelabel].push_back(category);
    }
    stringstream max_instance_str;
    max_instance_str<<maxI;
    max_instances[namelabel]=maxI;
    string nameobs = namelabel;
    string prefix = nameprefix;

    for (;itF != variablesFloat.end();++itF){
      
      string name=itF->instance()+"_"+itF->label();
      string nameinstance=itF->instance();
      string nameshort=itF->instance();
      
      string nametobranch = makeBranchName(namelabel,prefix,nameinstance);
      
      name = nametobranch;
      nameshort = nametobranch;
    
      if(saveNoCat && (saveBaseVariables|| isInVector(toSave,itF->instance()))) trees["noSyst"]->Branch(nameshort.c_str(), &vfloats_values[name],(nameshort+"["+max_instance_str.str()+"]/F").c_str());
      names.push_back(name);
      obj_to_floats[namelabel].push_back(name);
      obs_to_obj[name] = nameobs;
      obj_to_pref[nameobs] = prefix;

      t_floats[ name ] = consumes< std::vector<float> >( *itF );
      
      for(size_t sc = 0; sc< categories.size() ;++sc){
	string category = categories.at(sc);
	string nametobranchcat = makeBranchNameCat(namelabel,category,prefix,nameinstance);
	string namecat = nametobranchcat;
	nameshort = nametobranch;
	if(saveBaseVariables|| isInVector(toSave,itF->instance())) trees["noSyst"]->Branch(namecat.c_str(), &vfloats_values[namecat],(namecat+"["+max_instance_str.str()+"]/F").c_str());
      }
    }
  
    for (;itI != variablesInt.end();++itI){
      string name=itI->instance()+"_"+itI->label();
      string nameshort=itF->instance();
      string nametobranch = makeBranchName(namelabel,prefix,nameshort);
      name = nametobranch;
      nameshort = nametobranch;

      if(saveNoCat && (saveBaseVariables|| isInVector(toSave,itI->instance())) ) trees["noSyst"]->Branch(nameshort.c_str(), &vints_values[name],(nameshort+"["+max_instance_str.str()+"]/I").c_str());
      for(size_t sc = 0; sc< categories.size() ;++sc){
	string category = categories.at(sc);
	string nametobranchcat = makeBranchNameCat(namelabel,category,prefix,nameshort);
	string namecat = nametobranchcat;
	nameshort = nametobranch;
	if(saveBaseVariables|| isInVector(toSave,itF->instance())) trees["noSyst"]->Branch(namecat.c_str(), &vfloats_values[namecat],(namecat+"["+max_instance_str.str()+"]/I").c_str());
      }

      names.push_back(name);
      obj_to_ints[namelabel].push_back(name);
      obs_to_obj[name] = nameobs;
      obj_to_pref[nameobs] = prefix;

      t_ints[ name ] = consumes< std::vector<int> >( *itI );

    }
    
    if (variablesFloat.size()>0){
      string nameshortv = namelabel;
      vector<string> extravars = additionalVariables(nameshortv);
      for(size_t addv = 0; addv < extravars.size();++addv){
	string name = nameshortv+"_"+extravars.at(addv);
	if (saveNoCat && (saveBaseVariables || isInVector(toSave, extravars.at(addv)) || isInVector(toSave, "allExtra") ) )trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+max_instance_str.str()+"]/F").c_str());
	for(size_t sc = 0; sc< categories.size() ;++sc){
	  string category = categories.at(sc);
	  string nametobranchcat = nameshortv+category+"_"+extravars.at(addv);
	  string namecat = nametobranchcat;
	  cout << "extra var "<< extravars.at(addv)<<endl;
	  cout << " namecat "<< namecat<< endl;
	  if(saveBaseVariables|| isInVector(toSave,extravars.at(addv)) || isInVector(toSave,"allExtra")) trees["noSyst"]->Branch(namecat.c_str(), &vfloats_values[namecat],(namecat+"["+max_instance_str.str()+"]/F").c_str());
	}

	obj_to_floats[namelabel].push_back(name);
	obs_to_obj[name] = nameobs;
	obj_to_pref[nameobs] = prefix;
      }
    }
    names.push_back(nameobs);
    cout << "size part: nameobs is  "<< nameobs<<endl;
    if(saveNoCat) trees["noSyst"]->Branch((nameobs+"_size").c_str(), &sizes[nameobs]);
    for(size_t sc = 0; sc< categories.size() ;++sc){
      string category = categories.at(sc);
      trees["noSyst"]->Branch((nameobs+category+"_size").c_str(), &sizes[nameobs+category]);
    }
    
    //Initialize single pset objects
     for (;itsF != singleFloat.end();++itsF){
      string name=itsF->instance()+itsF->label();
      string nameshort=itsF->instance();
      string nametobranch = makeBranchName(namelabel,prefix,nameshort);
      name = nametobranch;
      nameshort = nametobranch;
      t_float[ name ] = consumes< float >( *itsF );
      if(saveBaseVariables|| isInVector(toSave,itsF->instance()))trees["noSyst"]->Branch(nameshort.c_str(), &float_values[name]);
    }
 
    for (;itsD != singleDouble.end();++itsD){
      string name=itsD->instance()+itsD->label();
      string nameshort=itsD->instance();
      string nametobranch = makeBranchName(namelabel,prefix,nameshort);
      name = nametobranch;
      nameshort = nametobranch;
      t_double[ name ] = consumes< double >( *itsD );
      if(saveBaseVariables|| isInVector(toSave,itsD->instance()))trees["noSyst"]->Branch(nameshort.c_str(), &double_values[name]);
    }
    for (;itsI != singleInt.end();++itsI){
      string name=itsI->instance()+itsI->label();
      string nameshort=itsI->instance();
      string nametobranch = makeBranchName(namelabel,prefix,nameshort);
      name = nametobranch;
      nameshort = nametobranch;
      t_int[ name ] = consumes< int >( *itsI );
      if(saveBaseVariables|| isInVector(toSave,itsI->instance()))trees["noSyst"]->Branch(nameshort.c_str(), &int_values[name]);
    }
  }
  if(doResolvedTopSemiLep){
    string nameshortv= "resolvedTopSemiLep";
    vector<string> extravarstop = additionalVariables(nameshortv);
    double max_instances_top = max_bjets_for_top*max((int)(max_instances[ele_label]),(int)(max_instances[mu_label]) );
    max_instances[nameshortv]=max_instances_top;
    stringstream mtop;
    mtop << max_instances_top;
    cout << " max instances top is "<< max_instances_top << " max_leading_jets_for_top "<< max_leading_jets_for_top << " max_instances[jets_label]  " <<max_instances[jets_label]<<endl;
    for(size_t addv = 0; addv < extravarstop.size();++addv){
      string name = nameshortv+"_"+extravarstop.at(addv);
      trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+mtop.str()+"]/F").c_str());
    }
    trees["noSyst"]->Branch((nameshortv+"_size").c_str(), &sizes[nameshortv]);
  }
  if(doResolvedTopHad){
    string nameshortv= "resolvedTopHad";
    vector<string> extravarstop = additionalVariables(nameshortv);
    double max_instances_top = TMath::Binomial(min((int)(max_instances[jets_label]),max_leading_jets_for_top),4);
    max_instances[nameshortv]=max_instances_top;
    stringstream mtop;
    mtop << max_instances_top;
    cout << " max instances top is "<< max_instances_top<< " max_leading_jets_for_top "<< max_leading_jets_for_top << " max_instances[jets_label]  " <<max_instances[jets_label]<<endl;
    for(size_t addv = 0; addv < extravarstop.size();++addv){
      string name = nameshortv+"_"+extravarstop.at(addv);
      trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+mtop.str()+"]/F").c_str());
    }
    trees["noSyst"]->Branch((nameshortv+"_size").c_str(), &sizes[nameshortv]);
  }
  
  if(!isData){
    string nameshortv= "genPart";
    vector<string> extravarstop = additionalVariables(nameshortv);
    //double max_instances_gen = max_genparticles*max((int)(max_instances[gen_label]),(int)(max_instances[gen_label]) );
    double max_instances_gen = 30*max((int)0,(int)1);
    max_instances[nameshortv]=max_instances_gen;
    stringstream mtop;
    mtop << max_instances_gen;
    //cout << " max instances top is "<< max_instances_top<< " max_leading_jets_for_top "<< max_leading_jets_for_top << " max_instances[jets_label]  " <<max_instances[jets_label]<<endl;
    for(size_t addv = 0; addv < extravarstop.size();++addv){
      string name = nameshortv+"_"+extravarstop.at(addv);
      trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+mtop.str()+"]/F").c_str());
    }
    trees["noSyst"]->Branch((nameshortv+"_size").c_str(), &sizes[nameshortv]);
  }

  string nameshortv= "Event";
  vector<string> extravars = additionalVariables(nameshortv);
  for(size_t addv = 0; addv < extravars.size();++addv){
    string name = nameshortv+"_"+extravars.at(addv);

    if (name.find("EventNumber")!=std::string::npos){
      std::cout<<"=====================sto riempendo il branch event number"<<std::endl;
      trees["noSyst"]->Branch(name.c_str(), &double_values[name],(name+"/D").c_str());
    }
    else trees["noSyst"]->Branch(name.c_str(), &float_values[name],(name+"/F").c_str());
  }
  
  //Prepare the trees cloning all branches and setting the correct names/titles:
  if(!addNominal){
    DMTrees = fs->mkdir( "systematics_trees" );
  }
  
  trees["EventHistory"] =  new TTree("EventHistory","EventHistory");
  trees["WeightHistory"] =  new TTree("WeightHistory","WeightHistory");
  trees["EventHistory"]->Branch("initialEvents",&nInitEvents);

  for(size_t s=0;s< systematics.size();++s){
    std::string syst  = systematics.at(s);
    if(syst=="noSyst")continue;
    trees[syst]= trees["noSyst"]->CloneTree();
    trees[syst]->SetName((channel+"__"+syst).c_str());
    trees[syst]->SetTitle((channel+"__"+syst).c_str());
  }

  initTreeWeightHistory(useLHEWeights);

  string L1Name ="Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt"; 
  string L1RCName = "Summer16_23Sep2016V3_MC_L1RC_AK4PFchs.txt"; 
  string L2Name = "Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt";
  string L3Name = "Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt";
  string L2L3ResName = "Summer16_23Sep2016V3_MC_L2L3Residual_AK4PFchs.txt";

  if(isData){
    L1Name   = ("Summer16_23Sep2016"+EraLabel+"V3_DATA_L1FastJet_AK4PFchs.txt").c_str() ;
    L1RCName = ("Summer16_23Sep2016"+EraLabel+"V3_DATA_L1RC_AK4PFchs.txt").c_str() ;  
    L2Name   = ("Summer16_23Sep2016"+EraLabel+"V3_DATA_L2Relative_AK4PFchs.txt").c_str() ;
    L3Name   = ("Summer16_23Sep2016"+EraLabel+"V3_DATA_L3Absolute_AK4PFchs.txt").c_str() ;
    L2L3ResName = ("Summer16_23Sep2016"+EraLabel+"V3_DATA_L2L3Residual_AK4PFchs.txt").c_str() ;
  }

  if(isV2){
    L1Name = "Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt";
    L2Name = "Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt";
    L3Name = "Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt";
    L2L3ResName = "Summer16_23Sep2016V3_MC_L2L3Residual_AK4PFchs.txt"; 
  }

  jecParsL1  = new JetCorrectorParameters(L1Name);
  jecParsL2  = new JetCorrectorParameters(L2Name);
  jecParsL3  = new JetCorrectorParameters(L3Name);
  jecParsL2L3Residuals  = new JetCorrectorParameters(L2L3ResName);
  jecPars.push_back(*jecParsL1);
  jecPars.push_back(*jecParsL2);
  jecPars.push_back(*jecParsL3);
  if(isData)jecPars.push_back(*jecParsL2L3Residuals);

  jecParsL1_vect.push_back(*jecParsL1);
  
  jecCorr = new FactorizedJetCorrector(jecPars);
  jecCorr_L1 = new FactorizedJetCorrector(jecParsL1_vect);
  jecUnc  = new JetCorrectionUncertainty(*(new JetCorrectorParameters(("Summer16_23Sep2016"+EraLabel+"V3_DATA_UncertaintySources_AK4PFchs.txt").c_str() , "Total")));

  //corrections on AK8

  string L1Name8 = "Summer16_23Sep2016V3_MC_L1FastJet_AK8PFchs.txt"; 
  string L1RCName8 = "Summer16_23Sep2016V3_MC_L1RC_AK8PFchs.txt"; 
  string L2Name8 = "Summer16_23Sep2016V3_MC_L2Relative_AK8PFchs.txt";
  string L3Name8 = "Summer16_23Sep2016V3_MC_L3Absolute_AK8PFchs.txt";
  string L2L3ResName8 = "Summer16_23Sep2016V3_MC_L2L3Residual_AK8PFchs.txt";

  if(isData){
    L1Name8   = ("Summer16_23Sep2016"+EraLabel+"V3_DATA_L1FastJet_AK8PFchs.txt").c_str() ;
    L1RCName8 = ("Summer16_23Sep2016"+EraLabel+"V3_DATA_L1RC_AK8PFchs.txt").c_str() ;  
    L2Name8   = ("Summer16_23Sep2016"+EraLabel+"V3_DATA_L2Relative_AK8PFchs.txt").c_str() ;
    L3Name8   = ("Summer16_23Sep2016"+EraLabel+"V3_DATA_L3Absolute_AK8PFchs.txt").c_str() ;
    L2L3ResName8 = ("Summer16_23Sep2016"+EraLabel+"V3_DATA_L2L3Residual_AK8PFchs.txt").c_str() ;
  }

  if(isV2){
    L1Name8 = "Summer16_23Sep2016V3_MC_L1FastJet_AK8PFchs.txt";
    L2Name8 = "Summer16_23Sep2016V3_MC_L2Relative_AK8PFchs.txt";
    L3Name8 = "Summer16_23Sep2016V3_MC_L3Absolute_AK8PFchs.txt";
    L2L3ResName8 = "Summer16_23Sep2016V3_MC_L2L3Residual_AK8PFchs.txt"; 
  }

  jecParsL18  = new JetCorrectorParameters(L1Name8);
  jecParsL28  = new JetCorrectorParameters(L2Name8);
  jecParsL38  = new JetCorrectorParameters(L3Name8);
  jecParsL2L3Residuals8  = new JetCorrectorParameters(L2L3ResName8);
  jecPars8.push_back(*jecParsL18);
  jecPars8.push_back(*jecParsL28);
  jecPars8.push_back(*jecParsL38);
  if(isData)jecPars8.push_back(*jecParsL2L3Residuals8);

  jecParsL1_vect8.push_back(*jecParsL18);
  jecParsNoL1_vect8.push_back(*jecParsL28);
  jecParsNoL1_vect8.push_back(*jecParsL38);
  if(isData)jecParsNoL1_vect8.push_back(*jecParsL2L3Residuals8);
    
  jecCorr8 = new FactorizedJetCorrector(jecPars8);
  jecCorr_L18 = new FactorizedJetCorrector(jecParsL1_vect8);
  jecCorr_NoL18 = new FactorizedJetCorrector(jecParsNoL1_vect8);
  jecUnc8  = new JetCorrectionUncertainty(*(new JetCorrectorParameters(("Summer16_23Sep2016"+EraLabel+"V3_DATA_UncertaintySources_AK8PFchs.txt").c_str() , "Total")));

  isFirstEvent = true;
  doBTagSF= true;
  if(isData)doPU= false;
  
  season = "Summer11";
  distr = "pileUpDistr" + season + ".root";
}


void DMAnalysisTreeMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
    iRun.getByLabel(edm::InputTag("TriggerUserData","triggerNameTree"), triggerNamesR);
    iRun.getByLabel(metNames_, metNames);
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      //cout << "trigger test tname "<< tname <<endl; 
    }
}


void DMAnalysisTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::vector<edm::ParameterSet >::const_iterator itPsets = physObjects.begin();

  nInitEventsHisto->Fill(0.1);
  nInitEvents+=1;

  // event info
  iEvent.getByToken(t_lumiBlock_,lumiBlock );
  iEvent.getByToken(t_runNumber_,runNumber );
  iEvent.getByToken(t_eventNumber_,eventNumber );

  if(useLHE){
    iEvent.getByToken(t_lhes_, lhes);
  }
  if(addLHAPDFWeights){
    iEvent.getByToken(t_genprod_, genprod);
  }

  iEvent.getByToken(t_jetKeys_, jetKeys);
  iEvent.getByToken(t_muKeys_, muKeys);
  
  vector<TLorentzVector> genlep,gentop,genantitop,genw,gennu,genz,gena;
  vector<TLorentzVector> pare,parebar,parmu,parmubar,parnumu,parnumubar,partau,partaubar,parnue,parnuebar,parnutau,parnutaubar,parz,parw;

  //Parton-level info                                                                                                                               
  if(isData==true){
    getPartonW=false;
    getParticleWZ=false;
    getPartonTop=false;
    doWReweighting=false;
    doTopReweighting=false;
    useLHE=false;
    useLHEWeights=false;
  }
  
  if(getPartonW || getPartonTop || doWReweighting || doTopReweighting){
    if(!useLHE)return;
    genlep.clear();
    gentop.clear();
    genantitop.clear();
    genw.clear();
    genz.clear();
    gena.clear();
    gennu.clear();
    pare.clear();parebar.clear();parmu.clear();parmubar.clear();parnumu.clear();parnumubar.clear();parnue.clear();parnuebar.clear();parnutau.clear();parnutaubar.clear();parz.clear();parw.clear();

    float_values["Event_Z_EW_Weight"]= 1.0;
    float_values["Event_W_EW_Weight"]= 1.0;
    float_values["Event_Z_QCD_Weight"]= 1.0;
    float_values["Event_W_QCD_Weight"]= 1.0;
    float_values["Event_Z_Weight"]= 1.0;
    float_values["Event_W_Weight"]= 1.0;
    float_values["Event_T_Weight"]= 1.0;
    float_values["Event_T_Ext_Weight"]= 1.0;
    size_t nup=lhes->hepeup().NUP;

    for( size_t i=0;i<nup;++i){
      //      cout << " particle number " << i << endl;
      int id = lhes->hepeup().IDUP[i];
      float px = lhes->hepeup().PUP[i][0];
      float py = lhes->hepeup().PUP[i][1];
      float pz = lhes->hepeup().PUP[i][2];
      float energy = lhes->hepeup().PUP[i][3];
    
      TLorentzVector vec;
      math::XYZTLorentzVector part = math::XYZTLorentzVector(px, py, pz, energy);
      float pt = part.pt();
      float phi = part.phi();
      float eta = part.eta();
      
      if(pt>0){
	vec.SetPtEtaPhiE(pt, eta, phi, energy);
	
	if(abs (id) == 11 || abs (id) == 13 || abs(id) == 15){
	  genlep.push_back(vec);
	}
	if(abs (id) == 12 || abs (id) == 14 || abs(id) == 16){
	  gennu.push_back(vec);
	}
	if(id == 6 ){
	  gentop.push_back(vec);
	}
	if(id == -6 ){
	  genantitop.push_back(vec);
	}
	if(abs (id) == 24 ){
	  genw.push_back(vec);
	}
	if(abs (id) == 23 ){
	  genz.push_back(vec);
	}
	if(abs (id) == 22 ){
	  gena.push_back(vec);
	}
      }      
    }
    
    if(getPartonTop && gentop.size()==1){
      float_values["Event_T_Pt"]= gentop.at(0).Pt();
      float_values["Event_T_Eta"]= gentop.at(0).Eta();
      float_values["Event_T_Phi"]= gentop.at(0).Phi();
      float_values["Event_T_E"]= gentop.at(0).Energy();
      float_values["Event_T_Mass"]= gentop.at(0).M();
    }
    if(getPartonTop && genantitop.size()==1){
      float_values["Event_Tbar_Pt"]= genantitop.at(0).Pt();
      float_values["Event_Tbar_Eta"]= genantitop.at(0).Eta();
      float_values["Event_Tbar_Phi"]= genantitop.at(0).Phi();
      float_values["Event_Tbar_E"]= genantitop.at(0).Energy();
      float_values["Event_Tbar_Mass"]= genantitop.at(0).M();
      
    }
    if((getPartonW || doWReweighting )) {
      if(genw.size()==1){
	float_values["Event_W_Pt"]= genw.at(0).Pt();
	float_values["Event_W_Eta"]= genw.at(0).Eta();
	float_values["Event_W_Phi"]= genw.at(0).Phi();
	float_values["Event_W_E"]= genw.at(0).Energy();
	float_values["Event_W_Mass"]= genw.at(0).M();	
	
	double ptW = genw.at(0).Pt();
	double wweight = getWPtWeight(ptW);			
	float_values["Event_W_QCD_Weight"]= wweight;
      }
      else (float_values["Event_W_QCD_Weight"]=1.0);
    }
    if((getPartonW || doWReweighting )){ 
      if(genz.size()==1){
	float_values["Event_Z_Pt"]= genz.at(0).Pt();
	float_values["Event_Z_Eta"]= genz.at(0).Eta();
	float_values["Event_Z_Phi"]= genz.at(0).Phi();
	float_values["Event_Z_E"]= genz.at(0).Energy();
	float_values["Event_Z_Mass"]= genz.at(0).M();	
	
	double ptW = genz.at(0).Pt();
	double wweight = getZPtWeight(ptW);			
	cout << wweight << endl;
	float_values["Event_Z_QCD_Weight"]= wweight;
      }
      else (float_values["Event_Z_QCD_Weight"]=1.0);
    }
    //cout << "getPartonW is \t" << getPartonW << "and doWReweighting is\t" << doWReweighting << endl;
    if((getPartonW || doWReweighting ) ) {
      if(gena.size()==1){       
      float_values["Event_a_Pt"]= gena.at(0).Pt();
      float_values["Event_a_Eta"]= gena.at(0).Eta();
      float_values["Event_a_Phi"]= gena.at(0).Phi();
      float_values["Event_a_E"]= gena.at(0).Energy();
      float_values["Event_a_Mass"]= gena.at(0).M();	
      
      double ptW = gena.at(0).Pt();
      double wweight = getAPtWeight(ptW);			
      //cout << wweight << endl;
      float_values["Event_a_Weight"]= wweight;
      }
      else (float_values["Event_a_Weight"]=1.0);
    }
    if( (getPartonTop || doTopReweighting)) {
      if (gentop.size()==1 && genantitop.size()==1 && getPartonTop){
	double ptT = gentop.at(0).Pt();
	double ptTbar = genantitop.at(0).Pt();
	double tweight = getTopPtWeight(ptT,ptTbar);			
	double tweightext = getTopPtWeight(ptT,ptTbar,true);			
	float_values["Event_T_Weight"]= tweight;
	float_values["Event_T_Ext_Weight"]= tweightext;
	//cout << tweight << endl;
      }
      else {(float_values["Event_T_Weight"]=1.0);
	(float_values["Event_T_Ext_Weight"]=1.0);}
    }
    if(useLHEWeights){
      getEventLHEWeights();
    }

    //if(!isData) {//GEN PARTICLES HAVE TO BE FIXED
    if(0>1){ 
    iEvent.getByToken(t_genParticleCollection_ , genParticles);
    
      //const bool print = false;
      int momId=-999;
      int momStatus=-999;
  
      std::vector< int > MomId;
      std::vector< int > MomStatus;
      for(size_t i=0; i<genParticles->size(); ++i) {
	  const reco::GenParticle& p = (*genParticles)[i];
	  momId = p.numberOfMothers() ? p.mother()->pdgId() : 0;                                                                                                                            
	  momStatus = p.numberOfMothers() ? p.mother()->status() : 0;    
	  MomId.push_back(momId);
	  MomStatus.push_back(momStatus);
      }

      for(size_t i=0;/*i<genParticles->size()*/(size_t)max_instances[gen_label]/*max((int)genParticles->size(),(int)30)*/;++i){
	string nameshortv= "genPart";
	string pref = obj_to_pref[gen_label];

	const reco::GenParticle& p = (*genParticles)[i];

	/*vfloats_values[makeName(gen_label,pref,"Pt")][i]=(float)((p.p4()).Pt());
	vfloats_values[makeName(gen_label,pref,"Eta")][i]=(float)((p.p4()).Eta());
	vfloats_values[makeName(gen_label,pref,"Phi")][i]=(float)((p.p4()).Phi());
	vfloats_values[makeName(gen_label,pref,"E")][i]=(float)((p.p4()).E());
	vfloats_values[makeName(gen_label,pref,"Status")][i]=(int)(p.status());
	vfloats_values[makeName(gen_label,pref,"Id")][i]=(int)(p.pdgId());
	*/

	int momId=-999;
	//int momStatus =-999;
	/*momId = p.numberOfMothers() ? p.mother()->pdgId() : 0;
	momStatus = p.numberOfMothers() ? p.mother()->status() : 0;

	vfloats_values[makeName(gen_label,pref,"Mom0Id")][i]=(float)(momId);
	vfloats_values[makeName(gen_label,pref,"Mom0Status")][i]=(float)(momStatus);
	*/
	/*for(size_t j = 0, n=p.numberOfDaughters(); j<n; ++j) {
	  vfloats_values[makeName(gen_label,pref,"dauId1")][i]=(float)(p.daughter(j)->pdgId());
	  vfloats_values[makeName(gen_label,pref,"dauStatus1")][i]=(float)(p.daughter(j)->status());
	  vfloats_values[makeName(gen_label,pref,"dauId2")][i]=(float)(p.daughter(j+1)->pdgId());
	  vfloats_values[makeName(gen_label,pref,"dauStatus2")][i]=(float)(p.daughter(j+1)->status());
	  vfloats_values[makeName(gen_label,pref,"dauId3")][i]=(float)(p.daughter(j+1)->pdgId());
	  vfloats_values[makeName(gen_label,pref,"dauStatus3")][i]=(float)(p.daughter(j+1)->status());
	  }*/

	if(p.pdgId()==momId && isEWKID(p.pdgId()) && getParticleWZ){
	  //cout << "inside gen loop" << endl;
	  TLorentzVector vec;
	  vec.SetPtEtaPhiE((p.p4()).Pt(),(p.p4()).Eta(),(p.p4()).Phi(), (p.p4()).E());
	  if(p.pdgId()==11  && p.status()==1){
	    pare.push_back(vec);	  }
	  if(p.pdgId()==-11  && p.status()==1){
	    parebar.push_back(vec);	  }
	  if(p.pdgId()==13 && p.status()==1){
	    parmu.push_back(vec);	  }
	  if(p.pdgId()==-13 && p.status()==1){
	    parmubar.push_back(vec);	  }

	  if(p.pdgId()==15  && p.status()==2){
	    partau.push_back(vec);  }
	  if(p.pdgId()==-15  && p.status()==2){
	    partaubar.push_back(vec);  }
	  
	  if(p.pdgId()==12  && p.status()==1){
	    parnue.push_back(vec);	  }
	  if(p.pdgId()==14  && p.status()==1){
	    parnumu.push_back(vec);	  }
	  if(p.pdgId()==16  && p.status()==1){
	    parnutau.push_back(vec);	  }

	  if(p.pdgId()==-12  && p.status()==1){
	    parnuebar.push_back(vec);	  }
	  if(p.pdgId()==-14  && p.status()==1){
	    parnumubar.push_back(vec);	  }
	  if(p.pdgId()==-16  && p.status()==1){
	    parnutaubar.push_back(vec);	  }
	  
	  if(abs(p.pdgId())==23 && p.status()==62){
	    double wweight = getZEWKPtWeight((p.p4()).Pt());
	    float_values["Event_Z_EW_Weight"]= wweight;

	  }
	  if(abs(p.pdgId())==24 && p.status()==62){
	    double wweight = getWEWKPtWeight((p.p4()).Pt());
	    cout << "third if ok"<< endl;
	    float_values["Event_W_EW_Weight"]= wweight;
	     cout << "z weight:\t" << wweight << endl;
	     }   
	 
	}

      }
  
      //Z
      if(parmu.size()>0&& parmubar.size()>0){parz.push_back(parmu.at(0)+parmubar.at(0)) ;}
      if(pare.size()>0&& parebar.size()>0){parz.push_back(pare.at(0)+parebar.at(0)) ;}
      if(partau.size()>0&& partaubar.size()>0){parz.push_back(partau.at(0)+partaubar.at(0)) ;}
      if(parnumu.size()>0&& parnumubar.size()>0){parz.push_back(parnumu.at(0)+parnumubar.at(0)) ;}
      if(parnue.size()>0&& parnuebar.size()>0){parz.push_back(parnue.at(0)+parnuebar.at(0)) ;}
      if(parnutau.size()>0&& parnutaubar.size()>0){parz.push_back(parnutau.at(0)+parnutaubar.at(0)) ;}
      if(   float_values["Event_Z_EW_Weight"] ==1 &&parz.size()>0 )    float_values["Event_Z_EW_Weight"]= getZEWKPtWeight(parz.at(0).Pt());
      if(   float_values["Event_Z_QCD_Weight"] ==1 &&parz.size()>0 )    float_values["Event_Z_QCD_Weight"]= getWPtWeight(parz.at(0).Pt());

      //W
      if(parmu.size()>0&& parnumubar.size()>0){parw.push_back(parmu.at(0)+parnumubar.at(0)) ;}

      if(pare.size()>0&& parnuebar.size()>0){parw.push_back(pare.at(0)+parnuebar.at(0)) ;}

      if(partau.size()>0&& parnutaubar.size()>0){parw.push_back(partau.at(0)+parnutaubar.at(0)) ;}
      
      if(parnumu.size()>0&& parmubar.size()>0){parw.push_back(parnumu.at(0)+parmubar.at(0)) ;}

      if(parnue.size()>0&& parebar.size()>0){parw.push_back(parnue.at(0)+parebar.at(0)) ;}


      if(parnutau.size()>0&& partaubar.size()>0){parw.push_back(parnutau.at(0)+partaubar.at(0)) ;}

      if(   float_values["Event_W_EW_Weight"] ==1 &&parw.size()>0 )    {
	float_values["Event_W_EW_Weight"]= getWEWKPtWeight(parw.at(0).Pt());
      }
      if(   float_values["Event_W_QCD_Weight"] ==1 &&parw.size()>0 )    float_values["Event_W_QCD_Weight"]= getWPtWeight(parw.at(0).Pt());
    }
    float_values["Event_W_Weight"]= float_values["Event_W_EW_Weight"]*float_values["Event_W_QCD_Weight"];
    float_values["Event_Z_Weight"]= float_values["Event_Z_EW_Weight"]*float_values["Event_Z_QCD_Weight"];   
    //cout << "ending mega if" << endl;
  }
  trees["WeightHistory"]->Fill();

  //Part 3: filling the additional variables

  //boosted part
  iEvent.getByToken(jetAK8vSubjetIndex0, ak8jetvSubjetIndex0);
  iEvent.getByToken(jetAK8vSubjetIndex1, ak8jetvSubjetIndex1);
 
  //Part 0: trigger preselection
 if(useTriggers){
    //Trigger names are retrieved from the run tree
    iEvent.getByToken(t_triggerBits_,triggerBits );
    iEvent.getByToken(t_triggerPrescales_,triggerPrescales );
    bool triggerOr = getEventTriggers();
    
    if(isFirstEvent){
      for(size_t bt = 0; bt < triggerNamesR->size();++bt){
	std::string tname = triggerNamesR->at(bt);
	cout << "trigger test tname "<< tname << " passes "<< triggerBits->at(bt)<< endl;
      }
    }
    
    if(cutOnTriggers && !triggerOr) return;
  }

  if(useMETFilters){
   iEvent.getByToken(t_metBits_,metBits );
   iEvent.getByToken(t_BadChargedCandidateFilter_, BadChargedCandidateFilter);
   iEvent.getByToken(t_BadPFMuonFilter_, BadPFMuonFilter);
   if(isFirstEvent){
      for(size_t bt = 0; bt < metNames->size();++bt){
	std::string tname = metNames->at(bt);
      }
    }
    getMETFilters();
  }
  
  if(changeJECs || recalculateEA){
    iEvent.getByToken(t_Rho_ ,rho);
    Rho = *rho; 
  }
  float_values["Event_Rho"] = (double)Rho;
  if(isFirstEvent){
    isFirstEvent = false;
  }
    
  if( addPV || changeJECs){
    iEvent.getByToken(t_pvZ_,pvZ);
    iEvent.getByToken(t_pvChi2_,pvChi2);
    iEvent.getByToken(t_pvNdof_,pvNdof);
    iEvent.getByToken(t_pvRho_,pvRho);
    nPV = pvZ->size();
  }
  
  //Part 1 taking the obs values from the edm file
  for (;itPsets!=physObjects.end();++itPsets){ 
    variablesFloat = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesF"); 
    variablesInt = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesI"); 
    singleFloat = itPsets->template getParameter<std::vector<edm::InputTag> >("singleF"); 
    singleInt = itPsets->template getParameter<std::vector<edm::InputTag> >("singleI"); 
    std::vector<edm::InputTag >::const_iterator itF = variablesFloat.begin();
    std::vector<edm::InputTag >::const_iterator itI = variablesInt.begin();
    std::vector<edm::InputTag >::const_iterator itsF = singleFloat.begin();
    std::vector<edm::InputTag >::const_iterator itsI = singleInt.begin();
    string namelabel = itPsets->getParameter< string >("label");
    string nameprefix = itPsets->getParameter< string >("prefix");
    size_t maxInstance=(size_t)max_instances[namelabel];


    variablesDouble = itPsets->template getParameter<std::vector<edm::InputTag> >("variablesD"); 
    singleDouble = itPsets->template getParameter<std::vector<edm::InputTag> >("singleD"); 

    std::vector<edm::InputTag >::const_iterator itD = variablesDouble.begin();
    std::vector<edm::InputTag >::const_iterator itsD = singleDouble.begin();
    
  
    //Vectors of floats
    for (;itF != variablesFloat.end();++itF){
      string varname=itF->instance();
      
      string name = makeBranchName(namelabel,nameprefix,varname);
      
      //string namelabel;
      float tmp =1.0;
      iEvent.getByToken(t_floats[name] ,h_floats[name]);
      //iEvent.getByLabel(*(itF),h_floats[name]);
      //      cout << "name "<< name <<endl;
      for (size_t fi = 0;fi < maxInstance ;++fi){
	if(fi <h_floats[name]->size()){tmp = h_floats[name]->at(fi);}
	else { tmp = -9999.; }
	//	cout << " setting name "<< name<< " at instance "<< fi <<" to value "<< tmp <<endl;
	vfloats_values[name][fi]=tmp;
	for (size_t sc=0;sc< obj_cats[namelabel].size();++sc){
	  string category = obj_cats[namelabel].at(sc);
	  string namecat = makeBranchNameCat(namelabel,category,nameprefix,varname);
	  vfloats_values[namecat][fi]=-9999.;
	}
      }
      sizes[namelabel]=h_floats[name]->size();
    }

    for (;itD != variablesDouble.end();++itD){
      string varname=itD->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
      //string namelabel;
      float tmp =1.0;
      iEvent.getByToken(t_doubles[name] ,h_doubles[name]);
      for (size_t fi = 0;fi < maxInstance ;++fi){
	if(fi <h_doubles[name]->size()){tmp = h_doubles[name]->at(fi);}
	else { tmp = -9999.; }
	vdoubles_values[name][fi]=tmp;
      }
      sizes[namelabel]=h_doubles[name]->size();
    }
    
    //Vectors of ints
    for (;itI != variablesInt.end();++itI){
      string varname=itI->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
      int tmp = 1;
      iEvent.getByToken(t_ints[name] ,h_ints[name]);
      //iEvent.getByLabel(*(itI),h_ints[name]);
      for (size_t fi = 0;fi < maxInstance;++fi){
	if(fi <h_ints[name]->size()){tmp = h_ints[name]->at(fi);}
	else { tmp = -9999.; }
	vints_values[name][fi]=tmp;
      }
    }  
    
    //    std::cout << " checkpoint ints"<<endl;
    //Single floats/ints
    for (;itsF != singleFloat.end();++itsF){
      string varname=itsF->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
      iEvent.getByToken(t_float[name],h_float[name]);
      //iEvent.getByLabel(*(itsF),h_float[name]);
      float_values[name]=*h_float[name];
    }

    for (;itsD != singleDouble.end();++itsD){
      string varname=itsD->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
      iEvent.getByToken(t_double[name] ,h_double[name]);
      //iEvent.getByLabel(*(itsD),h_double[name]);
      double_values[name]=*h_double[name];
    }
    for (;itsI != singleInt.end();++itsI){
      string varname=itsI->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
      iEvent.getByToken(t_int[name],h_int[name]);
      //iEvent.getByLabel(*(itsI),h_int[name]);
      int_values[name]=*h_int[name];
    }
    //    std::cout << " checkpoint singles"<<endl;
  }

  //  std::cout << " checkpoint part 1"<<endl;


  //Part 2: selection and analysis-level changes
  //This might change for each particular systematics, 
  //e.g. for each jet energy scale variation, for MET or even lepton energy scale variations


  vector<TLorentzVector> photons;
  vector<TLorentzVector> electrons;
  vector<TLorentzVector> muons;
  vector<TLorentzVector> leptons;
  vector<TLorentzVector> loosemuons;
  vector<TLorentzVector> jets;
  vector<TLorentzVector> jetsnob;
  vector<TLorentzVector> bjets;
  vector<TLorentzVector> subjets;
  vector<TLorentzVector> topjets;
  vector<TLorentzVector> type1topjets;
  vector<TLorentzVector> type2topjets;
  vector<TLorentzVector> resolvedtops;

  vector<float> leptonsCharge;

  vector<int> flavors;

  for (size_t s = 0; s< systematics.size();++s){
    
    int nb=0,nc=0,nudsg=0;

    int ncsvl_tags=0,ncsvt_tags=0,ncsvm_tags=0;
    int ncsvl_subj_tags=0,ncsvm_subj_tags=0;
    getEventTriggers();

    photons.clear();
    leptons.clear();
    electrons.clear();
    muons.clear();
    loosemuons.clear();
    jets.clear();
    jetsnob.clear();
    bjets.clear();
    type2topjets.clear();
    type1topjets.clear();
    topjets.clear();
    resolvedtops.clear();
    string syst = systematics.at(s);
    nTightJets=0;

    jsfscsvt.clear();
    jsfscsvt_b_tag_up.clear(); 
    jsfscsvt_b_tag_down.clear(); 
    jsfscsvt_mistag_up.clear(); 
    jsfscsvt_mistag_down.clear();

    jsfscsvm.clear(); 
    jsfscsvm_b_tag_up.clear(); 
    jsfscsvm_b_tag_down.clear(); 
    jsfscsvm_mistag_up.clear(); 
    jsfscsvm_mistag_down.clear();
    
    jsfscsvl.clear(); 
    jsfscsvl_b_tag_up.clear(); 
    jsfscsvl_b_tag_down.clear(); 
    jsfscsvl_mistag_up.clear();
    jsfscsvl_mistag_down.clear();

    jsfscsvm_subj.clear(); 
    jsfscsvm_subj_b_tag_up.clear(); 
    jsfscsvm_subj_b_tag_down.clear(); 
    jsfscsvm_subj_mistag_up.clear(); 
    jsfscsvm_subj_mistag_down.clear();
    
    jsfscsvl_subj.clear(); 
    jsfscsvl_subj_b_tag_up.clear(); 
    jsfscsvl_subj_b_tag_down.clear(); 
    jsfscsvl_subj_mistag_up.clear();
    jsfscsvl_subj_mistag_down.clear();


    //---------------- Soureek Adding PU Info ------------------------------
    //if(doPU_){
    //  iEvent.getByToken(t_ntrpu_,ntrpu);
    //  nTruePU=*ntrpu;
    //  getPUSF();
    //}
    
    int mapBJets[20], mapEle[20], mapMu[20];
    for(int i = 0; i<20;++i){
      mapBJets[i]=-1; mapEle[i]=-1; mapMu[i]=-1;
    } 
    int lepidx=0;
    int bjetidx=0;

    //Photons
    for(int ph = 0;ph < max_instances[photon_label] ;++ph){
      string pref = obj_to_pref[photon_label];
      
      
      float pt = vfloats_values[makeName(photon_label,pref,"Pt")][ph];
      float eta = vfloats_values[makeName(photon_label,pref,"Eta")][ph];      
      
      float sieie = vfloats_values[makeName(photon_label,pref,"SigmaIEtaIEta")][ph];      
      float hoe = vfloats_values[makeName(photon_label,pref,"HoverE")][ph];      
      
      float abseta = fabs(eta);

      float pho_isoC  = vfloats_values[makeName(photon_label,pref,"ChargedHadronIso")][ph];      
      float pho_isoP  = vfloats_values[makeName(photon_label,pref,"NeutralHadronIso")][ph];      
      float pho_isoN     =  vfloats_values[makeName(photon_label,pref,"PhotonIso")][ph];      

      float pho_isoCea  = vfloats_values[makeName(photon_label,pref,"ChargedHadronIsoEAcorrected")][ph];      
      float pho_isoPea  = vfloats_values[makeName(photon_label,pref,"PhotonIsoEAcorrected")][ph];      
      float pho_isoNea     =  vfloats_values[makeName(photon_label,pref,"NeutralHadronIsoEAcorrected")][ph];      

      if(recalculateEA){
	pho_isoCea     = std::max( double(0.0) ,(pho_isoC - Rho*getEffectiveArea("ch_hadrons",abseta)));
	pho_isoPea     = std::max( double(0.0) ,(pho_isoP - Rho*getEffectiveArea("photons",abseta)));
	pho_isoNea     = std::max( double(0.0) ,(pho_isoN - Rho*getEffectiveArea("neu_hadrons",abseta)));
      }
      
      bool isBarrel = (abseta<1.479);
      bool isEndcap = (abseta>1.479 && abseta < 2.5);

      vfloats_values[photon_label+"_isLooseSpring15"][ph]=0.0;
      vfloats_values[photon_label+"_isMediumSpring15"][ph]=0.0;
      vfloats_values[photon_label+"_isTightSpring15"][ph]=0.0;
 
      if(isBarrel){

	if( sieie < 0.0103 &&   hoe < 0.05 &&   pho_isoCea < 2.44 &&   pho_isoNea < (2.57+exp(0.0044*pt +0.5809) ) &&   pho_isoPea < (1.92+0.0043*pt ) )vfloats_values[photon_label+"_isLooseSpring15"][ph]=1.0;
	if( sieie < 0.01 &&   hoe < 0.05 &&   pho_isoCea < 1.31 &&   pho_isoNea < (0.60+exp(0.0044*pt +0.5809) ) &&   pho_isoPea < (1.33+0.0043*pt ) )vfloats_values[photon_label+"_isMediumSpring15"][ph]=1.0;
	if( sieie < 0.01 &&   hoe < 0.05 &&   pho_isoCea < 0.91 &&   pho_isoNea < (0.33+exp(0.0044*pt +0.5809) ) &&   pho_isoPea < (0.61+0.0043*pt ) )vfloats_values[photon_label+"_isTightSpring15"][ph]=1.0;
      }
      if(isEndcap){
	if( sieie < 0.0277 &&   hoe < 0.05 &&   pho_isoCea < 1.84 &&   pho_isoNea < (4.00+exp(0.0040*pt +0.9402) ) &&   pho_isoPea < (1.92+0.0043*pt ) )vfloats_values[photon_label+"_isLooseSpring15"][ph]=1.0;
	if( sieie < 0.0267 &&   hoe < 0.05 &&   pho_isoCea < 1.25 &&   pho_isoNea < (1.65+exp(0.0040*pt +0.9402) ) &&   pho_isoPea < (1.33+0.0043*pt ) )vfloats_values[photon_label+"_isMediumSpring15"][ph]=1.0;
	if( sieie < 0.0267 &&   hoe < 0.05 &&   pho_isoCea < 0.65 &&   pho_isoNea < (0.93+exp(0.0040*pt +0.9402) ) &&   pho_isoPea < (0.61+0.0043*pt ) )vfloats_values[photon_label+"_isTightSpring15"][ph]=1.0;
      }
    
    }

    //Muons
    for(int mu = 0;mu < max_instances[mu_label] ;++mu){
      string pref = obj_to_pref[mu_label];
      float isTight = vfloats_values[makeName(mu_label,pref,"IsTightMuon")][mu];
      float isLoose = vfloats_values[makeName(mu_label,pref,"IsLooseMuon")][mu];
      float isMedium = vfloats_values[makeName(mu_label,pref,"IsMediumMuon")][mu];
      float isSoft = vfloats_values[makeName(mu_label,pref,"IsSoftMuon")][mu];

      float pt = vfloats_values[makeName(mu_label,pref,"Pt")][mu];
      float eta = vfloats_values[makeName(mu_label,pref,"Eta")][mu];
      float phi = vfloats_values[makeName(mu_label,pref,"Phi")][mu];
      float energy = vfloats_values[makeName(mu_label,pref,"E")][mu];
      float iso = vfloats_values[makeName(mu_label,pref,"Iso04")][mu];
      
      float muCharge = vfloats_values[makeName(mu_label,pref,"Charge")][mu];

      if(isMedium>0 && pt> 30 && abs(eta) < 2.1 && iso <0.25){ 
	++float_values["Event_nMediumMuons"];
	TLorentzVector muon;
	muon.SetPtEtaPhiE(pt, eta, phi, energy);
	muons.push_back(muon);
	leptons.push_back(muon);
	flavors.push_back(13);
	
	leptonsCharge.push_back(muCharge);
	
	mapMu[lepidx]=mu; 
	if(isInVector(obj_cats[mu_label],"Medium")){
	  fillCategory(mu_label,"Medium",mu,float_values["Event_nMediumMuons"]-1);
	}
	++lepidx;
      }
      
      if(isInVector(obj_cats[mu_label],"Medium")){
	sizes[mu_label+"Medium"]=(int)float_values["Event_nMediumMuons"];
      }
      
      if(isLoose>0 && pt> 10 && abs(eta) < 2.4 && iso<0.25){
	if(isInVector(obj_cats[mu_label],"Loose")){
	  ++float_values["Event_nLooseMuons"];
	  if(isInVector(obj_cats[mu_label],"Loose")){
	    fillCategory(mu_label,"Loose",mu,float_values["Event_nLooseMuons"]-1);
	  }
	}
      }
      if(isInVector(obj_cats[mu_label],"Loose")){
	sizes[mu_label+"Loose"]=(int)float_values["Event_nLooseMuons"];
      }


      if(isTight>0 && pt> 10 && abs(eta) < 2.4 && iso<0.25){
        if(isInVector(obj_cats[mu_label],"Tight")){
          ++float_values["Event_nTightMuons"];
          if(isInVector(obj_cats[mu_label],"Tight")){
            fillCategory(mu_label,"Tight",mu,float_values["Event_nTightMuons"]-1);
          }
        }
      }
      if(isInVector(obj_cats[mu_label],"Tight")){
        sizes[mu_label+"Tight"]=(int)float_values["Event_nTightMuons"];
      }
      
      if(isSoft>0 && pt> 30 && abs(eta) < 2.4){
	++float_values["Event_nSoftMuons"]; 
      }
      if(isLoose>0 && pt > 15){
	TLorentzVector muon;
	muon.SetPtEtaPhiE(pt, eta, phi, energy);
	loosemuons.push_back(muon);
      }
    }

    //Electrons:
    for(int el = 0;el < max_instances[ele_label] ;++el){
      string pref = obj_to_pref[ele_label];
      float pt = vfloats_values[makeName(ele_label,pref,"Pt")][el];
      float isTight = vfloats_values[makeName(ele_label,pref,"isTight")][el];
      float isLoose = vfloats_values[makeName(ele_label,pref,"isLoose")][el];
      float isMedium = vfloats_values[makeName(ele_label,pref,"isMedium")][el];
      float isVeto = vfloats_values[makeName(ele_label,pref,"isVeto")][el];

      isTight = vfloats_values[makeName(ele_label,pref,"vidTight")][el];
      isLoose = vfloats_values[makeName(ele_label,pref,"vidLoose")][el];
      isMedium = vfloats_values[makeName(ele_label,pref,"vidMedium")][el];
      isVeto = vfloats_values[makeName(ele_label,pref,"vidVeto")][el];
      float eta = vfloats_values[makeName(ele_label,pref,"Eta")][el];
      float scEta = vfloats_values[makeName(ele_label,pref,"scEta")][el];
      float phi = vfloats_values[makeName(ele_label,pref,"Phi")][el];
      float energy = vfloats_values[makeName(ele_label,pref,"E")][el];      
      float iso = vfloats_values[makeName(ele_label,pref,"Iso03")][el];

      float elCharge = vfloats_values[makeName(ele_label,pref,"Charge")][el];

      bool passesDRmu = true;
      bool passesTightCuts = false;
      if(fabs(scEta)<=1.479){
	passesTightCuts = isTight >0.0 && iso < 0.0588 ;
      } //is barrel electron
      if (fabs(scEta)>1.479){
	passesTightCuts = isTight >0.0 && iso < 0.0571 ;
      }

      if(pt> 30 && fabs(eta) < 2.1 && passesTightCuts){
	TLorentzVector ele;
	ele.SetPtEtaPhiE(pt, eta, phi, energy);	
	double minDR=999;
	double deltaR = 999;
	for (size_t m = 0; m < (size_t)loosemuons.size(); ++m){
	  deltaR = ele.DeltaR(loosemuons[m]);
	  minDR = min(minDR, deltaR);
	}
	if(!loosemuons.size()) minDR=999;
	if(minDR>0.1){ 
	  electrons.push_back(ele); 
	  flavors.push_back(11);
	  leptons.push_back(ele);
	  
	  leptonsCharge.push_back(elCharge);

	  ++float_values["Event_nTightElectrons"];
	  mapEle[lepidx]=el;
	  ++lepidx;
	  if(isInVector(obj_cats[ele_label],"Tight")){
	    fillCategory(ele_label,"Tight",el,float_values["Event_nTightElectrons"]-1);
	  }
	}
	else {passesDRmu = false;}
      }
      if(isInVector(obj_cats[ele_label],"Tight")){
	sizes[ele_label+"Tight"]=(int)float_values["Event_nTightElectrons"];
      }

      if(isLoose>0 && pt> 30 && fabs(eta) < 2.5){
	++float_values["Event_nLooseElectrons"];

      }

      if(isMedium>0 && pt> 30 && fabs(eta) < 2.5){
	++float_values["Event_nMediumElectrons"]; 
      }
      
      if(isVeto>0 && pt> 10 && fabs(eta) < 2.5 ){
	if((fabs(scEta)<=1.479 && (iso<0.175)) 
	   || ((fabs(scEta)>1.479) && (iso<0.159))){
	  ++float_values["Event_nVetoElectrons"]; 
	  if(isInVector(obj_cats[ele_label],"Veto")){
	    fillCategory(ele_label,"Veto",el,float_values["Event_nVetoElectrons"]-1);
	  }
	}
      }
      if(isInVector(obj_cats[ele_label],"Veto")){
	sizes[ele_label+"Veto"]=(int)float_values["Event_nVetoElectrons"];
      }
      
      vfloats_values[ele_label+"_PassesDRmu"][el]=(float)passesDRmu;
    } 
    int firstidx=-1, secondidx=-1;
    double maxpt=0.0;

    for(size_t l =0; l< leptons.size();++l){
      //double maxpt= 0.0;
      double lpt= leptons.at(l).Pt();
      if(lpt>maxpt){maxpt = lpt;firstidx=l;}
      
    }

    maxpt=0.0;
    for(size_t l =0; l< leptons.size();++l){
      //double maxpt= 0.0;
      double lpt= leptons.at(l).Pt();
      if(lpt>maxpt&&firstidx!=(int)l){maxpt = lpt;secondidx=l;}
    }
    if(firstidx>-1){
      float_values["Event_Lepton1_Pt"]=leptons.at(firstidx).Pt(); 
      float_values["Event_Lepton1_Phi"]=leptons.at(firstidx).Phi(); 
      float_values["Event_Lepton1_Eta"]=leptons.at(firstidx).Eta(); 
      float_values["Event_Lepton1_E"]=leptons.at(firstidx).E(); 
      float_values["Event_Lepton1_Flavour"]=flavors.at(firstidx);

      float_values["Event_Lepton1_Charge"]=leptonsCharge.at(firstidx);

    }
    if(secondidx>-1){
      float_values["Event_Lepton2_Pt"]=leptons.at(secondidx).Pt(); 
      float_values["Event_Lepton2_Phi"]=leptons.at(secondidx).Phi(); 
      float_values["Event_Lepton2_Eta"]=leptons.at(secondidx).Eta(); 
      float_values["Event_Lepton2_E"]=leptons.at(secondidx).E(); 
      float_values["Event_Lepton2_Flavour"]=flavors.at(secondidx);

      float_values["Event_Lepton2_Charge"]=leptonsCharge.at(secondidx);

    }

    //Jets:
    double Ht=0;
    double corrMetPx =0;
    double corrMetPy =0;
    double corrBaseMetPx =0;
    double corrBaseMetPy =0;
    double corrMetT1Px =0;
    double corrMetT1Py =0;
    double corrMetPxNoHF =0;
    double corrMetPyNoHF =0;
    double DUnclusteredMETPx=0.0;
    double DUnclusteredMETPy=0.0;

    string prefm = obj_to_pref[met_label];
    float metZeroCorrPt = vfloats_values[makeName(met_label,prefm,"UncorrPt")][0];
    float metZeroCorrPhi = vfloats_values[makeName(met_label,prefm,"UncorrPhi")][0];
    float metZeroCorrY = metZeroCorrPt*sin(metZeroCorrPhi);
    float metZeroCorrX = metZeroCorrPt*cos(metZeroCorrPhi);

    for(int j = 0;j < max_instances[jets_label] ;++j){
      string pref = obj_to_pref[jets_label];
      float pt = vfloats_values[makeName(jets_label,pref,"Pt")][j];
      float ptnomu = pt;
      float ptzero = vfloats_values[makeName(jets_label,pref,"Pt")][j];
      float genpt = vfloats_values[makeName(jets_label,pref,"GenJetPt")][j];
      float eta = vfloats_values[makeName(jets_label,pref,"Eta")][j];
      float phi = vfloats_values[makeName(jets_label,pref,"Phi")][j];
      float energy = vfloats_values[makeName(jets_label,pref,"E")][j];
     
      float ptCorr = -9999;
      float ptCorrSmearZero = -9999;
      float ptCorrSmearZeroNoMu = -9999;
      float energyCorr = -9999;
      float smearfact = -9999;

      float jecscale = vfloats_values[makeName(jets_label,pref,"jecFactor0")][j];
      float area = vfloats_values[makeName(jets_label,pref,"jetArea")][j];
      
      float juncpt=0.;
      float junce=0.;
      float ptCorr_mL1 = 0;

      float chEmEnFrac = vfloats_values[makeName(jets_label,pref,"chargedEmEnergyFrac")][j];
      float neuEmEnFrac = vfloats_values[makeName(jets_label,pref,"neutralEmEnergyFrac")][j];
      
      if(pt>0){

	TLorentzVector jetUncorr_, jetCorr, jetUncorrNoMu_,jetCorrNoMu, jetL1Corr;
	TLorentzVector T1Corr;
	T1Corr.SetPtEtaPhiE(0,0,0,0);	
	
	jetUncorr_.SetPtEtaPhiE(pt,eta,phi,energy);
	jetUncorrNoMu_.SetPtEtaPhiE(pt,eta,phi,energy);
	
	jetUncorr_= jetUncorr_*jecscale;
	jetUncorrNoMu_= jetUncorrNoMu_*jecscale;
	
	DUnclusteredMETPx+=jetUncorr_.Pt()*cos(phi);
	DUnclusteredMETPy+=jetUncorr_.Pt()*sin(phi);

	juncpt=jetUncorr_.Perp();
	junce=jetUncorr_.E();
	
	if(changeJECs){
	   
	  jecCorr->setJetPhi(jetUncorr_.Phi());
	  jecCorr->setJetEta(jetUncorr_.Eta());
	  jecCorr->setJetE(jetUncorr_.E());
	  jecCorr->setJetPt(jetUncorr_.Perp());
	  jecCorr->setJetA(area);
	  jecCorr->setRho(Rho);
	  jecCorr->setNPV(nPV);
	  
	  double recorr =  jecCorr->getCorrection();
	  jetCorr = jetUncorr_ *recorr;
	  
	  pt = jetCorr.Pt();
	  eta = jetCorr.Eta();
	  energy = jetCorr.Energy();
	  phi = jetCorr.Phi();
	 

	  for( size_t c=0;c<jetKeys->at(j).size();++c){
	    for( size_t mk=0;mk<muKeys->size();++mk){
	      if(muKeys->at(mk).size()>0){
		if(muKeys->at(mk).at(0)  == jetKeys->at(j).at(c)){
		  
		  string prefmu = obj_to_pref[mu_label];
		  float mupt = vfloats_values[makeName(mu_label,prefmu,"Pt")][mk];
		  float mueta = vfloats_values[makeName(mu_label,prefmu,"Eta")][mk];
		  float muphi = vfloats_values[makeName(mu_label,prefmu,"Phi")][mk];
		  float mue = vfloats_values[makeName(mu_label,prefmu,"E")][mk];
		  float muIsGlobal = vfloats_values[makeName(mu_label,prefmu,"IsGlobalMuon")][mk];
		  float muIsTK = vfloats_values[makeName(mu_label,prefmu,"IsTrackerMuon")][mk];
		  float muISSAOnly = ((!muIsGlobal && !muIsTK));
		  
		  if(muIsGlobal || muISSAOnly){
		    TLorentzVector muP4_;
		    muP4_.SetPtEtaPhiE(mupt,mueta,muphi,mue);
		    jetUncorrNoMu_ -=muP4_;
		    
		  } 
		}
	      }
	    }	    
	    
	    jecCorr->setJetPhi(jetUncorrNoMu_.Phi());
	    jecCorr->setJetEta(jetUncorrNoMu_.Eta());
	    jecCorr->setJetE(jetUncorrNoMu_.E());
	    jecCorr->setJetPt(jetUncorrNoMu_.Perp());
	    jecCorr->setJetA(area);
	    jecCorr->setRho(Rho);
	    jecCorr->setNPV(nPV);
	    
	    double recorrMu =  jecCorr->getCorrection();
	    jetCorrNoMu = jetUncorrNoMu_ * recorrMu;
	    
	    //// Jet corrections for level 1
	    jecCorr_L1->setJetPhi(jetUncorrNoMu_.Phi()); /// deve essere raw
	    jecCorr_L1->setJetEta(jetUncorrNoMu_.Eta());
	    jecCorr_L1->setJetE(jetUncorrNoMu_.E());
	    jecCorr_L1->setJetPt(jetUncorrNoMu_.Perp());
	    jecCorr_L1->setJetA(area);
	    jecCorr_L1->setRho(Rho);
	    jecCorr_L1->setNPV(nPV);
	    
	    double recorr_L1 =  jecCorr_L1->getCorrection();
	    jetL1Corr = jetUncorrNoMu_ * recorr_L1;
	      
	    //cout << "AK4 pt Corr: " << recorr << " MET Corr noL1: " << recorr_L1 << endl;
	    ptnomu = jetCorrNoMu.Pt();
	    if(pt>15.0 && ( chEmEnFrac + neuEmEnFrac <0.9)){ 
	      T1Corr += jetCorrNoMu - jetL1Corr;
	      ptCorr_mL1 = T1Corr.Pt();
	    }
	  }
	}
	
	smearfact = smear(pt, genpt, eta, syst);
	
	ptCorr = pt * smearfact;
	energyCorr = energy * smearfact;
	
	float unc = jetUncertainty(ptCorr,eta,syst);
	float unc_nosmear = jetUncertainty(pt,eta,syst);
	
	ptCorr = ptCorr * (1 + unc);
	ptCorrSmearZero = pt * (1 + unc);//For full correction, including new JECs and MET
	ptCorrSmearZeroNoMu = ptnomu * (1 + unc);//For full correction, including new JECs and MET
	ptCorr_mL1 = ptCorr_mL1 * (1 + unc);//For T1 correction

	energyCorr = energyCorr * (1 + unc);
	corrMetPx -=(cos(phi)*(pt*unc_nosmear));//Difference between jesUp/down and nominal
        corrMetPy -=(sin(phi)*(pt*unc_nosmear));

	if(ptCorrSmearZeroNoMu>15.0 && ( chEmEnFrac+  neuEmEnFrac <0.9)){ // should be == T1 corrections minus the L1 difference.
	  corrBaseMetPx -=(cos(phi)*(ptCorrSmearZero-ptzero));
	  corrBaseMetPy -=(sin(phi)*(ptCorrSmearZero-ptzero));
	}
	
	if( ptCorrSmearZeroNoMu>15.0 && jetCorrNoMu.Pt()>0.0){ 
	  corrMetT1Px -= (T1Corr.Px()) + (cos(phi) * (ptCorr_mL1 -  ptCorrSmearZero));
	  corrMetT1Py -= (T1Corr.Py()) + (sin(phi) * (ptCorr_mL1 -  ptCorrSmearZero));
	}
	
	if(fabs(eta)<3.0){
	  corrMetPxNoHF -=(cos(phi)*(ptCorr-ptzero));
	  corrMetPyNoHF -=(sin(phi)*(ptCorr-ptzero));
	}
	
      }//closing pt>0
          
      float csv = vfloats_values[makeName(jets_label,pref,"CSVv2")][j];
      float partonFlavour = vfloats_values[makeName(jets_label,pref,"PartonFlavour")][j];
      int flavor = int(partonFlavour);
  
      if(getWZFlavour){
	if(abs(flavor)==5){++nb;}
	else{ 
	  if(abs(flavor)==4){++nc;}
	  else {++nudsg;}
	}
      }
 
      vfloats_values[jets_label+"_NoCorrPt"][j]=juncpt;
      vfloats_values[jets_label+"_NoCorrE"][j]=junce;

      vfloats_values[jets_label+"_CorrPt"][j]=ptCorr;
      vfloats_values[jets_label+"_CorrE"][j]=energyCorr;
      vfloats_values[jets_label+"_CorrEta"][j]=eta;
      vfloats_values[jets_label+"_CorrPhi"][j]=phi;

      bool isCSVT = csv  > 0.9535;
      bool isCSVM = csv  > 0.8484;
      bool isCSVL = csv  > 0.5426;
      vfloats_values[jets_label+"_IsCSVT"][j]=isCSVT;
      vfloats_values[jets_label+"_IsCSVM"][j]=isCSVM;
      vfloats_values[jets_label+"_IsCSVL"][j]=isCSVL;
      
      float bsf = getScaleFactor(ptCorr,eta,partonFlavour,"noSyst");
      float bsfup = getScaleFactor(ptCorr,eta,partonFlavour,"up");
      float bsfdown = getScaleFactor(ptCorr,eta,partonFlavour,"down");
      
      vfloats_values[jets_label+"_BSF"][j]=bsf;
      vfloats_values[jets_label+"_BSFUp"][j]=bsfup;
      vfloats_values[jets_label+"_BSFDown"][j]=bsfdown;
      
      bool passesID = true;
      
      if(!(jecscale*energy > 0))passesID = false;
      else{
        float neuMulti = vfloats_values[makeName(jets_label,pref,"neutralMultiplicity")][j];
        float chMulti = vfloats_values[makeName(jets_label,pref,"chargedMultiplicity")][j];
        float chHadEnFrac = vfloats_values[makeName(jets_label,pref,"chargedHadronEnergyFrac")][j];
        float chEmEnFrac = vfloats_values[makeName(jets_label,pref,"chargedEmEnergyFrac")][j];
        float neuEmEnFrac = vfloats_values[makeName(jets_label,pref,"neutralEmEnergyFrac")][j];
        float neuHadEnFrac = vfloats_values[makeName(jets_label,pref,"neutralHadronEnergyFrac")][j];
	float numConst = chMulti + neuMulti;

        if(fabs(eta)<=2.7){
          passesID =  (neuHadEnFrac<0.90 && neuEmEnFrac<0.90 && numConst>1) && ( (abs(eta)<=2.4 && chHadEnFrac>0 && chMulti>0 && chEmEnFrac<0.99) || abs(eta)>2.4);
        }
	else if(fabs(eta)>2.7 && fabs(eta)<=3.0){
	  passesID = neuMulti > 2 && neuHadEnFrac < 0.98 && neuEmEnFrac > 0.01;
	}
        else if(fabs(eta)>3){
          passesID = neuEmEnFrac<0.90 && neuMulti>10 ;
        }
      }
      
      vfloats_values[jets_label+"_PassesID"][j]=(float)passesID;
      
      //Remove overlap with tight electrons/muons
      double minDR=9999;
      double minDRThrEl=0.3;
      double minDRThrMu=0.4;
      bool passesDR=true;
 
      for (size_t e = 0; e < (size_t)electrons.size(); ++e){
	minDR = min(minDR,deltaR(math::PtEtaPhiELorentzVector(electrons.at(e).Pt(),electrons.at(e).Eta(),electrons.at(e).Phi(),electrons.at(e).Energy() ) ,math::PtEtaPhiELorentzVector(ptCorr, eta, phi, energyCorr)));
	if(minDR<minDRThrEl)passesDR = false;
      }
      for (size_t m = 0; m < (size_t)muons.size(); ++m){
	minDR = min(minDR,deltaR(math::PtEtaPhiELorentzVector(muons.at(m).Pt(),muons.at(m).Eta(),muons.at(m).Phi(),muons.at(m).Energy() ) ,math::PtEtaPhiELorentzVector(ptCorr, eta, phi, energyCorr)));
	//minDR = min(minDR,deltaR(muons.at(m) ,math::PtEtaPhiELorentzVector(ptCorr, eta, phi, energyCorr)));
	if(minDR<minDRThrMu)passesDR = false;
      }
      
      vfloats_values[jets_label+"_MinDR"][j]=minDR;
      vfloats_values[jets_label+"_PassesDR"][j]=(float)passesDR;
      
      vfloats_values[jets_label+"_IsTight"][j]=0.0;
      vfloats_values[jets_label+"_IsLoose"][j]=0.0;
	
      if(passesID && passesDR) vfloats_values[jets_label+"_IsLoose"][j]=1.0;

      if(passesID && passesDR && pt>50 && abs(eta)<2.4){
	Ht+=pt;
      }

      float_values["Event_Ht"] = (float)Ht;
      
      for (size_t ji = 0; ji < (size_t)jetScanCuts.size(); ++ji){
	stringstream j_n;
	double jetval = jetScanCuts.at(ji);
	j_n << "Cut" <<jetval;
	bool passesCut = ( ptCorr > jetval && fabs(eta) < 4.);

	if(!passesID || !passesCut || !passesDR) continue;
	if(ji==0){
	  vfloats_values[jets_label+"_IsTight"][j]=1.0;
	  TLorentzVector jet;
	  jet.SetPtEtaPhiE(ptCorr, eta, phi, energyCorr);
	  jets.push_back(jet);

	  if(!isCSVM)     jetsnob.push_back(jet);

	}

	if(passesCut &&  passesID && passesDR){	
	  float_values["Event_nJets"+j_n.str()]+=1;if(ji==0){
	    nTightJets+=1;
	    if(isInVector(obj_cats[jets_label],"Tight")){
	      fillCategory(jets_label,"Tight",j,float_values["Event_nJets"+j_n.str()]-1);
	    }
	  }
	}

	if(passesCut &&  passesID && passesDR){
	  double csvteff = MCTagEfficiency("csvt",flavor, ptCorr, eta);
	  double sfcsvt = TagScaleFactor("csvt", flavor, "noSyst", ptCorr);
	  
	  double csvleff = MCTagEfficiency("csvl",flavor,ptCorr, eta);
	  double sfcsvl = TagScaleFactor("csvl", flavor, "noSyst", ptCorr);

	  double csvmeff = MCTagEfficiency("csvm",flavor,ptCorr, eta);
	  double sfcsvm = TagScaleFactor("csvm", flavor, "noSyst", ptCorr);

	  double sfcsvt_mistag_up = TagScaleFactor("csvt", flavor, "mistag_up", ptCorr);
	  double sfcsvl_mistag_up = TagScaleFactor("csvl", flavor, "mistag_up", ptCorr);
	  double sfcsvm_mistag_up = TagScaleFactor("csvm", flavor, "mistag_up", ptCorr);

	  double sfcsvt_mistag_down = TagScaleFactor("csvt", flavor, "mistag_down", ptCorr);
	  double sfcsvl_mistag_down = TagScaleFactor("csvl", flavor, "mistag_down", ptCorr);
	  double sfcsvm_mistag_down = TagScaleFactor("csvm", flavor, "mistag_down", ptCorr);

	  double sfcsvt_b_tag_down = TagScaleFactor("csvt", flavor, "b_tag_down", ptCorr);
	  double sfcsvl_b_tag_down = TagScaleFactor("csvl", flavor, "b_tag_down", ptCorr);
	  double sfcsvm_b_tag_down = TagScaleFactor("csvm", flavor, "b_tag_down", ptCorr);
	  
	  double sfcsvt_b_tag_up = TagScaleFactor("csvt", flavor, "b_tag_up", ptCorr);
	  double sfcsvl_b_tag_up = TagScaleFactor("csvl", flavor, "b_tag_up", ptCorr);
	  double sfcsvm_b_tag_up = TagScaleFactor("csvm", flavor, "b_tag_up", ptCorr);

	  jsfscsvt.push_back(BTagWeight::JetInfo(csvteff, sfcsvt));
	  jsfscsvl.push_back(BTagWeight::JetInfo(csvleff, sfcsvl));
	  jsfscsvm.push_back(BTagWeight::JetInfo(csvmeff, sfcsvm));

	  jsfscsvt_mistag_up.push_back(BTagWeight::JetInfo(csvteff, sfcsvt_mistag_up));
	  jsfscsvl_mistag_up.push_back(BTagWeight::JetInfo(csvleff, sfcsvl_mistag_up));
	  jsfscsvm_mistag_up.push_back(BTagWeight::JetInfo(csvmeff, sfcsvm_mistag_up));

	  jsfscsvt_b_tag_up.push_back(BTagWeight::JetInfo(csvteff, sfcsvt_b_tag_up));
	  jsfscsvl_b_tag_up.push_back(BTagWeight::JetInfo(csvleff, sfcsvl_b_tag_up));
	  jsfscsvm_b_tag_up.push_back(BTagWeight::JetInfo(csvmeff, sfcsvm_b_tag_up));

	  jsfscsvt_mistag_down.push_back(BTagWeight::JetInfo(csvteff, sfcsvt_mistag_down));
	  jsfscsvl_mistag_down.push_back(BTagWeight::JetInfo(csvleff, sfcsvl_mistag_down));
	  jsfscsvm_mistag_down.push_back(BTagWeight::JetInfo(csvmeff, sfcsvm_mistag_down));

	  jsfscsvt_b_tag_down.push_back(BTagWeight::JetInfo(csvteff, sfcsvt_b_tag_down));
	  jsfscsvl_b_tag_down.push_back(BTagWeight::JetInfo(csvleff, sfcsvl_b_tag_down));
	  jsfscsvm_b_tag_down.push_back(BTagWeight::JetInfo(csvmeff, sfcsvm_b_tag_down));

	}
	
	if(isCSVT && passesCut &&  passesID && passesDR && fabs(eta) < 2.4) {
	  float_values["Event_nCSVTJets"+j_n.str()]+=1.0;
	  ncsvt_tags +=1;
	}

	if(isCSVL && passesCut &&  passesID && passesDR && fabs(eta) < 2.4) { 
	  ncsvl_tags +=1;
	}
	if(isCSVM && passesCut &&  passesID && passesDR && fabs(eta) < 2.4) { 
	  float_values["Event_nCSVMJets"+j_n.str()]+=1.0;
	  if(ji==0){
	    ncsvm_tags +=1;
	    TLorentzVector bjet;
	    bjet.SetPtEtaPhiE(ptCorr, eta, phi, energyCorr);
	    bjets.push_back(bjet);
	    mapBJets[bjetidx]=j;
	    ++bjetidx;
	  }
	}
	
	if(isCSVL && passesCut &&  passesID && passesDR && abs(eta) < 2.4) float_values["Event_nCSVLJets"+j_n.str()]+=1;
	
      }
    }
 
    if(isInVector(obj_cats[jets_label],"Tight")){
      sizes[jets_label+"Tight"]=(int)nTightJets;
    }
    
    float_values["Event_eventFlavour"]=eventFlavour(getWZFlavour, nb, nc, nudsg);
    if(syst.find("unclusteredMet")!= std::string::npos ){
      
      DUnclusteredMETPx=metZeroCorrX+DUnclusteredMETPx;
      DUnclusteredMETPy=metZeroCorrY+DUnclusteredMETPy;
      
      double signmet = 1.0; 
      if(syst.find("down")!=std::string::npos) signmet=-1.0;
      corrMetPx -=signmet*DUnclusteredMETPx*0.1;
      corrMetPy -=signmet*DUnclusteredMETPy*0.1;
    }
 
    string pref = obj_to_pref[met_label];
    float metpt = vfloats_values[makeName(met_label,pref,"Pt")][0];
    float metphi = vfloats_values[makeName(met_label,pref,"Phi")][0];
    
    float metPyCorrBase = metpt*sin(metphi);
    float metPxCorrBase = metpt*cos(metphi);
    metPxCorrBase+=corrBaseMetPx; metPyCorrBase+=corrBaseMetPy; // add JEC/JER contribution

    float metPyCorr = metpt*sin(metphi);
    float metPxCorr = metpt*cos(metphi);
    metPxCorr+=corrMetPx; metPyCorr+=corrMetPy; // add JEC/JER contribution

    float metPx = metPxCorr;
    float metPy = metPyCorr;

    float metptunc = vfloats_values[makeName(met_label,pref,"uncorPt")][0];
    float metphiunc = vfloats_values[makeName(met_label,pref,"uncorPhi")][0];

    float metT1Py = metptunc*sin(metphiunc);
    float metT1Px = metptunc*cos(metphiunc);

    //Correcting the pt
    metT1Px+=corrMetT1Px; metT1Py+=corrMetT1Py; // add JEC/JER contribution

    float metptT1Corr = sqrt(metT1Px*metT1Px + metT1Py*metT1Py);
    vfloats_values[met_label+"_CorrT1Pt"][0]=metptT1Corr;
    
    float metptCorr = sqrt(metPxCorr*metPxCorr + metPyCorr*metPyCorr);
    vfloats_values[met_label+"_CorrPt"][0]=metptCorr;

    float metptCorrBase = sqrt(metPxCorrBase*metPxCorrBase + metPyCorrBase*metPyCorrBase);
    vfloats_values[met_label+"_CorrBasePt"][0]=metptCorrBase;
    
    //Correcting the phi
    float metphiCorr = metphi;
    if(metPxCorr<0){
      if(metPyCorr>0)metphiCorr = atan(metPyCorr/metPxCorr)+3.141592;
      if(metPyCorr<0)metphiCorr = atan(metPyCorr/metPxCorr)-3.141592;
    }
    else  metphiCorr = (atan(metPyCorr/metPxCorr));
    
    float metphiCorrBase = metphi;
    if(metPxCorrBase<0){
      if(metPyCorrBase>0)metphiCorrBase = atan(metPyCorrBase/metPxCorrBase)+3.141592;
      if(metPyCorrBase<0)metphiCorrBase = atan(metPyCorrBase/metPxCorrBase)-3.141592;
    }
    else  metphiCorrBase = (atan(metPyCorrBase/metPxCorrBase));
    
    float metphiCorrT1 = metphi;
    if(metT1Px<0){
      if(metT1Py>0)metphiCorrT1 = atan(metT1Py/metT1Px)+3.141592;
      if(metT1Py<0)metphiCorrT1 = atan(metT1Py/metT1Px)-3.141592;
    }
    else  metphiCorr = (atan(metPyCorr/metPxCorr));
    
    vfloats_values[met_label+"_CorrPhi"][0]=metphiCorr;
    vfloats_values[met_label+"_CorrBasePhi"][0]=metphiCorrBase;
    vfloats_values[met_label+"_CorrT1Phi"][0]=metphiCorrT1;

    //Preselection part

    if(doPreselection){
      bool passes = true;
      bool metCondition = (metptCorr > 100.0 || Ht > 400.);

      float lep1phi = float_values["Event_Lepton1_Phi"];
      float lep1pt = float_values["Event_Lepton1_Pt"];

      float lep2phi = float_values["Event_Lepton2_Phi"];
      float lep2pt = float_values["Event_Lepton2_Pt"];

      float lep1px = 0.0;
      float lep1py = 0.0;
      float lep2px = 0.0;
      float lep2py = 0.0;
      float metPxLep=metPx;
      float metPyLep=metPy;
      
      if(lep1pt>0.0){
	lep1px = lep1pt*cos(lep1phi);
	lep1py = lep1pt*sin(lep1phi);

	if(lep2pt>0.0){
	  lep2px = lep2pt*cos(lep2phi);
	  lep2py = lep2pt*sin(lep2phi);
	}	
      } 
    
      metPxLep+=lep1px;
      metPyLep+=lep1py;

      metPxLep+=lep2px;
      metPyLep+=lep2py;
      
      double metLep=sqrt(metPxLep*metPxLep+metPyLep*metPyLep);
      
      metPxLep+=lep1px;

      metCondition = metCondition || (metLep>100.0);
	
      passes = passes && metCondition;

      if (!passes ) {
	//Reset event weights/#objects
	string nameshortv= "Event";
	vector<string> extravars = additionalVariables(nameshortv);
	for(size_t addv = 0; addv < extravars.size();++addv){
	  string name = nameshortv+"_"+extravars.at(addv);
	  if(!isMCWeightName(extravars.at(addv))) float_values[name]=0.0;

	}
	continue;
      }
    }
    
    TLorentzVector lepton;
    
    if( ( (electrons.size()==1 && muons.size()==0 ) || (muons.size()==1 && electrons.size()==0) ) && bjets.size()>0 ){
      if(electrons.size()==1) lepton = electrons[0];
      else if(muons.size()==1) lepton = muons[0];
      
      TVector2 met( metptCorr*cos(metphiCorr), metptCorr*sin(metphiCorr));
      float phi_lmet = fabs(deltaPhi(lepton.Phi(), metphiCorr) );
      float mt = sqrt(2* lepton.Pt() * metptCorr * ( 1- cos(phi_lmet)));
      float_values["Event_mt"] = (float)mt;
      Mt2Com_bisect *Mt2cal = new Mt2Com_bisect();
      double Mt2w = Mt2cal->calculateMT2w(jetsnob,bjets,lepton, met,"MT2w");
      float_values["Event_Mt2w"] = (float)Mt2w;    
    }

    for(int s = 0;s < min(max_instances[boosted_tops_subjets_label],sizes[boosted_tops_subjets_label]) ;++s){
      string pref = obj_to_pref[boosted_tops_subjets_label];
      float pt  = vfloats_values[makeName(boosted_tops_subjets_label,pref,"Pt")][s];
      float eta = vfloats_values[makeName(boosted_tops_subjets_label,pref,"Eta")][s];
      float phi = vfloats_values[makeName(boosted_tops_subjets_label,pref,"Phi")][s];
      float e   = vfloats_values[makeName(boosted_tops_subjets_label,pref,"E")][s];

      float partonFlavourSubjet = vfloats_values[makeName(boosted_tops_subjets_label,pref,"PartonFlavour")][s];
      int flavorSubjet = int(partonFlavourSubjet);

      float bsfsubj = getScaleFactor(pt,eta,partonFlavourSubjet,"noSyst");
      float bsfupsubj = getScaleFactor(pt,eta,partonFlavourSubjet,"up");
      float bsfdownsubj = getScaleFactor(pt,eta,partonFlavourSubjet,"down");
      
      vfloats_values[boosted_tops_subjets_label+"_BSF"][s]=bsfsubj;
      vfloats_values[boosted_tops_subjets_label+"_BSFUp"][s]=bsfupsubj;
      vfloats_values[boosted_tops_subjets_label+"_BSFDown"][s]=bsfdownsubj;
      
      TLorentzVector subjet;
      subjet.SetPtEtaPhiE(pt, eta, phi, e);       
      double minDR=999;
      float subjcsv = vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][s];
     
      bool isCSVM = (subjcsv>0.8484);

      vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][s] = (float)subjcsv;
      
      if(subjcsv>0.5426 && fabs(eta) < 2.4) {
	ncsvl_subj_tags +=1;
      }
      
      if(subjcsv>0.8484 && fabs(eta) < 2.4) {
	ncsvm_subj_tags +=1;
      }
      
      double csvleff_subj = MCTagEfficiencySubjet("csvl",flavorSubjet,pt, eta);
      double sfcsvl_subj = TagScaleFactorSubjet("csvl", flavorSubjet, "noSyst", pt);
      
      double csvmeff_subj = MCTagEfficiencySubjet("csvm",flavorSubjet,pt, eta);
      double sfcsvm_subj = TagScaleFactorSubjet("csvm", flavorSubjet, "noSyst", pt);

      double sfcsvl_mistag_up_subj = TagScaleFactorSubjet("csvl", flavorSubjet, "mistag_up", pt);
      double sfcsvm_mistag_up_subj = TagScaleFactorSubjet("csvm", flavorSubjet, "mistag_up", pt);

      double sfcsvl_mistag_down_subj = TagScaleFactorSubjet("csvl", flavorSubjet, "mistag_down", pt);
      double sfcsvm_mistag_down_subj = TagScaleFactorSubjet("csvm", flavorSubjet, "mistag_down", pt);

      double sfcsvl_b_tag_down_subj = TagScaleFactorSubjet("csvl", flavorSubjet, "b_tag_down", pt);
      double sfcsvm_b_tag_down_subj = TagScaleFactorSubjet("csvm", flavorSubjet, "b_tag_down", pt);
      
      double sfcsvl_b_tag_up_subj = TagScaleFactorSubjet("csvl", flavorSubjet, "b_tag_up", pt);
      double sfcsvm_b_tag_up_subj = TagScaleFactorSubjet("csvm", flavorSubjet, "b_tag_up", pt);     

      jsfscsvl_subj.push_back(BTagWeight::JetInfo(csvleff_subj, sfcsvl_subj));
      jsfscsvm_subj.push_back(BTagWeight::JetInfo(csvmeff_subj, sfcsvm_subj));
      jsfscsvl_subj_mistag_up.push_back(BTagWeight::JetInfo(csvleff_subj, sfcsvl_mistag_up_subj));
      jsfscsvm_subj_mistag_up.push_back(BTagWeight::JetInfo(csvmeff_subj, sfcsvm_mistag_up_subj));
      
      jsfscsvl_subj_b_tag_up.push_back(BTagWeight::JetInfo(csvleff_subj, sfcsvl_b_tag_up_subj));
      jsfscsvm_subj_b_tag_up.push_back(BTagWeight::JetInfo(csvmeff_subj, sfcsvm_b_tag_up_subj));
      
      jsfscsvl_subj_mistag_down.push_back(BTagWeight::JetInfo(csvleff_subj, sfcsvl_mistag_down_subj));
      jsfscsvm_subj_mistag_down.push_back(BTagWeight::JetInfo(csvmeff_subj, sfcsvm_mistag_down_subj));
      
      jsfscsvl_subj_b_tag_down.push_back(BTagWeight::JetInfo(csvleff_subj, sfcsvl_b_tag_down_subj));
      jsfscsvm_subj_b_tag_down.push_back(BTagWeight::JetInfo(csvmeff_subj, sfcsvm_b_tag_down_subj));
           
      
      for(int t = 0;t < min(max_instances[boosted_tops_label],sizes[boosted_tops_label]) ;++t){
	  
	  float ptj = vfloats_values[makeName(boosted_tops_label,pref,"Pt")][t];
	  if (ptj<0.0)continue;
	  float etaj = vfloats_values[makeName(boosted_tops_label,pref,"Eta")][t];
	  float phij = vfloats_values[makeName(boosted_tops_label,pref,"Phi")][t];
	  float ej = vfloats_values[makeName(boosted_tops_label,pref,"E")][t];
	  TLorentzVector topjet;
	  topjet.SetPtEtaPhiE(ptj, etaj, phij, ej);       
	  
	  float DR = subjet.DeltaR(topjet); 
	    if(DR < minDR){
	    minDR = DR;
	    subj_jet_map[s]=t;
	    }
      }
	size_t tm = subj_jet_map[s];
	if(isCSVM)vfloats_values[boosted_tops_label+"_nCSVM"][tm]+=1;
	vfloats_values[boosted_tops_label+"_nJ"][tm]+=1;
    }
    
    for(int t = 0;t < max_instances[boosted_tops_label] ;++t){
      string pref = obj_to_pref[boosted_tops_label];
      float prunedMass   = vfloats_values[makeName(boosted_tops_label,pref,"prunedMassCHS")][t];
      float softDropMass = vfloats_values[makeName(boosted_tops_label,pref,"softDropMass")][t];

      float topPt        = vfloats_values[makeName(boosted_tops_label,pref,"Pt")][t];
      float topEta = vfloats_values[makeName(boosted_tops_label,pref,"Eta")][t];
      float topPhi = vfloats_values[makeName(boosted_tops_label,pref,"Phi")][t];
      float topE = vfloats_values[makeName(boosted_tops_label,pref,"E")][t];
      float tau1         = vfloats_values[makeName(boosted_tops_label,pref,"tau1")][t];
      float tau2         = vfloats_values[makeName(boosted_tops_label,pref,"tau2")][t];
      float tau3         = vfloats_values[makeName(boosted_tops_label,pref,"tau3")][t];

      float genpt8 = vfloats_values[makeName(boosted_tops_label,pref,"GenJetPt")][t];

      float jecscale8 =  vfloats_values[makeName(boosted_tops_label,pref,"jecFactor0")][t];
      float area8 =  vfloats_values[makeName(boosted_tops_label,pref,"jetArea")][t];
      
      float ptCorr8 = -9999;
      float energyCorr8 = -9999;
      float smearfact8 = -9999;
      float prunedMassCorr=-9999;
      
      float juncpt8=0.;
      float junce8=0.;
      
      if(topPt>0){	
      	TLorentzVector jetUncorr8_, jetCorr8,  jetUncorrNoMu8_, jetNoL1Corr8;
	
	jetUncorr8_.SetPtEtaPhiE(topPt, topEta, topPhi, topE);
	jetUncorrNoMu8_.SetPtEtaPhiE(topPt, topEta, topPhi, topE);
	
	jetUncorr8_= jetUncorr8_*jecscale8;
	jetUncorrNoMu8_= jetUncorrNoMu8_*jecscale8;
	
	juncpt8=jetUncorr8_.Perp();
	junce8=jetUncorr8_.E();

	if(changeJECs){
	  
	  jecCorr8->setJetPhi(jetUncorr8_.Phi());
	  jecCorr8->setJetEta(jetUncorr8_.Eta());
	  jecCorr8->setJetE(jetUncorr8_.E());
	  jecCorr8->setJetPt(jetUncorr8_.Perp());
	  jecCorr8->setJetA(area8);
	  jecCorr8->setRho(Rho);
	  jecCorr8->setNPV(nPV);
	  
	  double recorr8 =  jecCorr8->getCorrection();
	  jetCorr8 = jetUncorr8_ *recorr8;
	  
	  topPt = jetCorr8.Pt();
	  topEta = jetCorr8.Eta();
	  topE = jetCorr8.Energy();
	  topPhi = jetCorr8.Phi();
	  
	  //// Jet corrections without level 1
	  jecCorr_NoL18->setJetPhi(jetUncorr8_.Phi()); /// deve essere raw
	  jecCorr_NoL18->setJetEta(jetUncorr8_.Eta());
	  jecCorr_NoL18->setJetE(jetUncorr8_.E());
	  jecCorr_NoL18->setJetPt(jetUncorr8_.Perp());
	  jecCorr_NoL18->setJetA(area8);
	  jecCorr_NoL18->setRho(Rho);
	  jecCorr_NoL18->setNPV(nPV);
	  
	  double recorr_NoL18 =  jecCorr_NoL18->getCorrection();

	  cout << "AK8 pt Corr: " << recorr8 << " AK8 pt Corr NoL1: " << recorr_NoL18 << endl;
	  prunedMassCorr = recorr_NoL18 * prunedMass;
	  
	}//close changeJECs

	smearfact8 = smear(topPt, genpt8, topEta, syst); 
	
	ptCorr8 = topPt * smearfact8;
	energyCorr8 = topE * smearfact8;
	float unc8 = jetUncertainty8(ptCorr8,topEta,syst);
	cout << "topPt: " << topPt << " smearfact8: " << smearfact8 << " genPt: " << genpt8 << " topEta: " << topEta<<  endl; 
	ptCorr8 = ptCorr8 * (1 + unc8);
	energyCorr8 = energyCorr8 * (1 + unc8);

      }//close pt>0

      vfloats_values[boosted_tops_label+"_NoCorrPt"][t]=juncpt8;
      vfloats_values[boosted_tops_label+"_NoCorrE"][t]=junce8;
      
      vfloats_values[boosted_tops_label+"_CorrPrunedMassCHS"][t]=prunedMassCorr;
      vfloats_values[boosted_tops_label+"_CorrPt"][t]=ptCorr8;
      vfloats_values[boosted_tops_label+"_CorrE"][t]=energyCorr8;
      
      float tau3OVERtau2 = (tau2!=0. ? tau3/tau2 : 9999.);
      float tau2OVERtau1 = (tau1!=0. ? tau2/tau1 : 9999.);

      vfloats_values[makeName(boosted_tops_label,pref,"tau3OVERtau2")][t]=(float)tau3OVERtau2;
      vfloats_values[makeName(boosted_tops_label,pref,"tau2OVERtau1")][t]=(float)tau2OVERtau1;
      
      math::PtEtaPhiELorentzVector p4bestTop;
      math::PtEtaPhiELorentzVector p4bestB;
      
      int indexv0 = vfloats_values[makeName(boosted_tops_label,pref,"vSubjetIndex0")][t];
      int indexv1 = vfloats_values[makeName(boosted_tops_label,pref,"vSubjetIndex1")][t];

      int nCSVsubj = 0;
      if( vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv0] > 0.8484) ++nCSVsubj;
      if( vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv1] > 0.8484) ++nCSVsubj;

      int nCSVsubj_tm = 0;

      if( vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv0] > 0.5426 && vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv0]<0.8484) ++nCSVsubj_tm;
      if( vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv1] > 0.5426 && vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv1] < 0.8484) ++nCSVsubj_tm;

      int nCSVsubj_t = 0;
      if(vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv0]<0.5426) ++nCSVsubj_t;
      if(vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][indexv1]< 0.5426) ++nCSVsubj_t;
      
      vfloats_values[makeName(boosted_tops_label,pref,"nCSVsubj")][t]=(float)nCSVsubj;
      vfloats_values[makeName(boosted_tops_label,pref,"nCSVsubj_tm")][t]=(float)nCSVsubj_tm;
    
      bool isTop = ( ( softDropMass <= 220 && softDropMass >=105 )
		     and (tau3OVERtau2 <= 0.81 )
		     and ( topPt > 500)
		     );

       bool isW = ( ( prunedMass <= 95 && prunedMass >=65 )
		    and (tau2OVERtau1 <= 0.6) //0.5
		    and ( topPt > 200)
		    );
    
      math::PtEtaPhiELorentzVector p4ak8;
   
       if (isW) {
	p4ak8 = math::PtEtaPhiELorentzVector(vfloats_values[makeName(boosted_tops_label,pref,"Pt")][t],
					     vfloats_values[makeName(boosted_tops_label,pref,"Eta")][t],
					     vfloats_values[makeName(boosted_tops_label,pref,"Phi")][t],
					     vfloats_values[makeName(boosted_tops_label,pref,"E")][t]);

	float bestTopMass = 0.;
	float dRmin = 2.6;
	
	for(int i = 0; i < sizes[jets_label]; ++i){
	 
	  math::PtEtaPhiELorentzVector p4ak4 = math::PtEtaPhiELorentzVector(vfloats_values[jets_label+"_CorrPt"][i], 
									    vfloats_values[jets_label+"_CorrEta"][i], 
									    vfloats_values[jets_label+"_CorrPhi"][i], 
									    vfloats_values[jets_label+"_CorrE"][i] );
	  double dR = ROOT::Math::VectorUtil::DeltaR(p4ak8,p4ak4);

	  if (dR <= 0.8 or dR > 2.5) continue;

	  math::PtEtaPhiELorentzVector p4top = (p4ak8+p4ak4);
	  float topMass = p4top.mass();
	  if ( dR < dRmin ) {
	    dRmin = dR;
	    bestTopMass = topMass;
	    p4bestB= p4ak4;
	  }
	}

	if (bestTopMass > 250 or bestTopMass < 140) isW=false; 
      }

      if(isW){
	p4bestTop = p4bestB + p4ak8;
      }
      if(isTop) p4bestTop = p4ak8;
      
      vfloats_values[makeName(boosted_tops_label,pref,"isType2")][t]=(float)isW;
      vfloats_values[makeName(boosted_tops_label,pref,"isType1")][t]=(float)isTop;

      if(isW || isTop){
	TLorentzVector topjet;
	float ptj  = vfloats_values[makeName(boosted_tops_label,pref,"Pt")][t];
	float etaj = vfloats_values[makeName(boosted_tops_label,pref,"Eta")][t];
	float phij = vfloats_values[makeName(boosted_tops_label,pref,"Phi")][t];
	float ej   = vfloats_values[makeName(boosted_tops_label,pref,"E")][t];
	
	vfloats_values[makeName(boosted_tops_label,pref,"TopPt")][t]   = p4bestTop.pt();
	vfloats_values[makeName(boosted_tops_label,pref,"TopEta")][t]  = p4bestTop.eta();
	vfloats_values[makeName(boosted_tops_label,pref,"TopPhi")][t]  = p4bestTop.phi();
	vfloats_values[makeName(boosted_tops_label,pref,"TopE")][t]    = p4bestTop.energy();
	vfloats_values[makeName(boosted_tops_label,pref,"TopMass")][t] = p4bestTop.mass();
	if (isW)
	  vfloats_values[makeName(boosted_tops_label,pref,"TopWMass")][t] = vfloats_values[makeName(boosted_tops_label,pref,"prunedMass")][t];


	topjet.SetPtEtaPhiE(ptj, etaj, phij, ej);       
	if(isW){
	  float_values["Event_nType2TopJets"]+=1;
	  type2topjets.push_back(topjet);
	}
	if(isTop){
	  float_values["Event_nType1TopJets"]+=1;
	  type1topjets.push_back(topjet);
	}
      }
    }
    
    //    int nTightLeptons = electrons.size()+muons.size();
    size_t cat = 0;

    size_t ni = 9;
    cat+= 100000*(min(ni,muons.size()));
    cat+= 10000*(min(ni,electrons.size()));
    cat+= 1000*(min(ni,jets.size()));
    cat+= 100*(min(ni,bjets.size()));
    cat+= 10*(min(ni,type2topjets.size()));
    cat+= 1*(min(ni,type1topjets.size()));

    float_values["Event_category"]=(float)cat;
    string namelabel= "resolvedTopSemiLep";
    //Resolved tops Semileptonic:
    if(doResolvedTopSemiLep){
      sizes[namelabel]=0;
      if(((electrons.size()==1 && muons.size()==0 ) || (muons.size()==1 && electrons.size()==0)) &&  bjets.size()>0 &&   jets.size()>0){
	//	if(jets.size()==2 && bjets.size()==1) cout << " check this one "<<endl;
	TopUtilities topUtils;
	size_t t = 0;
	//	cout << " size b "<< bjets.size()<< " size l  "<< leptons.size() << " size 0 "<< sizes[namelabel]<<endl ;
	for(size_t b =0; b<bjets.size();++b){
	  for(size_t l =0; l<leptons.size();++l){
	    //	  double metPx= 1.0, metPy =1.0;
	    math::PtEtaPhiELorentzVector topSemiLep = topUtils.top4Momentum(leptons.at(l), bjets.at(b),metPx , metPy);
	    if(t > (size_t)max_instances[namelabel])continue;
	    vfloats_values[namelabel+"_Pt"][t]=topSemiLep.pt();
	    vfloats_values[namelabel+"_Eta"][t]=topSemiLep.eta();
	    vfloats_values[namelabel+"_Phi"][t]=topSemiLep.phi();
	    vfloats_values[namelabel+"_E"][t]=topSemiLep.energy();
	    vfloats_values[namelabel+"_Mass"][t]=topSemiLep.mass();
	    vfloats_values[namelabel+"_MT"][t]= topUtils.topMtw(leptons.at(l),bjets.at(b),metPx,metPy);
	    vfloats_values[namelabel+"_LBMPhi"][t]=deltaPhi((leptons.at(l)+bjets.at(b)).Phi(), metphiCorr);
	    vfloats_values[namelabel+"_LMPhi"][t]= deltaPhi(leptons.at(l).Phi(), metphiCorr);
	    vfloats_values[namelabel+"_BMPhi"][t]= deltaPhi(bjets.at(b).Phi(), metphiCorr);
	    vfloats_values[namelabel+"_TMPhi"][t]= deltaPhi(topSemiLep.Phi(), metphiCorr);
	    vfloats_values[namelabel+"_LBPhi"][t]= deltaPhi(leptons.at(l).Phi(), bjets.at(b).Phi());

	    vfloats_values[namelabel+"_IndexB"][t]= mapBJets[b];
	    float lidx = -9999;
	    float flav = -9999;
	    if(mapMu[l]!=-1){lidx=mapMu[l]; flav = 11;}
	    if(mapEle[l]!=-1){lidx=mapEle[l]; flav = 13;}
	    if(mapMu[l]!=-1 && mapEle[l]!=-1) cout<<" what is going on? " <<endl;
	    
	    vfloats_values[namelabel+"_IndexL"][t]= lidx;
	    vfloats_values[namelabel+"_LeptonFlavour"][t]= flav;
	    
	    ++t;
	    ++sizes[namelabel];
	  }
	}
      }
      if((bjets.size()==0 || leptons.size()==0)){
	vector<string> extravars = additionalVariables(namelabel);
	for(size_t t =0;t<(size_t)max_instances[namelabel];++t){
	  for(size_t addv = 0; addv < extravars.size();++addv){
	    vfloats_values[namelabel+"_"+extravars.at(addv)][t]=-9999;
	    //vfloats_values[namelabel+"_"+extravars.at(addv)][t]=0.0;
	  }
	}
      }
      //      cout << "after all loop size is "<<  sizes[namelabel]<<endl;
    }
    //Resolved tops Fullhadronic:
   
    if(doResolvedTopHad || doResolvedTopSemiLep) {
      //      if(getwobjets) cout << " test 1 "<<endl;
      string namelabel= "resolvedTopHad";
      sizes[namelabel]=0;
      
      // ==================== Implementing kinematic fitter ==========================

      if(jets.size()>2 && bjets.size()>0   ){
      
	int maxJetLoop = min((int)(max_instances[jets_label]),max_leading_jets_for_top);
	maxJetLoop = min(maxJetLoop,sizes[jets_label]);
	string pref = obj_to_pref[jets_label];
	size_t t = 0;

	TLorentzVector jet1, jet2, jet3;
	
	
	for(int i = 0; i < maxJetLoop ;++i){
	  for(int j = i+1; j < maxJetLoop ;++j){
	    for(int k = j+1; k < maxJetLoop ;++k){
		      
	      if(vfloats_values[jets_label+"_CorrPt"][i]> 0. && vfloats_values[jets_label+"_CorrPt"][j]> 0. && vfloats_values[jets_label+"_CorrPt"][k]> 0.){
	      
	     
	      //
	      // Kinematic fitter
	      //
		//	      KinematicFitter fitter;
	      
	      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	      // Stuff that goes into event loop
	      //
	      //*********************************
	      //
	      // Example tri-jet combination
	      // The MVA training is done such that the "b-jet" is the jet in the triplet
	      // that has the highest CSVv2+IVF value. I've called this one "jet3" in the example.
	      //

	      float  csv_i= (float)vfloats_values[jets_label+"_CSVv2"][i];
	      float  csv_j= (float)vfloats_values[jets_label+"_CSVv2"][j];
	      float  csv_k= (float)vfloats_values[jets_label+"_CSVv2"][k];
	      
	      vector<float> csv = {csv_i, csv_j, csv_k};
	      vector<int> idx = {i, j, k};
	   
	      ++t;
      
	      size_t tt = 0;
	      //	cout << " test 3 "<<endl;
	      string pref = obj_to_pref[jets_label];
	      
	      bool isIBJet= (bool)vfloats_values[jets_label+"_IsCSVM"][i] && (fabs(vfloats_values[jets_label+"_CorrEta"][i]) < 2.4);
	      bool isITight= (bool)vfloats_values[jets_label+"_IsTight"][i];
	      
	      bool isJBJet= (bool)vfloats_values[jets_label+"_IsCSVM"][j] && (fabs(vfloats_values[jets_label+"_CorrEta"][j]) < 2.4);
	      bool isJTight=(bool)vfloats_values[jets_label+"_IsTight"][j];
		
		//for(int k = j+1; k < maxJetLoop ;++k){
		int nBJets =0;
	
		bool isKBJet= (bool)vfloats_values[jets_label+"_IsCSVM"][k] && (fabs(vfloats_values[jets_label+"_CorrEta"][k]) < 2.4);
		bool isKTight=(bool)vfloats_values[jets_label+"_IsTight"][k];
		
		if(isIBJet)nBJets++;
		if(isJBJet)nBJets++;
		if(isKBJet)nBJets++;
	
		if(nBJets !=1  || !(isITight && isJTight && isKTight) ) continue;
		math::PtEtaPhiELorentzVector p4i = math::PtEtaPhiELorentzVector(vfloats_values[jets_label+"_CorrPt"][i], vfloats_values[jets_label+"_CorrEta"][i], vfloats_values[jets_label+"_CorrPhi"][i], vfloats_values[jets_label+"_CorrE"][i] );
		math::PtEtaPhiELorentzVector p4j = math::PtEtaPhiELorentzVector(vfloats_values[jets_label+"_CorrPt"][j], vfloats_values[jets_label+"_CorrEta"][j], vfloats_values[jets_label+"_CorrPhi"][j], vfloats_values[jets_label+"_CorrE"][j] );
		math::PtEtaPhiELorentzVector p4k = math::PtEtaPhiELorentzVector(vfloats_values[jets_label+"_CorrPt"][k], vfloats_values[jets_label+"_CorrEta"][k], vfloats_values[jets_label+"_CorrPhi"][k], vfloats_values[jets_label+"_CorrE"][k] );
		math::PtEtaPhiELorentzVector topHad = (p4i+p4j)+p4k;
		vfloats_values[namelabel+"_Pt"][tt]=topHad.pt(); 
		vfloats_values[namelabel+"_Eta"][tt]=topHad.eta();
		vfloats_values[namelabel+"_Phi"][tt]=topHad.phi();
		vfloats_values[namelabel+"_E"][tt]=topHad.e();
		vfloats_values[namelabel+"_Mass"][tt]=topHad.mass();
		if(isKBJet){
		  vfloats_values[namelabel+"_IndexB"][tt]= k; vfloats_values[namelabel+"_IndexJ1"][tt]= i;  vfloats_values[namelabel+"_IndexJ2"][tt]= j; 
		}
		if(isJBJet){
		  vfloats_values[namelabel+"_IndexB"][tt]= j; vfloats_values[namelabel+"_IndexJ1"][tt]= i;  vfloats_values[namelabel+"_IndexJ2"][tt]= k; }
		if(isIBJet){
		  //		cout << " isibjet"<<endl;
		  vfloats_values[namelabel+"_IndexB"][tt]= i; vfloats_values[namelabel+"_IndexJ1"][tt]= j;  vfloats_values[namelabel+"_IndexJ2"][tt]= k; }

		math::PtEtaPhiELorentzVector p4w, p4b ;
		float massdrop_=0., deltaRjets=0.;
		if(isIBJet){ p4w = p4j+p4k; p4b= p4i;
		  massdrop_ = max( p4j.mass(),  p4k.mass() )/p4w.mass() ;
		  deltaRjets = deltaR( p4j,  p4k);}
		if(isJBJet){p4w = p4i+p4k; p4b= p4j;
		  massdrop_ = max( p4i.mass(),  p4k.mass() )/p4w.mass() ;
		  deltaRjets = deltaR( p4i,  p4k);}
		if(isKBJet){p4w = p4j+p4i; p4b= p4k;
		  massdrop_ = max( p4i.mass(),  p4j.mass() )/p4w.mass() ;
		  deltaRjets = deltaR( p4i,  p4j);}
		vfloats_values[namelabel+"_Pt"][tt]=topHad.pt(); 
		vfloats_values[namelabel+"_Eta"][tt]=topHad.eta();
		vfloats_values[namelabel+"_Phi"][tt]=topHad.phi();
		vfloats_values[namelabel+"_E"][tt]=topHad.e();
		vfloats_values[namelabel+"_Mass"][tt]=topHad.mass();
		vfloats_values[namelabel+"_WMass"][tt]=p4w.mass();
		vfloats_values[namelabel+"_massDrop"][tt]=massdrop_*deltaRjets;
		//cout<<"Mass drop: "<<massdrop_*deltaRjets<<endl;
		
		if(topHad.mass()<0. || p4w.mass()<0){
		  float genpti = vfloats_values[jets_label+"_GenJetPt"][i];
		  float genptj = vfloats_values[jets_label+"_GenJetPt"][j];
		  float genptk = vfloats_values[jets_label+"_GenJetPt"][k];
		  
		  float pti = vfloats_values[jets_label+"_Pt"][i];
		  float ptj = vfloats_values[jets_label+"_Pt"][j];
		  float ptk = vfloats_values[jets_label+"_Pt"][k];
		  
		  std::cout<<"Top Mass: "<<topHad.mass()<<std::endl;
		  std::cout<<"W Mass: "<<p4w.mass()<<std::endl;
		  std::cout<<"Pt1: "<<p4i.pt()<<" gen "<< genpti<<" pt0 "<<pti <<" eta "<< p4i.eta()<< " phi "<< p4i.phi()<< " e "<< p4i.e()<<std::endl;
		  std::cout<<"Pt2: "<<p4j.pt()<< " gen "<< genptj<< " pt0 "<<ptj <<" eta "<< p4j.eta()<< " phi "<< p4j.phi()<< " e "<< p4j.e()<<std::endl;
		  std::cout<<"Pt3> "<<p4k.pt()<< " gen "<< genptk<<" pt0 "<<ptk<< " eta "<< p4k.eta()<< " phi "<< p4k.phi()<< " e "<< p4k.e()<<std::endl;
		}
		if(t > (size_t)max_instances[namelabel])continue;
		vfloats_values[namelabel+"_WMPhi"][tt]=deltaPhi(p4w.Phi(), metphiCorr);
		vfloats_values[namelabel+"_TMPhi"][tt]= deltaPhi(topHad.Phi(), metphiCorr);
		vfloats_values[namelabel+"_BMPhi"][tt]= deltaPhi(p4b.Phi(), metphiCorr);
		vfloats_values[namelabel+"_WBPhi"][tt]= deltaPhi(p4w.Phi(), p4b.Phi());

		++tt;
		++sizes[namelabel];
		
		//  } //end of if statement on the number of bjets
	      }//end if statement on jets
	    }// end 3rd loop on jets
	  }// end 2st loop on jets
	}// end 1st loop on jets
	
      }// end if statement (at least 3 jets)
      // =============================================================================
      
      //	    }
      //	  }	
      //	}
      //}//end if on number of jets and bjets
      
      //      if(getwobjets)cout << " nresolvedtophad "<< sizes["resolvedTopHad"]<<endl;
      //      cout << " namelabel? "<< namelabel<< endl;
      if(sizes[namelabel]==0){
	vector<string> extravars = additionalVariables(namelabel);
	for(size_t t =0;t<(size_t)max_instances[namelabel];++t){
	  //	  if(t> 45)continue;
	  for(size_t addv = 0; addv < extravars.size();++addv){
	    //	    cout << "resetting variable "<< extravars.at(addv)<< " name "<< namelabel+"_"+extravars.at(addv)<<" t is "<< t << endl;
	    vfloats_values[namelabel+"_"+extravars.at(addv)][t]=-9999;
	  }
	}
      }
    }

    //BTagging part
    if(doBTagSF){
    //CSVT
      //0 tags
      b_weight_csvt_0_tags = b_csvt_0_tags.weight(jsfscsvt, ncsvt_tags);  
      b_weight_csvt_0_tags_mistag_up = b_csvt_0_tags.weight(jsfscsvt_mistag_up, ncsvt_tags);  
      b_weight_csvt_0_tags_mistag_down = b_csvt_0_tags.weight(jsfscsvt_mistag_down, ncsvt_tags);  
      b_weight_csvt_0_tags_b_tag_up = b_csvt_0_tags.weight(jsfscsvt_b_tag_up, ncsvt_tags);  
      b_weight_csvt_0_tags_b_tag_down = b_csvt_0_tags.weight(jsfscsvt_b_tag_down, ncsvt_tags);

      //1
      b_weight_csvt_1_tag = b_csvt_1_tag.weight(jsfscsvt, ncsvt_tags);  
      b_weight_csvt_1_tag_mistag_up = b_csvt_1_tag.weight(jsfscsvt_mistag_up, ncsvt_tags);  
      b_weight_csvt_1_tag_mistag_down = b_csvt_1_tag.weight(jsfscsvt_mistag_down, ncsvt_tags);  
      b_weight_csvt_1_tag_b_tag_up = b_csvt_1_tag.weight(jsfscsvt_b_tag_up, ncsvt_tags);  
      b_weight_csvt_1_tag_b_tag_down = b_csvt_1_tag.weight(jsfscsvt_b_tag_down, ncsvt_tags);

    //CSVM
      //0 tags
      b_weight_csvm_0_tags = b_csvm_0_tags.weight(jsfscsvm, ncsvm_tags);  
      b_weight_csvm_0_tags_mistag_up = b_csvm_0_tags.weight(jsfscsvm_mistag_up, ncsvm_tags);  
      b_weight_csvm_0_tags_mistag_down = b_csvm_0_tags.weight(jsfscsvm_mistag_down, ncsvm_tags);  
      b_weight_csvm_0_tags_b_tag_up = b_csvm_0_tags.weight(jsfscsvm_b_tag_up, ncsvm_tags);  
      b_weight_csvm_0_tags_b_tag_down = b_csvm_0_tags.weight(jsfscsvm_b_tag_down, ncsvm_tags);  
      
      b_weight_subj_csvm_0_tags = b_subj_csvm_0_tags.weight(jsfscsvm_subj, ncsvm_subj_tags);  
      b_weight_subj_csvm_0_tags_mistag_up = b_subj_csvm_0_tags.weight(jsfscsvm_subj_mistag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_0_tags_mistag_down = b_subj_csvm_0_tags.weight(jsfscsvm_subj_mistag_down, ncsvm_subj_tags);  
      b_weight_subj_csvm_0_tags_b_tag_up = b_subj_csvm_0_tags.weight(jsfscsvm_subj_b_tag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_0_tags_b_tag_down = b_subj_csvm_0_tags.weight(jsfscsvm_subj_b_tag_down, ncsvm_subj_tags);
   
      //1 tag
      b_weight_csvm_1_tag = b_csvm_1_tag.weight(jsfscsvm, ncsvm_tags);  
      b_weight_csvm_1_tag_mistag_up = b_csvm_1_tag.weight(jsfscsvm_mistag_up, ncsvm_tags);  
      b_weight_csvm_1_tag_mistag_down = b_csvm_1_tag.weight(jsfscsvm_mistag_down, ncsvm_tags);  
      b_weight_csvm_1_tag_b_tag_up = b_csvm_1_tag.weight(jsfscsvm_b_tag_up, ncsvm_tags);  
      b_weight_csvm_1_tag_b_tag_down = b_csvm_1_tag.weight(jsfscsvm_b_tag_down, ncsvm_tags);  

      b_weight_subj_csvm_1_tag = b_subj_csvm_1_tag.weight(jsfscsvm_subj, ncsvm_subj_tags);  
      b_weight_subj_csvm_1_tag_mistag_up = b_subj_csvm_1_tag.weight(jsfscsvm_subj_mistag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_1_tag_mistag_down = b_subj_csvm_1_tag.weight(jsfscsvm_subj_mistag_down, ncsvm_subj_tags);  
      b_weight_subj_csvm_1_tag_b_tag_up = b_subj_csvm_1_tag.weight(jsfscsvm_subj_b_tag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_1_tag_b_tag_down = b_subj_csvm_1_tag.weight(jsfscsvm_subj_b_tag_down, ncsvm_subj_tags);
      
      //2 tags

      b_weight_subj_csvm_2_tags = b_subj_csvm_2_tags.weight(jsfscsvm_subj, ncsvm_subj_tags);  
      b_weight_subj_csvm_2_tags_mistag_up = b_subj_csvm_2_tags.weight(jsfscsvm_subj_mistag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_2_tags_mistag_down = b_subj_csvm_2_tags.weight(jsfscsvm_subj_mistag_down, ncsvm_subj_tags);  
      b_weight_subj_csvm_2_tags_b_tag_up = b_subj_csvm_2_tags.weight(jsfscsvm_subj_b_tag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_2_tags_b_tag_down = b_subj_csvm_2_tags.weight(jsfscsvm_subj_b_tag_down, ncsvm_subj_tags);
      
      //0-1 tags  

      b_weight_subj_csvm_0_1_tags = b_subj_csvm_0_1_tags.weight(jsfscsvm_subj, ncsvm_subj_tags);  
      b_weight_subj_csvm_0_1_tags_mistag_up = b_subj_csvm_0_1_tags.weight(jsfscsvm_subj_mistag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_0_1_tags_mistag_down = b_subj_csvm_0_1_tags.weight(jsfscsvm_subj_mistag_down, ncsvm_subj_tags);  
      b_weight_subj_csvm_0_1_tags_b_tag_up = b_subj_csvm_0_1_tags.weight(jsfscsvm_subj_b_tag_up, ncsvm_subj_tags);  
      b_weight_subj_csvm_0_1_tags_b_tag_down = b_subj_csvm_0_1_tags.weight(jsfscsvm_subj_b_tag_down, ncsvm_subj_tags);
      
    //CSVL
      //0 tags
      b_weight_csvl_0_tags = b_csvl_0_tags.weight(jsfscsvl, ncsvl_tags);  
      b_weight_csvl_0_tags_mistag_up = b_csvl_0_tags.weight(jsfscsvl_mistag_up, ncsvl_tags);  
      b_weight_csvl_0_tags_mistag_down = b_csvl_0_tags.weight(jsfscsvl_mistag_down, ncsvl_tags);  
      b_weight_csvl_0_tags_b_tag_up = b_csvl_0_tags.weight(jsfscsvl_b_tag_up, ncsvl_tags);  
      b_weight_csvl_0_tags_b_tag_down = b_csvl_0_tags.weight(jsfscsvl_b_tag_down, ncsvl_tags);  

      b_weight_subj_csvl_0_tags = b_subj_csvl_0_tags.weight(jsfscsvl_subj, ncsvl_subj_tags);  
      b_weight_subj_csvl_0_tags_mistag_up = b_subj_csvl_0_tags.weight(jsfscsvl_subj_mistag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_0_tags_mistag_down = b_subj_csvl_0_tags.weight(jsfscsvl_subj_mistag_down, ncsvl_subj_tags);  
      b_weight_subj_csvl_0_tags_b_tag_up = b_subj_csvl_0_tags.weight(jsfscsvl_subj_b_tag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_0_tags_b_tag_down = b_subj_csvl_0_tags.weight(jsfscsvl_subj_b_tag_down, ncsvl_subj_tags);
      
      //1 tag
      b_weight_csvl_1_tag = b_csvl_1_tag.weight(jsfscsvl, ncsvl_tags);  
      b_weight_csvl_1_tag_mistag_up = b_csvl_1_tag.weight(jsfscsvl_mistag_up, ncsvl_tags);  
      b_weight_csvl_1_tag_mistag_down = b_csvl_1_tag.weight(jsfscsvl_mistag_down, ncsvl_tags);  
      b_weight_csvl_1_tag_b_tag_up = b_csvl_1_tag.weight(jsfscsvl_b_tag_up, ncsvl_tags);  
      b_weight_csvl_1_tag_b_tag_down = b_csvl_1_tag.weight(jsfscsvl_b_tag_down, ncsvl_tags);  

      b_weight_subj_csvl_1_tag = b_subj_csvl_1_tag.weight(jsfscsvl_subj, ncsvl_subj_tags);  
      b_weight_subj_csvl_1_tag_mistag_up = b_subj_csvl_1_tag.weight(jsfscsvl_subj_mistag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_1_tag_mistag_down = b_subj_csvl_1_tag.weight(jsfscsvl_subj_mistag_down, ncsvl_subj_tags);  
      b_weight_subj_csvl_1_tag_b_tag_up = b_subj_csvl_1_tag.weight(jsfscsvl_subj_b_tag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_1_tag_b_tag_down = b_subj_csvl_1_tag.weight(jsfscsvl_subj_b_tag_down, ncsvl_subj_tags);

      //2 tags  

      b_weight_subj_csvl_2_tags = b_subj_csvl_2_tags.weight(jsfscsvl_subj, ncsvl_subj_tags);  
      b_weight_subj_csvl_2_tags_mistag_up = b_subj_csvl_2_tags.weight(jsfscsvl_subj_mistag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_2_tags_mistag_down = b_subj_csvl_2_tags.weight(jsfscsvl_subj_mistag_down, ncsvl_subj_tags);  
      b_weight_subj_csvl_2_tags_b_tag_up = b_subj_csvl_2_tags.weight(jsfscsvl_subj_b_tag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_2_tags_b_tag_down = b_subj_csvl_2_tags.weight(jsfscsvl_subj_b_tag_down, ncsvl_subj_tags);
      
      b_weight_csvl_2_tag = b_csvl_2_tag.weight(jsfscsvl, ncsvl_tags);
      b_weight_csvl_2_tag_mistag_up = b_csvl_2_tag.weight(jsfscsvl_mistag_up, ncsvl_tags);
      b_weight_csvl_2_tag_mistag_down = b_csvl_2_tag.weight(jsfscsvl_mistag_down, ncsvl_tags);
      b_weight_csvl_2_tag_b_tag_up = b_csvl_2_tag.weight(jsfscsvl_b_tag_up, ncsvl_tags);
      b_weight_csvl_2_tag_b_tag_down = b_csvl_2_tag.weight(jsfscsvl_b_tag_down, ncsvl_tags);

      //0-1 tags
      
      b_weight_subj_csvl_0_1_tags = b_subj_csvl_0_1_tags.weight(jsfscsvl_subj, ncsvl_subj_tags);  
      b_weight_subj_csvl_0_1_tags_mistag_up = b_subj_csvl_0_1_tags.weight(jsfscsvl_subj_mistag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_0_1_tags_mistag_down = b_subj_csvl_0_1_tags.weight(jsfscsvl_subj_mistag_down, ncsvl_subj_tags);  
      b_weight_subj_csvl_0_1_tags_b_tag_up = b_subj_csvl_0_1_tags.weight(jsfscsvl_subj_b_tag_up, ncsvl_subj_tags);  
      b_weight_subj_csvl_0_1_tags_b_tag_down = b_subj_csvl_0_1_tags.weight(jsfscsvl_subj_b_tag_down, ncsvl_subj_tags);
      
      float_values["Event_bWeight0CSVL_subj"]=b_weight_subj_csvl_0_tags;
      float_values["Event_bWeight1CSVL_subj"]=b_weight_subj_csvl_1_tag;
      float_values["Event_bWeight2CSVL_subj"]=b_weight_subj_csvl_2_tags;
      float_values["Event_bWeight0_1CSVL_subj"]=b_weight_subj_csvl_0_1_tags;
      
      float_values["Event_bWeight0CSVM_subj"]=b_weight_subj_csvm_0_tags;
      float_values["Event_bWeight1CSVM_subj"]=b_weight_subj_csvm_1_tag;
      float_values["Event_bWeight2CSVM_subj"]=b_weight_subj_csvm_2_tags;
      float_values["Event_bWeight0_1CSVM_subj"]=b_weight_subj_csvm_0_1_tags;

      //Mistag
      float_values["Event_bWeightMisTagUp0CSVL_subj"]=b_weight_subj_csvl_0_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1CSVL_subj"]=b_weight_subj_csvl_1_tag_mistag_up;
      float_values["Event_bWeightMisTagUp2CSVL_subj"]=b_weight_subj_csvl_2_tags_mistag_up;
      float_values["Event_bWeightMisTagUp0_1CSVL_subj"]=b_weight_subj_csvl_0_1_tags_mistag_up;
      
      float_values["Event_bWeightMisTagUp0CSVM_subj"]=b_weight_subj_csvm_0_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1CSVM_subj"]=b_weight_subj_csvm_1_tag_mistag_up;
      float_values["Event_bWeightMisTagUp2CSVM_subj"]=b_weight_subj_csvm_2_tags_mistag_up;
      float_values["Event_bWeightMisTagUp0_1CSVM_subj"]=b_weight_subj_csvm_0_1_tags_mistag_up;
    
      float_values["Event_bWeightMisTagDown0CSVL_subj"]=b_weight_subj_csvl_0_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1CSVL_subj"]=b_weight_subj_csvl_1_tag_mistag_down;
      float_values["Event_bWeightMisTagDown2CSVL_subj"]=b_weight_subj_csvl_2_tags_mistag_down;
      float_values["Event_bWeightMisTagDown0_1CSVL_subj"]=b_weight_subj_csvl_0_1_tags_mistag_down;
      
      float_values["Event_bWeightMisTagDown0CSVM_subj"]=b_weight_subj_csvm_0_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1CSVM_subj"]=b_weight_subj_csvm_1_tag_mistag_down;
      float_values["Event_bWeightMisTagDown2CSVM_subj"]=b_weight_subj_csvm_2_tags_mistag_down;
      float_values["Event_bWeightMisTagDown0_1CSVM_subj"]=b_weight_subj_csvm_0_1_tags_mistag_down;
      
      //Btag
      float_values["Event_bWeightBTagUp0CSVL_subj"]=b_weight_subj_csvl_0_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1CSVL_subj"]=b_weight_subj_csvl_1_tag_b_tag_up;
      float_values["Event_bWeightBTagUp2CSVL_subj"]=b_weight_subj_csvl_2_tags_b_tag_up;
      float_values["Event_bWeightBTagUp0_1CSVL_subj"]=b_weight_subj_csvl_0_1_tags_b_tag_up;
      
      float_values["Event_bWeightBTagUp0CSVM_subj"]=b_weight_subj_csvm_0_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1CSVM_subj"]=b_weight_subj_csvm_1_tag_b_tag_up;
      float_values["Event_bWeightBTagUp2CSVM_subj"]=b_weight_subj_csvm_2_tags_b_tag_up;
      float_values["Event_bWeightBTagUp0_1CSVM_subj"]=b_weight_subj_csvm_0_1_tags_b_tag_up;

      float_values["Event_bWeightBTagDown0CSVL_subj"]=b_weight_subj_csvl_0_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1CSVL_subj"]=b_weight_subj_csvl_1_tag_b_tag_down;
      float_values["Event_bWeightBTagDown2CSVL_subj"]=b_weight_subj_csvl_2_tags_b_tag_down;
      float_values["Event_bWeightBTagDown0_1CSVL_subj"]=b_weight_subj_csvl_0_1_tags_b_tag_down;
      
      float_values["Event_bWeightBTagDown0CSVM_subj"]=b_weight_subj_csvm_0_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1CSVM_subj"]=b_weight_subj_csvm_1_tag_b_tag_down;
      float_values["Event_bWeightBTagDown2CSVM_subj"]=b_weight_subj_csvm_2_tags_b_tag_down;
      float_values["Event_bWeightBTagDown0_1CSVM_subj"]=b_weight_subj_csvm_0_1_tags_b_tag_down;
      
      /////AK4
      
      float_values["Event_bWeight0CSVL"]=b_weight_csvl_0_tags;
      float_values["Event_bWeight1CSVL"]=b_weight_csvl_1_tag;
      float_values["Event_bWeight2CSVL"]=b_weight_csvl_2_tag;

      float_values["Event_bWeight0CSVM"]=b_weight_csvm_0_tags;
      float_values["Event_bWeight1CSVM"]=b_weight_csvm_1_tag;
      float_values["Event_bWeight2CSVM"]=b_weight_csvm_2_tag;

      float_values["Event_bWeight0CSVT"]=b_weight_csvt_0_tags;
      float_values["Event_bWeight1CSVT"]=b_weight_csvt_1_tag;
      float_values["Event_bWeight2CSVT"]=b_weight_csvt_2_tag;
      //Mistag
      float_values["Event_bWeightMisTagUp0CSVL"]=b_weight_csvl_0_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1CSVL"]=b_weight_csvl_1_tag_mistag_up;
      float_values["Event_bWeightMisTagUp2CSVL"]=b_weight_csvl_2_tag_mistag_up;

      float_values["Event_bWeightMisTagUp0CSVM"]=b_weight_csvm_0_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1CSVM"]=b_weight_csvm_1_tag_mistag_up;
      float_values["Event_bWeightMisTagUp2CSVM"]=b_weight_csvm_2_tag_mistag_up;

      float_values["Event_bWeightMisTagUp0CSVT"]=b_weight_csvt_0_tags_mistag_up;
      float_values["Event_bWeightMisTagUp1CSVT"]=b_weight_csvt_1_tag_mistag_up;   
      float_values["Event_bWeightMisTagUp2CSVT"]=b_weight_csvt_2_tag_mistag_up;

      float_values["Event_bWeightMisTagDown0CSVL"]=b_weight_csvl_0_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1CSVL"]=b_weight_csvl_1_tag_mistag_down;
      float_values["Event_bWeightMisTagDown2CSVL"]=b_weight_csvl_2_tag_mistag_down;

      float_values["Event_bWeightMisTagDown0CSVM"]=b_weight_csvm_0_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1CSVM"]=b_weight_csvm_1_tag_mistag_down;
      float_values["Event_bWeightMisTagDown2CSVM"]=b_weight_csvm_2_tag_mistag_down;

      float_values["Event_bWeightMisTagDown0CSVT"]=b_weight_csvt_0_tags_mistag_down;
      float_values["Event_bWeightMisTagDown1CSVT"]=b_weight_csvt_1_tag_mistag_down;
      float_values["Event_bWeightMisTagDown2CSVT"]=b_weight_csvt_2_tag_mistag_down;

      //Btag
      float_values["Event_bWeightBTagUp0CSVL"]=b_weight_csvl_0_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1CSVL"]=b_weight_csvl_1_tag_b_tag_up;
      float_values["Event_bWeightBTagUp2CSVL"]=b_weight_csvl_2_tag_b_tag_up;

      float_values["Event_bWeightBTagUp0CSVM"]=b_weight_csvm_0_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1CSVM"]=b_weight_csvm_1_tag_b_tag_up;
      float_values["Event_bWeightBTagUp2CSVM"]=b_weight_csvm_2_tag_b_tag_up;

      float_values["Event_bWeightBTagUp0CSVT"]=b_weight_csvt_0_tags_b_tag_up;
      float_values["Event_bWeightBTagUp1CSVT"]=b_weight_csvt_1_tag_b_tag_up;
      float_values["Event_bWeightBTagUp2CSVT"]=b_weight_csvt_2_tag_b_tag_up;

      float_values["Event_bWeightBTagDown0CSVL"]=b_weight_csvl_0_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1CSVL"]=b_weight_csvl_1_tag_b_tag_down;
      float_values["Event_bWeightBTagDown2CSVL"]=b_weight_csvl_2_tag_b_tag_down;

      float_values["Event_bWeightBTagDown0CSVM"]=b_weight_csvm_0_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1CSVM"]=b_weight_csvm_1_tag_b_tag_down;
      float_values["Event_bWeightBTagDown2CSVM"]=b_weight_csvm_2_tag_b_tag_down;

      float_values["Event_bWeightBTagDown0CSVT"]=b_weight_csvt_0_tags_b_tag_down;
      float_values["Event_bWeightBTagDown1CSVT"]=b_weight_csvt_1_tag_b_tag_down;
      float_values["Event_bWeightBTagDown2CSVT"]=b_weight_csvt_2_tag_b_tag_down;
    }
    
    
    float LHEWeightSign=1.0;
    if(useLHE){
      //LHE and luminosity weights:
      float weightsign = lhes->hepeup().XWGTUP;
      float_values["Event_LHEWeight"]=weightsign;
      LHEWeightSign = weightsign/fabs(weightsign);
      float_values["Event_LHEWeightSign"]=LHEWeightSign;
     }
    float weightLumi = crossSection/originalEvents;
    float_values["Event_weight"]=weightLumi*LHEWeightSign;
    
    //Part 3: filling the additional variables

    if(useLHEWeights){
      getEventLHEWeights();
    }
    if(addLHAPDFWeights){
      getEventPdf();
    }
    
    if(doPU){
      iEvent.getByToken(t_ntrpu_,ntrpu);
      int nTruePV=*ntrpu;
      float_values["Event_nTruePV"]=(float)(nTruePV);
    }

    if(addPV){
      float nGoodPV = 0.0;
      for (size_t v = 0; v < pvZ->size();++v){
	bool isGoodPV = (
			 fabs(pvZ->at(v)) < 24.0 &&
			 pvNdof->at(v) > 4.0 &&
			 pvRho->at(v) <2.0
			 );
	if (isGoodPV)nGoodPV+=1.0;
      }	
      float_values["Event_nGoodPV"]=(float)(nGoodPV);
     float_values["Event_nPV"]=(float)(nPV);
    }

    float_values["Event_passesBadChargedCandidateFilter"] = (float)(*BadChargedCandidateFilter);
    float_values["Event_passesBadPFMuonFilter"] = (float)(*BadPFMuonFilter);
      
    //technical event informationx
    double_values["Event_EventNumber"]=*eventNumber;
    float_values["Event_LumiBlock"]=*lumiBlock;
    float_values["Event_RunNumber"]=*runNumber;
 
    trees[syst]->Fill();
    
    //Reset event weights/#objects
    string nameshortv= "Event";
    vector<string> extravars = additionalVariables(nameshortv);
    for(size_t addv = 0; addv < extravars.size();++addv){
      string name = nameshortv+"_"+extravars.at(addv);

      if(!isMCWeightName(extravars.at(addv))) float_values[name]=0.0;
    }
  }
  for(int t = 0;t < max_instances[boosted_tops_label] ;++t){
    vfloats_values[boosted_tops_label+"_nCSVM"][t]=0;
    vfloats_values[boosted_tops_label+"_nJ"][t]=0;
  }
   
}

bool DMAnalysisTreeMaker::flavourFilter(string ch, int nb, int nc, int nl)
{

  if (ch == "WJets_wbb" || ch == "ZJets_wbb") return (nb > 0 );
  if (ch == "WJets_wcc" || ch == "ZJets_wcc") return (nb == 0 && nc > 0);
  if (ch == "WJets_wlight" || ch == "ZJets_wlight") return (nb == 0 && nc == 0);
  return true;
}



int DMAnalysisTreeMaker::eventFlavour(bool getFlavour, int nb, int nc, int nl)
{
  if (!getFlavour) return 0;
  else
    {
      if ( flavourFilter("WJets_wlight", nb, nc, nl) ) return 1;
      if ( flavourFilter("WJets_wcc", nb, nc, nl) ) return 2;
      if ( flavourFilter("WJets_wbb", nb, nc, nl) ) return 3;
    }
  return 0;
}


string DMAnalysisTreeMaker::makeBranchName(string label, string pref, string var){ //substitutes the "pref" word with "label+'_'"
  string outVar = var;
  size_t prefPos=outVar.find(pref);
  size_t prefLength = pref.length();
  outVar.replace(prefPos,prefLength,label+"_");
  return outVar;
}

string DMAnalysisTreeMaker::makeBranchNameCat(string label, string cat, string pref, string var){
  return makeBranchName (label+cat,pref,var);
}

string DMAnalysisTreeMaker::makeName(string label,string pref,string var){
  return label+"_"+var;
}

void DMAnalysisTreeMaker::initCategoriesSize(string label){
  for(size_t sc = 0; sc< obj_cats[label].size() ;++sc){
    string category = obj_cats[label].at(sc);
    sizes[label+category]=0;
  }
}

void DMAnalysisTreeMaker::setCategorySize(string label, string category, size_t size){
    sizes[label+category]=size;
}
void DMAnalysisTreeMaker::fillCategory(string label, string category, int pos_nocat, int pos_cat){
  //  for(size_t sc = 0; sc< obj_cats[label].size() ;++sc){
  //    string category = obj_cats.at(sc);
  for (size_t obj =0; obj< obj_to_floats[label].size(); ++obj){
    string var = obj_to_floats[label].at(obj);
    string varCat = makeBranchNameCat(label,category,label+"_",var);
    //    cout << " var "<< var << " varcat "<< varCat<<endl;
    vfloats_values[varCat][pos_cat]= vfloats_values[var][pos_nocat];
  }
}


vector<string> DMAnalysisTreeMaker::additionalVariables(string object){
  vector<string> addvar;
  bool ismuon=object.find("muon")!=std::string::npos;
  bool isphoton=object.find("photon")!=std::string::npos;
  bool iselectron=object.find("electron")!=std::string::npos;
  bool ismet=object.find("met")!=std::string::npos;
  bool isjet=object.find("jet")!=std::string::npos && object.find("AK4")!=std::string::npos;
  bool isak8=object.find("jet")!=std::string::npos && object.find("AK8")!=std::string::npos && object.find("sub")==std::string::npos;
  bool isak8subjet=object.find("jet")!=std::string::npos && object.find("AK8")!=std::string::npos && object.find("sub")!=std::string::npos;
  bool isevent=object.find("Event")!=std::string::npos;
  bool isResolvedTopHad=object.find("resolvedTopHad")!=std::string::npos;
  bool isResolvedTopSemiLep=object.find("resolvedTopSemiLep")!=std::string::npos;
  //  bool isgen = object.find("genPart")!=std::string::npos;

  if(ismuon || iselectron){
    addvar.push_back("SFTrigger");
  }
  
  if(isphoton){
    addvar.push_back("isLooseSpring15");
    addvar.push_back("isMediumSpring15");
    addvar.push_back("isTightSpring15");
  }
  if(iselectron){
    addvar.push_back("PassesDRmu");
  }
  if(ismet){
    addvar.push_back("CorrPt");
    addvar.push_back("CorrPhi");
    addvar.push_back("CorrBasePt");
    addvar.push_back("CorrBasePhi");
    addvar.push_back("CorrT1Pt");
    addvar.push_back("CorrT1Phi");
    addvar.push_back("CorrPtNoHF");
    addvar.push_back("CorrPhiNoHF");
  }
  if(isjet){
    addvar.push_back("CorrPt");
    addvar.push_back("CorrE");
    addvar.push_back("NoCorrPt");
    addvar.push_back("NoCorrE");
    addvar.push_back("MinDR");
    addvar.push_back("IsCSVT");
    addvar.push_back("IsCSVM");
    addvar.push_back("IsCSVL");
    addvar.push_back("PassesID");
    addvar.push_back("PassesDR");
    addvar.push_back("CorrMass");
    addvar.push_back("IsTight");
    addvar.push_back("IsLoose");
  }
  if(isak8){
    addvar.push_back("CorrPrunedMassCHS");
    addvar.push_back("CorrPt");
    addvar.push_back("CorrE");
    addvar.push_back("NoCorrPt");
    addvar.push_back("NoCorrE");
    addvar.push_back("isType1");
    addvar.push_back("isType2");
    addvar.push_back("TopPt");
    addvar.push_back("TopEta");
    addvar.push_back("TopPhi");
    addvar.push_back("TopE");
    addvar.push_back("TopMass");
    addvar.push_back("TopWMass");
    addvar.push_back("nJ");
    addvar.push_back("nCSVM");
    addvar.push_back("CSVv2");
    addvar.push_back("nCSVsubj");
    addvar.push_back("nCSVsubj_tm");
    addvar.push_back("tau3OVERtau2");
    addvar.push_back("tau2OVERtau1");

  }  
  if(1>0){
    addvar.push_back("Pt");
    addvar.push_back("Eta");
    addvar.push_back("Phi");
    addvar.push_back("E");
    addvar.push_back("Id");
    addvar.push_back("Status");
    addvar.push_back("MomId");
    addvar.push_back("Mom0Status");
    /*addvar.push_back("dauId1");
    addvar.push_back("dauStatus1");
    addvar.push_back("dauId2");
    addvar.push_back("dauStatus2");
    addvar.push_back("dauId3");
    addvar.push_back("dauStatus3");*/
  }
  if(isak8subjet){
    ;//    addvar.push_back("CorrPt");
  }

  if(isResolvedTopHad ){
    addvar.push_back("Pt");    addvar.push_back("Eta");    addvar.push_back("Phi");    addvar.push_back("E"); addvar.push_back("Mass");  addvar.push_back("massDrop");
    addvar.push_back("WMass"); addvar.push_back("BMPhi");  addvar.push_back("WMPhi");  addvar.push_back("TMPhi");  addvar.push_back("WBPhi");
    addvar.push_back("IndexB");    addvar.push_back("IndexJ1");    addvar.push_back("IndexJ2");  addvar.push_back("IndexB_MVA");    addvar.push_back("IndexJ1_MVA");    addvar.push_back("IndexJ2_MVA");  addvar.push_back("MVA");  addvar.push_back("WMassPreFit"); addvar.push_back("MassPreFit"); addvar.push_back("PtPreFit"); addvar.push_back("EtaPreFit"); addvar.push_back("PhiPreFit"); addvar.push_back("BMPhiPreFit");  addvar.push_back("WMPhiPreFit");  addvar.push_back("TMPhiPreFit");  addvar.push_back("WBPhiPreFit"); addvar.push_back("WMassPostFit"); addvar.push_back("MassPostFit"); addvar.push_back("PtPostFit"); addvar.push_back("EtaPostFit"); addvar.push_back("PhiPostFit"); addvar.push_back("BMPhiPostFit");  addvar.push_back("WMPhiPostFit");  addvar.push_back("TMPhiPostFit");  addvar.push_back("WBPhiPostFit");  addvar.push_back("FitProb");  addvar.push_back("DPhiJet1b"); addvar.push_back("DPhiJet2b"); addvar.push_back("DRJet1b"); addvar.push_back("DRJet2b"); 

  }

  if(isResolvedTopSemiLep ){
    addvar.push_back("Pt");    addvar.push_back("Eta");    addvar.push_back("Phi");    addvar.push_back("E"); addvar.push_back("Mass");   
    addvar.push_back("MT");    addvar.push_back("LBMPhi");    addvar.push_back("LMPhi");    addvar.push_back("LBPhi");     addvar.push_back("BMPhi");  addvar.push_back("TMPhi"); 
    addvar.push_back("IndexL");    addvar.push_back("LeptonFlavour");    addvar.push_back("IndexB");
  }
  
  if(isevent){
    addvar.push_back("weight");
    addvar.push_back("nTightMuons");
    addvar.push_back("nSoftMuons");
    addvar.push_back("nLooseMuons");
    addvar.push_back("nMediumMuons");    
    addvar.push_back("nTightElectrons");
    addvar.push_back("nMediumElectrons");
    addvar.push_back("nLooseElectrons");
    addvar.push_back("nVetoElectrons");
    addvar.push_back("nElectronsSF");
    addvar.push_back("mt");
    addvar.push_back("Mt2w");
    addvar.push_back("category");
    addvar.push_back("nMuonsSF");
    addvar.push_back("nCSVTJets");
    addvar.push_back("nCSVMJets");
    addvar.push_back("nCSVLJets");
    addvar.push_back("nTightJets");
    addvar.push_back("nLooseJets");
    addvar.push_back("nType1TopJets");
    addvar.push_back("nType2TopJets");
    addvar.push_back("Ht");
    addvar.push_back("nGoodPV");
    addvar.push_back("nPV");
    addvar.push_back("nTruePV");
    addvar.push_back("Rho");

    addvar.push_back("bWeight0CSVT");
    addvar.push_back("bWeight1CSVT");
    addvar.push_back("bWeight2CSVT");

    addvar.push_back("bWeight0CSVM");
    addvar.push_back("bWeight1CSVM");
    addvar.push_back("bWeight2CSVM");

    addvar.push_back("bWeight0CSVL");
    addvar.push_back("bWeight1CSVL");
    addvar.push_back("bWeight2CSVL");

    addvar.push_back("bWeightMisTagDown0CSVT");
    addvar.push_back("bWeightMisTagDown1CSVT");
    addvar.push_back("bWeightMisTagDown2CSVT");

    addvar.push_back("bWeightMisTagDown0CSVM");
    addvar.push_back("bWeightMisTagDown1CSVM");
    addvar.push_back("bWeightMisTagDown2CSVM");

    addvar.push_back("bWeightMisTagDown0CSVL");
    addvar.push_back("bWeightMisTagDown1CSVL");
    addvar.push_back("bWeightMisTagDown2CSVL");

    addvar.push_back("bWeightMisTagUp0CSVT");
    addvar.push_back("bWeightMisTagUp1CSVT");
    addvar.push_back("bWeightMisTagUp2CSVT");

    addvar.push_back("bWeightMisTagUp0CSVM");
    addvar.push_back("bWeightMisTagUp1CSVM");
    addvar.push_back("bWeightMisTagUp2CSVM");

    addvar.push_back("bWeightMisTagUp0CSVL");
    addvar.push_back("bWeightMisTagUp1CSVL");
    addvar.push_back("bWeightMisTagUp2CSVL");

    addvar.push_back("bWeightBTagUp0CSVT");
    addvar.push_back("bWeightBTagUp1CSVT");
    addvar.push_back("bWeightBTagUp2CSVT");
    
    addvar.push_back("bWeightBTagUp0CSVM");
    addvar.push_back("bWeightBTagUp1CSVM");
    addvar.push_back("bWeightBTagUp2CSVM");

    addvar.push_back("bWeightBTagUp0CSVL");
    addvar.push_back("bWeightBTagUp1CSVL");
    addvar.push_back("bWeightBTagUp2CSVL");

    addvar.push_back("bWeightBTagDown0CSVT");
    addvar.push_back("bWeightBTagDown1CSVT");
    addvar.push_back("bWeightBTagDown2CSVT");

    addvar.push_back("bWeightBTagDown0CSVM");
    addvar.push_back("bWeightBTagDown1CSVM");
    addvar.push_back("bWeightBTagDown2CSVM");

    addvar.push_back("bWeightBTagDown0CSVL");
    addvar.push_back("bWeightBTagDown1CSVL");
    addvar.push_back("bWeightBTagDown2CSVL");

    //subjets

    addvar.push_back("bWeight0CSVM_subj");
    addvar.push_back("bWeight1CSVM_subj");
    addvar.push_back("bWeight2CSVM_subj");
    addvar.push_back("bWeight0_1CSVM_subj");

    addvar.push_back("bWeight0CSVL_subj");
    addvar.push_back("bWeight1CSVL_subj");
    addvar.push_back("bWeight2CSVL_subj");
    addvar.push_back("bWeight0_1CSVL_subj");

    addvar.push_back("bWeightMisTagDown0CSVM_subj");
    addvar.push_back("bWeightMisTagDown1CSVM_subj");
    addvar.push_back("bWeightMisTagDown2CSVM_subj");
    addvar.push_back("bWeightMisTagDown0_1CSVM_subj");

    addvar.push_back("bWeightMisTagDown0CSVL_subj");
    addvar.push_back("bWeightMisTagDown1CSVL_subj");
    addvar.push_back("bWeightMisTagDown2CSVL_subj");
    addvar.push_back("bWeightMisTagDown0_1CSVL_subj");

    addvar.push_back("bWeightMisTagUp0CSVM_subj");
    addvar.push_back("bWeightMisTagUp1CSVM_subj");
    addvar.push_back("bWeightMisTagUp2CSVM_subj");
    addvar.push_back("bWeightMisTagUp0_1CSVM_subj");

    addvar.push_back("bWeightMisTagUp0CSVL_subj");
    addvar.push_back("bWeightMisTagUp1CSVL_subj");
    addvar.push_back("bWeightMisTagUp2CSVL_subj");
    addvar.push_back("bWeightMisTagUp0_1CSVL_subj");

    addvar.push_back("bWeightBTagUp0CSVM_subj");
    addvar.push_back("bWeightBTagUp1CSVM_subj");
    addvar.push_back("bWeightBTagUp2CSVM_subj");
    addvar.push_back("bWeightBTagUp0_1CSVM_subj");

    addvar.push_back("bWeightBTagUp0CSVL_subj");
    addvar.push_back("bWeightBTagUp1CSVL_subj");
    addvar.push_back("bWeightBTagUp2CSVL_subj");
    addvar.push_back("bWeightBTagUp0_1CSVL_subj");

    addvar.push_back("bWeightBTagDown0CSVM_subj");
    addvar.push_back("bWeightBTagDown1CSVM_subj");
    addvar.push_back("bWeightBTagDown2CSVM_subj");
    addvar.push_back("bWeightBTagDown0_1CSVM_subj");

    addvar.push_back("bWeightBTagDown0CSVL_subj");
    addvar.push_back("bWeightBTagDown1CSVL_subj");
    addvar.push_back("bWeightBTagDown2CSVL_subj");
    addvar.push_back("bWeightBTagDown0_1CSVL_subj");

    addvar.push_back("T_Pt");
    addvar.push_back("T_Eta");
    addvar.push_back("T_Phi");
    addvar.push_back("T_E");

    addvar.push_back("Tbar_Pt");
    addvar.push_back("Tbar_Eta");
    addvar.push_back("Tbar_Phi");
    addvar.push_back("Tbar_E");

    addvar.push_back("W_Pt");
    addvar.push_back("W_Eta");
    addvar.push_back("W_Phi");
    addvar.push_back("W_E");

    addvar.push_back("Z_Pt");
    addvar.push_back("Z_Eta");
    addvar.push_back("Z_Phi");
    addvar.push_back("Z_E");

    addvar.push_back("Z_QCD_Weight");
    addvar.push_back("W_QCD_Weight");
    
    addvar.push_back("Z_Weight");
    addvar.push_back("W_Weight");

    addvar.push_back("Z_EW_Weight");
    addvar.push_back("W_EW_Weight");

    addvar.push_back("T_Weight");
    addvar.push_back("T_Ext_Weight");
    addvar.push_back("eventFlavour");

    addvar.push_back("Lepton1_Pt");
    addvar.push_back("Lepton1_Eta");
    addvar.push_back("Lepton1_Phi");
    addvar.push_back("Lepton1_E");
    addvar.push_back("Lepton1_Charge");
    addvar.push_back("Lepton1_Flavour");

    addvar.push_back("Lepton2_Pt");
    addvar.push_back("Lepton2_Eta");
    addvar.push_back("Lepton2_Phi");
    addvar.push_back("Lepton2_E");
    addvar.push_back("Lepton2_Charge");
    addvar.push_back("Lepton2_Flavour");
    
    addvar.push_back("LHEWeightSign");
    //addvar.push_back("LHEWeightAVG");
    addvar.push_back("LHEWeight");
    addvar.push_back("EventNumber");
    addvar.push_back("LumiBlock");
    addvar.push_back("RunNumber");
    
    for (size_t j = 0; j < (size_t)jetScanCuts.size(); ++j){
      stringstream j_n;
      double jetval = jetScanCuts.at(j);
      j_n << "Cut" <<jetval;
      
      addvar.push_back("nJets"+j_n.str());
      addvar.push_back("nCSVTJets"+j_n.str());
      addvar.push_back("nCSVMJets"+j_n.str());
      addvar.push_back("nCSVLJets"+j_n.str());
      addvar.push_back("bWeight1CSVTWeight"+j_n.str());
      addvar.push_back("bWeight2CSVTWeight"+j_n.str());
      addvar.push_back("bWeight1CSVMWeight"+j_n.str());
      addvar.push_back("bWeight2CSVMWeight"+j_n.str());
      addvar.push_back("bWeight1CSVLWeight"+j_n.str());
      addvar.push_back("bWeight2CSVLWeight"+j_n.str());
      addvar.push_back("bWeight0CSVLWeight"+j_n.str());
      addvar.push_back("bWeight0CSVLWeight"+j_n.str());
    }

    if(useLHEWeights){
      for (size_t w = 0; w <= (size_t)maxWeights; ++w)  {
	stringstream w_n;
	w_n << w;
	addvar.push_back("LHEWeight"+w_n.str());
      }
    }
    if(addLHAPDFWeights){
      for (size_t p = 1; p <= (size_t)maxPdf; ++p)  {
	//cout << " pdf # " << pdf_weights[p - 1] << " test "<< endl; 
	stringstream w_n;
	w_n << p;
	addvar.push_back("PDFWeight" + w_n.str());
      }
    }
    if(useMETFilters){
      for (size_t lt = 0; lt < metFilters.size(); ++lt)  {
	string trig = metFilters.at(lt);
	addvar.push_back("passes"+trig);
      }
      addvar.push_back("passesMETFilters");
      addvar.push_back("passesBadChargedCandidateFilter");
      addvar.push_back("passesBadPFMuonFilter");
    }
    if(useTriggers){
      for (size_t lt = 0; lt < SingleElTriggers.size(); ++lt)  {
	string trig = SingleElTriggers.at(lt);
	addvar.push_back("passes"+trig);
	addvar.push_back("prescale"+trig);
      }
      for (size_t lt = 0; lt < SingleMuTriggers.size(); ++lt)  {
	string trig = SingleMuTriggers.at(lt);
	addvar.push_back("passes"+trig);
	addvar.push_back("prescale"+trig);
      }
      for (size_t lt = 0; lt < PhotonTriggers.size(); ++lt)  {
	string trig = PhotonTriggers.at(lt);
	addvar.push_back("passes"+trig);
	addvar.push_back("prescale"+trig);
      }
      for (size_t ht = 0; ht < hadronicTriggers.size(); ++ht)  {
	string trig = hadronicTriggers.at(ht);
	addvar.push_back("passes"+trig);
	addvar.push_back("prescale"+trig);
      }
      for (size_t lt = 0; lt < HadronPFHT900Triggers.size(); ++lt)  {
        string trig = HadronPFHT900Triggers.at(lt);
        addvar.push_back("passes"+trig);
        addvar.push_back("prescale"+trig);
      }
      for (size_t lt = 0; lt < HadronPFHT800Triggers.size(); ++lt)  {
        string trig = HadronPFHT800Triggers.at(lt);
        addvar.push_back("passes"+trig);
        addvar.push_back("prescale"+trig);
      }
      for (size_t lt = 0; lt < HadronPFJet450Triggers.size(); ++lt)  {
        string trig = HadronPFJet450Triggers.at(lt);
        addvar.push_back("passes"+trig);
        addvar.push_back("prescale"+trig);
      }
      
      addvar.push_back("passesSingleElTriggers");
      addvar.push_back("passesSingleMuTriggers");
      addvar.push_back("passesPhotonTriggers");
      addvar.push_back("passesHadronicTriggers");
      addvar.push_back("passesHadronPFHT900Triggers");
      addvar.push_back("passesHadronPFHT800Triggers");
      addvar.push_back("passesHadronPFJet450Triggers");
    }
  }
  
  return addvar;
}

void DMAnalysisTreeMaker::initializePdf(string central, string varied){

    if(central == "CT") {  LHAPDF::initPDFSet(1, "cteq66.LHgrid"); }
    if(central == "CT10") {  LHAPDF::initPDFSet(1, "CT10.LHgrid"); }
    if(central == "CT10f4") {  LHAPDF::initPDFSet(1, "CT10f4.LHgrid"); }
    if(central == "NNPDF") { LHAPDF::initPDFSet(1, "NNPDF21_100.LHgrid");  }
    if(central == "MSTW") { LHAPDF::initPDFSet(1, "MSTW2008nlo68cl.LHgrid");  }

    if(varied == "CT") {  LHAPDF::initPDFSet(2, "cteq66.LHgrid"); maxPdf = 44; }
    if(varied == "CT10") {  LHAPDF::initPDFSet(2, "CT10.LHgrid"); maxPdf = 52; }
    if(varied == "CT10f4") {  LHAPDF::initPDFSet(2, "CT10f4.LHgrid"); maxPdf = 52; }
    if(varied == "NNPDF") { LHAPDF::initPDFSet(2, "NNPDF21_100.LHgrid");  maxPdf = 100; }
    if(varied == "MSTW") { LHAPDF::initPDFSet(2, "MSTW2008nlo68cl.LHgrid"); maxPdf = 40; }

}

double DMAnalysisTreeMaker::getWEWKPtWeight(double ptW){ 

  if(ptW<150.)return 0.980859;
  if(ptW>=150. && ptW <200.)return 0.962119;
  if(ptW>=200. && ptW <250.)return 0.944429;
  if(ptW>=250. && ptW <300.)return 0.927686;
  if(ptW>=300. && ptW <350.)return 0.911802;
  if(ptW>=350. && ptW <400.)return 0.8967;
  if(ptW>=400. && ptW <500.)return 0.875368;
  if(ptW>=500. && ptW <600.)return 0.849097;
  if(ptW>=600. && ptW <1000)return 0.792159;
  return 1.0;
}

double DMAnalysisTreeMaker::getZEWKPtWeight(double ptW){ 
  if(ptW<150.)return 0.984525;
  if(ptW>=150. && ptW <200.)return 0.969079;
  if(ptW>=200. && ptW <250.)return 0.954627;
  if(ptW>=250. && ptW <300.)return 0.941059;
  if(ptW>=300. && ptW <350.)return 0.928284;
  if(ptW>=350. && ptW <400.)return 0.91622;
  if(ptW>=400. && ptW <500.)return 0.899312;
  if(ptW>=500. && ptW <600.)return 0.878693;
  if(ptW>=600. && ptW <1000)return 0.834718;
  return 1.0;

}

double DMAnalysisTreeMaker::getWPtWeight(double ptW){ 
  //QCD
  
  if(ptW<150.)return 1.89123;
  if(ptW>=150. && ptW <200.)return 1.70414;
  if(ptW>=200. && ptW <250.)return 1.60726;
  if(ptW>=250. && ptW <300.)return 1.57206;
  if(ptW>=300. && ptW <350.)return 1.51689;
  if(ptW>=350. && ptW <400.)return 1.4109;
  if(ptW>=400. && ptW <500.)return 1.30758;
  if(ptW>=500. && ptW <600.)return 1.32046;
  if(ptW>=600. && ptW <1000)return 1.26853;
  return 1.0;
}

double DMAnalysisTreeMaker::getAPtWeight(double ptW){ 
  if(ptW<150.)return 1.24087;
  if(ptW>=150. && ptW <200.)return 1.55807;
  if(ptW>=200. && ptW <250.)return 1.51043;
  if(ptW>=250. && ptW <300.)return 1.47333;
  if(ptW>=300. && ptW <350.)return 1.43497;
  if(ptW>=350. && ptW <400.)return 1.37846;
  if(ptW>=400. && ptW <500.)return 1.29202;
  if(ptW>=500. && ptW <600.)return 1.31414;
  if(ptW>=600.)return 1.20454;
  return 1.0;
}

double DMAnalysisTreeMaker::getZPtWeight(double ptW){
  if(ptW<150.)return 1.685005;
  if(ptW>=150. && ptW <200.)return 1.552560;
  if(ptW>=200. && ptW <250.)return 1.522595;
  if(ptW>=250. && ptW <300.)return 1.520624;
  if(ptW>=300. && ptW <350.)return 1.432282;
  if(ptW>=350. && ptW <400.)return 1.457417;
  if(ptW>=400. && ptW <500.)return 1.368499;
  if(ptW>=500. && ptW <600.)return 1.358024;
  if(ptW>=600.)return 1.164847;;
  return 1.0;
}

double DMAnalysisTreeMaker::getTopPtWeight(double ptT, double ptTbar, bool extrap){ 
  if((ptT>0.0 && ptTbar>0.0) ){
    if (extrap || (ptT<=400.0 && ptTbar <=400.0)){
      double a = 0.0615;
      double b = -0.0005;
      double sfT = exp(a+b*ptT);
      double sfTbar = exp(a+b*ptTbar);
    return sqrt(sfT*sfTbar); 
    }
  }
  return 1.0;
}

bool DMAnalysisTreeMaker::getMETFilters(){
  bool METFilterAND=true;
  for(size_t mf =0; mf< metFilters.size();++mf){
    string fname = metFilters.at(mf);
    for(size_t bt = 0; bt < metNames->size();++bt){
      std::string tname = metNames->at(bt);
      //      cout << "test tname "<<endl;
      if(tname.find(fname)!=std::string::npos){
	METFilterAND = METFilterAND && (metBits->at(bt)>0);
	float_values["Event_passes"+fname]=metBits->at(bt);
      }
    }
  }
  float_values["Event_passesMETFilters"]=(float)METFilterAND;
  return METFilterAND;
}

bool DMAnalysisTreeMaker::getEventTriggers(){
  bool eleOR=false, muOR=false, hadronOR=false, phOR=false, HHThOR=false, HHTlOR=false, HJetOR=false;
  for(size_t lt =0; lt< SingleElTriggers.size();++lt){
    string lname = SingleElTriggers.at(lt);
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      if(tname.find(lname)!=std::string::npos){
	eleOR = muOR || (triggerBits->at(bt)>0);
	float_values["Event_passes"+lname]=triggerBits->at(bt);
	float_values["Event_prescale"+lname]=triggerPrescales->at(bt);
      }
    }
  }

  for(size_t bt = 0; bt < triggerNamesR->size();++bt){
    std::string tname = triggerNamesR->at(bt);
    //    std::cout << " tname is " << tname << " passes "<< triggerBits->at(bt)<< std::endl;
  }
  
  for(size_t lt =0; lt< SingleMuTriggers.size();++lt){
    string lname = SingleMuTriggers.at(lt);
    //    std::cout << lname << std::endl;
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      //      std::cout << " tname is " << tname << " passes "<< triggerBits->at(bt)<< std::endl;
      if(tname.find(lname)!=std::string::npos){
	//	cout << " matches "<<endl;
	muOR = muOR || (triggerBits->at(bt)>0);
	float_values["Event_passes"+lname]=triggerBits->at(bt);
	float_values["Event_prescale"+lname]=triggerPrescales->at(bt);
      }
    }
  }
  for(size_t pt =0; pt< PhotonTriggers.size();++pt){
    string pname = PhotonTriggers.at(pt);
    //std::cout << pname << std::endl;
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      if(tname.find(pname)!=std::string::npos){
	phOR = phOR || (triggerBits->at(bt)>0);
	float_values["Event_passes"+pname]=triggerBits->at(bt);
	float_values["Event_prescale"+pname]=triggerPrescales->at(bt);
      }
    }
  }

  for(size_t ht =0; ht<hadronicTriggers.size();++ht){
    string hname = hadronicTriggers.at(ht);
    //std::cout << hname << std::endl;
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      if(tname.find(hname)!=std::string::npos){
	//bool before = hadronOR;
	hadronOR = hadronOR || (triggerBits->at(bt)>0);
	//bool after = hadronOR;
	//if(before != after){
	//std::cout<< "hadr name "<< hname << std::endl;
	//std::cout<< "trig name "<< tname << std::endl;
	//std::cout << "hadronOR before " << before << std::endl;
	//std::cout << "hadronOR after " << after << std::endl;
	//}
	//hadronOR = hadronOR || (triggerBits->at(bt)>0);
	float_values["Event_passes"+hname]=triggerBits->at(bt);
	float_values["Event_prescale"+hname]=triggerPrescales->at(bt);
      }
    }
  }

  for(size_t lt =0; lt< HadronPFHT900Triggers.size();++lt){
    string lname = HadronPFHT900Triggers.at(lt);
    //std::cout << "==>lname: " <<  lname << std::endl;                                                                                                                                                
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      //      std::cout << " tname is " << tname << " passes "<< triggerBits->at(bt)<< std::endl;                                                                                         
      if(tname.find(lname)!=std::string::npos){
        //      cout << " matches "<<endl;                                                               
	//std::cout << "==>trigger bit: " << triggerBits->at(bt) << " HHThOR: " << HHThOR << std::endl;
        HHThOR = HHThOR || (triggerBits->at(bt)>0);
        float_values["Event_passes"+lname]=triggerBits->at(bt);
        float_values["Event_prescale"+lname]=triggerPrescales->at(bt);
      }
    }
  }
  for(size_t lt =0; lt< HadronPFHT800Triggers.size();++lt){
    string lname = HadronPFHT800Triggers.at(lt);
    //    std::cout << lname << std::endl;                                                                                                                                         
                                                                                                                                                                                    
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      //      std::cout << " tname is " << tname << " passes "<< triggerBits->at(bt)<< std::endl;                                                                                   
                                                                                                                                                                                       
      if(tname.find(lname)!=std::string::npos){
	//      cout << " matches "<<endl;                                                                                                                                        
                                                                                                                                                                                   
	HHTlOR = HHTlOR || (triggerBits->at(bt)>0);
        float_values["Event_passes"+lname]=triggerBits->at(bt);
	float_values["Event_prescale"+lname]=triggerPrescales->at(bt);
      }
    }
  }

  for(size_t lt =0; lt< HadronPFJet450Triggers.size();++lt){
    string lname = HadronPFJet450Triggers.at(lt);
    //    std::cout << lname << std::endl;                                                                                                                                            
                                                                                                                                                                                     
    for(size_t bt = 0; bt < triggerNamesR->size();++bt){
      std::string tname = triggerNamesR->at(bt);
      //      std::cout << " tname is " << tname << " passes "<< triggerBits->at(bt)<< std::endl;                                                                                     
                                                                                                                                                                                      
      if(tname.find(lname)!=std::string::npos){
	//      cout << " matches "<<endl;                                                                                                                                                                                                                                                                                                                                     
	HJetOR = HJetOR || (triggerBits->at(bt)>0);
        float_values["Event_passes"+lname]=triggerBits->at(bt);
	float_values["Event_prescale"+lname]=triggerPrescales->at(bt);
      }
    }
  }

  float_values["Event_passesSingleElTriggers"]=(float)eleOR;
  float_values["Event_passesSingleMuTriggers"]=(float)muOR;
  float_values["Event_passesPhotonTriggers"]=(float)phOR;
  float_values["Event_passesHadronicTriggers"]=(float)hadronOR;
  //std::cout << "............" << HHThOR << ".............." << std::endl;
  float_values["Event_passesHadronPFHT900Triggers"]=(float)HHThOR;
  float_values["Event_passesHadronPFHT800Triggers"]=(float)HHTlOR;
  float_values["Event_passesHadronPFJet450Triggers"]=(float)HJetOR;
  return (eleOR || muOR || hadronOR || phOR || HHThOR || HHTlOR || HJetOR);
}


void DMAnalysisTreeMaker::getEventPdf(){

  //  std::cout << " getting pdf "<<endl;

  double scalePDF = genprod->pdf()->scalePDF;
  double x1 =  genprod->pdf()->x.first;
  double x2 =  genprod->pdf()->x.second;
  int id1 =  genprod->pdf()->id.first;
  int id2 =  genprod->pdf()->id.second;

  //  std::cout << " maxpdf "<< maxPdf << " accessing x1 " << x1<< id1<<std::endl;


  LHAPDF::usePDFMember(1, 0);
  double xpdf1 = LHAPDF::xfx(1, x1, scalePDF, id1);
  double xpdf2 = LHAPDF::xfx(1, x2, scalePDF, id2);
  double w0 = xpdf1 * xpdf2;
  int maxPDFCount = maxPdf;

  std::cout << "entering pdf loop" <<std::endl;
  for (int p = 1; p <= maxPdf; ++p)
    {
      
      if ( p > maxPDFCount ) continue;
      LHAPDF::usePDFMember(2, p);
      double xpdf1_new = LHAPDF::xfx(2, x1, scalePDF, id1);
      double xpdf2_new = LHAPDF::xfx(2, x2, scalePDF, id2);
      double pweight = xpdf1_new * xpdf2_new / w0;
      stringstream w_n;
      w_n << p;
      float_values["PDFWeight"+w_n.str()]= pweight;
    }
  
}


void DMAnalysisTreeMaker::getEventLHEWeights(){
  //  std::cout << " in weight "<<endl;
  size_t wgtsize=  lhes->weights().size();
  //  std::cout << "weight size "<< wgtsize<<endl;
  for (size_t i = 0; i <  wgtsize; ++i)  {
    if (i<= (size_t)maxWeights){ 
      stringstream w_n;
      w_n << i;

      float ww = (float)lhes->weights().at(i).wgt;
      
      //      cout << "ww # " << i<< "is "<<ww <<endl;
      //      cout << "id  is "<< std::string(lhes->weights().at(i).id.data()) <<endl;
      //      cout <<" floatval before "<< float_values["Event_LHEWeight"+w_n.str()]<<endl;

      float_values["Event_LHEWeight"+w_n.str()]= ww;
      //if(i>=11)float_values["Event_LHEWeightAVG"]+= ww;

      //      cout <<" floatval after "<< float_values["Event_LHEWeight"+w_n.str()]<<endl;

    }
    //    float_values["Event_LHEWeightAVG"]+= ww;

    //    else cout << "WARNING! there are " << wgtsize << " weights, and you accomodated for only "<< maxWeights << " weights, check your configuration file/your lhe!!!"<<endl;
  }
  
}

void DMAnalysisTreeMaker::initTreeWeightHistory(bool useLHEW){
  cout << " preBranch "<<endl;

  trees["WeightHistory"]->Branch("Event_Z_EW_Weight",&float_values["Event_Z_EW_Weight"]);
  trees["WeightHistory"]->Branch("Event_W_EW_Weight",&float_values["Event_W_EW_Weight"]);
  trees["WeightHistory"]->Branch("Event_Z_QCD_Weight",&float_values["Event_Z_QCD_Weight"]);
  trees["WeightHistory"]->Branch("Event_W_QCD_Weight",&float_values["Event_W_QCD_Weight"]);
  trees["WeightHistory"]->Branch("Event_Z_Weight",&float_values["Event_Z_Weight"]);
  trees["WeightHistory"]->Branch("Event_W_Weight",&float_values["Event_W_Weight"]);

  trees["WeightHistory"]->Branch("Event_T_Weight",&float_values["Event_T_Weight"]);
  trees["WeightHistory"]->Branch("Event_T_Ext_Weight",&float_values["Event_T_Ext_Weight"]);
  trees["WeightHistory"]->Branch("Event_T_Pt",&float_values["Event_T_Pt"]);
  trees["WeightHistory"]->Branch("Event_Tbar_Pt",&float_values["Event_Tbar_Pt"]);

  trees["WeightHistory"]->Branch("Event_W_Pt",&float_values["Event_W_Pt"]);
  trees["WeightHistory"]->Branch("Event_Z_Pt",&float_values["Event_Z_Pt"]);
  
  cout << " preBranch Weight "<<endl;

  //  size_t wgtsize=  lhes->weights().size();
  //  std::cout << "weight size "<< wgtsize<<endl;
  if(useLHEW){
    for (size_t w = 1; w <=  (size_t)maxWeights; ++w)  {
      stringstream w_n;
      w_n << w;
      string name = "Event_LHEWeight"+w_n.str();
      cout << " pre single w # "<< w <<endl;
      trees["WeightHistory"]->Branch(name.c_str(),&float_values[name],(name+"/F").c_str());
      //trees["noSyst"]->Branch(name.c_str(), &float_values[name],(name+"/F").c_str());
    }
  }
}


double DMAnalysisTreeMaker::smear(double pt, double genpt, double eta, string syst){
  double resolScale = resolSF(fabs(eta), syst);
  double smear =1.0;
  if(genpt>0) smear = std::max((double)(0.0), (double)(pt + (pt - genpt) * resolScale) / pt);
  return  smear;
}

double DMAnalysisTreeMaker::resolSF(double eta, string syst)
{//from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
  double fac = 0.;
  if (syst == "jer__up")fac = 1.;
  if (syst == "jer__down")fac = -1.;
  if (eta <= 0.5)                       return 0.109 + (0.008 * fac);
  else if ( eta > 0.5 && eta <= 0.8 )   return 0.138 + (0.013 * fac);
  else if ( eta > 0.8 && eta <= 1.1 )   return 0.114 + (0.013 * fac);
  else if ( eta > 1.1 && eta <= 1.3 )   return 0.123 + (0.024 * fac);
  else if ( eta > 1.3 && eta <= 1.7 )   return 0.084 + (0.011 * fac);
  else if ( eta > 1.7 && eta <= 1.9 )   return 0.082 + (0.035 * fac);
  else if ( eta > 1.9 && eta <= 2.1 )   return 0.140 + (0.047 * fac);
  else if ( eta > 2.1 && eta <= 2.3 )   return 0.067 + (0.053 *fac);
  else if ( eta > 2.3 && eta <= 2.5 )   return 0.177 + (0.041 *fac);
  else if ( eta > 2.5 && eta <= 2.8 )   return 0.364 + (0.039 *fac);
  else if ( eta > 2.8 && eta <= 3.0 )   return 0.857 + (0.071 *fac);
  else if ( eta > 3.0 && eta <= 3.2 )   return 0.328 + (0.022 *fac);
  else if ( eta > 3.2 && eta <= 5.0 )   return 0.160 + (0.029 *fac);
  return 0.1;
 }

double DMAnalysisTreeMaker::getEffectiveArea(string particle, double eta){ 
  double aeta = fabs(eta);
  if(particle=="photon"){
    if(aeta<1.0)return 0.0725;
    if(aeta<1.479 && aeta >1.0)return 0.0604;
    if(aeta<2.0 && aeta >1.479)return 0.0320;
    if(aeta<2.2 && aeta >2.0)return 0.0512;
    if(aeta<2.3 && aeta >2.2)return 0.0766;
    if(aeta<2.4 && aeta >2.3)return 0.0949;
    if(aeta>2.4)return 0.1160;
  }
  if(particle=="ch_hadron"){
    if(aeta<1.0)return 0.0157;
    if(aeta<1.479 && aeta >1.0)return 0.0143;
    if(aeta<2.0 && aeta >1.479)return 0.0115;
    if(aeta<2.2 && aeta >2.0)return 0.0094;
    if(aeta<2.3 && aeta >2.2)return 0.0095;
    if(aeta<2.4 && aeta >2.3)return 0.0068;
    if(aeta>2.4)return 0.0053;
  }
  if(particle=="neu_hadron"){
    if(aeta<1.0)return 0.0143;
    if(aeta<1.479 && aeta >1.0)return 0.0210;
    if(aeta<2.0 && aeta >1.479)return 0.0147;
    if(aeta<2.2 && aeta >2.0)return 0.0082;
    if(aeta<2.3 && aeta >2.2)return 0.0124;
    if(aeta<2.4 && aeta >2.3)return 0.0186;
    if(aeta>2.4)return 0.0320;
  }
  return 0.0;


};
double DMAnalysisTreeMaker::jetUncertainty(double ptCorr, double eta, string syst)
{
  if(ptCorr<0)return ptCorr;
  if(syst == "jes__up" || syst == "jes__down"){
    double fac = 1.;
    if (syst == "jes__down")fac = -1.;
    jecUnc->setJetEta(eta);
    jecUnc->setJetPt(ptCorr);
    double JetCorrection = jecUnc->getUncertainty(true);
    return JetCorrection*fac;
  }
  return 0.0;
}

double DMAnalysisTreeMaker::jetUncertainty8(double ptCorr, double eta, string syst)
{
  if(ptCorr<0)return ptCorr;
  if(syst == "jes__up" || syst == "jes__down"){
    double fac = 1.;
    if (syst == "jes__down")fac = -1.;
    jecUnc8->setJetEta(eta);
    jecUnc8->setJetPt(ptCorr);
    double JetCorrection = jecUnc8->getUncertainty(true);
    return JetCorrection*fac;
  }
  return 0.0;
}

double DMAnalysisTreeMaker::getScaleFactor(double ptCorr,double etaCorr,double partonFlavour, string syst){
  return 1.0;
}


bool DMAnalysisTreeMaker::isMCWeightName(string s){
  
  if(s=="Z_Weight")return true;
  if(s=="W_Weight")return true;

  if(s=="Z_QCD_Weight")return true;
  if(s=="W_QCD_Weight")return true;

  if(s=="Z_EW_Weight")return true;
  if(s=="W_EW_Weight")return true;
  
  if(s=="T_Weight")return true;
  if(s=="T_Ext_Weight")return true;
    
  return false;

}

bool DMAnalysisTreeMaker::isInVector(std::vector<std::string> v, std::string s){
  for(size_t i = 0;i<v.size();++i){
    if(v.at(i)==s)return true;
  }
  return false;
}


//BTag weighter
bool DMAnalysisTreeMaker::BTagWeight::filter(int t)
{
    return (t >= minTags && t <= maxTags);
}

float DMAnalysisTreeMaker::BTagWeight::weight(vector<JetInfo> jetTags, int tags)
{
    if (!filter(tags))
    {
      //std::cout << "nThis event should not pass the selection, what is it doing here?" << std::endl;
        return 0;
    }
    int njetTags = jetTags.size();
    int comb = 1 << njetTags;
    float pMC = 0;
    float pData = 0;
    for (int i = 0; i < comb; i++)
    {
        float mc = 1.;
        float data = 1.;
        int ntagged = 0;
        for (int j = 0; j < njetTags; j++)
        {
            bool tagged = ((i >> j) & 0x1) == 1;
            if (tagged)
            {
                ntagged++;
                mc *= jetTags[j].eff;
                data *= jetTags[j].eff * jetTags[j].sf;
            }
            else
            {
                mc *= (1. - jetTags[j].eff);
                data *= (1. - jetTags[j].eff * jetTags[j].sf);
            }
        }

        if (filter(ntagged))
        {
	  //	  std::cout << mc << " " << data << endl;
            pMC += mc;
            pData += data;
        }
    }

    if (pMC == 0) return 0;
    return pData / pMC;
}


double DMAnalysisTreeMaker::MCTagEfficiencySubjet(string algo, int flavor, double pt, double eta){
  if ((abs(flavor) ==5)){ 
    if (algo=="csvm") {
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.284091; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.361111; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.283784; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.476489; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.454861; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.376923; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.502451; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.480892; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.439024; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.532209; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.527831; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.474074; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.54755; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.557831; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.467187; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.552198; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.573951; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.492498; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.535885; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.531685; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.453575; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.449269; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.44352; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.400552; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.390766; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.394737; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.317647; 
    }
    if (algo=="csvl") {
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.477273; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.555556; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.472973; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.658307; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.649306; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.588462; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.64951; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.707006; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.61324; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.691718; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.738964; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.654321; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.728146; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.736145; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.690625; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.730769; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.741722; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.732568; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.720694; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.70711; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.694771; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.677206; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.664588; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.680939; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.656532; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.642857; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.711765; 
    }
  }
  if ((abs(flavor) ==4)){ 
    if (algo=="csvm") {
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.0704225; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.0432692; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.0432692; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.0807834; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.0852941; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.0872675; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.0861298; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.0989957; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.0802469; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.101236; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.106457; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.087146; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.107405; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.112926; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.112971; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.110882; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.109333; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.092999; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.0974194; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.0953184; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.0986214; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.0934178; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.109236; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.113351; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.0968069; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.119436; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.1; 
    }
    if (algo=="csvl") {
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.232394; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.206731; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.192308; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.277846; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.305882; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.253219; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.289709; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.315638; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.257716; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.304482; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.316754; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.305011; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.307518; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.308949; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.320502; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.312049; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.332; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.338036; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.296452; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.31337; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.344115; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.304004; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.330776; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.385831; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.329954; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.397626; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.426087; 
    }
  }
  if ((abs(flavor)!= 5 && (abs(flavor)!= 4 ))){ 
    if (algo=="csvm") {
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.0127099; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.0142752; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.0150551; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.0149497; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.0148449; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.025; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.0102449; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.0128765; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.016835; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.010187; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.0139552; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.0218688; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.00983413; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.0115432; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.0188073; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.00940008; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.0106913; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.019997; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.01057; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.0139478; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.022888; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.0148112; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.0187542; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.0307771; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.022779; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.0282619; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.0396389; 
    }
    if (algo=="csvl") {
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.0755788; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.0950037; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.101357; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.10048; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.108998; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.137057; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.0830857; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.100368; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.129032; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.0839288; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.10876; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.137448; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.0784764; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.100403; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.150217; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.0774176; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.0937824; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.155002; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.0850308; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.10305; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.171304; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.109577; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.133381; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.215083; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.157024; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.179517; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.243524; 
    }
  }
  
  return 1.0;
}

double DMAnalysisTreeMaker::MCTagEfficiency(string algo, int flavor, double pt, double eta){

  if ((abs(flavor) ==5)){ 
    if (algo=="csvt") {
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.342571; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.288474; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.216505; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.437115; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.4047; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.338496; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.486694; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.422745; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.361853; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.481402; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.475984; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.374789; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.491533; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.460526; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.369161; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.448853; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.470953; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.347471; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.431072; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.421002; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.315734; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.361362; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.34923; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.295325; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.330469; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.348276; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.265306; 
    }
    if (algo=="csvm") {
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.534198; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.499001; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.414671; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.647904; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.627937; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.56167; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.67714; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.651373; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.576598; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.68872; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.699531; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.592855; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.695597; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.679656; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.576595; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.6739; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.680093; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.556281; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.644785; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.635838; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.55755; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.589675; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.590609; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.523375; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.53125; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.55977; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.483965; 
    }
    if (algo=="csvl") {
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.732311; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.714191; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.677025; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.810026; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.807441; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.779314; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.836807; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.826275; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.787961; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.848171; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.853738; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.799461; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.850127; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.857794; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.804119; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.83695; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.846631; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.797716; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.81984; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.806358; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.795143; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.804683; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.796038; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.777651; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.794531; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.791954; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.787172; 
    }
  }
  if ((abs(flavor) ==4)){ 
    if (algo=="csvt") {
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.0262541; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.0185376; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.0136986; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.0319807; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.0197183; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.0279493; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.0255741; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.01772; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.0205867; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.0232438; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.0279521; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.0255624; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.0260586; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.0359205; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.0274143; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.0343137; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.0388926; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.0204082; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.0273909; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.0384615; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.0315594; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.0278928; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.0291363; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.030622; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.0356135; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.0491315; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.0286369; 
    }
    if (algo=="csvm") { 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.112049; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.101442; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.0826069; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.138239; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.128773; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.132922; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.129436; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.113408; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.121462; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.137913; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.12664; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.121166; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.123779; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.140475; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.141433; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.143791; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.151615; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.115412; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.129062; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.142534; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.133045; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.118476; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.141172; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.146411; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.124824; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.144417; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.153494; 
    }
    if (algo=="csvl") {
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.320206; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.34243; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.288501; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.358666; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.4; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.375691; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.368476; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.384525; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.384457; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.370351; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.396463; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.373211; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.357763; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.381655; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.384424; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.383987; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.423863; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.399719; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.370474; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.397059; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.429455; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.360398; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.417967; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.465072; 
      if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.386812; 
      if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.456079; 
      if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.510882; 
    }
  }
  if ((abs(flavor)!= 5 && (abs(flavor)!= 4))){ 
      if (algo=="csvt") {
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.00259623; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.00213082; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.00150227; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.0014928; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.0017178; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.00184959; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.000814536; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.00163144; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.000850871; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.00123497; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.00141717; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.00166174; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.00107009; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.00175009; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.00151395; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.00120192; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.00189036; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.00296673; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.00119175; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.00238975; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.00293906; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.00236889; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.00339697; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.00545897; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.00495133; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.00795053; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.00723691; 
      }
      if (algo=="csvm") {
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.0128081; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.0132035; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.0141569; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.0117292; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.0147411; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.0181379; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.00733083; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.0127678; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.0150888; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.00708482; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.0106661; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.0147095; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.00825499; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.0119881; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.0188162; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.00773738; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.0116859; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.0221798; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.00991031; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.0137591; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.028191; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.013274; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.0217214; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.0386715; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.0226748; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.0373743; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.0433209; 
      }
      if (algo=="csvl") {
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  20) && (pt < 30))) return 0.0864027; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  20) && (pt < 30))) return 0.115483; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  20) && (pt < 30))) return 0.182538; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  30) && (pt < 50))) return 0.0854452; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  30) && (pt < 50))) return 0.116731; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  30) && (pt < 50))) return 0.164494; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  50) && (pt < 70))) return 0.0652256; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  50) && (pt < 70))) return 0.0963967; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  50) && (pt < 70))) return 0.130183; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  70) && (pt < 100))) return 0.0654534; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  70) && (pt < 100))) return 0.0936824; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  70) && (pt < 100))) return 0.137001; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  100) && (pt < 140))) return 0.0640526; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  100) && (pt < 140))) return 0.0920546; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  100) && (pt < 140))) return 0.145051; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  140) && (pt < 200))) return 0.0706881; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  140) && (pt < 200))) return 0.100962; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  140) && (pt < 200))) return 0.169457; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  200) && (pt < 300))) return 0.0755818; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  200) && (pt < 300))) return 0.107032; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  200) && (pt < 300))) return 0.183781; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  300) && (pt < 600))) return 0.0990443; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  300) && (pt < 600))) return 0.143629; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  300) && (pt < 600))) return 0.222533; 
	if(((fabs(eta)>=0) && (fabs(eta) <0.6)) && ((pt >=  600) && (pt < 1000))) return 0.152535; 
	if(((fabs(eta)>=0.6) && (fabs(eta) <1.2)) && ((pt >=  600) && (pt < 1000))) return 0.207257; 
	if(((fabs(eta)>=1.2) && (fabs(eta) <2.4)) && ((pt >=  600) && (pt < 1000))) return 0.287868; 
      }
    }
    
    return 1.0;
}
  
double DMAnalysisTreeMaker::TagScaleFactor(string algo, int flavor, string syst, double pt){
  double x = pt;
if(algo == "csvt"){
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){
        if (pt >= 20  && pt < 1000) return 0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x)));
      }
      if(abs(flavor)==4){
        if (pt >= 20  && pt < 1000) return 0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x)));
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
        if( pt>=20 && pt < 1000) return 0.971945+163.215/(x*x)+0.000517836*x;
      }
    }
     if(syst ==  "mistag_up") {
      if(abs(flavor)==5){
	if (pt >= 20  && pt < 1000) return 0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x)));
      }
      if(abs(flavor)==4){
	if (pt >= 20  && pt < 1000) return 0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x)));
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if (pt >= 20 && pt < 1000) return (0.971945+163.215/(x*x)+0.000517836*x)*(1+(0.291298+-0.000222983*x+1.69699e-07*x*x));
      }
    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){
	if (pt >= 20  && pt < 1000) return 0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x)));

      }
      if(abs(flavor)==4){
	if (pt >= 20  && pt < 1000) return 0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x)));
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if (pt >= 20 && pt < 1000) return (0.971945+163.215/(x*x)+0.000517836*x)*(1-(0.291298+-0.000222983*x+1.69699e-07*x*x));
      }
    }
  if(syst ==  "b_tag_up") {
      if(abs(flavor)==4){
	if (pt >= 20  && pt < 30 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.11806446313858032;
        if (pt >= 30  && pt < 50 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.054699532687664032;
        if (pt >= 50  && pt < 70 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.047356218099594116;
        if (pt >= 70  && pt < 100) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.047634456306695938;
        if (pt >= 100 && pt < 140) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.04632849246263504;
        if (pt >= 140 && pt < 200) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.048323042690753937;
        if (pt >= 200 && pt < 300) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.068715795874595642;
        if (pt >= 300 && pt < 600) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.10824859887361526;
	if (pt >= 600  && pt < 1000 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.18500012159347534;
      }
      if(abs(flavor)==5){
	if (pt >= 20  && pt < 30 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.033732704818248749;
        if (pt >= 30  && pt < 50 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.01562843844294548 ;
        if (pt >= 50  && pt < 70 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.013530348427593708;
        if (pt >= 70  && pt < 100) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.013609844259917736;
        if (pt >= 100 && pt < 140) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.013236711733043194;
        if (pt >= 140 && pt < 200) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.013806583359837532;
        if (pt >= 200 && pt < 300) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.019633084535598755;
        if (pt >= 300 && pt < 600) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.030928170308470726;
	if (pt >= 600  && pt < 1000 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.052857179194688797;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
        if( pt>= 20 && pt < 1000) return (0.971945+163.215/(x*x)+0.000517836*x)*(1+(0.291298+-0.000222983*x+1.69699e-07*x*x));
      }
    }
if(syst ==  "b_tag_down") {
      if(abs(flavor)==4){
	if (pt >= 20  && pt < 30 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.11806446313858032;
        if (pt >= 30  && pt < 50 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.054699532687664032;
        if (pt >= 50  && pt < 70 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.047356218099594116;
        if (pt >= 70  && pt < 100) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.047634456306695938;
        if (pt >= 100 && pt < 140) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.04632849246263504;
        if (pt >= 140 && pt < 200) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.048323042690753937;
        if (pt >= 200 && pt < 300) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.068715795874595642;
        if (pt >= 300 && pt < 600) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.10824859887361526;
	if (pt >= 600  && pt < 1000 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.18500012159347534;
      }
      if(abs(flavor)==5){
	if (pt >= 20  && pt < 30 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.033732704818248749;
        if (pt >= 30  && pt < 50 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.01562843844294548;
        if (pt >= 50  && pt < 70 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.013530348427593708;
        if (pt >= 70  && pt < 100) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.013609844259917736;
        if (pt >= 100 && pt < 140) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.013236711733043194;
        if (pt >= 140 && pt < 200) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.013806583359837532;
        if (pt >= 200 && pt < 300) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.019633084535598755;
        if (pt >= 300 && pt < 600) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.03092817030847072 ;
	if (pt >= 600  && pt < 1000 ) return (0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.052857179194688797;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
        if ( pt>=20 && pt<1000) return (0.971945+163.215/(x*x)+0.000517836*x)*(1-(0.291298+-0.000222983*x+1.69699e-07*x*x));
      }
    }
 }                 
  //Medium WP
  if(algo == "csvm"){
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){
	if (pt >= 20  && pt < 1000) return 0.561694*((1.+(0.31439*x))/(1.+(0.17756*x)));
      }
      if(abs(flavor)==4){
	if (pt >= 20  && pt < 1000) return 0.561694*((1.+(0.31439*x))/(1.+(0.17756*x)));
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if (pt >= 20 && pt < 1000) return 1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x;
      }
    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 670) return 0.561694*((1.+(0.31439*x))/(1.+(0.17756*x)));
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 670) return 0.561694*((1.+(0.31439*x))/(1.+(0.17756*x)));
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if (pt >= 20 && pt < 1000) return (1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x)*(1+(0.100485+3.95509e-05*x+-4.90326e-08*x*x));
      }
    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 670)  return 0.561694*((1.+(0.31439*x))/(1.+(0.17756*x)));

      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 670) return 0.561694*((1.+(0.31439*x))/(1.+(0.17756*x)));
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if (pt >= 20 && pt < 1000) return (1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x)*(1-(0.100485+3.95509e-05*x+-4.90326e-08*x*x));
      }
    }

    if(syst ==  "b_tag_up") {
      if(abs(flavor)==4){
	if (pt >= 20  && pt < 30 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.12064050137996674;
	if (pt >= 30  && pt < 50 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.042138919234275818;
	if (pt >= 50  && pt < 70 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.03711806982755661;
	if (pt >= 70  && pt < 100) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.036822021007537842;
	if (pt >= 100 && pt < 140) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.034397732466459274;
	if (pt >= 140 && pt < 200) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.0362386554479599;
	if (pt >= 200 && pt < 300) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.044985830783843994;
	if (pt >= 300 && pt < 600) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.064243391156196594;
	if (pt >= 600  && pt < 1000 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.097131341695785522;
      }
      if(abs(flavor)==5){
	if (pt >= 20  && pt < 30 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.040213499218225479;
	if (pt >= 30  && pt < 50 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.014046305790543556;
	if (pt >= 50  && pt < 70 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.012372690252959728;
	if (pt >= 70  && pt < 100) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.012274007312953472;
	if (pt >= 100 && pt < 140) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.011465910822153091;
	if (pt >= 140 && pt < 200) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.012079551815986633;
	if (pt >= 200 && pt < 300) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.014995276927947998;
	if (pt >= 300 && pt < 600) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.021414462476968765;
	if (pt >= 600  && pt < 1000 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.032377112656831741;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if( pt>=20 && pt<1000) return (1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x)*(1+(0.100485+3.95509e-05*x+-4.90326e-08*x*x));
      }
    }

    if(syst ==  "b_tag_down") {
      if(abs(flavor)==4){
	if (pt >= 20  && pt < 30 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.12064050137996674;
	if (pt >= 30  && pt < 50 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.042138919234275818;
	if (pt >= 50  && pt < 70 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.03711806982755661;
	if (pt >= 70  && pt < 100) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.036822021007537842;
	if (pt >= 100 && pt < 140) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.034397732466459274;
	if (pt >= 140 && pt < 200) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.0362386554479599;
	if (pt >= 200 && pt < 300) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.044985830783843994;
	if (pt >= 300 && pt < 670) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.064243391156196594;
	if (pt >= 600  && pt < 1000 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.097131341695785522;
      }
      if(abs(flavor)==5){
	if (pt >= 20  && pt < 30 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.040213499218225479 ;
	if (pt >= 30  && pt < 50 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.014046305790543556;
	if (pt >= 50  && pt < 70 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.012372690252959728;
	if (pt >= 70  && pt < 100) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.012274007312953472;
	if (pt >= 100 && pt < 140) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.011465910822153091;
	if (pt >= 140 && pt < 200) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.012079551815986633;
	if (pt >= 200 && pt < 300) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.014995276927947998;
	if (pt >= 300 && pt < 600) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.021414462476968765;
	if (pt >= 600  && pt < 1000 ) return (0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.032377112656831741;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return ((1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x)*(1-(0.100485+3.95509e-05*x+-4.90326e-08*x*x)));
      }
    }
  }

  //Loose WP
  if(algo == "csvl"){
    double x=pt;
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){
	if (pt >= 20  && pt < 1000) return 0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x)));
      }
      if(abs(flavor)==4){
	if (pt >= 20  && pt < 1000) return 0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x)));
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return 1.13904+-0.000594946*x+1.97303e-06*x*x+-1.38194e-09*x*x*x;
      }
    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 670) return 0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x)));
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 670) return 0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x)));
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return (1.13904+-0.000594946*x+1.97303e-06*x*x+-1.38194e-09*x*x*x)*(1+(0.0996438+-8.33354e-05*x+4.74359e-08*x*x));
      }
    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){
	if (pt >= 30  && pt < 670) return 0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x)));
      }
      if(abs(flavor)==4){
	if (pt >= 30  && pt < 670) return 0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x)));
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return (1.13904+-0.000594946*x+1.97303e-06*x*x+-1.38194e-09*x*x*x)*(1-(0.0996438+-8.33354e-05*x+4.74359e-08*x*x));
      }
    }
    
    if(syst ==  "b_tag_up") {
      if(abs(flavor)==4){
	if (pt >= 20  && pt < 30 ) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.063454590737819672; 
	if (pt >= 30  && pt < 50 ) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.031410016119480133;
	if (pt >= 50  && pt < 70 ) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.02891194075345993;
	if (pt >= 70  && pt < 100) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.028121808543801308;
	if (pt >= 100 && pt < 140) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.027028990909457207;
	if (pt >= 140 && pt < 200) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.027206243947148323;
	if (pt >= 200 && pt < 300) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.033642303198575974;
	if (pt >= 300 && pt < 600) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.04273652657866478;
	if (pt >= 600 && pt< 1000) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.054665762931108475;
      }
      if(abs(flavor)==5){
	if (pt >= 20  && pt < 30 ) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.025381835177540779;
      	if (pt >= 30  && pt < 50 ) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.012564006261527538;
	if (pt >= 50  && pt < 70 ) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.011564776301383972;
	if (pt >= 70  && pt < 100) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.011248723603785038;
	if (pt >= 100 && pt < 140) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.010811596177518368;
	if (pt >= 140 && pt < 200) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.010882497765123844;
	if (pt >= 200 && pt < 300) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.013456921093165874;
	if (pt >= 300 && pt < 600) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.017094610258936882;
	if (pt >= 600 && pt < 1000) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.02186630479991436;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return (1.13904+-0.000594946*x+1.97303e-06*x*x+-1.38194e-09*x*x*x)*(1+(0.0996438+-8.33354e-05*x+4.74359e-08*x*x));
      }
    }

    if(syst ==  "b_tag_down") {
      if(abs(flavor)==4){
	if (pt >= 20  && pt < 30 ) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.063454590737819672;
	if (pt >= 30  && pt < 50 ) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.031410016119480133;
	if (pt >= 50  && pt < 70 ) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.02891194075345993;
	if (pt >= 70  && pt < 100) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.028121808543801308;
	if (pt >= 100 && pt < 140) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.027028990909457207;
	if (pt >= 140 && pt < 200) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.027206243947148323;
	if (pt >= 200 && pt < 300) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.033642303198575974;
	if (pt >= 300 && pt < 600) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.04273652657866478;
	if (pt >= 600 && pt < 1000) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.054665762931108475;
      }
      if(abs(flavor)==5){
	if (pt >= 20 && pt < 30) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.025381835177540779;
      	if (pt >= 30  && pt < 50 ) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.012564006261527538;
	if (pt >= 50  && pt < 70 ) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.011564776301383972;
	if (pt >= 70  && pt < 100) return(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.011248723603785038;
	if (pt >= 100 && pt < 140) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.010811596177518368;
	if (pt >= 140 && pt < 200) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.010882497765123844;
	if (pt >= 200 && pt < 300) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.013456921093165874;
	if (pt >= 300 && pt < 600) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.017094610258936882;
	if (pt >= 600 && pt < 1000) return (0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.02186630479991436;
      }      
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return (1.13904+-0.000594946*x+1.97303e-06*x*x+-1.38194e-09*x*x*x)*(1-(0.0996438+-8.33354e-05*x+4.74359e-08*x*x));
      }
    }
  }

  return 1.0;

}

double DMAnalysisTreeMaker::TagScaleFactorSubjet(string algo, int flavor, string syst, double pt){
  double x = pt;

  //Medium WP
  if(algo == "csvm"){
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){
	if (pt >= 30   && pt < 120 ) return 0.97841;
	if (pt >= 120  && pt < 180 ) return 1.00499;
	if (pt >= 180  && pt < 240 ) return 1.01106;
	if (pt >= 240  && pt < 450 ) return 1.02189;
      }
      if(abs(flavor)==4){
	if (pt >= 30   && pt < 120 ) return 0.97841;
	if (pt >= 120  && pt < 180 ) return 1.00499;
	if (pt >= 180  && pt < 240 ) return 1.01106;
	if (pt >= 240  && pt < 450 ) return 1.02189;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if (pt >= 20 && pt < 1000) return 0.629961+0.00245187*x+-3.64539e-06*x*x+2.04999e-09*x*x*x;
      }
    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){
	if (pt >= 30   && pt < 120 ) return 0.97841;
	if (pt >= 120  && pt < 180 ) return 1.00499;
	if (pt >= 180  && pt < 240 ) return 1.01106;
	if (pt >= 240  && pt < 450 ) return 1.02189;
      }
      if(abs(flavor)==4){
	if (pt >= 30   && pt < 120 ) return 0.97841;
	if (pt >= 120  && pt < 180 ) return 1.00499;
	if (pt >= 180  && pt < 240 ) return 1.01106;
	if (pt >= 240  && pt < 450 ) return 1.02189;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if (pt >= 20 && pt < 1000) return 0.676736+0.00286128*x+-4.34618e-06*x*x+2.44485e-09*x*x*x;
      }
    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){
	if (pt >= 30   && pt < 120 ) return 0.97841;
	if (pt >= 120  && pt < 180 ) return 1.00499;
	if (pt >= 180  && pt < 240 ) return 1.01106;
	if (pt >= 240  && pt < 450 ) return 1.02189;

      }
      if(abs(flavor)==4){
	if (pt >= 30   && pt < 120 ) return 0.97841;
	if (pt >= 120  && pt < 180 ) return 1.00499;
	if (pt >= 180  && pt < 240 ) return 1.01106;
	if (pt >= 240  && pt < 450 ) return 1.02189;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if (pt >= 20 && pt < 1000) return 0.582177+0.00204348*x+-2.94226e-06*x*x+1.6536e-09*x*x*x;
      }
    }

    if(syst ==  "b_tag_up") {
      if(abs(flavor)==5){
	if (pt >= 30   && pt < 120 ) return 1.01236;
	if (pt >= 120  && pt < 180 ) return 1.02287;
	if (pt >= 180  && pt < 240 ) return 1.0347;
	if (pt >= 240  && pt < 450 ) return 1.05521;
      }
      if(abs(flavor)==4){
	if (pt >= 30   && pt < 120 ) return 1.04631;
	if (pt >= 120  && pt < 180 ) return 1.04074;
	if (pt >= 180  && pt < 240 ) return 1.05834;
	if (pt >= 240  && pt < 450 ) return 1.08854;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if( pt>=20 && pt<1000) return 0.676736+0.00286128*x+-4.34618e-06*x*x+2.44485e-09*x*x*x;
      }
    }

    if(syst ==  "b_tag_down") {
      if(abs(flavor)==5){
	if (pt >= 30   && pt < 120 ) return 0.94446;
	if (pt >= 120  && pt < 180 ) return 0.98711;
	if (pt >= 180  && pt < 240 ) return 0.98742;
	if (pt >= 240  && pt < 450 ) return 0.98857;
      }
      if(abs(flavor)==4){
	if (pt >= 30   && pt < 120 ) return 0.9105;
	if (pt >= 120  && pt < 180 ) return 0.96924;
	if (pt >= 180  && pt < 240 ) return 0.96378;
	if (pt >= 240  && pt < 450 ) return 0.95524;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return 0.582177+0.00204348*x+-2.94226e-06*x*x+1.6536e-09*x*x*x;
      }
    }
  }

  //Loose WP
  if(algo == "csvl"){
    double x=pt;
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){
	if (pt >= 30   && pt < 120 ) return 0.99839;
	if (pt >= 120  && pt < 180 ) return 1.0022;
	if (pt >= 180  && pt < 240 ) return 1.00468;
	if (pt >= 240  && pt < 450 ) return 1.01764;
      }
      if(abs(flavor)==4){
	if (pt >= 30   && pt < 120 ) return 0.99839;
	if (pt >= 120  && pt < 180 ) return 1.0022;
	if (pt >= 180  && pt < 240 ) return 1.00468;
	if (pt >= 240  && pt < 450 ) return 1.01764;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return 0.954689+0.000316059*x+3.22024e-07*x*x+-4.06201e-10*x*x*x;
      }
    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){
	if (pt >= 30   && pt < 120 ) return 0.99839;
	if (pt >= 120  && pt < 180 ) return 1.0022;
	if (pt >= 180  && pt < 240 ) return 1.00468;
	if (pt >= 240  && pt < 450 ) return 1.01764;
      }
      if(abs(flavor)==4){
	if (pt >= 30   && pt < 120 ) return 0.99839;
	if (pt >= 120  && pt < 180 ) return 1.0022;
	if (pt >= 180  && pt < 240 ) return 1.00468;
	if (pt >= 240  && pt < 450 ) return 1.01764;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return 1.0358+0.000107516*x+9.58049e-07*x*x+-8.59906e-10*x*x*x;
      }
    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){
	if (pt >= 30   && pt < 120 ) return 0.99839;
	if (pt >= 120  && pt < 180 ) return 1.0022;
	if (pt >= 180  && pt < 240 ) return 1.00468;
	if (pt >= 240  && pt < 450 ) return 1.01764;
      }
      if(abs(flavor)==4){
	if (pt >= 30   && pt < 120 ) return 0.99839;
	if (pt >= 120  && pt < 180 ) return 1.0022;
	if (pt >= 180  && pt < 240 ) return 1.00468;
	if (pt >= 240  && pt < 450 ) return 1.01764;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return 0.873638+0.0005247*x+-3.15316e-07*x*x+4.83633e-11*x*x*x;
      }
    }
    
    if(syst ==  "b_tag_up") {
      if(abs(flavor)==5){
	if (pt >= 30   && pt < 120 ) return 1.0091;
	if (pt >= 120  && pt < 180 ) return 1.01303;
	if (pt >= 180  && pt < 240 ) return 1.01659;
	if (pt >= 240  && pt < 450 ) return 1.03703;
      }
      if(abs(flavor)==4){
	if (pt >= 30   && pt < 120 ) return 1.01981;
	if (pt >= 120  && pt < 180 ) return 1.02386;
	if (pt >= 180  && pt < 240 ) return 1.02849;
	if (pt >= 240  && pt < 450 ) return 1.05642;
      }
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return 1.0358+0.000107516*x+9.58049e-07*x*x+-8.59906e-10*x*x*x;
      }
    }

    if(syst ==  "b_tag_down") {
      if(abs(flavor)==5){
	if (pt >= 30   && pt < 120 ) return 0.98769;
	if (pt >= 120  && pt < 180 ) return 0.99137;
	if (pt >= 180  && pt < 240 ) return 0.99277;
	if (pt >= 240  && pt < 450 ) return 0.99825;
      }
      if(abs(flavor)==4){
	if (pt >= 30   && pt < 120 ) return 0.97698;
	if (pt >= 120  && pt < 180 ) return 0.98054;
	if (pt >= 180  && pt < 240 ) return 0.98087;
	if (pt >= 240  && pt < 450 ) return 0.97886;
      }      
      if(abs(flavor)!=5 && abs(flavor)!=4){
	if(pt>=20 && pt<1000) return 0.873638+0.0005247*x+-3.15316e-07*x*x+4.83633e-11*x*x*x;
      }
    }
  }

  return 1.0;

}

float DMAnalysisTreeMaker::BTagWeight::weightWithVeto(vector<JetInfo> jetsTags, int tags, vector<JetInfo> jetsVetoes, int vetoes)
{//This function takes into account cases where you have n b-tags and m vetoes, but they have different thresholds. 
    if (!filter(tags))
    {
        //   std::cout << "nThis event should not pass the selection, what is it doing here?" << std::endl;
        return 0;
    }
    int njets = jetsTags.size();
    if(njets != (int)(jetsVetoes.size()))return 0;//jets tags and vetoes must have same size!
    int comb = 1 << njets;
    float pMC = 0;
    float pData = 0;
    for (int i = 0; i < comb; i++)
    {
        float mc = 1.;
        float data = 1.;
        int ntagged = 0;
        for (int j = 0; j < njets; j++)
        {
            bool tagged = ((i >> j) & 0x1) == 1;
            if (tagged)
            {
                ntagged++;
                mc *= jetsTags[j].eff;
                data *= jetsTags[j].eff * jetsTags[j].sf;
            }
            else
            {
                mc *= (1. - jetsVetoes[j].eff);
                data *= (1. - jetsVetoes[j].eff * jetsVetoes[j].sf);
            }
        }

        if (filter(ntagged))
        {
            //  std::cout << mc << " " << data << endl;
            pMC += mc;
            pData += data;
        }
    }

    if (pMC == 0) return 0;
    return pData / pMC;
}

bool DMAnalysisTreeMaker::isEWKID(int id){
  bool isewk=false;
  int aid = abs(id);
  if((aid>10 && aid <17) || (aid==23 || aid==24)){isewk=true;}

  return isewk;
}

double DMAnalysisTreeMaker::pileUpSF(string syst)
{
  // if (syst == "PUUp" )return LumiWeightsUp_.weight( n0);
  // if (syst == "PUDown" )return LumiWeightsDown_.weight( n0);
  // return LumiWeights_.weight( n0);

  if (syst == "PUUp" )return LumiWeightsUp_.weight(n0);
  if (syst == "PUDown" )return LumiWeightsDown_.weight(n0);
  return LumiWeights_.weight( n0);

}

void DMAnalysisTreeMaker::endJob(){
  //  for(size_t s=0;s< systematics.size();++s){
  //    std::string syst  = systematics.at(s);
  /*  cout <<" init events are "<< nInitEvents <<endl;
  trees["EventHistory"]->SetBranchAddress("initialEvents",&nInitEvents);
  for(size_t i = 0; i < (size_t)nInitEvents;++i){
      //      trees["EventHistory"]->GetBranch("initialEvents")->Fill();
      
      //      trees["EventHistory"]->GetBranch("initialEvents")->Fill();
      trees["EventHistory"]->Fill();
      //      cout <<" i is "<< i << " entry is now "<< trees["EventHistory"]->GetBranch("initialEvents")->GetEntry()<<endl;
      
      }*/
  ;
}



//DMAnalysisTreeMaker::~DMAnalysisTreeMaker(const edm::ParameterSet& iConfig)
// ------------ method called once each job just after ending the event loop  ------------


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(DMAnalysisTreeMaker);

//  LocalWords:  firstidx
