/**
 * DMAnalysisTreeMaker
 * 
 * Produces analysis trees from edm-ntuples adding extra variables for resolved and unresolved tops
 * For custom systematics scenarios
 * 
 * \Author A. Orso M. Iorio
 * 
 * 
 *\version  $Id:
 * 
 * 
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
//#include "./EquationSolver.h"
//#include "./DMTopVariables.h"
#include "./Weights.h"


#include "TFile.h"
#include "TSpline.h"
#include "TTree.h"
#include "TMath.h"
#include <vector>
#include <algorithm>
#include <TLorentzVector.h>
#include <TMVA/Reader.h>
#include <string>
#include <iostream>
//#include "TopTagger/Resolved/interface/KinematicFitter.hh"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
//#include "DataFormats/Common/interface/
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
//using namespace reco;
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
  void fillCategory(string label, string category, int pos_nocat,int pos_cat);
  void fillSystCategory(string label, string category, string syst, int pos_nocat,int pos_cat, string scancut="");
  void setCatCategoryValue(string label,string category, int pos_cat,string var, float value);

  void fillScanCuts(string label, string category, int pos_nocat);
  void fillScanSystsCuts(string label, string category, int pos_nocat);
  void fillSysts(string label, string category, int pos_nocat, string var, string cut, string prefix="");
  bool passesScanCut(string label, string category, string scanCut, int pos_nocat);
  string parseScanCut(string category, int obj);
  bool isScanCut(string label, string category);
  bool isSysCat(string label, string category);
  void resetVariables(string cat, int sizecat);


  void setEventBTagSF(string label, string category, string algo);
  void setEventLeptonSF(string label, string category);


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
  
  float getReshapedBTagValue(float flavor, float btag, float pt, float eta, string algo, string syst);
  
  //  helpers
  double getPhi(double px, double py);
  double getMCTagEfficiencyFuncParam(float flavor, float ptCorr, float eta, string algo,string syst, string wp, string region);
  float avgRatioFunc(float x0,float x1,float a1, float b1 , float a2, float b2);
  float getEffRatioFunc(float flavor,float btag,float pt, float eta, string algo ,string syst,bool normalize);
  //
  // Set up MVA reader
  //
  // spectator variables, not used for MVA evaluation
  int isSig, b_mis, w_mis, wb_mis;
  float mtop,leadingLeptonCharge,topCharge;
  // MVA input variables
  //  float bdt_qgid1, bdt_qgid2;
  //  float bdt_dphij1b, bdt_dphij2b, bdt_drj1b, bdt_drj2b;
  //  float bdt_bjcsv, bdt_jet1csv, bdt_jet2csv;
  //  float bdt_prob;
  //  TMVA::Reader res_topmvaReader;

  double jetUncertainty(double pt, double eta, string syst);
  // double smearPt(double pt, double genpt, double eta, string syst);
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
  edm::EDGetTokenT< GenEventInfoProduct > t_genprod_;
  edm::EDGetTokenT< std::vector<string> > t_triggerNames_;
  edm::EDGetTokenT< std::vector<float> > t_triggerBits_;
  edm::EDGetTokenT< std::vector<int> > t_triggerPrescales_;

  edm::EDGetTokenT<std::vector<string> > t_triggerNamesR_;

  edm::EDGetTokenT< unsigned int > t_lumiBlock_;
  edm::EDGetTokenT< unsigned int > t_runNumber_;
  edm::EDGetTokenT< ULong64_t > t_eventNumber_;
  edm::EDGetTokenT< bool > t_HBHEFilter_;
  edm::EDGetTokenT< bool > t_HBHEIsoFilter_;
  edm::EDGetTokenT< std::vector<string> > t_metNames_;
  edm::EDGetTokenT< std::vector<float> > t_metBits_;

  edm::EDGetTokenT< std::vector<float> > t_pvZ_,t_pvChi2_,t_pvRho_;
  edm::EDGetTokenT< std::vector<int> > t_pvNdof_;

  edm::EDGetTokenT< double > t_Rho_;
  edm::EDGetTokenT<int> t_ntrpu_;
  
  edm::EDGetTokenT< std::vector<float> > jetAK8topSubjetIndex0;
  edm::EDGetTokenT< std::vector<float> > jetAK8topSubjetIndex1;
  edm::EDGetTokenT< std::vector<float> > jetAK8topSubjetIndex2;
  edm::EDGetTokenT< std::vector<float> > jetAK8topSubjetIndex3;

  edm::EDGetTokenT <std::vector<std::vector<int>>> t_jetKeys_, t_muKeys_;
  
  //  edm::EDGetTokenT<edm::Association<reco::GenParticleCollection> > t_genParticles;
  edm::EDGetTokenT< reco::GenParticleCollection> t_genParticles_;
  edm::EDGetTokenT< reco::GenJetCollection> t_genJets_;
  

  edm::EDGetTokenT< std::vector<float> > genPartID ;
  edm::EDGetTokenT< std::vector<float> > genPartStatus;
  edm::EDGetTokenT< std::vector<float> > genPartMom0ID; 
  edm::EDGetTokenT< std::vector<float> > genPartPt;
  edm::EDGetTokenT< std::vector<float> > genPartPhi;
  edm::EDGetTokenT< std::vector<float> > genPartEta;
  edm::EDGetTokenT< std::vector<float> > genPartE;

  //--------------------------------------------------------------------------------------

  edm::LumiReWeighting LumiWeights_, LumiWeightsUp_, LumiWeightsDown_;
  
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
  map< string , std::vector<string> > obj_cats, obj_scanCuts, obj_systCats;
  map< string , bool > passesJecCuts;

  map< string , double[100] > vdoubles_values;
  map< string , double[100] > vdouble_values;
  map< string , double > double_values;


  map< string , float > float_values;
  map< string , int > int_values;
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
  

  string mu_label, ele_label, jets_label, boosted_tops_label, boosted_tops_subjets_label, met_label, photon_label;//metNoHF_label

  bool getPartonW, getPartonTop, doWReweighting, doTopReweighting;
  bool getParticleWZ, getWZFlavour;
  bool usePrunedGenParticles,useGenJets;

  //Do resolved top measurement:
  bool doResolvedTopHad,doResolvedTopSemiLep;
  int max_leading_jets_for_top;
  int max_bjets_for_top;
  int n0;

  //EffectiveAreas
  bool recalculateEA;

  //MC info:
  edm::ParameterSet channelInfo;
  std::string channel;
  double crossSection, originalEvents;
  bool useLHEWeights, useLHE, useTriggers,cutOnTriggers, useMETFilters, addPV;
  bool addLHAPDFWeights;
  string centralPdfSet,variationPdfSet;
  std::vector<string> SingleElTriggers, SingleMuTriggers, PhotonTriggers, hadronicTriggers,metFilters;
  int maxPdf, maxWeights;
  edm::Handle<LHEEventProduct > lhes;
  edm::Handle<GenEventInfoProduct> genprod;

  // Do top decay reshaping to reduce b contamination/ 
  bool doTopDecayReshaping;

  //Trigger info
  edm::Handle<std::vector<float> > triggerBits;
  edm::Handle<std::vector<string> > triggerNames;
  edm::Handle<std::vector<int> > triggerPrescales;
  edm::Handle<std::vector<string> > triggerNamesR ;  
  edm::Handle<bool> HBHE;
  edm::Handle<bool> HBHEIso;
  
  edm::Handle<unsigned int> lumiBlock;
  edm::Handle<unsigned int> runNumber;
  edm::Handle<ULong64_t> eventNumber;
  //  edm::Handle<double> eventNumber;
  edm::Handle<reco::GenParticleCollection> genParts;
  edm::Handle<reco::GenJetCollection> genJets;


  edm::Handle<std::vector<float> > metBits;
  edm::Handle<std::vector<string> > metNames;
  
  edm::Handle<std::vector<float> > pvZ,pvChi2,pvRho;
  edm::Handle<std::vector<int> > pvNdof;

  edm::InputTag metNames_;

  edm::Handle<std::vector<float> > ak8jetSubjetIndex0;
  edm::Handle<std::vector<float> > ak8jetSubjetIndex1;
  edm::Handle<std::vector<float> > ak8jetSubjetIndex2;
  edm::Handle<std::vector<float> > ak8jetSubjetIndex3;
  


edm::Handle<std::vector<std::vector<int>>> jetKeys;
edm::Handle<std::vector<std::vector<int>>> muKeys;
  
  //  edm::InputTag partID_,partStatus_,partMomID_,partPt_;
  edm::Handle<std::vector<float> > partID;
  edm::Handle<std::vector<float> > partStatus;
  edm::Handle<std::vector<float> > partMom0ID;
  edm::Handle<std::vector<float> > partPt;
  edm::Handle<std::vector<float> > partEta;
  edm::Handle<std::vector<float> > partPhi;
  edm::Handle<std::vector<float> > partE;

  float nPV;
  edm::Handle<int> ntrpu;

  //JEC info
  bool changeJECs,doT1MET, doResol;
  bool isData, applyRes;
  edm::Handle<double> rho;
  double Rho;
  string JECVersion;
  
  //edm::Handle<double> Rho;
  std::vector<double> jetScanCuts;
  std::vector<JetCorrectorParameters> jecPars, jecParsL1_vect
;
  JetCorrectorParameters *jecParsL1, *jecParsL1RC, *jecParsL2, *jecParsL3, *jecParsL2L3Residuals;
  JetCorrectionUncertainty *jecUnc;
  FactorizedJetCorrector *jecCorr,*jecCorr_L1;

  bool isFirstEvent;
  //Do preselection
  bool doPreselection;
  //Reweight to get tWs/tWd decay verices:
  bool doTopBToLightQuarkReweight;
  //BTag
  BTagCalibration *calib;
  BTagCalibration *calib_cmvav2;
  BTagCalibrationReader *readerCSVLoose,*readerCSVLooseUDSG, *readerCSVMedium, *readerCSVTight, *readerCSVReshape;
  BTagCalibrationReader *readerCMVALoose,*readerCMVALooseUDSG, *readerCMVAMedium, *readerCMVATight, *readerCMVAReshape;

  string filename_cmva;
  TFile* file_cmva;
  Weights *cmvaeffbt,*cmvaeffbm,*cmvaeffbl;
  Weights *cmvaeffct,*cmvaeffcm,*cmvaeffcl;
  Weights *cmvaeffot,*cmvaeffom,*cmvaeffol;

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
  
  BTagWeight b_csvt_0_tags= BTagWeight(0,0),
    b_csvt_1_tag= BTagWeight(1,1),
    b_csvt_1_2_tags= BTagWeight(1,4),
    b_csvt_2_tags= BTagWeight(2,4);
  
  double b_weight_csvt_0_tags,
    b_weight_csvt_1_tag,
    b_weight_csvt_1_2_tags,
    b_weight_csvt_2_tags;
  double b_weight_csvt_0_tags_mistag_up,
    b_weight_csvt_1_tag_mistag_up,
    b_weight_csvt_1_2_tags_mistag_up,
    b_weight_csvt_2_tags_mistag_up;
  double b_weight_csvt_0_tags_mistag_down,
    b_weight_csvt_1_tag_mistag_down,
    b_weight_csvt_1_2_tags_mistag_down,
    b_weight_csvt_2_tags_mistag_down;
  double b_weight_csvt_0_tags_b_tag_down,
    b_weight_csvt_1_tag_b_tag_down,
    b_weight_csvt_1_2_tags_b_tag_down,
    b_weight_csvt_2_tags_b_tag_down;
  double b_weight_csvt_0_tags_b_tag_up,
    b_weight_csvt_1_tag_b_tag_up,
    b_weight_csvt_1_2_tags_b_tag_up,
    b_weight_csvt_2_tags_b_tag_up;

  BTagWeight b_csvm_0_tags= BTagWeight(0,0),
    b_csvm_1_tag= BTagWeight(1,1),
    b_csvm_1_2_tags= BTagWeight(1,4),
    b_csvm_2_tags= BTagWeight(2,4);
  
  double b_weight_csvm_0_tags,
    b_weight_csvm_1_tag,
    b_weight_csvm_1_2_tags,
    b_weight_csvm_2_tags;
  double b_weight_csvm_0_tags_mistag_up,
    b_weight_csvm_1_tag_mistag_up,
    b_weight_csvm_1_2_tags_mistag_up,
    b_weight_csvm_2_tags_mistag_up;
  double b_weight_csvm_0_tags_mistag_down,
    b_weight_csvm_1_tag_mistag_down,
    b_weight_csvm_1_2_tags_mistag_down,
    b_weight_csvm_2_tags_mistag_down;
  double b_weight_csvm_0_tags_b_tag_down,
    b_weight_csvm_1_tag_b_tag_down,
    b_weight_csvm_1_2_tags_b_tag_down,
    b_weight_csvm_2_tags_b_tag_down;
  double b_weight_csvm_0_tags_b_tag_up,
    b_weight_csvm_1_tag_b_tag_up,
    b_weight_csvm_1_2_tags_b_tag_up,
    b_weight_csvm_2_tags_b_tag_up;

  BTagWeight b_csvl_0_tags= BTagWeight(0,0),
    b_csvl_1_tag= BTagWeight(1,1),
    b_csvl_1_2_tags= BTagWeight(1,4),
    b_csvl_2_tags= BTagWeight(2,4);
  
  double b_weight_csvl_0_tags_mistag_up,
    b_weight_csvl_1_tag_mistag_up,
    b_weight_csvl_1_2_tags_mistag_up,
    b_weight_csvl_2_tags_mistag_up;
  double b_weight_csvl_0_tags_mistag_down,
    b_weight_csvl_1_tag_mistag_down,
    b_weight_csvl_1_2_tags_mistag_down,
    b_weight_csvl_2_tags_mistag_down;
  double b_weight_csvl_0_tags_b_tag_down,
    b_weight_csvl_1_tag_b_tag_down,
    b_weight_csvl_1_2_tags_b_tag_down,
    b_weight_csvl_2_tags_b_tag_down;
  double b_weight_csvl_0_tags_b_tag_up,
    b_weight_csvl_1_tag_b_tag_up,
    b_weight_csvl_1_2_tags_b_tag_up,
    b_weight_csvl_2_tags_b_tag_up;
  double b_weight_csvl_0_tags,
    b_weight_csvl_1_tag,
    b_weight_csvl_1_2_tags,
    b_weight_csvl_2_tags;

  double MCTagEfficiency(string algo, int flavor, double pt, double eta=0.); 
  double TagScaleFactor(string algo, int flavor, string syst,double pt, double eta=0.);
  double getMCTagEfficiencyFunc(float flavor, float btag, float pt, float eta, string algo,string syst, bool norm=false, bool spline=false);
  double getMCTagEfficiencyFuncInt(float flavor, float ptCorr, float eta, string algo, string syst);
  float getWPAlgo(string algo, string wp);
  bool doBTagSF;
  bool doPU;
  
  string season;
  
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
  physObjects = iConfig.template getParameter<std::vector<edm::ParameterSet> >("physicsObjects");
  
  channelInfo = iConfig.getParameter<edm::ParameterSet >("channelInfo"); // The physics of the channel, e.g. the cross section, #original events, etc.
  channel = channelInfo.getParameter<std::string>("channel");
  crossSection = channelInfo.getParameter<double>("crossSection");
  originalEvents = channelInfo.getParameter<double>("originalEvents");

  doPreselection = iConfig.getUntrackedParameter<bool>("doPreselection",false);
  doPU = iConfig.getUntrackedParameter<bool>("doPU",true);

  doTopBToLightQuarkReweight = iConfig.getUntrackedParameter<bool>("doTopBToQReweight",false);

  useLHEWeights = channelInfo.getUntrackedParameter<bool>("useLHEWeights",false);
  useLHE = channelInfo.getUntrackedParameter<bool>("useLHE",false);
  addLHAPDFWeights = channelInfo.getUntrackedParameter<bool>("addLHAPDFWeights",false);

  getPartonW = channelInfo.getUntrackedParameter<bool>("getPartonW",false);
  getParticleWZ = channelInfo.getUntrackedParameter<bool>("getParticleWZ",false);
  getPartonTop = channelInfo.getUntrackedParameter<bool>("getPartonTop",false);
  doWReweighting = channelInfo.getUntrackedParameter<bool>("doWReweighting",false);
  doTopReweighting = channelInfo.getUntrackedParameter<bool>("doTopReweighting",false);

  getWZFlavour = channelInfo.getUntrackedParameter<bool>("getWZFlavour",false);

  doTopDecayReshaping = channelInfo.getUntrackedParameter<bool>("doTopDecayReshaping",false);

  doResolvedTopSemiLep = iConfig.getUntrackedParameter<bool>("doResolvedTopSemiLep",false);
  doResolvedTopHad = iConfig.getUntrackedParameter<bool>("doResolvedTopHad",false);

  usePrunedGenParticles = iConfig.getUntrackedParameter<bool>("usePrunedGenParticles",true);
  useGenJets = iConfig.getUntrackedParameter<bool>("useGenJets",false);


  edm::InputTag genprod_ = iConfig.getParameter<edm::InputTag>( "genprod" );
  t_genprod_ = consumes<GenEventInfoProduct>( genprod_ );
  
  t_genParticles_ = consumes< reco::GenParticleCollection >(iConfig.getParameter<edm::InputTag>("genParticles"));
  if(useGenJets)  t_genJets_ = consumes< reco::GenJetCollection >(iConfig.getParameter<edm::InputTag>("genJets"));
  

  useTriggers = iConfig.getUntrackedParameter<bool>("useTriggers",true);
  cutOnTriggers = iConfig.getUntrackedParameter<bool>("cutOnTriggers",true);

  edm::InputTag ak8jetSubjetIndex0_ = iConfig.getParameter<edm::InputTag>("ak8jetSubjetIndex0");
  jetAK8topSubjetIndex0 = consumes< std::vector<float> >( ak8jetSubjetIndex0_);
  edm::InputTag ak8jetSubjetIndex1_ = iConfig.getParameter<edm::InputTag>("ak8jetSubjetIndex1");
  jetAK8topSubjetIndex1 = consumes< std::vector<float> >( ak8jetSubjetIndex1_);
  edm::InputTag ak8jetSubjetIndex2_ = iConfig.getParameter<edm::InputTag>("ak8jetSubjetIndex2");
  jetAK8topSubjetIndex2 = consumes< std::vector<float> >( ak8jetSubjetIndex2_);
  edm::InputTag ak8jetSubjetIndex3_ = iConfig.getParameter<edm::InputTag>("ak8jetSubjetIndex3");
  jetAK8topSubjetIndex3 = consumes< std::vector<float> >( ak8jetSubjetIndex3_);
  
  edm::InputTag PartID_ = iConfig.getParameter<edm::InputTag>("partID");
  genPartID = consumes< std::vector<float> >( PartID_ );
  edm::InputTag PartStatus_ = iConfig.getParameter<edm::InputTag>("partStatus");
  genPartStatus = consumes< std::vector<float> >( PartStatus_ );
  edm::InputTag PartMom0ID_ = iConfig.getParameter<edm::InputTag>("partMom0ID");
  genPartMom0ID = consumes< std::vector<float> >( PartMom0ID_ );
  edm::InputTag PartPt_ = iConfig.getParameter<edm::InputTag>("partPt");
  genPartPt = consumes< std::vector<float> >( PartPt_ );
  edm::InputTag PartPhi_ = iConfig.getParameter<edm::InputTag>("partPhi");
  genPartPhi = consumes< std::vector<float> >( PartPhi_);
  edm::InputTag PartEta_ = iConfig.getParameter<edm::InputTag>("partEta");
  genPartEta = consumes< std::vector<float> >( PartEta_ );
  edm::InputTag PartE_ = iConfig.getParameter<edm::InputTag>("partE");
  genPartE = consumes< std::vector<float> >( PartE_);
  
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

    //testing to read from Run
    edm::InputTag triggerNamesR_ = iConfig.getParameter<edm::InputTag>("triggerNames");
    t_triggerNamesR_ = mayConsume< std::vector<string>, edm::InRun>(edm::InputTag("TriggerUserData","triggerNameTree"));

    SingleElTriggers= channelInfo.getParameter<std::vector<string> >("SingleElTriggers");
    SingleMuTriggers= channelInfo.getParameter<std::vector<string> >("SingleMuTriggers");
    PhotonTriggers= channelInfo.getParameter<std::vector<string> >("PhotonTriggers");
    hadronicTriggers= channelInfo.getParameter<std::vector<string> >("hadronicTriggers");


  }
  useMETFilters = iConfig.getUntrackedParameter<bool>("useMETFilters",true);
  if(useMETFilters){
    metFilters = channelInfo.getParameter<std::vector<string> >("metFilters");
    edm::InputTag metBits_ = iConfig.getParameter<edm::InputTag>("metBits");
    t_metBits_ = consumes< std::vector<float> >( metBits_ );
    metNames_ = iConfig.getParameter<edm::InputTag>("metNames");
    t_metNames_ = consumes< std::vector<string>, edm::InRun >( metNames_ );
  }
  
  addPV = iConfig.getUntrackedParameter<bool>("addPV",true);
  changeJECs = iConfig.getUntrackedParameter<bool>("changeJECs",false);
  doResol = iConfig.getUntrackedParameter<bool>("doResol",true);
  JECVersion = iConfig.getParameter<std::string >("JECVersion");
  doT1MET = iConfig.getUntrackedParameter<bool>("doT1MET",true);
  recalculateEA = iConfig.getUntrackedParameter<bool>("recalculateEA",true);

  isData = iConfig.getUntrackedParameter<bool>("isData",false);
  applyRes = iConfig.getUntrackedParameter<bool>("applyRes",false);
  
  t_Rho_ = consumes<double>( edm::InputTag( "fixedGridRhoFastjetAll" ) ) ;
  
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

  if (doPU){ t_ntrpu_ = consumes< int >( edm::InputTag( "eventUserData","puNtrueInt" ) );
  }
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
  TFileDirectory DMTrees;// = fs->mkdir( "systematics_trees" );

  if(addNominal){
    DMTrees = fs->mkdir( "systematics_trees" );
  }
  trees["noSyst"] =  new TTree((channel+"__noSyst").c_str(),(channel+"__noSyst").c_str());

  nInitEventsHisto = new TH1D("initialEvents","initalEvents",10,0,10);
  
  if (useLHE){
  edm::InputTag lhes_ = iConfig.getParameter<edm::InputTag>( "lhes" );
  t_lhes_ = consumes< LHEEventProduct >( lhes_ );
  }
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
    std::vector<std::string > scanCuts = itPsets->getParameter<std::vector<std::string> >("scanCuts");
    std::vector<std::string > systCats = itPsets->getParameter<std::vector<std::string> >("systCats");
    std::vector<std::string > toSave= itPsets->getParameter<std::vector<std::string> >("toSave");
    
    std::vector<edm::InputTag >::const_iterator itF = variablesFloat.begin();
    std::vector<edm::InputTag >::const_iterator itI = variablesInt.begin();
    
    std::vector<edm::InputTag >::const_iterator itsF = singleFloat.begin();
    std::vector<edm::InputTag >::const_iterator itsD = singleDouble.begin();
    std::vector<edm::InputTag >::const_iterator itsI = singleInt.begin();
    
    if(systCats.size() >=1){
      for(size_t sy = 0; sy< systCats.size() ;++sy){
	string systCat = systCats.at(sy);
	obj_systCats[namelabel].push_back(systCat);
      }
    }
    for(size_t sc = 0; sc< categories.size() ;++sc){
      string category = categories.at(sc);
      if(scanCuts.size() >=1){
	for(size_t sccut = 0; sccut< scanCuts.size() ;++sccut){
	  string scanCut = scanCuts.at(sccut);
	  obj_cats[namelabel].push_back(category+"_"+scanCut);
	}
      }
      if(systCats.size() >=1){
      	for(size_t sy = 0; sy< systCats.size() ;++sy){
       	  string syCat = systCats.at(sy);
	  if(scanCuts.size() >=1){
	    for(size_t sccut = 0; sccut< scanCuts.size() ;++sccut){
	      string scanCut = scanCuts.at(sccut);
	      obj_cats[namelabel].push_back(category+syCat+"_"+scanCut);
	    }
	  }
	  obj_cats[namelabel].push_back(category+syCat);
	}
      }
      obj_cats[namelabel].push_back(category);
    }
    for(size_t sccut = 0; sccut< scanCuts.size() ;++sccut){
      string scanCut = scanCuts.at(sccut);
      obj_scanCuts[namelabel].push_back(scanCut);
    }    
    stringstream max_instance_str;
    max_instance_str<<maxI;
    max_instances[namelabel]=maxI;
    string nameobs = namelabel;
    string prefix = nameprefix;
    
    //    if(changeJECs){
      edm::InputTag jetKeys_ = iConfig.getParameter<edm::InputTag >("jetKeysAK4CHS");
      edm::InputTag muKeys_ = iConfig.getParameter<edm::InputTag >("muonKeys");
      t_jetKeys_ = consumes<std::vector<std::vector<int> > > (jetKeys_);
      t_muKeys_ = consumes<std::vector<std::vector<int> > > (muKeys_);
  
      //    }
    if(saveNoCat) trees["noSyst"]->Branch((nameobs+"_size").c_str(), &sizes[nameobs]);
    for(size_t sc = 0; sc< obj_cats[namelabel].size() ;++sc){
      string category = obj_cats[namelabel].at(sc);
      trees["noSyst"]->Branch((nameobs+category+"_size").c_str(), &sizes[nameobs+category]);
    }
    
    for (;itF != variablesFloat.end();++itF){
      
      string name=itF->instance()+"_"+itF->label();
      string nameinstance=itF->instance();
      string nameshort=itF->instance();
      
      string nametobranch = makeBranchName(namelabel,prefix,nameinstance);
      name = nametobranch;
      nameshort = nametobranch;
      //      if(saveBaseVariables|| isInVector(toSave,itF->instance())) trees["noSyst"]->Branch(nameshort.c_str(), &vfloats_values[name],(nameshort+"["+max_instance_str.str()+"]/F").c_str());
      //      if (saveNoCat && (saveBaseVariables || isInVector(toSave, extravars.at(addv)) || isInVector(toSave, "allExtra") ) )trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+nameobs+"_size"+"]/F").c_str());
      //if(saveBaseVariables|| isInVector(toSave,extravars.at(addv)) || isInVector(toSave,"allExtra")) trees["noSyst"]->Branch(namecat.c_str(), &vfloats_values[namecat],(namecat+"["+nameobs+category+"_size"+"]/F").c_str());	for(size_t sc = 0; sc< obj_cats[namelabel].size() ;++sc){
      if(saveNoCat && (saveBaseVariables|| isInVector(toSave,itF->instance()))) trees["noSyst"]->Branch(nameshort.c_str(), &vfloats_values[name],(nameshort+"["+nameobs+"_size"+"]/F").c_str());
      names.push_back(name);
      obj_to_floats[namelabel].push_back(name);
      obs_to_obj[name] = nameobs;
      obj_to_pref[nameobs] = prefix;
      //      cout << " branching name "<< name<< " for obs "<< nameobs << " instance "<< nameinstance << endl;

      t_floats[ name ] = consumes< std::vector<float> >( *itF );
      
      for(size_t sc = 0; sc< obj_cats[namelabel].size() ;++sc){
	string category = obj_cats[namelabel].at(sc);
	string nametobranchcat = makeBranchNameCat(namelabel,category,prefix,nameinstance);
	string namecat = nametobranchcat;
	nameshort = nametobranch;
	//	if(saveBaseVariables|| isInVector(toSave,itF->instance())) trees["noSyst"]->Branch(namecat.c_str(), &vfloats_values[namecat],(namecat+"["+max_instance_str.str()+"]/F").c_str());
	if(saveBaseVariables|| isInVector(toSave,itF->instance())) trees["noSyst"]->Branch(namecat.c_str(), &vfloats_values[namecat],(namecat+"["+nameobs+category+"_size"+"]/F").c_str());
	//	cout << " branching category "<< category<< " for label "<< namelabel<< " name "<<namecat<< endl;
      }
    }
  
    for (;itI != variablesInt.end();++itI){
      string name=itI->instance()+"_"+itI->label();
      string nameshort=itF->instance();
      string nametobranch = makeBranchName(namelabel,prefix,nameshort);
      name = nametobranch;
      nameshort = nametobranch;
      //      if(saveNoCat && (saveBaseVariables|| isInVector(toSave,itI->instance())) ) trees["noSyst"]->Branch(nameshort.c_str(), &vints_values[name],(nameshort+"["+max_instance_str.str()+"]/I").c_str());
      if(saveNoCat && (saveBaseVariables|| isInVector(toSave,itI->instance())) ) trees["noSyst"]->Branch(nameshort.c_str(), &vints_values[name],(nameshort+"["+nameobs+"_size"+"]/I").c_str());
      for(size_t sc = 0; sc< obj_cats[namelabel].size() ;++sc){
	string category = obj_cats[namelabel].at(sc);
	string nametobranchcat = makeBranchNameCat(namelabel,category,prefix,nameshort);
	string namecat = nametobranchcat;
	nameshort = nametobranch;
	if(saveBaseVariables|| isInVector(toSave,itF->instance())) trees["noSyst"]->Branch(namecat.c_str(), &vfloats_values[namecat],(namecat+"["+nameobs+category+"_size"+"]/I").c_str());
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
	//	if (saveNoCat && (saveBaseVariables || isInVector(toSave, extravars.at(addv)) || isInVector(toSave, "allExtra") ) )trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+max_instance_str.str()+"]/F").c_str());
	if (saveNoCat && (saveBaseVariables || isInVector(toSave, extravars.at(addv)) || isInVector(toSave, "allExtra") ) )trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+nameobs+"_size"+"]/F").c_str());
	//	  if(saveBaseVariables|| isInVector(toSave,extravars.at(addv)) || isInVector(toSave,"allExtra")) trees["noSyst"]->Branch(namecat.c_str(), &vfloats_values[namecat],(namecat+"["+nameobs+category+"_size"+"]/F").c_str());	
	for(size_t sc = 0; sc< obj_cats[namelabel].size() ;++sc){
	  string category = obj_cats[namelabel].at(sc);
	  string nametobranchcat = nameshortv+category+"_"+extravars.at(addv);
	  string namecat = nametobranchcat;
	  //	  cout << "extra var "<< extravars.at(addv)<<endl;
	  //	  cout << " namecat "<< namecat<< endl;
	  //	  if(saveBaseVariables|| isInVector(toSave,extravars.at(addv)) || isInVector(toSave,"allExtra")) trees["noSyst"]->Branch(namecat.c_str(), &vfloats_values[namecat],(namecat+"["+max_instance_str.str()+"]/F").c_str());
	  if(saveBaseVariables|| isInVector(toSave,extravars.at(addv)) || isInVector(toSave,"allExtra")) trees["noSyst"]->Branch(namecat.c_str(), &vfloats_values[namecat],(namecat+"["+nameobs+category+"_size"+"]/F").c_str());
	}


	obj_to_floats[namelabel].push_back(name);
	obs_to_obj[name] = nameobs;
	obj_to_pref[nameobs] = prefix;
      }
    }
    names.push_back(nameobs);
  

    //    cout << "size part: nameobs is  "<< nameobs<<endl;
    
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
    //    cout << " max instances top is "<< max_instances_top << " max_leading_jets_for_top "<< max_leading_jets_for_top << " max_instances[jets_label]  " <<max_instances[jets_label]<<endl;
    for(size_t addv = 0; addv < extravarstop.size();++addv){
      string name = nameshortv+"_"+extravarstop.at(addv);
      trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+mtop.str()+"]/F").c_str());
    }
    trees["noSyst"]->Branch((nameshortv+"_size").c_str(), &sizes[nameshortv]);
  }

  if(usePrunedGenParticles) {
    //gen muons
    string nameshortv= "genMuon";
    vector<string> extravarsgenmu = additionalVariables(nameshortv);
    double max_instances_genmu = max_instances[mu_label] ;
    //    double max_instances_genel = max_instances[mu_label] ;
    //    double max_instances_genb = max_instances[mu_label] ;
    //    trees["noSyst"]->Branch(, &vfloats_values[],(namecat+"["+max_instance_str.str()+"]/F").c_str());
    
    //    cout<<  "nextravarsmu" << extravarsgenmu.size()<<endl;
    stringstream mgenmu,mgenel,mgenb;

    mgenmu << max_instances_genmu;
    mgenel << max_instances_genmu;
    mgenb << max_instances_genmu;
    for(size_t addv = 0; addv < extravarsgenmu.size();++addv){
      string name = nameshortv+"_"+extravarsgenmu.at(addv);
      //      cout << " branching "<< name <<endl;
      trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+mgenmu.str()+"]/F").c_str());
    }
    trees["noSyst"]->Branch((nameshortv+"_size").c_str(), &sizes[nameshortv]);
    //gen electrons
    nameshortv= "genElectron";
    vector<string> extravarsgenel = additionalVariables(nameshortv);
    for(size_t addv = 0; addv < extravarsgenel.size();++addv){
      string name = nameshortv+"_"+extravarsgenel.at(addv);
      trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+mgenel.str()+"]/F").c_str());
    }
    trees["noSyst"]->Branch((nameshortv+"_size").c_str(), &sizes[nameshortv]);

    //gen b
    nameshortv= "genB";
    vector<string> extravarsgenb = additionalVariables(nameshortv);
    for(size_t addv = 0; addv < extravarsgenb.size();++addv){
      string name = nameshortv+"_"+extravarsgenb.at(addv);
      trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+mgenb.str()+"]/F").c_str());
    }
    trees["noSyst"]->Branch((nameshortv+"_size").c_str(), &sizes[nameshortv]);
  
  }
  
  if(useGenJets) {
    double max_instances_genjet = max_instances[jets_label] ;
    //    trees["noSyst"]->Branch(, &vfloats_values[],(namecat+"["+max_instance_str.str()+"]/F").c_str());

    stringstream mgenjet;
    mgenjet << max_instances_genjet;
                cout <<" using genJets: not implemented atm"<<endl;
    string nameshortv= "genJets";
    vector<string> extravarsgenjets = additionalVariables(nameshortv);
    for(size_t addv = 0; addv < extravarsgenjets.size();++addv){
      string name = nameshortv+"_"+extravarsgenjets.at(addv);
      trees["noSyst"]->Branch(name.c_str(), &vfloats_values[name],(name+"["+mgenjet.str()+"]/F").c_str());
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
    //    cout << " max instances top is "<< max_instances_top<< " max_leading_jets_for_top "<< max_leading_jets_for_top << " max_instances[jets_label]  " <<max_instances[jets_label]<<endl;
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
      //      std::cout<<"=====================sto riempendo il branch event number"<<std::endl;
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
    //trees[syst]= treesBase->CloneTree();
    trees[syst]->SetName((channel+"__"+syst).c_str());
    trees[syst]->SetTitle((channel+"__"+syst).c_str());
  }
  
  initTreeWeightHistory(useLHEWeights);
  
  
  if(JECVersion == "")JECVersion = "Spring16_25nsV10";
  string L1Name = JECVersion+"_MC_L1FastJet_AK4PFchs.txt"; //
  string L1RCName = JECVersion+"_MC_L1RC_AK4PFchs.txt"; 
  string L2Name = JECVersion+"_MC_L2Relative_AK4PFchs.txt";
  string L3Name = JECVersion+"_MC_L3Absolute_AK4PFchs.txt";
  string L2L3ResName = JECVersion+"_MC_L2L3Residual_AK4PFchs.txt";
  if(isData){
    L1Name   = JECVersion+"_DATA_L1FastJet_AK4PFchs.txt";
    L1RCName = JECVersion+"_DATA_L1RC_AK4PFchs.txt";  
    L2Name   = JECVersion+"_DATA_L2Relative_AK4PFchs.txt";
    L3Name   = JECVersion+"_DATA_L3Absolute_AK4PFchs.txt";
    L2L3ResName = JECVersion+"_DATA_L2L3Residual_AK4PFchs.txt";
  }
  if (isData)doResol=false;
  //  string jecDir="JEC/";
  string jecDir="./";
  jecParsL1  = new JetCorrectorParameters(jecDir+L1Name);
  //  jecParsL1RC  = new JetCorrectorParameters(L1RCName);
  jecParsL2  = new JetCorrectorParameters(jecDir+L2Name);
  jecParsL3  = new JetCorrectorParameters(jecDir+L3Name);
  jecParsL2L3Residuals  = new JetCorrectorParameters(jecDir+L2L3ResName);
  jecPars.push_back(*jecParsL1);
  jecPars.push_back(*jecParsL2);
  jecPars.push_back(*jecParsL3);
  if(isData)jecPars.push_back(*jecParsL2L3Residuals);

  jecParsL1_vect.push_back(*jecParsL1);
  
  jecCorr = new FactorizedJetCorrector(jecPars);
  jecCorr_L1 = new FactorizedJetCorrector(jecParsL1_vect);
  jecUnc  = new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecDir+"Spring16_23Sep2016BCDV2_DATA_UncertaintySources_AK4PFchs.txt", "Total")));
  filename_cmva="btagging_cmva.root";
  file_cmva= TFile::Open(filename_cmva.c_str());
  cmvaeffbt = new Weights(file_cmva,"b__tight");
  cmvaeffbl = new Weights(file_cmva,"b__loose");
  cmvaeffbm = new Weights(file_cmva,"b__medium");

  cmvaeffct = new Weights(file_cmva,"c__tight");
  cmvaeffcl = new Weights(file_cmva,"c__medium");
  cmvaeffcm = new Weights(file_cmva,"c__loose");

  cmvaeffot = new Weights(file_cmva,"other__tight");
  cmvaeffol = new Weights(file_cmva,"other__medium");
  cmvaeffom = new Weights(file_cmva,"other__loose");


  //  calib = new BTagCalibration("CSVv2", "CSVv2_ichep.csv");
  calib = new BTagCalibration("CSVv2", "CSVv2_Moriond17_B_H.csv");
  readerCSVLoose = new BTagCalibrationReader(BTagEntry::OP_LOOSE, "central", {"up", "down"});      
  readerCSVLoose->load(*calib, BTagEntry::FLAV_B,   "comb");
  readerCSVLoose->load(*calib, BTagEntry::FLAV_C,   "comb");
  readerCSVLoose->load(*calib, BTagEntry::FLAV_UDSG,   "incl");

  readerCSVMedium = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"});      
  readerCSVMedium->load(*calib, BTagEntry::FLAV_B,   "comb");
  readerCSVMedium->load(*calib, BTagEntry::FLAV_C,   "comb");
  readerCSVMedium->load(*calib, BTagEntry::FLAV_UDSG,   "incl");

  readerCSVTight = new BTagCalibrationReader(BTagEntry::OP_TIGHT, "central", {"up", "down"});      
  readerCSVTight->load(*calib, BTagEntry::FLAV_B,   "comb");
  readerCSVTight->load(*calib, BTagEntry::FLAV_C,   "comb");
  readerCSVTight->load(*calib, BTagEntry::FLAV_UDSG,   "incl");

  readerCSVReshape = new BTagCalibrationReader(BTagEntry::OP_RESHAPING, "central", {"up_jes", "down_jes"});      
  readerCSVReshape->load(*calib, BTagEntry::FLAV_B,   "iterativefit");
  readerCSVReshape->load(*calib, BTagEntry::FLAV_C,   "iterativefit");
  readerCSVReshape->load(*calib, BTagEntry::FLAV_UDSG,   "iterativefit");



  //new combined MVAv2 tagger: pay load taken from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80X
  //cMVAv2_ichep.csv
  
  calib_cmvav2 = new BTagCalibration("CMVAv2", "cMVAv2_Moriond17_B_H.csv");
  readerCMVALoose = new BTagCalibrationReader(BTagEntry::OP_LOOSE, "central", {"up", "down"});      
  readerCMVALoose->load(*calib_cmvav2, BTagEntry::FLAV_B,   "ttbar");
  readerCMVALoose->load(*calib_cmvav2, BTagEntry::FLAV_C,   "ttbar");
  readerCMVALoose->load(*calib_cmvav2, BTagEntry::FLAV_UDSG,   "incl");

  readerCMVAMedium = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"});      
  readerCMVAMedium->load(*calib_cmvav2, BTagEntry::FLAV_B,   "ttbar");
  readerCMVAMedium->load(*calib_cmvav2, BTagEntry::FLAV_C,   "ttbar");
  readerCMVAMedium->load(*calib_cmvav2, BTagEntry::FLAV_UDSG,   "incl");

  readerCMVATight = new BTagCalibrationReader(BTagEntry::OP_TIGHT, "central", {"up", "down"});      
  readerCMVATight->load(*calib_cmvav2, BTagEntry::FLAV_B,   "ttbar");
  readerCMVATight->load(*calib_cmvav2, BTagEntry::FLAV_C,   "ttbar");
  readerCMVATight->load(*calib_cmvav2, BTagEntry::FLAV_UDSG,   "incl");

  readerCMVAReshape = new BTagCalibrationReader(BTagEntry::OP_RESHAPING, "central", {"up_jes", "down_jes"});      
  readerCMVAReshape->load(*calib_cmvav2, BTagEntry::FLAV_B,   "iterativefit");
  readerCMVAReshape->load(*calib_cmvav2, BTagEntry::FLAV_C,   "iterativefit");
  readerCMVAReshape->load(*calib_cmvav2, BTagEntry::FLAV_UDSG,   "iterativefit");

  // reader.load(...)     // for FLAV_C
  // reader.load(...)     // for FLAV_UDSG
  
  isFirstEvent = true;
  doBTagSF= true;
  if(isData)doPU= false;
  
  season = "Summer11";
    
  distr = "pileUpDistr" + season + ".root";
  //  cout << " end build "<<endl;
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
  
  if(usePrunedGenParticles)iEvent.getByToken(t_genParticles_,genParts);
  if(useGenJets)iEvent.getByToken(t_genJets_,genJets);
  
  //  cout << "number of gen particles"<<genParts->size()<<endl;
  if(useLHE){
    iEvent.getByToken(t_lhes_, lhes);
  }
  if(addLHAPDFWeights){
    iEvent.getByToken(t_genprod_, genprod);
  }
  
  //  if(changeJECs){
    iEvent.getByToken(t_jetKeys_, jetKeys);
    iEvent.getByToken(t_muKeys_, muKeys);
    //  }
  
  vector<TLorentzVector> genlep,gentop,genantitop,genw,gennu,genz,gena;
  vector<TLorentzVector> pare,parebar,parmu,parmubar,parnumu,parnumubar,partau,partaubar,parnue,parnuebar,parnutau,parnutaubar,parz,parw;
  
  //Parton-level info
  getPartonW=false;
  getParticleWZ=false;
  //  getPartonTop=false;
  doWReweighting=false;
  //doTopReweighting=false;
  //  useLHE=false;
  //  useLHEWeights=false;
  //  useLHE=false;
  //  useLHEWeights=false;

  if(useLHEWeights){
      getEventLHEWeights();
    }
  
  
  if((getPartonW || getPartonTop || doWReweighting || doTopReweighting )&& useLHE){
    //    if(!useLHE)return;
    genlep.clear();
    gentop.clear();
    genantitop.clear();
    genw.clear();
    genz.clear();
    gena.clear();
    gennu.clear();
    pare.clear();parebar.clear();parmu.clear();parmubar.clear();parnumu.clear();parnumubar.clear();parnue.clear();parnuebar.clear();parnutau.clear();parnutaubar.clear();parz.clear();parw.clear();

    if(getParticleWZ){
      iEvent.getByToken(genPartID, partID);
      iEvent.getByToken(genPartStatus, partStatus);
      iEvent.getByToken(genPartMom0ID, partMom0ID);
      iEvent.getByToken(genPartPt, partPt);
      iEvent.getByToken(genPartPhi, partPhi);
      iEvent.getByToken(genPartEta, partEta);
      iEvent.getByToken(genPartE, partE);
    
    }
    
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
      //      float mass = lhes->hepeup().PUP[i][4];
      
      //      if(abs (id) == 24 )  cout << " px is"<< px << " py "<< py << " pz "<< pz << " e "<<energy<<endl;
    
      TLorentzVector vec;
      math::XYZTLorentzVector part = math::XYZTLorentzVector(px, py, pz, energy);
      float pt = part.pt();
      float phi = part.phi();
      float eta = part.eta();
      
      
      if(pt>0){
	vec.SetPtEtaPhiE(pt, eta, phi, energy);
	//	if(abs (id) == 24 ) cout << " pt is"<< pt << " phi "<< phi << " eta "<< eta <<endl;
	
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
    topCharge=0;
    //    cout <<  " top charge"<<topCharge<<" gentop size "<<gentop.size()<<" antitop size"<< genantitop.size()<<endl;
    if(gentop.size()==1){
      float_values["Event_T_size"]=1;
      topCharge=1;
    }
    if(gentop.size()==0&&genantitop.size()==1){
      float_values["Event_Tbar_size"]= 1;
      topCharge=-1; 
    }
    //    cout <<  " top charge after"<<topCharge<<endl;
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
	float_values["Event_Z_QCD_Weight"]= wweight;
      }
      else (float_values["Event_Z_QCD_Weight"]=1.0);
    }
    if((getPartonW || doWReweighting ) ) {
      if(gena.size()==1){       
      float_values["Event_a_Pt"]= gena.at(0).Pt();
      float_values["Event_a_Eta"]= gena.at(0).Eta();
      float_values["Event_a_Phi"]= gena.at(0).Phi();
      float_values["Event_a_E"]= gena.at(0).Energy();
      float_values["Event_a_Mass"]= gena.at(0).M();	
      
      double ptW = gena.at(0).Pt();
      double wweight = getAPtWeight(ptW);			
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
      }
      else {(float_values["Event_T_Weight"]=1.0);
	(float_values["Event_T_Ext_Weight"]=1.0);}
    }
    if( (getParticleWZ)) {
      for(size_t p=0;p<partID->size();++p){
 	if(partID->at(p)==partMom0ID->at(p) && isEWKID(partID->at(p)) ){
	  TLorentzVector vec;
	  vec.SetPtEtaPhiE(partPt->at(p),partEta->at(p),partPhi->at(p),partE->at(p));
	  if(partID->at(p)==11  && partStatus->at(p)==1){
	    pare.push_back(vec);	  }
	  if(partID->at(p)==-11  && partStatus->at(p)==1){
	    parebar.push_back(vec);	  }
	  if(partID->at(p)==13 && partStatus->at(p)==1){
	    parmu.push_back(vec);	  }
	  if(partID->at(p)==-13 && partStatus->at(p)==1){
	    parmubar.push_back(vec);	  }

	  if(partID->at(p)==15  && partStatus->at(p)==2){
	    partau.push_back(vec);  }
	  if(partID->at(p)==-15  && partStatus->at(p)==2){
	    partaubar.push_back(vec);  }
	  
	  if(partID->at(p)==12  && partStatus->at(p)==1){
	    parnue.push_back(vec);	  }
	  if(partID->at(p)==14  && partStatus->at(p)==1){
	    parnumu.push_back(vec);	  }
	  if(partID->at(p)==16  && partStatus->at(p)==1){
	    parnutau.push_back(vec);	  }

	  if(partID->at(p)==-12  && partStatus->at(p)==1){
	    parnuebar.push_back(vec);	  }
	  if(partID->at(p)==-14  && partStatus->at(p)==1){
	    parnumubar.push_back(vec);	  }
	  if(partID->at(p)==-16  && partStatus->at(p)==1){
	    parnutaubar.push_back(vec);	  }
	  
	  if(abs(partID->at(p))==23 && partStatus->at(p)==62){
	    //	    getPtWeight
	    double wweight = getZEWKPtWeight(partPt->at(p));
	    float_values["Event_Z_EW_Weight"]= wweight;
	  }
	  if(abs(partID->at(p))==24 && partStatus->at(p)==62){
	    double wweight = getWEWKPtWeight(partPt->at(p));
	    float_values["Event_W_EW_Weight"]= wweight;
	  }
	  //	  cout << " in loop 2: p= "<<p<<endl;      

	}
      }
      //      cout << " parmusize "<< parmu.size()<< " parnumubarsize"<<parnumubar.size()<<endl;
      //      cout << " paresize "<< pare.size()<< " parnuebarsize"<<parnuebar.size()<<endl;
      //      cout << " partausize "<< partau.size()<< " parnutaubarsize"<<parnutaubar.size()<<endl;
      //      cout << " parnumusize "<< parnumu.size()<< " parmubarsize"<<parmubar.size()<<endl;
      //      cout << " parnuesize "<< parnue.size()<< " parebarsize"<<parebar.size()<<endl;
      //      cout << " parnutausize "<< parnutau.size()<< " partaubarsize"<<partaubar.size()<<endl;
      //      cout << "parw size "<< parw.size()<<endl;
      //      cout << "parz size "<< parz.size()<<endl;
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
	//	cout << "parw size "<< parw.size()<<endl;
	//	cout << " w is one new val is "<<getWEWKPtWeight(parw.at(0).Pt())<<endl;
	float_values["Event_W_EW_Weight"]= getWEWKPtWeight(parw.at(0).Pt());
      }
      if(   float_values["Event_W_QCD_Weight"] ==1 &&parw.size()>0 )    float_values["Event_W_QCD_Weight"]= getWPtWeight(parw.at(0).Pt());
    }
    //    cout << " after loop "<<endl;
    //    float_values["Event_W_EW_Weight"]=1.0;//*float_values["Event_W_QCD_Weight"];
    //    float_values["Event_Z_EW_Weight"]=1.0;//*float_values["Event_W_QCD_Weight"];
    float_values["Event_W_Weight"]= float_values["Event_W_EW_Weight"]*float_values["Event_W_QCD_Weight"];
    float_values["Event_Z_Weight"]= float_values["Event_Z_EW_Weight"]*float_values["Event_Z_QCD_Weight"];
    //    cout << " after loop 2 "<<endl;


  }

  if(usePrunedGenParticles){
    int nmu=0;
    int nel=0;
    for (size_t g = 0; g<genParts->size();++g){
      bool isPromptFinalState=genParts->at(g).isPromptFinalState();
      bool isDirectPromptTauDecayProductFinalState=genParts->at(g).isDirectPromptTauDecayProductFinalState();
      if(!(isPromptFinalState||isDirectPromptTauDecayProductFinalState))continue;
      else {
	float id = genParts->at(g).pdgId();
	if(fabs(id) ==13 || fabs(id)==11){
	float pt = genParts->at(g).pt();
	float eta = genParts->at(g).eta();
	float phi = genParts->at(g).phi();
	float energy = genParts->at(g).energy();
	float charge = genParts->at(g).charge();
	
	  if(fabs(id)==13 ){
	    if(nmu< max_instances[mu_label]){
	      vfloats_values["genMuon_Pt"][nmu]= pt;
	      vfloats_values["genMuon_Eta"][nmu]= eta;
	      vfloats_values["genMuon_Phi"][nmu]= phi;
	      vfloats_values["genMuon_E"][nmu]= energy;
	      vfloats_values["genMuon_Charge"][nmu]= charge;
	      ++nmu;
	    }
	  }
	  if(fabs(id)==11){
	    if(nel< max_instances[mu_label]){
	      vfloats_values["genElectron_Pt"][nel]= pt;
	      vfloats_values["genElectron_Eta"][nel]= eta;
	      vfloats_values["genElectron_Phi"][nel]= phi;
	      vfloats_values["genElectron_E"][nel]= energy;
	      vfloats_values["genElectron_Charge"][nel]= charge;
	    ++nel;
	    }
	  }
	}
      }
    }
    sizes["genMuon"]=nmu;    
    sizes["genElectron"]=nel;
  }

  if(doPU){
    iEvent.getByToken(t_ntrpu_,ntrpu);
    int nTruePV=*ntrpu;
    //      cout<< "ntpv is "<< nTruePV<<endl;
    float_values["Event_nTruePV"]=(float)(nTruePV);
  }

      
  trees["WeightHistory"]->Fill();
  
  //Part 3: filling the additional variables


  //boosted part

  iEvent.getByToken(jetAK8topSubjetIndex0, ak8jetSubjetIndex0);
  iEvent.getByToken(jetAK8topSubjetIndex1, ak8jetSubjetIndex1);
  iEvent.getByToken(jetAK8topSubjetIndex2, ak8jetSubjetIndex2);
  iEvent.getByToken(jetAK8topSubjetIndex3, ak8jetSubjetIndex3);

  //  if(isFirstEvent)
  //  cout << "jec version is "<<JECVersion<<" Change JECs? " << changeJECs<<endl;
  //Part 0: trigger preselection
  
      
  if(useTriggers){
    //Trigger names are retrieved from the run tree
    iEvent.getByToken(t_triggerBits_,triggerBits );
    iEvent.getByToken(t_triggerPrescales_,triggerPrescales );
    bool triggerOr = getEventTriggers();
    
    if(isFirstEvent){

      for(size_t bt = 0; bt < triggerNamesR->size();++bt){
	std::string tname = triggerNamesR->at(bt);
	// cout << "trigger test tname "<< tname << " passes "<< triggerBits->at(bt)<< endl;
      }
    }
    
    if(cutOnTriggers && !triggerOr) return;
  }

  if(useMETFilters){
    iEvent.getByToken(t_metBits_,metBits );
    //    iEvent.getByToken(t_metNames_,metNames );
    if(isFirstEvent){
      for(size_t bt = 0; bt < metNames->size();++bt){
	std::string tname = metNames->at(bt);
	//cout << "test tname "<< tname << " passes "<< metBits->at(bt)<< endl;
      }
    }
    getMETFilters();
  }
  if(changeJECs || recalculateEA){
    iEvent.getByToken(t_Rho_ ,rho);
    Rho = *rho; 
  }
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
    // singleDouble = itPsets->template getParameter<std::vector<edm::InputTag> >("singleD"); 
    singleInt = itPsets->template getParameter<std::vector<edm::InputTag> >("singleI"); 
    std::vector<edm::InputTag >::const_iterator itF = variablesFloat.begin();
    std::vector<edm::InputTag >::const_iterator itI = variablesInt.begin();
    std::vector<edm::InputTag >::const_iterator itsF = singleFloat.begin();
    //std::vector<edm::InputTag >::const_iterator itsD = singleDouble.begin();
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
      //Getting the names for the floats of the object
      float tmp =1.0;
      iEvent.getByToken(t_floats[name] ,h_floats[name]);
      for (size_t fi = 0;fi < maxInstance ;++fi){
	if(fi <h_floats[name]->size()){tmp = h_floats[name]->at(fi);}
	else { tmp = -9999.; }
	vfloats_values[name][fi]=tmp;
	for (size_t sc=0;sc< obj_cats[namelabel].size();++sc){
	  string category = obj_cats[namelabel].at(sc);
	  string namecat = makeBranchNameCat(namelabel,category,nameprefix,varname);
	  vfloats_values[namecat][fi]=-9999.;
	}
      }
      //Getting the sizes for the object
      sizes[namelabel]=h_floats[name]->size();
      initCategoriesSize(namelabel);
      //      initSystsSize(namelabel);
    }
    //    std::cout << " checkpoint floats"<<endl;
    for (;itD != variablesDouble.end();++itD){
      string varname=itD->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
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
    
    //Single floats/ints
    for (;itsF != singleFloat.end();++itsF){
      string varname=itsF->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
      iEvent.getByToken(t_float[name],h_float[name]);
      float_values[name]=*h_float[name];
    }

    for (;itsD != singleDouble.end();++itsD){
      string varname=itsD->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
      iEvent.getByToken(t_double[name] ,h_double[name]);
      double_values[name]=*h_double[name];
    }
    for (;itsI != singleInt.end();++itsI){
      string varname=itsI->instance();
      string name = makeBranchName(namelabel,nameprefix,varname);
      iEvent.getByToken(t_int[name],h_int[name]);
      int_values[name]=*h_int[name];
    }
  }



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

  //  cout << "pres"<<endl;
  for (size_t s = 0; s< systematics.size();++s){

    //    cout << "syst"<< systematics.at(s)<<endl;

    int nb=0,nc=0,nudsg=0;
    //    int ncsvl_tags=0,ncsvt_tags=0 
    int ncsvm_tags=0;//, njets_tottag=0;
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

    int lepidx=0;
    int bjetidx=0;
    
    
    initCategoriesSize(jets_label);
    initCategoriesSize(mu_label);
    initCategoriesSize(ele_label);

    //Photons
    //    initCategoriesSize(photon_label);

    //    cout << " prephoton "<<endl;
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
      //      cout << " sieie " << sieie << " hoe "<< hoe << " pho_isoCea pho_isoCea "<< " pho_isoNea "<< pho_isoNea << " cut isoNea "<< (2.57+exp(0.0044*pt +0.5809) ) << " pho_isoPea "<< pho_isoPea << " cut isoPea "<< (1.92+0.0043) << " eta "<< eta << " abseta "<< abseta <<" pt "<< pt <<endl;
      if(isBarrel){
	//	cout << " isbarrel "<<endl;
	if( sieie < 0.0103 &&   hoe < 0.05 &&   pho_isoCea < 2.44 &&   pho_isoNea < (2.57+exp(0.0044*pt +0.5809) ) &&   pho_isoPea < (1.92+0.0043*pt ) )vfloats_values[photon_label+"_isLooseSpring15"][ph]=1.0;
	if( sieie < 0.01 &&   hoe < 0.05 &&   pho_isoCea < 1.31 &&   pho_isoNea < (0.60+exp(0.0044*pt +0.5809) ) &&   pho_isoPea < (1.33+0.0043*pt ) )vfloats_values[photon_label+"_isMediumSpring15"][ph]=1.0;
	if( sieie < 0.01 &&   hoe < 0.05 &&   pho_isoCea < 0.91 &&   pho_isoNea < (0.33+exp(0.0044*pt +0.5809) ) &&   pho_isoPea < (0.61+0.0043*pt ) )vfloats_values[photon_label+"_isTightSpring15"][ph]=1.0;
      }
      if(isEndcap){
	//	cout << " isendcap "<<endl;
	if( sieie < 0.0277 &&   hoe < 0.05 &&   pho_isoCea < 1.84 &&   pho_isoNea < (4.00+exp(0.0040*pt +0.9402) ) &&   pho_isoPea < (1.92+0.0043*pt ) )vfloats_values[photon_label+"_isLooseSpring15"][ph]=1.0;
	if( sieie < 0.0267 &&   hoe < 0.05 &&   pho_isoCea < 1.25 &&   pho_isoNea < (1.65+exp(0.0040*pt +0.9402) ) &&   pho_isoPea < (1.33+0.0043*pt ) )vfloats_values[photon_label+"_isMediumSpring15"][ph]=1.0;
	if( sieie < 0.0267 &&   hoe < 0.05 &&   pho_isoCea < 0.65 &&   pho_isoNea < (0.93+exp(0.0040*pt +0.9402) ) &&   pho_isoPea < (0.61+0.0043*pt ) )vfloats_values[photon_label+"_isTightSpring15"][ph]=1.0;
      }
    
    }
    
    //Muons
    //    cout << "testing muons now"<<endl;
    //    initCategoriesSize(mu_label);
    for(int mu = 0;mu < max_instances[mu_label] ;++mu){
      string pref = obj_to_pref[mu_label];
      float isTight = vfloats_values[makeName(mu_label,pref,"IsTightMuon")][mu];
      float isLoose = vfloats_values[makeName(mu_label,pref,"IsLooseMuon")][mu];
      float isSoft = vfloats_values[makeName(mu_label,pref,"IsSoftMuon")][mu];
      
      float pt = vfloats_values[makeName(mu_label,pref,"Pt")][mu];
      float eta = vfloats_values[makeName(mu_label,pref,"Eta")][mu];
      float phi = vfloats_values[makeName(mu_label,pref,"Phi")][mu];
      float energy = vfloats_values[makeName(mu_label,pref,"E")][mu];
      float iso = vfloats_values[makeName(mu_label,pref,"Iso04")][mu];
      
      float muCharge = vfloats_values[makeName(mu_label,pref,"Charge")][mu];
      
      
      //if(isTight>0 && pt> 24 && abs(eta) < 2.4 /*&& iso <0.15 Isolation added afterwards*/){//UCL Selection
      if(isTight>0 && pt> 26 && fabs(eta) < 2.4 /*&& iso <0.15 Isolation added afterwards*/){//NA Selection
 	

	if(iso<0.06){// 2015 Selection
	  ++float_values["Event_nTightMuons"];
	  TLorentzVector muon;
	  muon.SetPtEtaPhiE(pt, eta, phi, energy);
	  muons.push_back(muon);
	  leptons.push_back(muon);
	  flavors.push_back(13);
	  leptonsCharge.push_back(muCharge);
	  if(isInVector(obj_cats[mu_label],"Tight")){
	    fillCategory(mu_label,"Tight",mu,float_values["Event_nTightMuons"]-1);
	  }
	  ++lepidx;
	}
	if(obj_scanCuts[mu_label].size()>=1) {
	  fillScanCuts(mu_label,"Tight",mu);
	}
      }

      if(isInVector(obj_cats[mu_label],"Tight")){
	sizes[mu_label+"Tight"]=(int)float_values["Event_nTightMuons"];
      }
      
      if(isInVector(obj_cats[mu_label],"TightAntiIso")){
	sizes[mu_label+"TightAntiIso"]=(int)float_values["Event_nTightAntiIsoMuons"];
      }
      
      if(isLoose>0 && pt> 10 && fabs(eta) < 2.4 ){//NA Selection
	if(isInVector(obj_cats[mu_label],"Loose")){
	  if(iso<0.2){
	    ++float_values["Event_nLooseMuons"];
	      fillCategory(mu_label,"Loose",mu,float_values["Event_nLooseMuons"]-1);
	  }
	  if(obj_scanCuts[mu_label].size()>=1) {
	    fillScanCuts(mu_label,"Loose",mu);
	  }
	}
      }
      if(isInVector(obj_cats[mu_label],"Loose")){
	sizes[mu_label+"Loose"]=(int)float_values["Event_nLooseMuons"];
      }
      
      if(isSoft>0 && pt> 30 && abs(eta) < 2.4){
	++float_values["Event_nSoftMuons"]; 
      }
      if(isLoose>0 && pt > 10){
	TLorentzVector muon;
	muon.SetPtEtaPhiE(pt, eta, phi, energy);
	loosemuons.push_back(muon);
      }
    }
    //    cout << "testing electrons now"<<endl;
    //    cout << " test 1 "<<endl;
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
      //float scEta = vfloats_values[makeName(ele_label,pref,"scEta")][el];
      float scEta = vfloats_values[makeName(ele_label,pref,"SCEta")][el];
      float phi = vfloats_values[makeName(ele_label,pref,"Phi")][el];
      float energy = vfloats_values[makeName(ele_label,pref,"E")][el];      
      float iso = vfloats_values[makeName(ele_label,pref,"Iso03")][el];
      float eldz = vfloats_values[makeName(ele_label,pref,"Dz")][el];
      float eldxy = vfloats_values[makeName(ele_label,pref,"Dxy")][el];

      float elCharge = vfloats_values[makeName(ele_label,pref,"Charge")][el];

      bool passesDRmu = true;
      bool passesTightCuts = false;
      bool passesTightAntiIsoCuts = false;

      if(fabs(scEta)<=1.479){
	//passesTightCuts = ( isTight >0.0 /*&& iso < 0.0588 */) && (fabs(eldz) < 0.10) && (fabs(eldxy) <0.05 ) ;
	passesTightCuts = (isTight >0.0) && (fabs(eldz) < 0.10) && (fabs(eldxy) <0.05 ) && (fabs(scEta)<1.4442 || fabs(scEta)>1.5660);
	passesTightAntiIsoCuts = isTight >0.0 && iso > 0.0588 ;

      } //is barrel electron
      //if( ( fabs(scEta)>1.479 && fabs(scEta)<2.5 ) && ( (fabs(eldz) < 0.20) && (fabs(eldxy) < 0.10) ) ){
      if( ( fabs(scEta)>1.479 && fabs(scEta)<2.5 ) && ( (fabs(eldz) < 0.20) && (fabs(eldxy) < 0.10) ) && (fabs(scEta)<1.4442 || fabs(scEta)>1.5660) ){
	passesTightCuts = isTight >0.0 /*&& iso < 0.0571*/ ;
	passesTightAntiIsoCuts = isTight >0.0 && iso > 0.0571 ;
      }

      if(pt> 30 && fabs(eta) < 2.1 && ((fabs(scEta)<1.4442 || fabs(scEta)>1.5660))){
	TLorentzVector ele;
	ele.SetPtEtaPhiE(pt, eta, phi, energy);	
	double minDR=999;
	double deltaR = 999;
	for (size_t m = 0; m < (size_t)loosemuons.size(); ++m){
	  deltaR = ele.DeltaR(loosemuons[m]);
	  minDR = min(minDR, deltaR);
	  minDR=999;
	}
	if(!loosemuons.size()) minDR=999;
	if(minDR>0.1){ 
	  if(passesTightCuts){
	    electrons.push_back(ele); 
	    flavors.push_back(11);
	    leptons.push_back(ele);
	    leptonsCharge.push_back(elCharge);
	    ++float_values["Event_nTightElectrons"];
	    ++lepidx;
	    if(isInVector(obj_cats[ele_label],"Tight")){
	      fillCategory(ele_label,"Tight",el,float_values["Event_nTightElectrons"]-1);
	    }
	  }
	  if(passesTightAntiIsoCuts){
	    ++float_values["Event_nTightAntiIsoElectrons"];
	    if(isInVector(obj_cats[ele_label],"TightAntiIso")){
	      fillCategory(ele_label,"TightAntiIso",el,float_values["Event_nTightAntiIsoElectrons"]-1);
	    }
	  }
	}
	else {passesDRmu = false;}
      }
      if(isInVector(obj_cats[ele_label],"Tight")){
	sizes[ele_label+"Tight"]=(int)float_values["Event_nTightElectrons"];
      }
      if(isInVector(obj_cats[ele_label],"TightAntiIso")){
	sizes[ele_label+"TightAntiIso"]=(int)float_values["Event_nTightAntiIsoElectrons"];
      }

      if(isLoose>0 && pt> 30 && fabs(eta) < 2.5 && ((fabs(scEta)<1.4442 || fabs(scEta)>1.5660))){
	++float_values["Event_nLooseElectrons"];

      }

      if(isMedium>0 && pt> 30 && fabs(eta) < 2.5 && ((fabs(scEta)<1.4442 || fabs(scEta)>1.5660)) ){
	++float_values["Event_nMediumElectrons"]; 
      }


      if(isVeto>0 && pt> 15 && fabs(eta) < 2.5 && ((fabs(scEta)<1.4442 || fabs(scEta)>1.5660))){
	  
      //if((fabs(scEta)<=1.479 && (iso<0.175)) || ((fabs(scEta)>1.479 && fabs(scEta)<2.5) && (iso<0.159))){
	
	if((fabs(scEta)<=1.479 && (fabs(eldz) < 0.10) && (fabs(eldxy) <0.05 )) || ((fabs(scEta)>1.479 && fabs(scEta)<2.5) && (fabs(eldz) < 0.20) && (fabs(eldxy) < 0.10) )){
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
    
    int nTightLeptons = float_values["Event_nTightMuons"]+float_values["Event_nTightElectrons"];
    int nTightAntiIsoLeptons = 0;
    for(size_t mc =0; mc < obj_cats[mu_label].size();++mc){
      string cat= obj_cats[mu_label].at(mc);
      if(cat.find("_Iso04_")!=std::string::npos && cat.find("GE")!=std::string::npos && cat.find("Tight")!=std::string::npos){
	//	cout<< " cat " << 
	nTightAntiIsoLeptons+=sizes[mu_label+cat];
      }
    } 
    //float_values["Event_nTightAntiIsoMuons"]+float_values["Event_nTightAntiIsoElectrons"];
    
    for(size_t l =0; l< leptons.size();++l){
      double lpt= leptons.at(l).Pt();
      if(lpt>maxpt){maxpt = lpt;firstidx=l;}
      
    }

    maxpt=0.0;
    for(size_t l =0; l< leptons.size();++l){
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
    //Lepton charge to be used later on:
    leadingLeptonCharge=0;
    //    cout<< "ntight leptons "<<nTightLeptons << " charge bef "<< leadingLeptonCharge<<endl;
    if(nTightLeptons>0){leadingLeptonCharge=float_values["Event_Lepton1_Charge"];
      //      cout<< "ntight leptons "<<nTightLeptons << " charge "<< leadingLeptonCharge<<endl;
    }
    if(nTightLeptons==0 && nTightAntiIsoLeptons>=1.0){
      if(float_values["muonsTightAntiIso_size"])leadingLeptonCharge=vfloats_values["muonsTightAntiIso_Charge"][0];
      if(!float_values["muonsTightAntiIso_size"] && 
	 float_values["electronsTightAntiIso_size"]) leadingLeptonCharge=vfloats_values["electronsTightAntiIso_Charge"][0]; //Prioritize anti-iso muons over anti-iso electrons
    }

    //Jets:
    double Ht=0;
    double corrBaseMetPx =0;
    double corrBaseMetPy =0;
    double corrMetPx =0;
    double corrMetPy =0;
    double corrMetT1Px =0;
    double corrMetT1Py =0;
    //    double corrMetPxNoHF =0;
    //    double corrMetPyNoHF =0;
    double DUnclusteredMETPx=0.0;
    double DUnclusteredMETPy=0.0;

    string prefm = obj_to_pref[met_label];
    float metZeroCorrPt = vfloats_values[makeName(met_label,prefm,"UncorrPt")][0];
    float metZeroCorrPhi = vfloats_values[makeName(met_label,prefm,"UncorrPhi")][0];
    float metZeroCorrY = metZeroCorrPt*sin(metZeroCorrPhi);
    float metZeroCorrX = metZeroCorrPt*cos(metZeroCorrPhi);

    //    cout << " max jets "<<max_instances[jets_label]<<endl;
    //met uncorrected;
    float metptunc = vfloats_values[makeName(met_label,prefm,"uncorPt")][0];
    float metphiunc = vfloats_values[makeName(met_label,prefm,"uncorPhi")][0];
    
    float metT1Py = metptunc*sin(metphiunc);
    float metT1Px = metptunc*cos(metphiunc);

    
    for (size_t ju = 0; ju < obj_systCats[jets_label].size();++ju){
      string systjet=obj_systCats[jets_label].at(ju);
      if (systjet.find("JES")!=std::string::npos || 
	  systjet.find("JER")!=std::string::npos ){
	//	vfloats_values[makeName(met_label,prefm,"CorrT1Px"+systjet)][0]=metT1Px;
	//vfloats_values[makeName(met_label,prefm,"CorrT1Py"+systjet)][0]=metT1Py;
	vfloats_values[makeName(met_label,prefm,"CorrT1Px"+systjet)][0]=0;
	vfloats_values[makeName(met_label,prefm,"CorrT1Py"+systjet)][0]=0;
      }
    }    
    //    if(isInVector(obj_cats[met_label],"CorrT1")){
    //  fillCategory(met_label,"CorrT1",0,sizes[met_label+"CorrT1"]);
    //}	    

    
    for (size_t ju = 0; ju < obj_cats[jets_label].size();++ju){
      if (obj_cats[jets_label].at(ju).find("JES")!=std::string::npos || 
	  obj_cats[jets_label].at(ju).find("JER")!=std::string::npos){  
	string catjet=obj_cats[jets_label].at(ju);
	//	cout << " ju "<< ju <<" --> catname "<< catjet << endl;
	passesJecCuts[catjet]=false;
      }
    }
    //    for(int j = 0;j < max_instances[jets_label] ;++j){
    for(int j = 0;j < max_instances[jets_label] ;++j){
      if(j > sizes[jets_label])continue;
      string pref = obj_to_pref[jets_label];
      float pt = vfloats_values[makeName(jets_label,pref,"Pt")][j];
      float ptzero = vfloats_values[makeName(jets_label,pref,"Pt")][j];
      float genpt = vfloats_values[makeName(jets_label,pref,"GenJetPt")][j];
      float eta = vfloats_values[makeName(jets_label,pref,"Eta")][j];
      float phi = vfloats_values[makeName(jets_label,pref,"Phi")][j];
      float energy = vfloats_values[makeName(jets_label,pref,"E")][j];
      float ptnomu = pt;
      
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

      //cout << " jet "<<j << " pt "<< pt << " jecfactor0 "<< jecscale <<endl;
      float chEmEnFrac = vfloats_values[makeName(jets_label,pref,"chargedEmEnergyFrac")][j];
      float neuEmEnFrac = vfloats_values[makeName(jets_label,pref,"neutralEmEnergyFrac")][j];


      if(pt>0){
	
	TLorentzVector jetUncorr_, jetCorr, jetUncorrNoMu_,jetCorrNoMu, jetL1Corr, jetWithUncCorr_;
	TLorentzVector T1Corr;
	T1Corr.SetPtEtaPhiE(0,0,0,0);	
	
	jetUncorr_.SetPtEtaPhiE(pt,eta,phi,energy);
	jetUncorrNoMu_.SetPtEtaPhiE(pt,eta,phi,energy);

	jetUncorr_= jetUncorr_*jecscale;
	jetUncorrNoMu_= jetUncorrNoMu_*jecscale;
	
	//Insert no mu
	
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
	}  
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
	  if(doT1MET){
	      
	    jecCorr->setJetPhi(jetUncorrNoMu_.Phi());
	    jecCorr->setJetEta(jetUncorrNoMu_.Eta());
	    jecCorr->setJetE(jetUncorrNoMu_.E());
	    jecCorr->setJetPt(jetUncorrNoMu_.Perp());
	    jecCorr->setJetA(area);
	    jecCorr->setRho(Rho);
	    jecCorr->setNPV(nPV);
	    
	    double recorrMu =  jecCorr->getCorrection();
	    jetCorrNoMu = jetUncorrNoMu_ *recorrMu;
	    
	    //// Jet corrections for level 1
	    jecCorr_L1->setJetPhi(jetUncorrNoMu_.Phi()); /// deve essere raw
	    jecCorr_L1->setJetEta(jetUncorrNoMu_.Eta());
	    jecCorr_L1->setJetE(jetUncorrNoMu_.E());
	    jecCorr_L1->setJetPt(jetUncorrNoMu_.Perp());
	    jecCorr_L1->setJetA(area);
	    jecCorr_L1->setRho(Rho);
	    jecCorr_L1->setNPV(nPV);
	    
	    double recorr_L1 =  jecCorr_L1->getCorrection();
	    jetL1Corr = jetUncorrNoMu_ *recorr_L1;
	      
	    ptnomu = jetCorrNoMu.Pt();
	    if(ptnomu>15.0 && ( chEmEnFrac+  neuEmEnFrac <0.9)){ T1Corr += jetCorrNoMu - jetL1Corr;
	      ptCorr_mL1 = T1Corr.Pt();
	      
	    }
	    for (size_t ju = 0; ju < obj_systCats[jets_label].size();++ju){
	      string systjet=obj_systCats[jets_label].at(ju);
	      if (systjet.find("JES")!=std::string::npos || 
		  systjet.find("JER")!=std::string::npos ){

		float pxCorrJet=jetCorr.Px()*jetUncertainty(jetCorr.Pt(),jetCorr.Eta(),systjet);
		float pyCorrJet=jetCorr.Py()*jetUncertainty(jetCorr.Pt(),jetCorr.Eta(),systjet);
		
		vfloats_values[makeName(met_label,prefm,"CorrT1Px"+systjet)][0]-=pxCorrJet;
		vfloats_values[makeName(met_label,prefm,"CorrT1Py"+systjet)][0]-=pyCorrJet;


	      }
	    }
	    for (size_t ju = 0; ju < obj_systCats[jets_label].size();++ju){
	      continue;
	      string systjet=obj_systCats[jets_label].at(ju);
	      if (systjet.find("JES")!=std::string::npos || 
		  systjet.find("JER")!=std::string::npos ){
		float smearfactJet= smear(ptnomu, genpt, eta,systjet);
		
		float uncJet=jetUncertainty(ptnomu,jetCorrNoMu.Eta(),systjet);
		//		uncJet=0;
		float ptCorrJet = ptnomu * (1+uncJet);// * smearfactJet;
		float energyCorrJet = jetCorrNoMu.E() * (1+uncJet);// * smearfactJet;
		float phinomu = jetCorrNoMu.Phi();
		jetWithUncCorr_.SetPtEtaPhiE(ptCorrJet,jetCorrNoMu.Eta(),phinomu,energyCorrJet);

		float pxCorrJet=(jetWithUncCorr_-jetL1Corr).Px();
		float pyCorrJet=(jetWithUncCorr_-jetL1Corr).Py();

		//		float ptCorrJet=(jetWithUncCorr_-jetL1Corr).Perp();
		//float phiCorrJet=(jetWithUncCorr_-jetL1Corr).Phi();
		//		cout << " jet unc absolute value " << ptnomu * (uncJet)*cos(phinomu)  << " contribution to la "<< pxCorrJet-(jetCorrNoMu-jetL1Corr).Px() <<endl;

		if(ptnomu*(1+uncJet)>15.0 && ( chEmEnFrac+  neuEmEnFrac <0.9)){ 
		  //cout << ptnomu * (1+uncJet)<< endl;
		  vfloats_values[makeName(met_label,prefm,"CorrT1Px"+systjet)][0]-=pxCorrJet;
		  vfloats_values[makeName(met_label,prefm,"CorrT1Py"+systjet)][0]-=pyCorrJet;
		  float pxt=vfloats_values[makeName(met_label,prefm,"CorrT1Px"+systjet)][0];
		  float pyt=vfloats_values[makeName(met_label,prefm,"CorrT1Py"+systjet)][0];
		  float ptt=sqrt(pxt*pxt +pyt*pyt);
		  vfloats_values[makeName(met_label,prefm,"CorrT1Pt"+systjet)][0]=ptt;
		}
	      }
	    }
	  }
	  else ptnomu=pt;
	}
	
	smearfact = smear(pt, genpt, eta, syst);
	if(!doResol) smearfact =1.0;
	ptCorr = pt * smearfact;
	energyCorr = energy * smearfact;

	float unc = jetUncertainty(ptCorr,eta,syst);
	float unc_nosmear = jetUncertainty(pt,eta,syst);
	
	for (size_t ju = 0; ju < obj_systCats[jets_label].size();++ju){
	  string systjet=obj_systCats[jets_label].at(ju);
	  if (systjet.find("JES")!=std::string::npos || 
	      systjet.find("JER")!=std::string::npos ){
	    //	    cout << " ju "<< ju <<endl;
	    float smearfactJet= smear(pt, genpt, eta,systjet);
	    float uncJet=jetUncertainty(ptCorr,eta,systjet);
	    float ptCorrJet = pt * (1+uncJet) * smearfactJet;
	    float energyCorrJet = energy * (1+uncJet) * smearfactJet;
	    //	    cout << " cat "<< systjet<< " pt "<< ptCorr;
	    //	    cout << " catname " << makeName(jets_label,pref,"CorrPt"+systjet)<< " value " << ptCorrJet<<endl;
	      //makeBranchNameCat(jets_label,systjet,jets_label+"_","CorrPt")<<endl;
	    //	    fillCategory()
	    //	    vfloats_values[makeName(jets_label+systjet,pref,"CorrPt")][j]=ptCorrJet;
	    vfloats_values[makeName(jets_label,pref,"CorrPt"+systjet)][j]=ptCorrJet;
	    vfloats_values[makeName(jets_label,pref,"CorrE"+systjet)][j]=energyCorrJet;
	    //	    cout <<" value after "<<vfloats_values[makeName(jets_label,pref,"CorrPt"+systjet)][j]<<endl;
	    //setCatCategoryValue(jets_label,systjet,sizes[jets_label+systjet],"CorrPt",ptCorrJet);
	
	  }
	}
	
	ptCorr = ptCorr * (1 + unc);
	energyCorr = energyCorr * (1 + unc);

	ptCorrSmearZero = pt * (1 + unc);//For full correction, including new JECs and MET
	ptCorrSmearZeroNoMu = ptnomu * (1 + unc);//For full correction, including new JECs and MET
	ptCorr_mL1 = ptCorr_mL1 * (1 + unc);//For T1 correction
	
	corrMetPx -=(cos(phi)*(ptCorr-ptzero));
        corrMetPy -=(sin(phi)*(ptCorr-ptzero));

	if(ptCorrSmearZeroNoMu>15.0 && ( chEmEnFrac+  neuEmEnFrac <0.9)){ // should be == T1 corrections minus the L1 difference.
	  corrBaseMetPx -=(cos(phi)*(ptCorrSmearZero-ptzero));
	  corrBaseMetPy -=(sin(phi)*(ptCorrSmearZero-ptzero));
	}
	
	

      	if( ptCorrSmearZeroNoMu>15.0 && jetCorrNoMu.Pt()>0.0 &&doT1MET){ 
	  //Adding in T1 corrections
	  corrMetT1Px -=(T1Corr.Px());
	  corrMetT1Py -=(T1Corr.Py());
	
	  
  
	  //Correcting by jes uncertainty if available
	  //	  corrMetT1Px -=(cos(phi)*(pt*unc_nosmear));
	  //rrMetT1Py -=(sin(phi)*(pt*unc_nosmear));
	}

      }
      
      float csv = vfloats_values[makeName(jets_label,pref,"CSVv2")][j];
      float cmva = vfloats_values[makeName(jets_label,pref,"CMVAv2")][j];
      
      float partonFlavour = vfloats_values[makeName(jets_label,pref,"PartonFlavour")][j];
      float hadronFlavour = vfloats_values[makeName(jets_label,pref,"HadronFlavour")][j];
      int flavor = int(hadronFlavour);
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

      //WP taken from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco

      bool isCSVT = csv  > getWPAlgo("CSV","T");
//0.9535;
      bool isCSVM = csv  > getWPAlgo("CSV","M");//0.8484;
      bool isCSVL = csv  > getWPAlgo("CSV","L");//0.5426;
      
      bool isCMVAT = cmva  > getWPAlgo("CMVA","T");//0.9432;
      bool isCMVAM = cmva  > getWPAlgo("CMVA","M");//0.4432;
      bool isCMVAL = cmva  > getWPAlgo("CMVA","L");//-0.5884;



      vfloats_values[jets_label+"_IsCSVT"][j]=isCSVT;
      vfloats_values[jets_label+"_IsCSVM"][j]=isCSVM;
      vfloats_values[jets_label+"_IsCSVL"][j]=isCSVL;
      
      vfloats_values[jets_label+"_IsCMVAT"][j]=isCMVAT;
      vfloats_values[jets_label+"_IsCMVAM"][j]=isCMVAM;
      vfloats_values[jets_label+"_IsCMVAL"][j]=isCMVAL;

      
      //      float bsf = getScaleFactor(ptCorr,eta,hadronFlavour,"noSyst");
      //      float bsfup = getScaleFactor(ptCorr,eta,hadronFlavour,"up");
      //      float bsfdown = getScaleFactor(ptCorr,eta,hadronFlavour,"down");
      
      //      vfloats_values[jets_label+"_BSF"][j]=bsf;
      //      vfloats_values[jets_label+"_BSFUp"][j]=bsfup;
      //      vfloats_values[jets_label+"_BSFDown"][j]=bsfdown;

      float flavForShaping = hadronFlavour;
      if(doTopDecayReshaping ){
	//	cout <<  " top charge before b"<<topCharge<<endl;

	//double leadingLeptonCharge=+1;//leadingleptonCharge;
	double product = partonFlavour*topCharge;
	float reshapeF_SD_CSV=1.0;
	float reshapeF_SD_CMVA=1.0;
	if(abs(partonFlavour)==5 && product >0 ){
	  reshapeF_SD_CSV=getReshapedBTagValue(partonFlavour, csv, pt, eta,"CSV_sd",syst);
	  reshapeF_SD_CMVA=getReshapedBTagValue(partonFlavour, cmva, pt, eta,"CMVA_sd",syst);
	  //reshapeF_SD = getReshapedBTagValue(1, csv, pt, eta,"csv_sd",syst)/reshapeF_SD;
	  //	  reshapeF_SD = getReshapedBTagValue(1, csv, pt, eta,"csv_sd",syst)/reshapeF_SD;

	  flavForShaping=1;//CSV needs to be reshaped to lihgt-quark data scale factors rather than b-tag scale factors
	}
	vfloats_values[jets_label+"_reshapeFactorCSV_SD"][j]=reshapeF_SD_CSV;
	vfloats_values[jets_label+"_reshapeFactorCMVA_SD"][j]=reshapeF_SD_CMVA;
	
      }
      float reshapeF=getReshapedBTagValue(flavForShaping, csv, pt, eta,"CSV",syst);
      float reshapeF_CMVA=getReshapedBTagValue(flavForShaping, cmva, pt, eta,"CMVA",syst);
      //      cout << " jet csv"<< csv<<endl;
      //      cout << "reshape factor "<< reshapeF << " value "<< reshapeF*csv<<endl;
      vfloats_values[jets_label+"_reshapeFactorCSV"][j]=reshapeF;
      vfloats_values[jets_label+"_reshapedCSV"][j]=reshapeF*csv;
      vfloats_values[jets_label+"_reshapeFactorCMVA"][j]=reshapeF_CMVA;
      vfloats_values[jets_label+"_reshapedCMVA"][j]=reshapeF_CMVA*cmva;
      
      bool passesID = true;
  
      if(!(jecscale*energy > 0))passesID = false;
      else{
        float neuMulti = vfloats_values[makeName(jets_label,pref,"neutralMultiplicity")][j];
        float chMulti = vfloats_values[makeName(jets_label,pref,"chargedMultiplicity")][j];
        float chHadEnFrac = vfloats_values[makeName(jets_label,pref,"chargedHadronEnergyFrac")][j];
        float neuHadEnFrac = vfloats_values[makeName(jets_label,pref,"neutralHadronEnergyFrac")][j];
        float numConst = chMulti+neuMulti;//vfloats_values[makeName(jets_label,pref,"NumConstituents")][j];
        //passesID =  (nDau >1.0 && fabs(eta) < 4.7 && (fabs(eta)>=2.4 ||(chHadEnFrac>0.0 && chMulti>0 && chEmEnFrac<0.99)) && neuEmEnFrac<0.99 && neuHadEnFrac <0.99 && muEnFrac<0.8) ;

        if(fabs(eta)<=2.7){
          passesID =  (neuHadEnFrac<0.99 && neuEmEnFrac<0.99 && numConst>1) && 
	    ( (fabs(eta)<=2.4 && chHadEnFrac>0 && chMulti>0 && chEmEnFrac<0.99) || fabs(eta)>2.4);
        }
	else if( (fabs(eta) >2.7) && (fabs(eta)<=3.0)) {
          //passesID = neuEmEnFrac<0.90 && neuMulti>2. ;
	  // spotted during the sync exe. https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016: For 2.7<|eta|<= 3.0 Apply
	  passesID = neuEmEnFrac > 0.01 && neuHadEnFrac < 0.98 && neuMulti > 2.;
	}
        else if(fabs(eta)>3.0){
          passesID = neuEmEnFrac<0.90 && neuMulti>10. ;
        }
      }
      
      vfloats_values[jets_label+"_PassesID"][j]=(float)passesID;
      //Remove overlap with tight electrons/muons
      double minDR=9999;
      double minDRThrEl=0.3;
      //double minDRThrMu=0.3;
      double minDRThrMu=0.4;
      bool passesDR=true;
 
      
      for(size_t nmuc=0;nmuc<obj_cats[mu_label].size();nmuc++){
	string catmu = obj_cats[mu_label].at(nmuc);
	if (catmu.find("Tight")==std::string::npos)continue;
	double minDRCat = 9999;
	bool passesDRCat=true;
	for (int m =0;m<sizes[mu_label+catmu];++m){
	  minDRCat = min(minDRCat,deltaR(math::PtEtaPhiELorentzVector(vfloats_values[makeName(mu_label+catmu,"","Pt")][m],vfloats_values[makeName(mu_label+catmu,"","Eta")][m],vfloats_values[makeName(mu_label+catmu,"","Phi")][m],vfloats_values[makeName(mu_label+catmu,"","E")][m]),math::PtEtaPhiELorentzVector(ptCorr, eta, phi, energyCorr)));
	  //	  cout << " ptetaphie "<< ptCorr<<" "<< eta <<" " <<phi <<" " <<energyCorr<<endl;
	  //	  cout << "catmu is "<<catmu<<" min DR Cat "<<minDRCat<<endl;
	  if(minDRCat <minDRThrMu){
	    passesDRCat=false;break;
	  }
	}
	vfloats_values[jets_label+"_PassesDRmu"+catmu][j]=(float)passesDRCat;
      }

      for(size_t nelc=0;nelc<obj_cats[ele_label].size();nelc++){
	string catel = obj_cats[ele_label].at(nelc);
	if (catel.find("Tight")==std::string::npos)continue;
	double minDRCat = 9999;
	bool passesDRCat=true;
	for (int m =0;m<sizes[ele_label+catel];++m){
	  minDRCat = min(minDRCat,deltaR(math::PtEtaPhiELorentzVector(vfloats_values[makeName(ele_label+catel,"","Pt")][m],vfloats_values[makeName(ele_label+catel,"","Eta")][m],vfloats_values[makeName(ele_label+catel,"","Phi")][m],vfloats_values[makeName(ele_label+catel,"","E")][m]),math::PtEtaPhiELorentzVector(ptCorr, eta, phi, energyCorr)));
	  //	  cout << " ptetaphie "<< ptCorr<<" "<< eta <<" " <<phi <<" " <<energyCorr<<endl;
	  //      cout << "catel is "<<catel<<" min DR Cat "<<minDRCat<<endl;
	  if(minDRCat <minDRThrEl){
	    passesDRCat=false;break;
	  }
	}
	vfloats_values[jets_label+"_PassesDRel"+catel][j]=(float)passesDRCat;
      }


      for (size_t e = 0; e < (size_t)electrons.size(); ++e){
	minDR = min(minDR,deltaR(math::PtEtaPhiELorentzVector(electrons.at(e).Pt(),electrons.at(e).Eta(),electrons.at(e).Phi(),electrons.at(e).Energy() ) ,math::PtEtaPhiELorentzVector(ptCorr, eta, phi, energyCorr)));
	if(minDR<minDRThrEl)passesDR = false;
      }
      for (size_t m = 0; m < (size_t)muons.size(); ++m){
	minDR = min(minDR,deltaR(math::PtEtaPhiELorentzVector(muons.at(m).Pt(),muons.at(m).Eta(),muons.at(m).Phi(),muons.at(m).Energy() ) ,math::PtEtaPhiELorentzVector(ptCorr, eta, phi, energyCorr)));
	if(minDR<minDRThrMu)passesDR = false;
      }
      //      cout<<" jet "<<j<< " minDR is "<< minDR<<endl;
      
      vfloats_values[jets_label+"_MinDR"][j]=minDR;
      vfloats_values[jets_label+"_PassesDR"][j]=(float)passesDR;
      
      vfloats_values[jets_label+"_IsTight"][j]=0.0;
      vfloats_values[jets_label+"_IsLoose"][j]=0.0;
	
      if( passesID && passesDR) vfloats_values[jets_label+"_IsLoose"][j]=1.0;

      if(passesID && passesDR && pt>50 && abs(eta)<2.4){
	Ht+=pt;
      }
      
      float_values["Event_Ht"] = (float)Ht;
      

      double jetval = 40.0;
      double jetbaseval = 20.0;
      
      bool passesCut = ( ptCorr > jetval && fabs(eta) < 4.7);
      
      if(passesID && passesDR){
	//	fillSysts(jets_label,"Tight",)
	
	if(passesCut){
	  vfloats_values[jets_label+"_IsTight"][j]=1.0;
	  TLorentzVector jet;
	  jet.SetPtEtaPhiE(ptCorr, eta, phi, energyCorr);
	  jets.push_back(jet);
	  
	  if(!isCSVM) jetsnob.push_back(jet);
	  
	  nTightJets+=1;
	  if(isInVector(obj_cats[jets_label],"Tight")){
	    //	    cout<< " tight jet precat "<<j << " sizes tight cat "<<sizes[jets_label+"Tight"]<<endl;

	    fillCategory(jets_label,"Tight",j,sizes[jets_label+"Tight"]);
	    //	    cout << " sizes post cat "<< sizes[jets_label+"Tight"]<<endl;
	    //	  cout<< " tight jet poscat "<<j <<endl;
	    
	  }
	}
	fillSysts(jets_label,"Tight",j,"CorrPt","40_GE");
      }
      if(isCSVM && passesCut &&  passesID && passesDR && fabs(eta) < 2.4) { 
	float_values["Event_nCSVMJets"]+=1.0;
	ncsvm_tags +=1;
	TLorentzVector bjet;
	bjet.SetPtEtaPhiE(ptCorr, eta, phi, energyCorr);
	bjets.push_back(bjet);
	++bjetidx;
      
      }
      
      if(isCSVL && passesCut &&  passesID && passesDR && abs(eta) < 2.4) float_values["Event_nCSVLJets"]+=1;
      
      bool passesBaseCut = ( ptCorr > jetbaseval && fabs(eta) < 4.7);
      //      cout <<" in cat loop "<<endl;
      if(passesID && passesDR) if(obj_scanCuts[jets_label].size()>=1) {
	  if(passesBaseCut) fillScanCuts(jets_label,"Tight",j);
	  if(fabs(eta) < 4.7) fillScanSystsCuts(jets_label,"Tight",j);
      
	  /*if(passesBaseCut) fillScanCuts(jets_label,"Tight",j);
	  	    if(isInVector(obj_cats[jets_label],"Tight")){
	      string catjet=obj_cats[jets_label].at(ju);
	      if ( ( catjet.find("JES")!=std::string::npos || 
		     catjet.find("JER")!=std::string::npos ) && !isScanCut(jets_label,catjet) )  {
		cout << " cat "<< catjet << " jet "<< j << " pt "<< makeName(jets_label+catjet,pref,"Pt")<< " value "<< vfloats_values[makeName(jets_label+catjet,pref,"Pt")][j] <<endl;
		passesJecCuts[catjet]= fabs(eta) < 4.7 && vfloats_values[makeName(jets_label+catjet,pref,"Pt")][j] >jetbaseval;
		if (passesJecCuts[catjet]){
		  cout << " cat "<< catjet<< " size is? ";
		  cout << sizes[jets_label+catjet];
		  fillScanCuts(jets_label,catjet,j);
		  cout << " after "<<sizes[jets_label+catjet]<< endl;
		}
		//	cout << "now testing scan jets "<<endl;
	      }
	    }
	  }
	  */
      
    }
    
      //      }//ENDS HERE
      
    }
    fillSysts(met_label,"CorrT1",0,"","");

    if(isInVector(obj_cats[jets_label],"Tight")){
      sizes[jets_label+"Tight"]=(int)nTightJets;
      setEventBTagSF(jets_label,"Tight","CSV");
      setEventBTagSF(jets_label,"Tight","CMVA");
      
      for(size_t jc = 0; jc < obj_scanCuts[jets_label].size();++jc){
	string scanCut=obj_scanCuts[jets_label].at(jc);
	
	//	if(false)
	setEventBTagSF(jets_label,"Tight_"+scanCut,"CSV");
	setEventBTagSF(jets_label,"Tight_"+scanCut,"CMVA");
      }
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
    
    //Met and mt

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

    
    if( doT1MET){
      //Correcting the pt
      metT1Px+=corrMetT1Px; metT1Py+=corrMetT1Py; // add JEC/JER contribution
    }
    float metptT1Corr = sqrt(metT1Px*metT1Px + metT1Py*metT1Py);
    vfloats_values[met_label+"_CorrT1Pt"][0]=metptT1Corr;
    
    for (size_t ju = 0; ju < obj_systCats[jets_label].size();++ju){
      string systjet=obj_systCats[jets_label].at(ju);
      if (systjet.find("JES")!=std::string::npos || 
	  systjet.find("JER")!=std::string::npos ){
	vfloats_values[makeName(met_label,prefm,"CorrT1Px"+systjet)][0]+=metT1Px;
	vfloats_values[makeName(met_label,prefm,"CorrT1Py"+systjet)][0]+=metT1Py;
	float pxt=vfloats_values[makeName(met_label,prefm,"CorrT1Px"+systjet)][0];
	float pyt=vfloats_values[makeName(met_label,prefm,"CorrT1Py"+systjet)][0];
	float ptt=sqrt(pxt*pxt +pyt*pyt);
	vfloats_values[makeName(met_label,prefm,"CorrT1Pt"+systjet)][0]=ptt;
	vfloats_values[makeName(met_label,prefm,"CorrT1Phi"+systjet)][0] = getPhi(pxt,pyt);

      }
    }


    float metptCorr = sqrt(metPxCorr*metPxCorr + metPyCorr*metPyCorr);
    vfloats_values[met_label+"_CorrPt"][0]=metptCorr;

    float metptCorrBase = sqrt(metPxCorrBase*metPxCorrBase + metPyCorrBase*metPyCorrBase);
    vfloats_values[met_label+"_CorrBasePt"][0]=metptCorrBase;
    
    //Correcting the phi
    float metphiCorr = metphi;
    metphiCorr=getPhi(metPx,metPy);
    //    if(metPx<0){
    //      if(metPy>0)metphiCorr = atan(metPy/metPx)+3.141592;
    //      if(metPy<0)metphiCorr = atan(metPy/metPx)-3.141592;
    //    }
    //    else  metphiCorr = (atan(metPy/metPx));

    float metphiCorrBase = metphi;
    metphiCorrBase=getPhi(metPxCorrBase,metPyCorrBase);
    /*    if(metPxCorrBase<0){
      if(metPyCorrBase>0)metphiCorrBase = atan(metPyCorrBase/metPxCorrBase)+3.141592;
      if(metPyCorrBase<0)metphiCorrBase = atan(metPyCorrBase/metPxCorrBase)-3.141592;
    }
    else  metphiCorrBase = (atan(metPyCorrBase/metPxCorrBase));
    */

    float metphiCorrT1 = metphi;
    metphiCorrT1 = getPhi(metT1Px,metT1Py);

    /*    if(metT1Px<0){
      if(metT1Py>0)metphiCorrT1 = atan(metT1Py/metT1Px)+3.141592;
      if(metT1Py<0)metphiCorrT1 = atan(metT1Py/metT1Px)-3.141592;
    }
    else  metphiCorr = (atan(metPyCorr/metPxCorr));
    */
    vfloats_values[met_label+"_CorrPhi"][0]=metphiCorr;
    vfloats_values[met_label+"_CorrBasePhi"][0]=metphiCorrBase;
    vfloats_values[met_label+"_CorrT1Phi"][0]=metphiCorrT1;

    //Preselection part

    //    cout << " before preselection "<<endl;
		       
    if(doPreselection){
      bool passes = true;
      //bool metCondition = (metptCorr >100.0);
      passes = passes && nTightJets>=1.0;
      //      passes = passes && nTightLeptons>=1.0;
      passes = passes && (nTightLeptons>=1.0 || nTightAntiIsoLeptons>=1.0);
     
      if (!passes ) {
	//Reset eventmax weights/#objects
	string nameshortv= "Event";
	vector<string> extravars = additionalVariables(nameshortv);
	for(size_t addv = 0; addv < extravars.size();++addv){
	  string name = nameshortv+"_"+extravars.at(addv);
	  if(!isMCWeightName(extravars.at(addv))) float_values[name]=0.0;
	}
	resetVariables("genMuon",max_instances[mu_label]);
	resetVariables("genElectron",max_instances[mu_label]);
	continue;
      }
    }
    

    //======================================================
    //    cout << " before lepton "<<endl;

    TLorentzVector lepton;
    //    cout << " test 1 "<<endl;
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
      
      TLorentzVector subjet;
      subjet.SetPtEtaPhiE(pt, eta, phi, e);       
      double minDR=999;
      float subjcsv = vfloats_values[makeName(boosted_tops_subjets_label,pref,"CSVv2")][s];
     
      bool isCSVM = (subjcsv>0.800);
      
      
      for(int t = 0;t < min(max_instances[boosted_tops_label],sizes[boosted_tops_label]) ;++t){
	
	float ptj = vfloats_values[makeName(boosted_tops_label,pref,"Pt")][t];
	if (ptj<0.0)continue;
	float etaj = vfloats_values[makeName(boosted_tops_label,pref,"Eta")][t];
	float phij = vfloats_values[makeName(boosted_tops_label,pref,"Phi")][t];
	float ej = vfloats_values[makeName(boosted_tops_label,pref,"E")][t];
	TLorentzVector topjet;
	topjet.SetPtEtaPhiE(ptj, etaj, phij, ej);       
	
	float DR = subjet.DeltaR(topjet); 
	//	  cout <<"subjet# "<<s<< " jet "<< t << " DR is "<< DR << " minDR is "<<minDR<<endl;
	if(DR < minDR){
	  minDR = DR;
	  subj_jet_map[s]=t;
	}
	//	    cout << " sjmap is "<< subj_jet_map[s]<<endl;
      }
      size_t tm = subj_jet_map[s];
      if(isCSVM)vfloats_values[boosted_tops_label+"_nCSVM"][tm]+=1;
      vfloats_values[boosted_tops_label+"_nJ"][tm]+=1;
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
    //Reset event weights/#objects
    if(useLHEWeights){
      getEventLHEWeights();
    }
    if(addLHAPDFWeights){
      getEventPdf();
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
    
    //technical event information
    double_values["Event_EventNumber"]=*eventNumber;
    float_values["Event_LumiBlock"]=*lumiBlock;
    float_values["Event_RunNumber"]=*runNumber;
    
    trees[syst]->Fill();

    string nameshortv= "Event";
    vector<string> extravars = additionalVariables(nameshortv);
    for(size_t addv = 0; addv < extravars.size();++addv){
      string name = nameshortv+"_"+extravars.at(addv);
      if(!isMCWeightName(extravars.at(addv))) float_values[name]=0.0;
    }
    //Reset gen variables
    resetVariables("genMuon",max_instances[mu_label]);
    resetVariables("genElectron",max_instances[mu_label]);
  }
  for(int t = 0;t < max_instances[boosted_tops_label] ;++t){
    vfloats_values[boosted_tops_label+"_nCSVM"][t]=0;
    vfloats_values[boosted_tops_label+"_nJ"][t]=0;
  }
  
  //treesBase->Fill(); 
}

void DMAnalysisTreeMaker::resetVariables(string cat, int maxinstances){
  string nameshortv= cat;
  vector<string> extravars = additionalVariables(cat);
  size_t smem=(size_t)maxinstances;
  for(size_t addv = 0; addv < extravars.size();++addv){
    string name = nameshortv+"_"+extravars.at(addv);
    for (size_t nmem=0; nmem < smem; ++nmem){
      vfloats_values[name][nmem]=0.0;
    }
  }
  sizes[cat]=0;
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
  //  string outVar = var;
  //  size_t prefPos=outVar.find(pref);
  //  size_t prefLength = pref.length();
  //  std::cout << " outvar is "<< outVar<< " prefpos "<< prefPos<< " length "<< prefLength<< endl;
  //  std::cout << "it is " << label+"_"+var<<endl;
  //  outVar.replace(prefPos,prefLength,label+"_");

  //  outVar.replace(prefPos,prefLength,label+"_");
  //  return outVar;
  //  return makeBranchName(label,pref,var);
  //  return pref+var+"_"+label;
}

void DMAnalysisTreeMaker::initCategoriesSize(string label){
  for(size_t sc = 0; sc< obj_cats[label].size() ;++sc){
    string category = obj_cats[label].at(sc);
    //cout << " label "<<label<< " cat "<<category << " size bef "<<sizes[label+category]<<endl;
    sizes[label+category]=0;
    //    cout << " size after "<<sizes[label+category]<<endl;
  }
}

/*void DMAnalysisTreeMaker::initSystSize(string label){
  for(size_t ss = 0; sc< obj_systs[label].size() ;++ss){
    for(size_t sc = 0; sc< obj_cats[label].size() ;++sc){
      string category = obj_cats[label].at(sc);
      string category = obj_cats[label].at(sc);
      //cout << " label "<<label<< " cat "<<category << " size bef "<<sizes[label+category]<<endl;
      sizes[label+category]=0;
      //    cout << " size after "<<sizes[label+category]<<endl;
    }
  }
  }*/

void DMAnalysisTreeMaker::setCategorySize(string label, string category, size_t size){
    sizes[label+category]=size;
}

void  DMAnalysisTreeMaker::setCatCategoryValue(string label,string category, int pos_cat,string var, float value){
  string varCat = makeName(label+category,"",var);
  vfloats_values[varCat][pos_cat]=value;
}

void DMAnalysisTreeMaker::fillCategory(string label, string category, int pos_nocat, int pos_cat){
  //  cout << " filling cat  " << category << " for label "<< label<<endl;
  sizes[label+category]=pos_cat+1;//Update the size with the latest position filled
 
  for (size_t obj =0; obj< obj_to_floats[label].size(); ++obj){
    
    string var = obj_to_floats[label].at(obj);
    string varCat = makeBranchNameCat(label,category,label+"_",var);
    //    cout << " var "<< var << " varcat "<< varCat<<endl;
    //    cout << " poscat "<< pos_cat << " posnocat "<< pos_nocat<<endl;
    vfloats_values[varCat][pos_cat]= vfloats_values[var][pos_nocat];
  }
  //  cout << "sizes is there "<< sizes[label+category]<<endl;
}


void DMAnalysisTreeMaker::fillSystCategory(string label, string category, string sys, int pos_nocat, int pos_cat, string scancut){
  //  cout << " filling cat  " << category << " for label "<< label<<endl;
  string catsys = category+sys+scancut;
  sizes[label+catsys]=pos_cat+1;//Update the size with the latest position filled

  
  for (size_t obj =0; obj< obj_to_floats[label].size(); ++obj){
    
    string var = obj_to_floats[label].at(obj);
    
    string varCat = makeBranchNameCat(label,catsys,label+"_",var);
    //    if(label == met_label){  cout << " var "<< var << " varcat "<< varCat<<endl;
    //      cout << " poscat "<< pos_cat << " posnocat "<< pos_nocat<<endl;}
    vfloats_values[varCat][pos_cat]= vfloats_values[var][pos_nocat];
  }
  for (size_t obj =0; obj< obj_to_floats[label].size(); ++obj){
    string var = obj_to_floats[label].at(obj);    
    //    if(label == met_label)cout << " var is "<< var <<" sys is "<< sys<<endl;
    if(var.find(sys)!=std::string::npos){
      string varnew = var;
      varnew.replace(varnew.find(sys),sys.length(),"");
      string varCat = makeBranchNameCat(label,catsys,label+"_",varnew);
      vfloats_values[varCat][pos_cat]= vfloats_values[var][pos_nocat];
	//      piecetmp.replace(piecetmp.find("p"),1,".");
    }
  }
  //  cout << "sizes is there "<< sizes[label+category]<<endl;
}

bool DMAnalysisTreeMaker::isScanCut(string label, string category){
  bool isScCut=false;
  for(size_t sccut = 0; sccut< obj_scanCuts[label].size() ;++sccut){
    string scanCut = obj_scanCuts[label].at(sccut);
    if(category.find(scanCut)!=std::string::npos){
      isScCut=true;
      break;
    };
  }
  return isScCut;
}

bool DMAnalysisTreeMaker::isSysCat(string label, string category){
  bool isScCut=false;
  for(size_t sccut = 0; sccut< obj_systCats[label].size() ;++sccut){
    string scanCut = obj_systCats[label].at(sccut);
    if(category.find(scanCut)!=std::string::npos){
      isScCut=true;
      break;
    };
  }
  return isScCut;
}

void DMAnalysisTreeMaker::fillScanCuts(string label, string category, int pos_nocat){
  for(size_t sccut = 0; sccut< obj_scanCuts[label].size() ;++sccut){
    string scanCut = obj_scanCuts[label].at(sccut);
    //cout << " filling scan cut "<< scanCut <<" for label "<< label <<endl;
    if(passesScanCut(label,category,scanCut,pos_nocat)){
      int pos_cat=sizes[label+category+"_"+scanCut];// the size is the position of the last element
      fillCategory(label,category+"_"+scanCut,pos_nocat,pos_cat); // this will also increment the size of the category
    }
  }    
}

void DMAnalysisTreeMaker::fillScanSystsCuts(string label, string category, int pos_nocat){
  for(size_t sccut = 0; sccut< obj_scanCuts[label].size() ;++sccut){
    string scanCut = obj_scanCuts[label].at(sccut);
    string variable = parseScanCut(scanCut,0);
    string endcut= scanCut;
    endcut.replace(endcut.find(variable),variable.length()+1,"");
    fillSysts(label, category, pos_nocat, variable, endcut, "_"+scanCut);
  }
}

void DMAnalysisTreeMaker::fillSysts(string label,string category,int pos_nocat, string variable, string cut, string postfix){
  for(size_t syCat = 0; syCat< obj_systCats[label].size() ;++syCat){
    string sysCat = obj_systCats[label].at(syCat);
    string finalcut = (variable+sysCat+"_"+cut);
    //    cout << " sy is  "<< sysCat << " cut is "<< cut << " final cut is "<<finalcut<< " pf "<< postfix<<endl;
    //cout << " posnocat "<<pos_nocat<<" value "<< vfloats_values[label+category+"_"+variable+sysCat][pos_nocat]<<endl;
    bool passesCut = (cut== "" && variable == "");
    if(!passesCut){passesCut= passesScanCut(label,category,finalcut,pos_nocat);}
    if(passesCut){
      //      cout<< " passes! filling cat "<< label+category+sysCat<< " pos "<< sizes[label+category+sysCat]<<endl; 
      int pos_cat=sizes[label+category+sysCat+postfix];// the size is the position of the last element
      //      cout << " posnocat "<<pos_nocat<<endl;
      fillSystCategory(label,category,sysCat,pos_nocat,pos_cat,postfix); // this will also increment the size of the category
    }
  }    
}

string DMAnalysisTreeMaker::parseScanCut(string category, int obj){
  std::stringstream cat;
  cat << category;
  std::string segment;
  std::vector<std::string> pieces;
  //  cout << " caat "<<category<<endl;
  while(std::getline(cat, segment, '_'))    {
    //    cout << " cat split "<<segment <<" "; 
      pieces.push_back(segment);
    }
  //  cout <<endl;
  if(pieces.size()<2)return "";
  if(obj==0){
    return (pieces.at(0));
  }
  if(obj==1){
    string piecetmp = pieces.at(1);

    if( piecetmp.find("p")!=std::string::npos){
      piecetmp.replace(piecetmp.find("p"),1,".");
    }
    return piecetmp;
  }
  if(obj==2){
    if (pieces.size()<3) return "GE";
    else{
      return (pieces.at(2));
    }
  }
  return "";
}


bool DMAnalysisTreeMaker::passesScanCut(string label, string category, string scanCut, int pos_nocat){
  string lookAtLabel = label+category;  
  string variable = parseScanCut(scanCut,0);
  stringstream cutvalue;
  cutvalue << parseScanCut(scanCut,1);
  float cutfloat=0.;
  cutvalue >> cutfloat;
  string sign = parseScanCut(scanCut,2);
  //  cout << " obs " << variable<< " sign "<< sign << " threshold " <<  cutfloat << endl;
  string pref = obj_to_pref[label];
  float var = vfloats_values[makeName(label,pref,variable)][pos_nocat];
  //  if( isSystCat(label, category))var=vfloats_values[makeName(label+category,pref,variable);
  //  cout << " variable full name is is "<< makeName(label,pref,variable)<<" value "<<var<< " passes cut? "<< cutfloat <<endl;
  
  if(sign=="GE") return (bool)(var>cutfloat);
  if(sign=="LE") return (bool)(var<cutfloat);

  return true;
  //  return false;
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
  
  bool isgenpart=object.find("genMuon")!=std::string::npos 
    || object.find("genEle")!=std::string::npos 
    || object.find("genTop")!=std::string::npos
    ||  object.find("genB")!=std::string::npos;
  bool isgenjet=object.find("genJet")!=std::string::npos;
  if(isgenpart || isgenjet){
    addvar.push_back("Pt");
    addvar.push_back("Eta");
    addvar.push_back("Phi");
    addvar.push_back("E");
    if(isgenpart){
      addvar.push_back("Mass");
      addvar.push_back("Id");
      addvar.push_back("MomId");
    }
    if(isgenjet){
      addvar.push_back("hadronFlavour");
      addvar.push_back("partonFlavour");
    }
  }

  if(ismuon || iselectron){
    addvar.push_back("SFTrigger");
  }
  /*
  if(ismuon || iselectron){
    addvar.push_back("SFTrigger");
    addvar.push_back("SFReco");
    addvar.push_back("isQCd-");
    //    Addvar.push_back("isTightOffline");
    //    addvar.push_back("isLooseOffline");
    }*/
  
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

    addvar.push_back("CorrT1Px");
    addvar.push_back("CorrT1Py");
    addvar.push_back("CorrPx");
    addvar.push_back("CorrPy");
    for(size_t sccut = 0; sccut< obj_systCats[jets_label].size() ;++sccut){
      string syCat=obj_systCats[jets_label].at(sccut);
      addvar.push_back("CorrT1Pt"+syCat);
      addvar.push_back("CorrT1Phi"+syCat);
      addvar.push_back("CorrT1Px"+syCat);
      addvar.push_back("CorrT1Py"+syCat);
    }

    //    addvar.push_back("CorrPtNoHF");
    // addvar.push_back("CorrPhiNoHF");
  }
  if(isjet){
    addvar.push_back("CorrPt");
    //    addvar.push_back("CorrEta");
    //    addvar.push_back("CorrPhi");
    addvar.push_back("CorrE");
    addvar.push_back("NoCorrPt");
    addvar.push_back("NoCorrE");
    addvar.push_back("MinDR");
    addvar.push_back("IsCSVT");
    addvar.push_back("IsCSVM");
    addvar.push_back("IsCSVL");

     addvar.push_back("IsCMVAT");
     addvar.push_back("IsCMVAM");
     addvar.push_back("IsCMVAL");

     addvar.push_back("BSF");
     addvar.push_back("BSFUp");
     addvar.push_back("BSFDown");
    addvar.push_back("PassesID");
    addvar.push_back("PassesDR");
    addvar.push_back("CorrMass");
    addvar.push_back("IsTight");
    addvar.push_back("IsLoose");

    for(size_t nmuc=0;nmuc<obj_cats[mu_label].size();nmuc++){
      string catmu = obj_cats[mu_label].at(nmuc);
      if (catmu.find("Tight")==std::string::npos)continue;
      addvar.push_back("PassesDRmu"+catmu);
    }
    for(size_t nelc=0;nelc<obj_cats[ele_label].size();nelc++){
      string catel = obj_cats[ele_label].at(nelc);
      if (catel.find("Tight")==std::string::npos)continue;
      addvar.push_back("PassesDRel"+catel);
    }

    //    addvar.push_back("CorrNJets");
    //    addvar.push_back("CorrPartonFlavour");
    
    addvar.push_back("reshapeFactorCSV");
    addvar.push_back("reshapeFactorCSV_SD");
    addvar.push_back("reshapeFactorCMVA");
    addvar.push_back("reshapeFactorCMVA_SD");
    addvar.push_back("reshapedCSV");
    addvar.push_back("reshapedCMVA");
    addvar.push_back("reshapedCSVJESUp");
    addvar.push_back("reshapedCSVJESDown");
    //    addvar.push_back("nReshapedCSV");
    //    addvar.push_back("nReshapedCSVMisTagDown");

    for(size_t sccut = 0; sccut< obj_systCats[jets_label].size() ;++sccut){
      string syCat=obj_systCats[jets_label].at(sccut);
      addvar.push_back("CorrPt"+syCat);
      addvar.push_back("CorrE"+syCat);
    }
  }
  if(isak8){
    addvar.push_back("CorrPt");
    addvar.push_back("CorrE");
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
    addvar.push_back("nCSVsubj");
    addvar.push_back("nCSVsubj_tm");
    addvar.push_back("tau3OVERtau2");
    addvar.push_back("tau2OVERtau1");

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
    addvar.push_back("nTightAntiIsoMuons");
    addvar.push_back("nSoftMuons");
    addvar.push_back("nLooseMuons");
    addvar.push_back("nTightElectrons");
    addvar.push_back("nTightAntiIsoElectrons");
    addvar.push_back("nMediumElectrons");
    addvar.push_back("nLooseElectrons");
    addvar.push_back("nVetoElectrons");
    addvar.push_back("nElectronsSF");
    addvar.push_back("nJets");
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
    
    std::vector<std::string>algos;
    algos.push_back("CSV");
    algos.push_back("CMVA");
    for (size_t jc=0;jc<obj_cats[jets_label].size();jc++){
      string category = obj_cats[jets_label].at(jc);
      
      addvar.push_back("nJets"+category);
      
      for(size_t alg =0; alg< algos.size();++alg){
	string algo= algos.at(alg);
	addvar.push_back("nJets"+algo+"T"+category);
	addvar.push_back("nJets"+algo+"M"+category);
	addvar.push_back("nJets"+algo+"L"+category);
	
	addvar.push_back("bWeight0"+algo+"T"+category);
	addvar.push_back("bWeight1"+algo+"T"+category);
	addvar.push_back("bWeight2"+algo+"T"+category);
	addvar.push_back("bWeight1_2"+algo+"T"+category);
	addvar.push_back("bWeightGE3"+algo+"T"+category);
	
	addvar.push_back("bWeight0"+algo+"M"+category);
	addvar.push_back("bWeight1"+algo+"M"+category);
	addvar.push_back("bWeight2"+algo+"M"+category);
	addvar.push_back("bWeight1_2"+algo+"M"+category);
	addvar.push_back("bWeightGE3"+algo+"M"+category);
	
	addvar.push_back("bWeight0"+algo+"L"+category);
	addvar.push_back("bWeight1"+algo+"L"+category);
	addvar.push_back("bWeight2"+algo+"L"+category);
	addvar.push_back("bWeight1_2"+algo+"L"+category);
	addvar.push_back("bWeightGE3"+algo+"L"+category);
	

	addvar.push_back("bWeightBTagUp0"+algo+"T"+category);
	addvar.push_back("bWeightBTagUp1"+algo+"T"+category);
	addvar.push_back("bWeightBTagUp2"+algo+"T"+category);
	addvar.push_back("bWeightBTagUp1_2"+algo+"T"+category);
	addvar.push_back("bWeightBTagUpGE3"+algo+"T"+category);
	
	addvar.push_back("bWeightBTagUp0"+algo+"M"+category);
	addvar.push_back("bWeightBTagUp1"+algo+"M"+category);
	addvar.push_back("bWeightBTagUp2"+algo+"M"+category);
	addvar.push_back("bWeightBTagUp1_2"+algo+"M"+category);
	addvar.push_back("bWeightBTagUpGE3"+algo+"M"+category);

	addvar.push_back("bWeightBTagUp0"+algo+"L"+category);
	addvar.push_back("bWeightBTagUp1"+algo+"L"+category);
	addvar.push_back("bWeightBTagUp2"+algo+"L"+category);
	addvar.push_back("bWeightBTagUp1_2"+algo+"L"+category);
	addvar.push_back("bWeightBTagUpGE3"+algo+"L"+category);
	
	addvar.push_back("bWeightBTagDown0"+algo+"T"+category);
	addvar.push_back("bWeightBTagDown1"+algo+"T"+category);
	addvar.push_back("bWeightBTagDown2"+algo+"T"+category);
	addvar.push_back("bWeightBTagDown1_2"+algo+"T"+category);
	addvar.push_back("bWeightBTagDownGE3"+algo+"T"+category);
	
	addvar.push_back("bWeightBTagDown0"+algo+"M"+category);
	addvar.push_back("bWeightBTagDown1"+algo+"M"+category);
	addvar.push_back("bWeightBTagDown2"+algo+"M"+category);
	addvar.push_back("bWeightBTagDown1_2"+algo+"M"+category);
	addvar.push_back("bWeightBTagDownGE3"+algo+"M"+category);
	
	addvar.push_back("bWeightBTagDown0"+algo+"L"+category);
	addvar.push_back("bWeightBTagDown1"+algo+"L"+category);
	addvar.push_back("bWeightBTagDown2"+algo+"L"+category);
	addvar.push_back("bWeightBTagDown1_2"+algo+"L"+category);
	addvar.push_back("bWeightBTagDownGE3"+algo+"L"+category);
	
	addvar.push_back("bWeightMisTagUp0"+algo+"T"+category);
	addvar.push_back("bWeightMisTagUp1"+algo+"T"+category);
	addvar.push_back("bWeightMisTagUp2"+algo+"T"+category);
	addvar.push_back("bWeightMisTagUp1_2"+algo+"T"+category);
	addvar.push_back("bWeightMisTagUpGE3"+algo+"T"+category);
	
	addvar.push_back("bWeightMisTagUp0"+algo+"M"+category);
	addvar.push_back("bWeightMisTagUp1"+algo+"M"+category);
	addvar.push_back("bWeightMisTagUp2"+algo+"M"+category);
	addvar.push_back("bWeightMisTagUp1_2"+algo+"M"+category);
	addvar.push_back("bWeightMisTagUpGE3"+algo+"M"+category);
	
	addvar.push_back("bWeightMisTagUp0"+algo+"L"+category);
	addvar.push_back("bWeightMisTagUp1"+algo+"L"+category);
	addvar.push_back("bWeightMisTagUp2"+algo+"L"+category);
	addvar.push_back("bWeightMisTagUp1_2"+algo+"L"+category);
	addvar.push_back("bWeightMisTagUpGE3"+algo+"L"+category);
	
	addvar.push_back("bWeightMisTagDown0"+algo+"T"+category);
	addvar.push_back("bWeightMisTagDown1"+algo+"T"+category);
	addvar.push_back("bWeightMisTagDown2"+algo+"T"+category);
	addvar.push_back("bWeightMisTagDown1_2"+algo+"T"+category);
	addvar.push_back("bWeightMisTagDownGE3"+algo+"T"+category);
	
	addvar.push_back("bWeightMisTagDown0"+algo+"M"+category);
	addvar.push_back("bWeightMisTagDown1"+algo+"M"+category);
	addvar.push_back("bWeightMisTagDown2"+algo+"M"+category);
	addvar.push_back("bWeightMisTagDown1_2"+algo+"M"+category);
	addvar.push_back("bWeightMisTagDownGE3"+algo+"M"+category);
	
	addvar.push_back("bWeightMisTagDown0"+algo+"L"+category);
	addvar.push_back("bWeightMisTagDown1"+algo+"L"+category);
	addvar.push_back("bWeightMisTagDown2"+algo+"L"+category);
	addvar.push_back("bWeightMisTagDown1_2"+algo+"L"+category);
	addvar.push_back("bWeightMisTagDownGE3"+algo+"L"+category);
	
	if(doTopBToLightQuarkReweight){
	  addvar.push_back("reweightBToQ"+category);
	  addvar.push_back("reweightBToQBTagUp"+category);
	  addvar.push_back("reweightBToQBTagDown"+category);
	  addvar.push_back("reweightBToQMisTagUp"+category);
	  addvar.push_back("reweightBToQMisTagDown"+category);
	}
      }
    }


    addvar.push_back("T_Pt");
    addvar.push_back("T_Eta");
    addvar.push_back("T_Phi");
    addvar.push_back("T_E");

    addvar.push_back("T_size");

    addvar.push_back("Tbar_Pt");
    addvar.push_back("Tbar_Eta");
    addvar.push_back("Tbar_Phi");
    addvar.push_back("Tbar_E");

    addvar.push_back("Tbar_size");

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
    //    addvar.push_back("LHEWeightAVG");
    addvar.push_back("LHEWeight");
    addvar.push_back("EventNumber");
    addvar.push_back("LumiBlock");
    addvar.push_back("RunNumber");
    
    if(useLHEWeights){
      for (size_t w = 0; w <= (size_t)maxWeights; ++w)  {
	//	cout << " weight # " << lhe_weights[w - 1] << " test "<< endl; 
	stringstream w_n;
	w_n << w;
	addvar.push_back("LHEWeight"+w_n.str());
	//addvar.push_back(("LHEWeight"+w_n.str())+"ID");
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
    }
    if(useTriggers){
      for (size_t lt = 0; lt < SingleElTriggers.size(); ++lt)  {
	string trig = SingleElTriggers.at(lt);
	addvar.push_back("passes"+trig);
	//addvar.push_back("prescale"+trig);
      }
      for (size_t lt = 0; lt < SingleMuTriggers.size(); ++lt)  {
	string trig = SingleMuTriggers.at(lt);
	addvar.push_back("passes"+trig);
	//	addvar.push_back("prescale"+trig);
      }
      for (size_t lt = 0; lt < PhotonTriggers.size(); ++lt)  {
	string trig = PhotonTriggers.at(lt);
	addvar.push_back("passes"+trig);
	//	addvar.push_back("prescale"+trig);
      }
      for (size_t ht = 0; ht < hadronicTriggers.size(); ++ht)  {
	string trig = hadronicTriggers.at(ht);
	addvar.push_back("passes"+trig);
	//	addvar.push_back("prescale"+trig);
      }
      addvar.push_back("passesSingleElTriggers");
      addvar.push_back("passesSingleMuTriggers");
      addvar.push_back("passesPhotonTriggers");
      addvar.push_back("passesHadronicTriggers");
    }
  }

  //--- Soureek adding PU info -----------------    
  //if(doPU_){
  //   addvar.push_back("puWeight");
  //    addvar.push_back("puWeightUp");
  //    addvar.push_back("puWeightDown");
  //    addvar.push_back("nTruePU");
  //  }
  return addvar;
}

double DMAnalysisTreeMaker::getMCTagEfficiencyFunc(float flavor, float btag, float ptCorr, float eta, string algo,string syst, bool norm, bool spline ){
  bool doSpline = spline;
  double scale=1.0;
  if (algo == "CSV"|| algo=="CMVA"){
    double btagmin=getWPAlgo(algo,"MIN"),btagmax=getWPAlgo(algo,"MAX");
    double thrloose=getWPAlgo(algo,"L"),thrmedium=getWPAlgo(algo,"M"),thrtight=getWPAlgo(algo,"T");
    if(btag<thrtight)doSpline=false;
    if( doSpline){
      Double_t xs[5] = {btagmin,thrloose,thrmedium,thrtight,btagmax};
      Double_t ys[5] = {1,MCTagEfficiency(algo+"L",flavor, ptCorr,eta),MCTagEfficiency(algo+"M",flavor, ptCorr,eta),MCTagEfficiency(algo+"T",flavor, ptCorr,eta),0};
      Double_t endpointder=0.0;
      if(fabs(flavor)==5){
	endpointder= -(MCTagEfficiency(algo+"T",flavor, ptCorr,eta))/(btagmax-thrtight);
	//cout << "endpointder "<<endpointer<<endl;
	//	endpointder=-1;
      }
      //      cout <<" xs[5]={";
      //      for( int i = 0; i < 5; ++i){cout << xs[i]<<", ";}
      //      cout <<"}"<<endl<<" ys[5]={";
      //      for( int i = 0; i < 5; ++i){cout << ys[i]<<", ";}
      //      cout <<"}"<<endl;
      //      cout <<"endpointder="<<endpointder<<endl;
      
      TSpline3 sp("sp",xs,ys,5,"sp",0.0,endpointder);

      scale = sp.Eval(btag);

    }

    else{
      if (btag > btagmin && btag<=thrloose ){
	scale= (1. - (btag-btagmin)*(1-MCTagEfficiency(algo+"L",flavor, ptCorr,eta))/(thrloose-btagmin));
      }
      if (btag > thrloose && btag<=thrmedium){
	scale= (MCTagEfficiency(algo+"L",flavor, ptCorr,eta) - (btag-thrloose)*(MCTagEfficiency(algo+"L",flavor, ptCorr,eta)-MCTagEfficiency(algo+"M",flavor, ptCorr,eta))/(thrmedium-thrloose));
      }
      if (btag > thrmedium && btag<=thrtight){
	//      cout << "efffM "<< MCTagEfficiency(algo+"M",flavor, ptCorr,eta) << " btag "<< btag << " thrmed "<< thrmedium << " effT " <<MCTagEfficiency(algo+"T",flavor, ptCorr,eta)<< " pt "<<ptCorr << " eta "<< eta <<endl;
	scale = MCTagEfficiency(algo+"M",flavor, ptCorr,eta) - (btag-thrmedium)*(MCTagEfficiency(algo+"M",flavor, ptCorr,eta)-MCTagEfficiency(algo+"T",flavor, ptCorr,eta))/(thrtight-thrmedium);
      }
      if (btag > thrtight && btag <= btagmax){
	scale = MCTagEfficiency(algo+"T",flavor, ptCorr,eta) - (btag-thrtight)*(MCTagEfficiency(algo+"T",flavor, ptCorr,eta)-0)/(btagmax-thrtight);
      }
      if(btag<btagmin || btag > btagmax)scale= 1.;
    }
  }
  return scale;
}

//Helper to get a/b of the working point:
double DMAnalysisTreeMaker::getMCTagEfficiencyFuncParam(float flavor, float ptCorr, float eta, string algo,string syst, string param, string region){
  if (algo == "CSV"|| algo=="CMVA"){
    double btagmin=getWPAlgo(algo,"MIN"),btagmax=getWPAlgo(algo,"MAX");
    double thrloose=getWPAlgo(algo,"L"),thrmedium=getWPAlgo(algo,"M"),thrtight=getWPAlgo(algo,"T");
    if (param=="a"){
      if ( region == "0L") return (1. + (btagmin)*(1-MCTagEfficiency(algo+"L",flavor, ptCorr,eta))/(thrloose-btagmin));
      if ( region == "LM") return (MCTagEfficiency(algo+"L",flavor, ptCorr,eta) + (thrloose)*(MCTagEfficiency(algo+"L",flavor, ptCorr,eta)-MCTagEfficiency(algo+"M",flavor, ptCorr,eta))/(thrmedium-thrloose));
      if ( region == "MT")  return MCTagEfficiency(algo+"M",flavor, ptCorr,eta) + (thrmedium)*(MCTagEfficiency(algo+"M",flavor, ptCorr,eta)-MCTagEfficiency(algo+"T",flavor, ptCorr,eta))/(thrtight-thrmedium);
      if ( region == "T1") return MCTagEfficiency(algo+"T",flavor, ptCorr,eta) + (thrtight)*(MCTagEfficiency(algo+"T",flavor, ptCorr,eta)-0)/(btagmax-thrtight);
    }
    if (param=="b"){
      if ( region == "0L") return -((1-MCTagEfficiency(algo+"L",flavor, ptCorr,eta))/(thrloose-btagmin));
      if ( region == "LM") return -((MCTagEfficiency(algo+"L",flavor, ptCorr,eta)-MCTagEfficiency(algo+"M",flavor, ptCorr,eta))/(thrmedium-thrloose));
      if ( region == "MT") return -(MCTagEfficiency(algo+"M",flavor, ptCorr,eta)-MCTagEfficiency(algo+"T",flavor, ptCorr,eta))/(thrtight-thrmedium);
      if ( region == "T1") return -(MCTagEfficiency(algo+"T",flavor, ptCorr,eta)-0)/(btagmax-thrtight);
    }
  }
  return 0.;
}

double DMAnalysisTreeMaker::getMCTagEfficiencyFuncInt(float flavor, float ptCorr,float eta, string algo, string syst){
  double integral=0.;
  if (algo == "CSV"){
    integral =0.0;
    double btagmin=getWPAlgo(algo,"MIN"),btagmax=getWPAlgo(algo,"MAX");
    double thrloose=getWPAlgo(algo,"L"),thrmedium=getWPAlgo(algo,"M"),thrtight=getWPAlgo(algo,"T");
    integral+= (thrloose-btagmin)*(1+MCTagEfficiency(algo+"L",flavor, ptCorr,eta))/2.;
    integral+= (thrmedium-thrloose)*(MCTagEfficiency(algo+"M",flavor, ptCorr,eta)+MCTagEfficiency(algo+"L",flavor, ptCorr,eta))/2.;
    integral+= (thrtight-thrmedium)*(MCTagEfficiency(algo+"T",flavor, ptCorr,eta)+MCTagEfficiency(algo+"M",flavor, ptCorr,eta))/2.;
    integral+= (btagmax-thrtight)*(MCTagEfficiency(algo+"T",flavor, ptCorr,eta)+0)/2.;
  } 
  return integral;
}


double DMAnalysisTreeMaker::getPhi(double px, double py){
  double phi = -10;
  if(px<0){
      if(py>0)phi = atan(py/px)+3.141592;
      if(py<0)phi = atan(py/px)-3.141592;
  }
  else  phi= (atan(py/px));
  return phi;
}

float DMAnalysisTreeMaker::getReshapedBTagValue(float flavor, float btag, float pt, float eta, string algo, string syst){
  if (fabs(eta)>2.4)return 1.;
  if (syst=="noSyst"){
    if(algo == "CSV_sd") {
      if(fabs(flavor)!=5) return 1.;
      double effRatio = getEffRatioFunc(flavor,btag,pt,eta,"CSV",syst,false);
      //      cout<< effRatio <<" effRatio "<<flavor << " flavor "<<endl;
      return effRatio;
    }
    if(algo == "CMVA_sd") {
      if(fabs(flavor)!=5) return 1.;
      double effRatio = getEffRatioFunc(flavor,btag,pt,eta,"CMVA",syst,false);
      //      cout<< effRatio <<" effRatio "<<flavor << " flavor "<<endl;
      return effRatio;
    }
    if(algo=="CSV"){
      if(fabs(flavor) == 5) return readerCSVReshape->eval_auto_bounds("central",BTagEntry::FLAV_B, eta, pt, btag);
      if(fabs(flavor) == 4) return readerCSVReshape->eval_auto_bounds("central",BTagEntry::FLAV_C, eta, pt, btag);
      if(fabs(flavor) != 4 && fabs(flavor) != 5) return readerCSVReshape->eval_auto_bounds("central",BTagEntry::FLAV_UDSG, eta, pt, btag);
    }
    if(algo=="CMVA"){
      if(fabs(flavor) == 5) return readerCMVAReshape->eval_auto_bounds("central",BTagEntry::FLAV_B, eta, pt, btag);
      if(fabs(flavor) == 4) return readerCMVAReshape->eval_auto_bounds("central",BTagEntry::FLAV_C, eta, pt, btag);
      if(fabs(flavor) != 4 && fabs(flavor) != 5) return readerCMVAReshape->eval_auto_bounds("central",BTagEntry::FLAV_UDSG, eta, pt, btag);
    }
  }
  
  return 1.;
}

float DMAnalysisTreeMaker::getWPAlgo(string algo, string wp){
  if(algo=="CSV"){
    if(wp=="MIN")return 0.0;
    if(wp=="L")return 0.5426;
    if(wp=="M")return 0.8484;
    if(wp=="T")return 0.9535;
    if(wp=="MAX")return 0.999999;
  }
  if(algo=="CMVA"){
    if(wp=="MIN")return -1.;
    if(wp=="L")return -0.5884;
    if(wp=="M")return 0.4432;
    if(wp=="T")return 0.9432;
    if(wp=="MAX")return 0.999999;
  }
  return -1.0;
}
//double DMAnalysisTreeMaker::getMCTagEfficiencyFunc(float flavor, float ptCorr, string algo,string syst, string wp, string region){
//float avgRatioFunc(float x0,float x1,float a1, string b1 , float a2, float b2){

float DMAnalysisTreeMaker::getEffRatioFunc(float flavor,float btag,float pt,float eta, string algo ,string syst,bool normalize){

  float localvalue=-1;
  if(algo=="CSV" || algo == "CMVA"){
    double btagmin=getWPAlgo(algo,"MIN"),btagmax=getWPAlgo(algo,"MAX");
    double thrloose=getWPAlgo(algo,"L"),thrmedium=getWPAlgo(algo,"M"),thrtight=getWPAlgo(algo,"T");
    //    cout << " btagmin "<<btagmin<<" btag "<< btag <<endl;
    
    //    localvalue= getMCTagEfficiencyFunc(1,btag,pt, eta, algo,syst,false,false)/getMCTagEfficiencyFunc(flavor,btag,pt,eta,algo,syst,false,false);//Get the MC efficiency ratio
   
    //Spline version:
    localvalue= getMCTagEfficiencyFunc(1,btag,pt, eta, algo,syst,false,true)/getMCTagEfficiencyFunc(flavor,btag,pt,eta,algo,syst,false,true);//Get the MC efficiency ratio
    //    cout << " btag "<< btag<<endl;
    //    cout<< " eff l "<< getMCTagEfficiencyFunc(1,btag,pt, eta, algo,syst,false,false)<< " eff b "<< getMCTagEfficiencyFunc(flavor,btag,pt,eta,algo,syst,false,false) <<endl;
    //    cout<< " lvl "<< localvalue <<endl;
    if( normalize){ 
      //Use analytically evaluated function average to get the reweighting:
      float average=0.;
      //      float valMin= 1;

      float numA0L= getMCTagEfficiencyFuncParam(1, pt, eta, algo,syst, "a",  "0L");
      float numB0L= getMCTagEfficiencyFuncParam(1, pt, eta, algo,syst, "b",  "0L");

      float numALM= getMCTagEfficiencyFuncParam(1, pt, eta, algo,syst, "a",  "LM");
      float numBLM= getMCTagEfficiencyFuncParam(1, pt, eta, algo,syst, "b",  "LM");

      float numAMT= getMCTagEfficiencyFuncParam(1, pt, eta, algo,syst, "a",  "MT");
      float numBMT= getMCTagEfficiencyFuncParam(1, pt, eta, algo,syst, "b",  "MT");

      float numAT1= getMCTagEfficiencyFuncParam(1, pt, eta, algo,syst, "a",  "T1");
      float numBT1= getMCTagEfficiencyFuncParam(1, pt, eta, algo,syst, "b",  "T1");

      float denA0L= getMCTagEfficiencyFuncParam(flavor, pt, eta, algo,syst, "a",  "0L");
      float denB0L= getMCTagEfficiencyFuncParam(flavor, pt, eta, algo,syst, "b",  "0L");

      float denALM= getMCTagEfficiencyFuncParam(flavor, pt, eta, algo,syst, "a",  "LM");
      float denBLM= getMCTagEfficiencyFuncParam(flavor, pt, eta, algo,syst, "b",  "LM");

      float denAMT= getMCTagEfficiencyFuncParam(flavor, pt, eta, algo,syst, "a",  "MT");
      float denBMT= getMCTagEfficiencyFuncParam(flavor, pt, eta, algo,syst, "b",  "MT");

      float denAT1= getMCTagEfficiencyFuncParam(flavor, pt, eta, algo,syst, "a",  "T1");
      float denBT1= getMCTagEfficiencyFuncParam(flavor, pt, eta, algo,syst, "b",  "T1");
     
      float averagetot =average;
      //      cout << " localvalue before "<< localvalue << endl;
      //average = -denB0L*avgRatioFunc(btagmin,thrloose,numA0L,numB0L,denA0L,denB0L);
      average = avgRatioFunc(btagmin,thrloose,numA0L,numB0L,denA0L,denB0L);
      //            cout<< " average 0L " <<average<<endl;
      averagetot +=average;
	    
      average = -denBLM*avgRatioFunc(thrloose,thrmedium,numALM,numBLM,denALM,denBLM);
      average = avgRatioFunc(thrloose,thrmedium,numALM,numBLM,denALM,denBLM);
      //      cout<< " average LM " <<average<<endl;
      averagetot +=average;

      //            average = -denBMT*avgRatioFunc(thrmedium,thrtight,numAMT,numBMT,denAMT,denBMT);
      average = avgRatioFunc(thrmedium,thrtight,numAMT,numBMT,denAMT,denBMT);
      //      cout<< " average MT " <<average<<endl;
      averagetot +=average;
      
      //      average = -denBT1*avgRatioFunc(thrtight,btagmax,numAT1,numBT1,denAT1,denBT1);
      //average = avgRatioFunc(thrtight,btagmax,numAT1,numBT1,denAT1,denBT1); Not-well defined because AT1 == BT1, use simplified integral

      average = numAT1/denAT1*(btagmax-thrtight);

      //      cout << " localvalue before "<< localvalue << endl;
      //      cout<< " average is " <<average<<endl;
      //      cout << "tight "<< thrtight <<" max "<< btagmax<< " numaT1 "<< numAT1<<" numBT1 "<<  numBT1<< " denAT1 "<< denAT1<< " denBT1 "<< denBT1<<endl;
      averagetot +=average;

      //      cout<< " average tot is " <<averagetot<<endl;
      localvalue=localvalue/averagetot;
      //      cout << " localvalue after "<< localvalue << endl;
    }
  }
  
  return localvalue;
 
  //		  if normalize{}
}
 
float DMAnalysisTreeMaker::avgRatioFunc(float x0,float x1,float anum, float bnum , float aden, float bden){
  float p2=(bnum/(2*bden)) * (x1*x1 -x0*x0);
  float p1=((anum*bden-aden*bnum)/(bden*bden)) * (x1-x0);
  float p0=(-aden*(anum*bden-aden*bnum)/(bden*bden*bden))*(log(fabs(bden*x1+aden))-log(fabs(bden*x0+aden)));
  
  //float p2=(bnum/3)* (x1*x1*x1 -x0*x0*x0);
  //  float p1= (anum/2) * (x1*x1 -x0*x0);
  //  float p0=0;

  return p0+p1+p2;
}

void DMAnalysisTreeMaker::setEventBTagSF(string label, string category, string algo){
  
  int ncsv_tmp_t_tags =0;
  int ncsv_tmp_m_tags =0;
  int ncsv_tmp_l_tags =0;
  string lc=label+category;

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

  
  for(int j =0; j<sizes[label+category];++j){
    if(fabs(vfloats_values[lc+"_Eta"][j])>2.4)continue;
    if(vfloats_values[lc+"_Is"+algo+"T"][j])++ncsv_tmp_t_tags;
    if(vfloats_values[lc+"_Is"+algo+"M"][j])++ncsv_tmp_m_tags;
    if(vfloats_values[lc+"_Is"+algo+"L"][j])++ncsv_tmp_l_tags;

    float ptCorr = vfloats_values[lc+"_CorrPt"][j];
    float eta = vfloats_values[lc+"_Eta"][j];
    //    int flavor = vfloats_values[lc+"_PartonFlavour"][j];
    int flavor = vfloats_values[lc+"_HadronFlavour"][j];
    
    //    cout <<"jet "<< j<< " istagged l "<<vfloats_values[lc+"_IsCSVL"][j]<< " flavor "<< flavor<<endl;
    
    double csvteff = MCTagEfficiency(algo+"T",flavor, ptCorr,eta);
    double sfcsvt = TagScaleFactor(algo+"T", flavor, "noSyst", ptCorr,eta);
    
    double csvleff = MCTagEfficiency(algo+"L",flavor,ptCorr,eta);
    double sfcsvl = TagScaleFactor(algo+"L", flavor, "noSyst", ptCorr,eta);
    
    double csvmeff = MCTagEfficiency(algo+"M",flavor,ptCorr,eta);
    double sfcsvm = TagScaleFactor(algo+"M", flavor, "noSyst", ptCorr,eta);
    
    if(doTopDecayReshaping){
      //double leadingLeptonCharge=+1;//leadingleptonCharge;
      double pFlavour=vfloats_values[lc+"_PartonFlavour"][j];
      double product=topCharge*pFlavour;
      if(fabs(pFlavour)==5){
	//	cout << "jet j"<< j<< " pflav "<< pFlavour <<" sf bef "<<sfcsvt <<endl;
	if(fabs(pFlavour)==5 && product >0){
	  sfcsvt = sfcsvt*MCTagEfficiency(algo+"T",1,ptCorr,eta)/csvteff;
	  sfcsvm = sfcsvm*MCTagEfficiency(algo+"M",1,ptCorr,eta)/csvmeff;
	  sfcsvl = sfcsvl*MCTagEfficiency(algo+"L",1,ptCorr,eta)/csvleff;
	}
	//	cout << "jet j"<< j<< " pflav "<< pFlavour <<" sf aft "<<sfcsvt <<endl;
      }
    }
    
    double sfcsvt_mistag_up = TagScaleFactor(algo+"T", flavor, "mistag_up", ptCorr,eta);
    double sfcsvl_mistag_up = TagScaleFactor(algo+"L", flavor, "mistag_up", ptCorr,eta);
    double sfcsvm_mistag_up = TagScaleFactor(algo+"M", flavor, "mistag_up", ptCorr,eta);
    
    double sfcsvt_mistag_down = TagScaleFactor(algo+"T", flavor, "mistag_down", ptCorr,eta);
    double sfcsvl_mistag_down = TagScaleFactor(algo+"L", flavor, "mistag_down", ptCorr,eta);
    double sfcsvm_mistag_down = TagScaleFactor(algo+"M", flavor, "mistag_down", ptCorr,eta);
    
    double sfcsvt_b_tag_down = TagScaleFactor(algo+"T", flavor, "b_tag_down", ptCorr,eta);
    double sfcsvl_b_tag_down = TagScaleFactor(algo+"L", flavor, "b_tag_down", ptCorr,eta);
    double sfcsvm_b_tag_down = TagScaleFactor(algo+"M", flavor, "b_tag_down", ptCorr,eta);
    
    double sfcsvt_b_tag_up = TagScaleFactor(algo+"T", flavor, "b_tag_up", ptCorr,eta);
    double sfcsvl_b_tag_up = TagScaleFactor(algo+"L", flavor, "b_tag_up", ptCorr,eta);
    double sfcsvm_b_tag_up = TagScaleFactor(algo+"M", flavor, "b_tag_up", ptCorr,eta);
    
    
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
  
  //BTagging part
  //if(doBTagSF){
  if(doBTagSF){
    

    //CSVT
    //0 tags
    b_weight_csvt_0_tags = b_csvt_0_tags.weight(jsfscsvt, ncsv_tmp_t_tags);  
    b_weight_csvt_0_tags_mistag_up = b_csvt_0_tags.weight(jsfscsvt_mistag_up, ncsv_tmp_t_tags);  
    b_weight_csvt_0_tags_mistag_down = b_csvt_0_tags.weight(jsfscsvt_mistag_down, ncsv_tmp_t_tags);  
    b_weight_csvt_0_tags_b_tag_up = b_csvt_0_tags.weight(jsfscsvt_b_tag_up, ncsv_tmp_t_tags);  
    b_weight_csvt_0_tags_b_tag_down = b_csvt_0_tags.weight(jsfscsvt_b_tag_down, ncsv_tmp_t_tags);  
    
    //1 tag
    b_weight_csvt_1_tag = b_csvt_1_tag.weight(jsfscsvt, ncsv_tmp_t_tags);  
    b_weight_csvt_1_tag_mistag_up = b_csvt_1_tag.weight(jsfscsvt_mistag_up, ncsv_tmp_t_tags);  
    b_weight_csvt_1_tag_mistag_down = b_csvt_1_tag.weight(jsfscsvt_mistag_down, ncsv_tmp_t_tags);  
    b_weight_csvt_1_tag_b_tag_up = b_csvt_1_tag.weight(jsfscsvt_b_tag_up, ncsv_tmp_t_tags);  
    b_weight_csvt_1_tag_b_tag_down = b_csvt_1_tag.weight(jsfscsvt_b_tag_down, ncsv_tmp_t_tags);  
    //      cout <<"w1t check: is"<< b_weight_csvt_1_tag<<endl;
    
    //2 tags
    b_weight_csvt_2_tags = b_csvt_2_tags.weight(jsfscsvt, ncsv_tmp_t_tags);  
    b_weight_csvt_2_tags_mistag_up = b_csvt_2_tags.weight(jsfscsvt_mistag_up, ncsv_tmp_t_tags);  
    b_weight_csvt_2_tags_mistag_down = b_csvt_2_tags.weight(jsfscsvt_mistag_down, ncsv_tmp_t_tags);  
    b_weight_csvt_2_tags_b_tag_up = b_csvt_2_tags.weight(jsfscsvt_b_tag_up, ncsv_tmp_t_tags);  
    b_weight_csvt_2_tags_b_tag_down = b_csvt_2_tags.weight(jsfscsvt_b_tag_down, ncsv_tmp_t_tags);  
    
    //1-2 tags
    b_weight_csvt_1_2_tags = b_csvt_1_2_tags.weight(jsfscsvt, ncsv_tmp_t_tags);  
    b_weight_csvt_1_2_tags_b_tag_up = b_csvt_1_2_tags.weight(jsfscsvt_b_tag_up, ncsv_tmp_t_tags);  
    b_weight_csvt_1_2_tags_b_tag_down = b_csvt_1_2_tags.weight(jsfscsvt_b_tag_down, ncsv_tmp_t_tags);  
    b_weight_csvt_1_2_tags_mistag_up = b_csvt_1_2_tags.weight(jsfscsvt_mistag_up, ncsv_tmp_t_tags);  
    b_weight_csvt_1_2_tags_mistag_down = b_csvt_1_2_tags.weight(jsfscsvt_mistag_down, ncsv_tmp_t_tags);  
    
    
    //CSVM
    //0 tags
    b_weight_csvm_0_tags = b_csvm_0_tags.weight(jsfscsvm, ncsv_tmp_m_tags);  
    b_weight_csvm_0_tags_mistag_up = b_csvm_0_tags.weight(jsfscsvm_mistag_up, ncsv_tmp_m_tags);  
    b_weight_csvm_0_tags_mistag_down = b_csvm_0_tags.weight(jsfscsvm_mistag_down, ncsv_tmp_m_tags);  
    b_weight_csvm_0_tags_b_tag_up = b_csvm_0_tags.weight(jsfscsvm_b_tag_up, ncsv_tmp_m_tags);  
    b_weight_csvm_0_tags_b_tag_down = b_csvm_0_tags.weight(jsfscsvm_b_tag_down, ncsv_tmp_m_tags);  
    
    //1 tag
    b_weight_csvm_1_tag = b_csvm_1_tag.weight(jsfscsvm, ncsv_tmp_m_tags);  
    b_weight_csvm_1_tag_mistag_up = b_csvm_1_tag.weight(jsfscsvm_mistag_up, ncsv_tmp_m_tags);  
    b_weight_csvm_1_tag_mistag_down = b_csvm_1_tag.weight(jsfscsvm_mistag_down, ncsv_tmp_m_tags);  
    b_weight_csvm_1_tag_b_tag_up = b_csvm_1_tag.weight(jsfscsvm_b_tag_up, ncsv_tmp_m_tags);  
    b_weight_csvm_1_tag_b_tag_down = b_csvm_1_tag.weight(jsfscsvm_b_tag_down, ncsv_tmp_m_tags);  
    //      cout <<"w1t check: is"<< b_weight_csvm_1_tag<<endl;
    
    //2 tags
    b_weight_csvm_2_tags = b_csvm_2_tags.weight(jsfscsvm, ncsv_tmp_m_tags);  
    b_weight_csvm_2_tags_mistag_up = b_csvm_2_tags.weight(jsfscsvm_mistag_up, ncsv_tmp_m_tags);  
    b_weight_csvm_2_tags_mistag_down = b_csvm_2_tags.weight(jsfscsvm_mistag_down, ncsv_tmp_m_tags);  
    b_weight_csvm_2_tags_b_tag_up = b_csvm_2_tags.weight(jsfscsvm_b_tag_up, ncsv_tmp_m_tags);  
    b_weight_csvm_2_tags_b_tag_down = b_csvm_2_tags.weight(jsfscsvm_b_tag_down, ncsv_tmp_m_tags);  
    
    //1-2 tags
    b_weight_csvm_1_2_tags = b_csvm_1_2_tags.weight(jsfscsvm, ncsv_tmp_m_tags);  
    b_weight_csvm_1_2_tags_b_tag_up = b_csvm_1_2_tags.weight(jsfscsvm_b_tag_up, ncsv_tmp_m_tags);  
    b_weight_csvm_1_2_tags_b_tag_down = b_csvm_1_2_tags.weight(jsfscsvm_b_tag_down, ncsv_tmp_m_tags);  
    b_weight_csvm_1_2_tags_mistag_up = b_csvm_1_2_tags.weight(jsfscsvm_mistag_up, ncsv_tmp_m_tags);  
    b_weight_csvm_1_2_tags_mistag_down = b_csvm_1_2_tags.weight(jsfscsvm_mistag_down, ncsv_tmp_m_tags);  
    
    
    //CSVL
    //0 tags
    b_weight_csvl_0_tags = b_csvl_0_tags.weight(jsfscsvl, ncsv_tmp_l_tags);  
    b_weight_csvl_0_tags_mistag_up = b_csvl_0_tags.weight(jsfscsvl_mistag_up, ncsv_tmp_l_tags);  
    b_weight_csvl_0_tags_mistag_down = b_csvl_0_tags.weight(jsfscsvl_mistag_down, ncsv_tmp_l_tags);  
    b_weight_csvl_0_tags_b_tag_up = b_csvl_0_tags.weight(jsfscsvl_b_tag_up, ncsv_tmp_l_tags);  
    b_weight_csvl_0_tags_b_tag_down = b_csvl_0_tags.weight(jsfscsvl_b_tag_down, ncsv_tmp_l_tags);  
    
    //1 tag
    b_weight_csvl_1_tag = b_csvl_1_tag.weight(jsfscsvl, ncsv_tmp_l_tags);  
    b_weight_csvl_1_tag_mistag_up = b_csvl_1_tag.weight(jsfscsvl_mistag_up, ncsv_tmp_l_tags);  
    b_weight_csvl_1_tag_mistag_down = b_csvl_1_tag.weight(jsfscsvl_mistag_down, ncsv_tmp_l_tags);  
    b_weight_csvl_1_tag_b_tag_up = b_csvl_1_tag.weight(jsfscsvl_b_tag_up, ncsv_tmp_l_tags);  
    b_weight_csvl_1_tag_b_tag_down = b_csvl_1_tag.weight(jsfscsvl_b_tag_down, ncsv_tmp_l_tags);  
    //      cout <<"w1t check: is"<< b_weight_csvl_1_tag<<endl;
    
    //2 tags
    b_weight_csvl_2_tags = b_csvl_2_tags.weight(jsfscsvl, ncsv_tmp_l_tags);  
    b_weight_csvl_2_tags_mistag_up = b_csvl_2_tags.weight(jsfscsvl_mistag_up, ncsv_tmp_l_tags);  
    b_weight_csvl_2_tags_mistag_down = b_csvl_2_tags.weight(jsfscsvl_mistag_down, ncsv_tmp_l_tags);  
    b_weight_csvl_2_tags_b_tag_up = b_csvl_2_tags.weight(jsfscsvl_b_tag_up, ncsv_tmp_l_tags);  
    b_weight_csvl_2_tags_b_tag_down = b_csvl_2_tags.weight(jsfscsvl_b_tag_down, ncsv_tmp_l_tags);  
    
    //1-2 tags
    b_weight_csvl_1_2_tags = b_csvl_1_2_tags.weight(jsfscsvl, ncsv_tmp_l_tags);  
    b_weight_csvl_1_2_tags_b_tag_up = b_csvl_1_2_tags.weight(jsfscsvl_b_tag_up, ncsv_tmp_l_tags);  
    b_weight_csvl_1_2_tags_b_tag_down = b_csvl_1_2_tags.weight(jsfscsvl_b_tag_down, ncsv_tmp_l_tags);  
    b_weight_csvl_1_2_tags_mistag_up = b_csvl_1_2_tags.weight(jsfscsvl_mistag_up, ncsv_tmp_l_tags);  
    b_weight_csvl_1_2_tags_mistag_down = b_csvl_1_2_tags.weight(jsfscsvl_mistag_down, ncsv_tmp_l_tags);  
    
    //    cout << " n tight tags "<< ncsv_tmp_l_tags  << " w0tag "<< b_weight_csvl_0_tags<<" w1tag " << b_weight_csvl_1_tag <<" w2tags "<<b_weight_csvt_2_tags  <<endl;

    
    float_values["Event_bWeight0"+algo+"L"+category]=b_weight_csvl_0_tags;
    float_values["Event_bWeight1"+algo+"L"+category]=b_weight_csvl_1_tag;
    float_values["Event_bWeight2"+algo+"L"+category]=b_weight_csvl_2_tags;
    float_values["Event_bWeight1_2"+algo+"L"+category]=b_weight_csvl_1_2_tags;
    
    float_values["Event_bWeight0"+algo+"M"+category]=b_weight_csvm_0_tags;
    float_values["Event_bWeight1"+algo+"M"+category]=b_weight_csvm_1_tag;
    float_values["Event_bWeight2"+algo+"M"+category]=b_weight_csvm_2_tags;
    float_values["Event_bWeight1_2"+algo+"M"+category]=b_weight_csvm_1_2_tags;
    
    float_values["Event_bWeight0"+algo+"T"+category]=b_weight_csvt_0_tags;
    float_values["Event_bWeight1"+algo+"T"+category]=b_weight_csvt_1_tag;
    float_values["Event_bWeight2"+algo+"T"+category]=b_weight_csvt_2_tags;
    float_values["Event_bWeight1_2"+algo+"T"+category]=b_weight_csvt_1_2_tags;
    
    //Mistag
    float_values["Event_bWeightMisTagUp0"+algo+"L"+category]=b_weight_csvl_0_tags_mistag_up;
    float_values["Event_bWeightMisTagUp1"+algo+"L"+category]=b_weight_csvl_1_tag_mistag_up;
    float_values["Event_bWeightMisTagUp2"+algo+"L"+category]=b_weight_csvl_2_tags_mistag_up;
    float_values["Event_bWeightMisTagUp1_2"+algo+"L"+category]=b_weight_csvl_1_2_tags_mistag_up;
    
    float_values["Event_bWeightMisTagUp0"+algo+"M"+category]=b_weight_csvm_0_tags_mistag_up;
    float_values["Event_bWeightMisTagUp1"+algo+"M"+category]=b_weight_csvm_1_tag_mistag_up;
    float_values["Event_bWeightMisTagUp2"+algo+"M"+category]=b_weight_csvm_2_tags_mistag_up;
    float_values["Event_bWeightMisTagUp1_2"+algo+"M"+category]=b_weight_csvm_1_2_tags_mistag_up;
    
    float_values["Event_bWeightMisTagUp0"+algo+"T"+category]=b_weight_csvt_0_tags_mistag_up;
    float_values["Event_bWeightMisTagUp1"+algo+"T"+category]=b_weight_csvt_1_tag_mistag_up;
    float_values["Event_bWeightMisTagUp2"+algo+"T"+category]=b_weight_csvt_2_tags_mistag_up;
    float_values["Event_bWeightMisTagUp1_2"+algo+"T"+category]=b_weight_csvt_1_2_tags_mistag_up;
    
    float_values["Event_bWeightMisTagDown0"+algo+"L"+category]=b_weight_csvl_0_tags_mistag_down;
    float_values["Event_bWeightMisTagDown1"+algo+"L"+category]=b_weight_csvl_1_tag_mistag_down;
    float_values["Event_bWeightMisTagDown2"+algo+"L"+category]=b_weight_csvl_2_tags_mistag_down;
    float_values["Event_bWeightMisTagDown1_2"+algo+"L"+category]=b_weight_csvl_1_2_tags_mistag_down;
    
    float_values["Event_bWeightMisTagDown0"+algo+"M"+category]=b_weight_csvm_0_tags_mistag_down;
    float_values["Event_bWeightMisTagDown1"+algo+"M"+category]=b_weight_csvm_1_tag_mistag_down;
    float_values["Event_bWeightMisTagDown2"+algo+"M"+category]=b_weight_csvm_2_tags_mistag_down;
    float_values["Event_bWeightMisTagDown1_2"+algo+"M"+category]=b_weight_csvm_1_2_tags_mistag_down;
    
    float_values["Event_bWeightMisTagDown0"+algo+"T"+category]=b_weight_csvt_0_tags_mistag_down;
    float_values["Event_bWeightMisTagDown1"+algo+"T"+category]=b_weight_csvt_1_tag_mistag_down;
    float_values["Event_bWeightMisTagDown2"+algo+"T"+category]=b_weight_csvt_2_tags_mistag_down;
    float_values["Event_bWeightMisTagDown1_2"+algo+"T"+category]=b_weight_csvt_1_2_tags_mistag_down;
    
    //Btag
    float_values["Event_bWeightBTagUp0"+algo+"L"+category]=b_weight_csvl_0_tags_b_tag_up;
    float_values["Event_bWeightBTagUp1"+algo+"L"+category]=b_weight_csvl_1_tag_b_tag_up;
    float_values["Event_bWeightBTagUp2"+algo+"L"+category]=b_weight_csvl_2_tags_b_tag_up;
    float_values["Event_bWeightBTagUp1_2"+algo+"L"+category]=b_weight_csvl_1_2_tags_b_tag_up;
    
    float_values["Event_bWeightBTagUp0"+algo+"M"+category]=b_weight_csvm_0_tags_b_tag_up;
    float_values["Event_bWeightBTagUp1"+algo+"M"+category]=b_weight_csvm_1_tag_b_tag_up;
    float_values["Event_bWeightBTagUp2"+algo+"M"+category]=b_weight_csvm_2_tags_b_tag_up;
    float_values["Event_bWeightBTagUp1_2"+algo+"M"+category]=b_weight_csvm_1_2_tags_b_tag_up;
    
    float_values["Event_bWeightBTagUp0"+algo+"T"+category]=b_weight_csvt_0_tags_b_tag_up;
    float_values["Event_bWeightBTagUp1"+algo+"T"+category]=b_weight_csvt_1_tag_b_tag_up;
    float_values["Event_bWeightBTagUp2"+algo+"T"+category]=b_weight_csvt_2_tags_b_tag_up;
    float_values["Event_bWeightBTagUp1_2"+algo+"T"+category]=b_weight_csvt_1_2_tags_b_tag_up;
    
    float_values["Event_bWeightBTagDown0"+algo+"L"+category]=b_weight_csvl_0_tags_b_tag_down;
    float_values["Event_bWeightBTagDown1"+algo+"L"+category]=b_weight_csvl_1_tag_b_tag_down;
    float_values["Event_bWeightBTagDown2"+algo+"L"+category]=b_weight_csvl_2_tags_b_tag_down;
    float_values["Event_bWeightBTagDown1_2"+algo+"L"+category]=b_weight_csvl_1_2_tags_b_tag_down;
    
    float_values["Event_bWeightBTagDown0"+algo+"M"+category]=b_weight_csvm_0_tags_b_tag_down;
    float_values["Event_bWeightBTagDown1"+algo+"M"+category]=b_weight_csvm_1_tag_b_tag_down;
    float_values["Event_bWeightBTagDown2"+algo+"M"+category]=b_weight_csvm_2_tags_b_tag_down;
    float_values["Event_bWeightBTagDown1_2"+algo+"M"+category]=b_weight_csvm_1_2_tags_b_tag_down;
    
    float_values["Event_bWeightBTagDown0"+algo+"T"+category]=b_weight_csvt_0_tags_b_tag_down;
    float_values["Event_bWeightBTagDown1"+algo+"T"+category]=b_weight_csvt_1_tag_b_tag_down;
    float_values["Event_bWeightBTagDown2"+algo+"T"+category]=b_weight_csvt_2_tags_b_tag_down;
    float_values["Event_bWeightBTagDown1_2"+algo+"T"+category]=b_weight_csvt_1_2_tags_b_tag_down;
  }
}

void DMAnalysisTreeMaker::setEventLeptonSF(string label, string category){
  ;
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
  //EWK
  //  [100, 150]:  0.980859
  //  [150, 200]:  0.962119
  //  [200, 250]:  0.944429
  //  [250, 300]:  0.927686
  //  [300, 350]:  0.911802
  //  [350, 400]:  0.8967  
  //  [400, 500]:  0.875368
  //  [500, 600]:  0.849097
  //  [600, 1000]: 0.792159

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

//  [100,150]:  0.984525
//  [150,200]:  0.969079
//   [200,250]:  0.954627
//   [250,300]:  0.941059
//    [300,350]:  0.928284
//    [350,400]:  0.91622
//    [400,500]:  0.899312
//    [500,600]:  0.878693
//    [600,1000]: 0.834718

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
      //      double a = 0.156;// 7/8 TeV
      //      double b = -0.00137;//
      double a = 0.0615;//
      double b = -0.0005;//
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
  bool eleOR=false, muOR=false, hadronOR=false, phOR=false;
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

  for(size_t ht =0; ht< hadronicTriggers.size();++ht){
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
  
  float_values["Event_passesSingleElTriggers"]=(float)eleOR;
  float_values["Event_passesSingleMuTriggers"]=(float)muOR;
  float_values["Event_passesPhotonTriggers"]=(float)phOR;
  float_values["Event_passesHadronicTriggers"]=(float)hadronOR;
  return (eleOR || muOR || hadronOR || phOR);
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

  //  std::cout << "entering pdf loop" <<std::endl;
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
      
	    //cout << "ww # " << i<< "is "<<ww <<endl;
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
  //  cout << " preBranch "<<endl;
  
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

  trees["WeightHistory"]->Branch("Event_T_size",&float_values["Event_T_size"]);
  trees["WeightHistory"]->Branch("Event_Tbar_size",&float_values["Event_Tbar_size"]);


  trees["WeightHistory"]->Branch("Event_W_Pt",&float_values["Event_W_Pt"]);
  trees["WeightHistory"]->Branch("Event_Z_Pt",&float_values["Event_Z_Pt"]);
  trees["WeightHistory"]->Branch("Event_nTruePV",&float_values["Event_nTruePV"]);

  //  cout << " preBranch Weight "<<endl;

  //  size_t wgtsize=  lhes->weights().size();
  //  std::cout << "weight size "<< wgtsize<<endl;
  if(useLHEW){
    for (size_t w = 0; w <  (size_t)maxWeights; ++w)  {
      stringstream w_n;
      w_n << w;
      string name = "Event_LHEWeight"+w_n.str();
      cout << " pre single w # "<< w <<endl;
      trees["WeightHistory"]->Branch(name.c_str(),&float_values[name],(name+"/F").c_str());
      cout << " branched "<< float_values[name]<<endl;
      //trees["noSyst"]->Branch(name.c_str(), &float_values[name],(name+"/F").c_str());
    }
  }
}

// double DMAnalysisTreeMaker::smearPt(double ptCorr, double genpt, double eta, string syst){
//   double resolScale = resolSF(fabs(eta), syst);
//   double smear =1.0;
//   if(genpt>0) smear = std::max((double)(0.0), (double)(ptCorr + (ptCorr - genpt) * resolScale) / ptCorr);
//   return ptCorr * smear;
// }

double DMAnalysisTreeMaker::smear(double pt, double genpt, double eta, string syst){
  double resolScale = resolSF(fabs(eta), syst);
  double smear =1.0;
  if(genpt>0) smear = std::max((double)(0.0), (double)(pt + (pt - genpt) * resolScale) / pt);
  return  smear;
}



double DMAnalysisTreeMaker::resolSF(double eta, string syst)
{
  double fac = 0.;
  if (syst == "jer__up")fac = 1.;
  if (syst == "jer__down")fac = -1.;
  if (eta <= 0.5)                       return 0.122 + (0.026 * fac);
  else if ( eta > 0.5 && eta <= 0.8 )   return 0.167 + (0.048 * fac);
  else if ( eta > 0.8 && eta <= 1.1 )   return 0.168 + (0.046 * fac);
  else if ( eta > 1.1 && eta <= 1.3 )   return 0.029 + (0.066 * fac);
  else if ( eta > 1.3 && eta <= 1.7 )   return 0.115 + (0.030 * fac);
  else if ( eta > 1.7 && eta <= 1.9 )   return 0.141 + (0.062 * fac);
  else if ( eta > 1.9 && eta <= 2.1 )   return 0.167 + (0.086 * fac);
  else if ( eta > 2.1 && eta <= 2.3 )   return 0.094 + (0.093 * fac);
  else if ( eta > 2.3 && eta <= 2.5 )   return 0.168 + (0.120 * fac);
  else if ( eta > 2.5 && eta <= 2.8 )   return 0.266 + (0.132 * fac);
  else if ( eta > 2.8 && eta <= 3.0 )   return 0.595 + (0.175 * fac);
  else if ( eta > 3.0 && eta <= 3.2 )   return -0.002 + (0.066 * fac);
  else if ( eta > 3.2 && eta <= 5.0 )   return 0.226 + (0.145 * fac);
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
  if(syst == "JESUp" || syst == "JESDown"){
    double fac = 1.;
    if (syst == "JESDown")fac = -1.;
    jecUnc->setJetEta(eta);
    jecUnc->setJetPt(ptCorr);
    double JetCorrection = jecUnc->getUncertainty(true);
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
    //    std::cout << " label is " << s << " vector i-th element  "<< v.at(i)<<" are they equal? "<< (v.at(i)==s) << " is v in s? "<< (s.find(v.at(i))!=std::string::npos)<<endl;
    if(v.at(i)==s)return true;
    //    if(s.find(v.at(i))!=std::string::npos)return true;
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
        //   std::cout << "nThis event should not pass the selection, what is it doing here?" << std::endl;
        return 0;
    }
    int njetTags = jetTags.size();
    //    cout<< " njettags "<< njetTags<<endl;
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

double DMAnalysisTreeMaker::MCTagEfficiency(string algo, int flavor, double pt, double eta){
  if(pt < 40)pt = 40.1;
  if (abs(flavor) ==5){
    if(algo=="CSVT") return 0.51;
    if(algo=="CSVM") return 0.71;
    if(algo=="CSVL") return 0.86;
    if(algo=="CMVAT") return cmvaeffbt->getEff(fabs(eta),pt);
    if(algo=="CMVAM") return cmvaeffbm->getEff(fabs(eta),pt);
    if(algo=="CMVAL") return cmvaeffbl->getEff(fabs(eta),pt);

  }

  if (abs(flavor) ==4){
    if(algo=="CSVT") return 0.015;
    if(algo=="CSVM") return 0.08;
    if(algo=="CSVL") return 0.28;
    if(algo=="CMVAT") return cmvaeffct->getEff(fabs(eta),pt);
    if(algo=="CMVAM") return cmvaeffcm->getEff(fabs(eta),pt);
    if(algo=="CMVAL") return cmvaeffcl->getEff(fabs(eta),pt);
  }

  if (abs(flavor) !=4 && abs(flavor) !=5){
    if(algo=="CSVT") return 0.003;
    if(algo=="CSVM") return 0.02;
    if(algo=="CSVL") return 0.16;
    if(algo=="CMVAT") return 0.003; //cmvaeffot->getEff(fabs(eta),pt);
    if(algo=="CMVAM") return 0.02;//cmvaeffom->getEff(fabs(eta),pt);
    if(algo=="CMVAL") return 0.13; //cmvaeffol->getEff(fabs(eta),pt);
  }
  return 1.0;
}


double DMAnalysisTreeMaker::TagScaleFactor(string algo, int flavor, string syst, double pt, double eta){
  // source (02/11):
  // https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation76X/CSVv2_prelim.csv
  if(algo == "CSVL"){
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){return readerCSVLoose->eval_auto_bounds("central",BTagEntry::FLAV_B, eta, pt);    }
      if(abs(flavor)==4){return readerCSVLoose->eval_auto_bounds("central",BTagEntry::FLAV_C, eta, pt);    }
      if(abs(flavor)!=5 && abs(flavor)!=4){ return readerCSVLoose->eval_auto_bounds("central",BTagEntry::FLAV_UDSG, eta, pt);    }    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){	return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)==4){	return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){	return readerCSVLoose->eval_auto_bounds("up",BTagEntry::FLAV_UDSG, eta, pt);      }    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){	return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)==4){	return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){	return readerCSVLoose->eval_auto_bounds("down",BTagEntry::FLAV_UDSG, eta, pt);      }    }
    if(syst ==  "b_tag_up") {
      if(abs(flavor)==5){return readerCSVLoose->eval_auto_bounds("up",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCSVLoose->eval_auto_bounds("up",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return 1.0*TagScaleFactor(algo,flavor,"noSyst",pt);      }    }
    if(syst ==  "b_tag_down") {
      if(abs(flavor)==5){return readerCSVLoose->eval_auto_bounds("down",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCSVLoose->eval_auto_bounds("down",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return 1.0*TagScaleFactor(algo,flavor,"noSyst",pt);      }    }
  }
  if(algo == "CSVM"){
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){return readerCSVMedium->eval_auto_bounds("central",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCSVMedium->eval_auto_bounds("central",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return readerCSVMedium->eval_auto_bounds("central",BTagEntry::FLAV_UDSG, eta, pt);      }
    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)==4){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return readerCSVMedium->eval_auto_bounds("up",BTagEntry::FLAV_UDSG, eta, pt);      }    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)==4){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return readerCSVMedium->eval_auto_bounds("down",BTagEntry::FLAV_UDSG, eta, pt);      }   }
    if(syst ==  "b_tag_up") {
      if(abs(flavor)==5){return readerCSVMedium->eval_auto_bounds("up",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCSVMedium->eval_auto_bounds("up",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return 1.0*TagScaleFactor(algo,flavor,"noSyst",pt);      }     }
    if(syst ==  "b_tag_down") {
      if(abs(flavor)==5){return readerCSVMedium->eval_auto_bounds("down",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCSVMedium->eval_auto_bounds("down",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return 1.0*TagScaleFactor(algo,flavor,"noSyst",pt);      }    }  
  }
  if(algo == "CSVT"){
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){return readerCSVTight->eval_auto_bounds("central",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCSVTight->eval_auto_bounds("central",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return readerCSVTight->eval_auto_bounds("central",BTagEntry::FLAV_UDSG, eta, pt);      }    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)==4){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return readerCSVTight->eval_auto_bounds("up",BTagEntry::FLAV_UDSG, eta, pt);      }    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)==4){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return readerCSVTight->eval_auto_bounds("down",BTagEntry::FLAV_UDSG, eta, pt);      }    }
    if(syst ==  "b_tag_up") {
      if(abs(flavor)==5){return readerCSVTight->eval_auto_bounds("up",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCSVTight->eval_auto_bounds("up",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return 1.0*TagScaleFactor(algo,flavor,"noSyst",pt);      }    }
    if(syst ==  "b_tag_down") {
      if(abs(flavor)==5){return readerCSVTight->eval_auto_bounds("down",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCSVTight->eval_auto_bounds("down",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return 1.0*TagScaleFactor(algo,flavor,"noSyst",pt);      }    }
  }
if(algo == "CMVAL"){
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){return readerCMVALoose->eval_auto_bounds("central",BTagEntry::FLAV_B, eta, pt);    }
      if(abs(flavor)==4){return readerCMVALoose->eval_auto_bounds("central",BTagEntry::FLAV_C, eta, pt);    }
      if(abs(flavor)!=5 && abs(flavor)!=4){ return readerCMVALoose->eval_auto_bounds("central",BTagEntry::FLAV_UDSG, eta, pt);    }    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){	return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)==4){	return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){	return readerCMVALoose->eval_auto_bounds("up",BTagEntry::FLAV_UDSG, eta, pt);      }    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){	return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)==4){	return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){	return readerCMVALoose->eval_auto_bounds("down",BTagEntry::FLAV_UDSG, eta, pt);      }    }
    if(syst ==  "b_tag_up") {
      if(abs(flavor)==5){return readerCMVALoose->eval_auto_bounds("up",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCMVALoose->eval_auto_bounds("up",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return 1.0*TagScaleFactor(algo,flavor,"noSyst",pt);      }    }
    if(syst ==  "b_tag_down") {
      if(abs(flavor)==5){return readerCMVALoose->eval_auto_bounds("down",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCMVALoose->eval_auto_bounds("down",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return 1.0*TagScaleFactor(algo,flavor,"noSyst",pt);      }    }
  }
  if(algo == "CMVAM"){
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){return readerCMVAMedium->eval_auto_bounds("central",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCMVAMedium->eval_auto_bounds("central",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return readerCMVAMedium->eval_auto_bounds("central",BTagEntry::FLAV_UDSG, eta, pt);      }
    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)==4){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return readerCMVAMedium->eval_auto_bounds("up",BTagEntry::FLAV_UDSG, eta, pt);      }    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)==4){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return readerCMVAMedium->eval_auto_bounds("down",BTagEntry::FLAV_UDSG, eta, pt);      }   }
    if(syst ==  "b_tag_up") {
      if(abs(flavor)==5){return readerCMVAMedium->eval_auto_bounds("up",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCMVAMedium->eval_auto_bounds("up",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return 1.0*TagScaleFactor(algo,flavor,"noSyst",pt);      }     }
    if(syst ==  "b_tag_down") {
      if(abs(flavor)==5){return readerCMVAMedium->eval_auto_bounds("down",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCMVAMedium->eval_auto_bounds("down",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return 1.0*TagScaleFactor(algo,flavor,"noSyst",pt);      }    }  
  }
  if(algo == "CMVAT"){
    if(syst ==  "noSyst") {
      if(abs(flavor)==5){return readerCMVATight->eval_auto_bounds("central",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCMVATight->eval_auto_bounds("central",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return readerCMVATight->eval_auto_bounds("central",BTagEntry::FLAV_UDSG, eta, pt);      }    }
    if(syst ==  "mistag_up") {
      if(abs(flavor)==5){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)==4){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return readerCMVATight->eval_auto_bounds("up",BTagEntry::FLAV_UDSG, eta, pt);      }    }
    if(syst ==  "mistag_down") {
      if(abs(flavor)==5){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)==4){return 1.00*TagScaleFactor(algo,flavor,"noSyst",pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return readerCMVATight->eval_auto_bounds("down",BTagEntry::FLAV_UDSG, eta, pt);      }    }
    if(syst ==  "b_tag_up") {
      if(abs(flavor)==5){return readerCMVATight->eval_auto_bounds("up",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCMVATight->eval_auto_bounds("up",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return 1.0*TagScaleFactor(algo,flavor,"noSyst",pt);      }    }
    if(syst ==  "b_tag_down") {
      if(abs(flavor)==5){return readerCMVATight->eval_auto_bounds("down",BTagEntry::FLAV_B, eta, pt);      }
      if(abs(flavor)==4){return readerCMVATight->eval_auto_bounds("down",BTagEntry::FLAV_C, eta, pt);      }
      if(abs(flavor)!=5 && abs(flavor)!=4){return 1.0*TagScaleFactor(algo,flavor,"noSyst",pt);      }    }
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
