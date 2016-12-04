#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include <vector>
#include <assert.h>
#include <TMVA/Reader.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>
#include "TFileCollection.h"
#include "THashList.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Analysis/NAAnaFW/bin/Weights.h"
//#include "Tprime/TprimeAnalysis/interface/Weights.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Analysis/NAAnaFW/src/DMTopVariables.h"

using namespace std;

typedef vector<double> vdouble;
typedef vector<float> vfloat;
typedef vector<int> vint;
typedef vector<bool> vbool;
typedef vector<string> vstring;

enum weightedSysts { NOSYST=0, BTAGUP = 1,BTAGDOWN=2,MISTAGUP=3,MISTAGDOWN=4, PUUP=5, PUDOWN=6, LEPUP=7, LEPDOWN=8, MAXSYSTS=9};

struct systWeights{
  void initHistogramsSysts(TH1F** histo, TString name, TString, int, float, float);
  void createFilesSysts(TFile ** allFiles, TString basename, TString opt="RECREATE");
  void fillHistogramsSysts(TH1F** histo, float v, float W, double *wcats= NULL,bool verbose=false);
  void fillHistogramsSysts(TH1F** histo, float v, float W,  float *systWeights, int nFirstSysts=(int)MAXSYSTS, double *wcats=NULL, bool verbose=false);

  void closeFilesSysts(TFile ** allFiles);
  void writeHistogramsSysts(TH1F** histo, TFile ** allFiles );
  void writeSingleHistogramSysts(TH1F* histo,TFile ** allMyFiles);
  void setMax(int max);
  void setMaxNonPDF(int max);
  //  TFile** initFilesSysts();
  void setSystValue(string name, double value, bool mult=false);
  void setSystValue(int systPlace, double value, bool mult=false);
  void setPDFWeights(float * wpdfs, int numPDFs, float wzero=1.0,bool mult=true);
  void setQ2Weights(float q2up, float q2down, float wzero=1.0,bool mult=true);
  void setTWeight(float tweight, float wtotsample=1.0,bool mult=true);
  void setVHFWeight(int vhf,bool mult=true, double shiftval=0.65);
  void setPDFValue(int numPDF, double value);
  double getPDFValue(int numPDF);
  void setWeight(string name, double value, bool mult=false);
  void setWeight(int systPlace, double value, bool mult=false);
  void prepareDefault(bool addDefault, bool addPDF, bool addQ2, bool addTopPt, bool addVHF, bool addTTSplit, int numPDF=102);
  void addSyst(string name);
  void addSystNonPDF(string name);
  void setWCats(double *wcats);

  void addkFact(string name);
  void setkFact(string name,float kfact_nom, float kfact_up,float kfact_down,  bool mult=true);

  void copySysts(systWeights sys);
  void calcPDFHisto(TH1F** histos, TH1F* singleHisto, double scalefactor=1.0, int c = 0);
  void setOnlyNominal(bool useOnlyNominal=false);
  bool onlyNominal;
  bool addPDF, addQ2, addTopPt, addVHF, addTTSplit;
  int maxSysts, maxSystsNonPDF;
  int nPDF;
  int nCategories;
  float weightedSysts[150];
  double wCats[10];
  string weightedNames[150];
  string categoriesNames[10];
  };

struct lept{
  TLorentzVector p4;
  int flavour;
  };

struct by_pt{
    bool operator()(lept const &l1, lept const &l2){
    return l1.p4.Pt() > l2.p4.Pt();
    }
  };

bool pt_ordered (double i,double j){ 
  return (j<i);
  }

float calculate_mtw(float met[0], float metphi[0], vector<TLorentzVector> & mu) {
  TVector2 met_(met[0]*cos(metphi[0]),met[0]*sin(metphi[0]));
  float phi_lmet = fabs(deltaPhi(mu.at(0).Phi(),metphi[0])); 
  float mtw=sqrt(2*mu.at(0).Pt()*met[0]*(1-cos(phi_lmet)));
  return mtw;
  }

void callme(){
  std::cout<<" NaN value"<<std::endl;
  }

int main(int argc, char **argv) {

  std::cout<<"Info: Starting to run ..."<<endl;
  std::cout<<"--------------------"<<endl; 
  std::cout<<"Usage:\t "<<argv[0]<<" "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<" "<<argv[6]<<" "<<argv[7]<<endl;
  

  string sample(argv[1]) ;
  std::cout<<"\t Sample: "<<sample<<endl;

  string path(argv[2]);
  std::cout<<"\t File list to open: "<<path<<endl;

  string channel(argv[3]);
  std::cout<<"\t Channel: "<<channel<<endl;
  
  string cat(argv[4]);
  std::cout<<"\t Category: " <<cat<<endl;
  
  string sys(argv[5]);
  std::cout<<"\t Systematics: "<<sys<<endl;

  string sync(argv[6]);
  std::cout<<"\t Synchro: "<<sync<<endl;

  string isData(argv[7]);
  std::cout<<"\t isData: "<<isData<<endl;

  TString path_ = path ; 
  std::cout<<"\t File to open: "<<path_<<endl;
  std::cout<<"--------------------"<<endl; 
  
  std::cout << "Info: Loading file collection from " << path << std::endl;
  TFileCollection fc(sample.c_str(),sample.c_str(),path.c_str());
  std::cout << "Info: Files found : " << fc.GetNFiles() << std::endl;

  bool useCSVMSelection = true;
  bool useCSVTSelection = false;
  //  cout << "cat is "<<endl;
  if(cat=="CSVT"){
     useCSVMSelection = false;
     useCSVTSelection = true;
  }
  string reportName = "SelectedEvents_"+channel+"_"+cat+"_"+sample+".txt";
  //  string reportName = "SelectedEvents_"+channel+""cat"_"+sample+".txt";
  ofstream fileout;
  bool doSynch=sync.find("sync")!=std::string::npos;
  if(doSynch){
    fileout.open(reportName.c_str(),ios::in | ios::out | ios::trunc);
    fileout<<"RunNumber EvtNumber Lumi "<<std::endl;
  }
  TString outfile = "test/"+sample + "_"+channel+".root";
  
  TString treePath = "DMTreesDumper/ttDM__noSyst";
  TString treePathNEvents = "DMTreesDumper/WeightHistory";

  TChain chain(treePath);
  chain.AddFileInfoList(fc.GetList());
  TChain chainNEvents(treePathNEvents);
  chainNEvents.AddFileInfoList(fc.GetList());
  Int_t nEventsTot = (Int_t)chainNEvents.GetEntries();

  TH1D totalWeightTop("w_top_total","Top pt reweighting: overall sample weight",2000,0,2.0);
  chainNEvents.Project("w_top_total","Event_T_Weight","Event_T_Weight!=1.00");
  double topWeight=totalWeightTop.GetMean();
  cout << "totaltopweight is "<< topWeight<<endl;
  if(topWeight==0)topWeight=1;
  


  Int_t nEvents = (Int_t)chain.GetEntries();
  std::cout<<"Info: Number of Events: "<<nEvents<< endl;
  //nEvents = std::min(nEvents, 1000);
  
  TString weightedSystsNames (weightedSysts sy);
  
  systWeights systZero,syst0BM,syst1BM,syst2BM; 
  int maxSysts=0; 
  int sizeMax=50;
  int muSize, elSize, jetSize,jet20Size, muLooseSize,elLooseSize;
  //float passTrigHT(0.), Ht(0.);
  float Ht(0.), mt(0.);
  float runNumber(0.), lumiSec(0.);
  double evtNumber(0.);
  float jetE[sizeMax], jetPt[sizeMax], jetPhi[sizeMax], jetEta[sizeMax];
  float jet20E[sizeMax], jet20Pt[sizeMax], jet20Phi[sizeMax], jet20Eta[sizeMax];
  float jetIsCSVL[sizeMax], jetIsCSVM[sizeMax], jetIsCSVT[sizeMax],jetIsLoose[sizeMax],jetIsTight[sizeMax],jetak4chs_csvv2[sizeMax];
  float jet20IsCSVL[sizeMax], jet20IsCSVM[sizeMax], jet20IsCSVT[sizeMax],jet20IsLoose[sizeMax],jet20IsTight[sizeMax],jet20ak4chs_csvv2[sizeMax];
  float nTightMuons, nTightElectrons, nVetoElectrons, nLooseMuons, nJets,nCSVJets;//, nCSVLJets;
  float nGoodPV, nPV, numTrueInt, w_pu;
  float metPt[1],metPhi[1],metPx[1],metPy[1];
  
  //double met,metpx,metpy;
  float met,metpx,metpy;
  float muE[sizeMax], muPt[sizeMax], muEta[sizeMax], muIso[sizeMax], muIsTight[sizeMax],muIsLoose[sizeMax], muPhi[sizeMax],elPt[sizeMax];
  float muAntiIsoE[sizeMax], muAntiIsoPt[sizeMax], muAntiIsoEta[sizeMax], muAntiIsoIso[sizeMax], muAntiIsoIsTight[sizeMax],muAntiIsoIsLoose[sizeMax], muAntiIsoPhi[sizeMax];
  int   muAntiIsoSize;
  float muCharge[sizeMax], elCharge[sizeMax];
  float muAntiIsoCharge[sizeMax];

  float elEta[sizeMax], elIso[sizeMax], elIsTight[sizeMax], elIsVeto[sizeMax], elE[sizeMax], elPhi[sizeMax], elPassesDRmu[sizeMax];    
  float lep1Pt(0), lep1Eta(0), lep1Phi(0), lep1E(0), lep2Pt(0), lep2Eta(0), lep2Phi(0),  lep2E(0), lep1Flavour(0), lep1Charge(0), lep2Flavour(0), lep2Charge(0);
  int nPDF=102;

  float slTrigIsoMu20_v1(0.), slTrigIsoMu20_v2(0.), slTrigIsoMu20_v3(0.);
  float slTrigIsoMu22_v1(0.), slTrigIsoMu22_v2(0.), slTrigIsoMu22_v3(0.);
  float slTrigIsoMu24_v1(0.), slTrigIsoMu24_v2(0.), slTrigIsoMu24_v3(0.);

  float slTrigIsoTkMu20_v1(0.), slTrigIsoTkMu20_v2(0.), slTrigIsoTkMu20_v3(0.);
  float slTrigIsoTkMu22_v1(0.), slTrigIsoTkMu22_v2(0.), slTrigIsoTkMu22_v3(0.);
  float slTrigIsoTkMu24_v1(0.), slTrigIsoTkMu24_v2(0.), slTrigIsoTkMu24_v3(0.);

  float slTrigEle_v1(0.), slTrigEle_v2(0.);
  
  bool  TrigIsoMu20 = false;
  bool  TrigIsoMu22 = false;
  bool  TrigIsoMu24 = false;

  
  float n_trig(0.),n_lepton(0.),n_loose_veto(0.),n_lepton_cross_veto(0),n_2j(0), n_2j1t(0), n_2j1tmtw(0.);
  int nev_trig(0),nev_lepton(0),nev_loose_veto(0),nev_lepton_cross_veto(0),nev_2j(0),nev_2j1t(0), nev_2j1tmtw(0.);

  float LHEWeightSign[1] = {1.};
  float w(1.),w_q2up(1.),w_q2down(1.),w_zero(1.),w_top(1.);
  float w_pdfs[nPDF];
  float passElTrig(0.), passMuTrig(0.);
  float passHadTrig(0.) ;

  float bWeight0CSVM(1.), bWeight0CSVMBTagUp(1.),bWeight0CSVMBTagDown(1.),bWeight0CSVMMisTagUp(1.),bWeight0CSVMMisTagDown(1.);
  float bWeight1CSVM(1.), bWeight1CSVMBTagUp(1.),bWeight1CSVMBTagDown(1.),bWeight1CSVMMisTagUp(1.),bWeight1CSVMMisTagDown(1.);
  float bWeight2CSVM(1.), bWeight2CSVMBTagUp(1.),bWeight2CSVMBTagDown(1.),bWeight2CSVMMisTagUp(1.),bWeight2CSVMMisTagDown(1.);

  float lepWeight1Mu(1.), lepWeight1MuLepUp(1.), lepWeight1MuLepDown(1.);
  float lepWeight1Ele(1.), lepWeight1EleLepUp(1.), lepWeight1EleLepDown(1.);
  
  //  float bWeight0CSVT, bWeight0CSVTBTagUp,bWeight0CSVTBTagDown,bWeight0CSVTMisTagUp,bWeight0CSVTMisTagDown;
  //  float bWeight1CSVT, bWeight1CSVTBTagUp,bWeight1CSVTBTagDown,bWeight1CSVTMisTagUp,bWeight1CSVTMisTagDown;
  //  float bWeight2CSVT, bWeight2CSVTBTagUp,bWeight2CSVTBTagDown,bWeight2CSVTMisTagUp,bWeight2CSVTMisTagDown;
  
  bool onlyNominal=false;
  systZero.setOnlyNominal(onlyNominal);

  bool addPDF=false,addQ2=false,addTopPt=false,addVHF=false,addTTSplit=false;
  if(sample=="TT"){
    addTTSplit=false;
    }


  systZero.prepareDefault(true, addQ2, addPDF, addTopPt,addVHF,addTTSplit);
  //Here must add all systematics we want to put there so that the size of the syst vector is set

  //this maxSysts is to set the number of systematics we do want IMPORTANT FOR ALL INITIALIZATIONS
  maxSysts= systZero.maxSysts;
  TFile * allMyFiles[maxSysts];
  
  systZero.createFilesSysts(allMyFiles,"./res/"+sample + "_" +channel);
  //addWZNLO=false, 
  //addPDF=true;
  //addQ2=true;
  //addTopPt=true;
  //addVHF=false;
  //addWZNLO=false;
  //if(isData=="DATA"){addPDF=false, addQ2=false;addVHF=false;addTTSplit=false;} 

  if(isData=="MC"){chain.SetBranchAddress("Event_LHEWeight0", &w_zero); 
  }
  
  if(addQ2){
    chain.SetBranchAddress("Event_LHEWeight4", &w_q2up);
    chain.SetBranchAddress("Event_LHEWeight8", &w_q2down);
  }
  if(addPDF){
    for (int p = 1; p <= nPDF;++p){
      stringstream pdfss;
      pdfss<<(p+8);
      string pstr =(pdfss.str());
      chain.SetBranchAddress(("Event_LHEWeight"+pstr).c_str(), &w_pdfs[p-1]);
    //    cout << "Event_LHEWeight"+pstr<<w_pdfs<<endl;
    }
    chain.SetBranchAddress("Event_T_Weight",&w_top);
    //    chain.SetBranchAddress("Event_T_Weight",&w_top);
  }
  
  chain.SetBranchAddress("electronsTight_E", elE);
  chain.SetBranchAddress("electronsTight_Phi", elPhi);
  chain.SetBranchAddress("electronsTight_Eta", elEta);
  chain.SetBranchAddress("electronsTight_Pt", elPt);
  chain.SetBranchAddress("electronsTight_Charge", elCharge);
  chain.SetBranchAddress("electronsTight_Iso03", elIso);
  chain.SetBranchAddress("electronsTight_vidTight", elIsTight);
  chain.SetBranchAddress("electronsTight_vidVeto", elIsVeto);
  chain.SetBranchAddress("electronsTight_PassesDRmu", elPassesDRmu);
  chain.SetBranchAddress("electronsTight_size", &elSize);
  chain.SetBranchAddress("electronsVeto_size", &elLooseSize);

  chain.SetBranchAddress("muonsTight_E", muE);
  chain.SetBranchAddress("muonsTight_Phi", muPhi);
  chain.SetBranchAddress("muonsTight_Eta", muEta);
  chain.SetBranchAddress("muonsTight_Pt", muPt);
  chain.SetBranchAddress("muonsTight_Charge", muCharge);
  chain.SetBranchAddress("muonsTight_Iso04", muIso);
  chain.SetBranchAddress("muonsTight_IsTightMuon", muIsTight);
  chain.SetBranchAddress("muonsTight_IsLooseMuon", muIsLoose);
  chain.SetBranchAddress("muonsTight_size", &muSize);
  chain.SetBranchAddress("muonsLoose_size", &muLooseSize);

  chain.SetBranchAddress("muonsTightAntiIso_E", muAntiIsoE);
  chain.SetBranchAddress("muonsTightAntiIso_Phi", muAntiIsoPhi);
  chain.SetBranchAddress("muonsTightAntiIso_Eta", muAntiIsoEta);
  chain.SetBranchAddress("muonsTightAntiIso_Pt", muAntiIsoPt);
  chain.SetBranchAddress("muonsTightAntiIso_Charge", muAntiIsoCharge);
  chain.SetBranchAddress("muonsTightAntiIso_Iso04", muAntiIsoIso);
  chain.SetBranchAddress("muonsTightAntiIso_IsTightMuon", muAntiIsoIsTight);
  chain.SetBranchAddress("muonsTightAntiIso_IsLooseMuon", muAntiIsoIsLoose);
  chain.SetBranchAddress("muonsTightAntiIso_size", &muAntiIsoSize);


  chain.SetBranchAddress("jetsAK4CHSTight_CorrE",      &jetE);
  chain.SetBranchAddress("jetsAK4CHSTight_CorrPt",     &jetPt);
  chain.SetBranchAddress("jetsAK4CHSTight_Phi",    &jetPhi);
  chain.SetBranchAddress("jetsAK4CHSTight_Eta",    &jetEta);
  chain.SetBranchAddress("jetsAK4CHSTight_size",   &jetSize);
  chain.SetBranchAddress("jetsAK4CHSTight_IsCSVL", &jetIsCSVL);
  chain.SetBranchAddress("jetsAK4CHSTight_IsCSVM", &jetIsCSVM);
  chain.SetBranchAddress("jetsAK4CHSTight_IsCSVT", &jetIsCSVT);
  chain.SetBranchAddress("jetsAK4CHSTight_IsLoose",&jetIsLoose);
  chain.SetBranchAddress("jetsAK4CHSTight_IsTight",&jetIsTight);
  chain.SetBranchAddress("jetsAK4CHSTight_CSVv2",  &jetak4chs_csvv2);

  chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_CorrE",      &jet20E);
  chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_CorrPt",     &jet20Pt);
  chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_Phi",    &jet20Phi);
  chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_Eta",    &jet20Eta);
  chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_size",   &jet20Size);
  chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_IsCSVL", &jet20IsCSVL);
  chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_IsCSVM", &jet20IsCSVM);
  chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_IsCSVT", &jet20IsCSVT);
  chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_IsLoose",&jet20IsLoose);
  chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_IsTight",&jet20IsTight);
  chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_CSVv2",  &jet20ak4chs_csvv2);
  
  chain.SetBranchAddress("Event_RunNumber", &runNumber);
  chain.SetBranchAddress("Event_LumiBlock", &lumiSec);
  chain.SetBranchAddress("Event_EventNumber", &evtNumber);
  
  chain.SetBranchAddress("Event_Ht", &Ht);
  chain.SetBranchAddress("Event_mt", &mt);
  chain.SetBranchAddress("Event_nGoodPV",&nGoodPV);
  chain.SetBranchAddress("Event_nPV",&nPV);
  chain.SetBranchAddress("Event_nTruePV",&numTrueInt);
 
  chain.SetBranchAddress("Event_Lepton1_Pt", &lep1Pt);
  chain.SetBranchAddress("Event_Lepton1_Eta", &lep1Eta);
  chain.SetBranchAddress("Event_Lepton1_Phi", &lep1Phi);
  chain.SetBranchAddress("Event_Lepton1_E", &lep1E);
  chain.SetBranchAddress("Event_Lepton1_Charge", &lep1Charge);
  chain.SetBranchAddress("Event_Lepton1_Flavour", &lep1Flavour);

  chain.SetBranchAddress("Event_Lepton2_Pt", &lep2Pt);
  chain.SetBranchAddress("Event_Lepton2_Eta", &lep2Eta);
  chain.SetBranchAddress("Event_Lepton2_Phi", &lep2Phi);
  chain.SetBranchAddress("Event_Lepton2_E", &lep2E);
  chain.SetBranchAddress("Event_Lepton2_Charge", &lep2Charge);
  chain.SetBranchAddress("Event_Lepton2_Flavour", &lep2Flavour);

  chain.SetBranchAddress("Event_nTightElectrons",&nTightElectrons);
  chain.SetBranchAddress("Event_nVetoElectrons",&nVetoElectrons);
  chain.SetBranchAddress("Event_nLooseMuons",&nLooseMuons);
  chain.SetBranchAddress("Event_nTightMuons",&nTightMuons);


  //  float bWeight0CSVM, bWeight0CSVMBTagUp,bWeight0CSVMBTagDown,bWeight0CSVMMisTagUp,bWeight0CSVMMisTagDown;
  //  float bWeight1CSVM, bWeight1CSVMBTagUp,bWeight1CSVMBTagDown,bWeight1CSVMMisTagUp,bWeight1CSVMMisTagDown;
  //  float bWeight2CSVM, bWeight2CSVMBTagUp,bWeight2CSVMBTagDown,bWeight2CSVMMisTagUp,bWeight2CSVMMisTagDown;
  //b-weight settings
  if(isData!="DATA"){
    
    if(useCSVMSelection){
      chain.SetBranchAddress("Event_bWeight2CSVMTight", &bWeight2CSVM);
      chain.SetBranchAddress("Event_bWeightMisTagUp2CSVMTight", &bWeight2CSVMMisTagUp);
      chain.SetBranchAddress("Event_bWeightMisTagDown2CSVMTight", &bWeight2CSVMMisTagDown);
      chain.SetBranchAddress("Event_bWeightBTagUp2CSVMTight", &bWeight2CSVMBTagUp);
      chain.SetBranchAddress("Event_bWeightBTagDown2CSVMTight", &bWeight2CSVMBTagDown);
      
      chain.SetBranchAddress("Event_bWeight0CSVMTight", &bWeight0CSVM);
      
      chain.SetBranchAddress("Event_bWeightMisTagUp0CSVMTight", &bWeight0CSVMMisTagUp);
      chain.SetBranchAddress("Event_bWeightMisTagDown0CSVMTight", &bWeight0CSVMMisTagDown);
      chain.SetBranchAddress("Event_bWeightBTagUp0CSVMTight", &bWeight0CSVMBTagUp);
      chain.SetBranchAddress("Event_bWeightBTagDown0CSVMTight", &bWeight0CSVMBTagDown);
      
      chain.SetBranchAddress("Event_bWeight1CSVMTight", &bWeight1CSVM);
      chain.SetBranchAddress("Event_bWeightMisTagUp1CSVMTight", &bWeight1CSVMMisTagUp);
      chain.SetBranchAddress("Event_bWeightMisTagDown1CSVMTight", &bWeight1CSVMMisTagDown);
      chain.SetBranchAddress("Event_bWeightBTagUp1CSVMTight", &bWeight1CSVMBTagUp);
      chain.SetBranchAddress("Event_bWeightBTagDown1CSVMTight", &bWeight1CSVMBTagDown);
    }
    if(useCSVTSelection){
      chain.SetBranchAddress("Event_bWeight2CSVTTight", &bWeight2CSVM);
      chain.SetBranchAddress("Event_bWeightMisTagUp2CSVTTight", &bWeight2CSVMMisTagUp);
      chain.SetBranchAddress("Event_bWeightMisTagDown2CSVTTight", &bWeight2CSVMMisTagDown);
      chain.SetBranchAddress("Event_bWeightBTagUp2CSVTTight", &bWeight2CSVMBTagUp);
      chain.SetBranchAddress("Event_bWeightBTagDown2CSVTTight", &bWeight2CSVMBTagDown);
      
      chain.SetBranchAddress("Event_bWeight0CSVTTight", &bWeight0CSVM);
      
      chain.SetBranchAddress("Event_bWeightMisTagUp0CSVTTight", &bWeight0CSVMMisTagUp);
      chain.SetBranchAddress("Event_bWeightMisTagDown0CSVTTight", &bWeight0CSVMMisTagDown);
      chain.SetBranchAddress("Event_bWeightBTagUp0CSVTTight", &bWeight0CSVMBTagUp);
      chain.SetBranchAddress("Event_bWeightBTagDown0CSVTTight", &bWeight0CSVMBTagDown);
      
      chain.SetBranchAddress("Event_bWeight1CSVTTight", &bWeight1CSVM);
      chain.SetBranchAddress("Event_bWeightMisTagUp1CSVTTight", &bWeight1CSVMMisTagUp);
      chain.SetBranchAddress("Event_bWeightMisTagDown1CSVTTight", &bWeight1CSVMMisTagDown);
      chain.SetBranchAddress("Event_bWeightBTagUp1CSVTTight", &bWeight1CSVMBTagUp);
      chain.SetBranchAddress("Event_bWeightBTagDown1CSVTTight", &bWeight1CSVMBTagDown);
    }
  }

  if(isData=="MC"){
    chain.SetBranchAddress("metFull_CorrPt",metPt);
    chain.SetBranchAddress("metFull_CorrPhi",metPhi);}
  else {
    chain.SetBranchAddress("metFull_CorrPt",metPt);
    chain.SetBranchAddress("metFull_CorrPhi",metPhi);}

  chain.SetBranchAddress("metFull_Px",metPx);
  chain.SetBranchAddress("metFull_Py",metPy);
  
  //Muon trigger
  if(isData=="MC" && false){
    chain.SetBranchAddress("Event_passesHLT_IsoMu22_v1", &slTrigIsoMu22_v1);
    chain.SetBranchAddress("Event_passesHLT_IsoMu22_v2", &slTrigIsoMu22_v2);
    chain.SetBranchAddress("Event_passesHLT_IsoMu22_v3", &slTrigIsoMu22_v3);
    }
  else {
    //string ptmu="20";
    //if(runNumber>=272023 && runNumber<=274443) ptmu="20";
    //if(runNumber>=274954 && runNumber<=276811) ptmu="22";
    //if(runNumber>=276824 && runNumber<=99999999) ptmu="24";
    //chain.SetBranchAddress(("Event_passesHLT_IsoMu"+ptmu+"_v1").c_str(), &slTrigIsoMu22_v1);
    //chain.SetBranchAddress(("Event_passesHLT_IsoMu"+ptmu+"_v2").c_str(), &slTrigIsoMu22_v2);
    //chain.SetBranchAddress(("Event_passesHLT_IsoMu"+ptmu+"_v3").c_str(), &slTrigIsoMu22_v3);

    //Event_passesHLT_IsoTkMu22_v3

    chain.SetBranchAddress("Event_passesHLT_IsoMu20_v1", &slTrigIsoMu20_v1);
    chain.SetBranchAddress("Event_passesHLT_IsoMu20_v2", &slTrigIsoMu20_v2);
    chain.SetBranchAddress("Event_passesHLT_IsoMu20_v3", &slTrigIsoMu20_v3);
    
    chain.SetBranchAddress("Event_passesHLT_IsoMu22_v1", &slTrigIsoMu22_v1);
    chain.SetBranchAddress("Event_passesHLT_IsoMu22_v2", &slTrigIsoMu22_v2);
    chain.SetBranchAddress("Event_passesHLT_IsoMu22_v3", &slTrigIsoMu22_v3);
    
    chain.SetBranchAddress("Event_passesHLT_IsoMu24_v1", &slTrigIsoMu24_v1);
    chain.SetBranchAddress("Event_passesHLT_IsoMu24_v2", &slTrigIsoMu24_v2);
    chain.SetBranchAddress("Event_passesHLT_IsoMu24_v3", &slTrigIsoMu24_v3);
    //
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu20_v1", &slTrigIsoTkMu20_v1);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu20_v2", &slTrigIsoTkMu20_v2);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu20_v3", &slTrigIsoTkMu20_v3);
    
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu22_v1", &slTrigIsoTkMu22_v1);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu22_v2", &slTrigIsoTkMu22_v2);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu22_v3", &slTrigIsoTkMu22_v3);
    
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu24_v1", &slTrigIsoTkMu24_v1);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu24_v2", &slTrigIsoTkMu24_v2);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu24_v3", &slTrigIsoTkMu24_v3);
    }

  
  //Electron Triger
  if(isData=="MC"){
    chain.SetBranchAddress("Event_passesHLT_Ele27_eta2p1_WP75_Gsf_v1", &slTrigEle_v1);
    chain.SetBranchAddress("Event_passesHLT_Ele27_eta2p1_WP75_Gsf_v2", &slTrigEle_v2);
    }
  
  else{
    chain.SetBranchAddress("Event_passesHLT_Ele27_eta2p1_WPLoose_Gsf_v1", &slTrigEle_v1);
    chain.SetBranchAddress("Event_passesHLT_Ele27_eta2p1_WPLoose_Gsf_v2", &slTrigEle_v2);
    }

  chain.SetBranchAddress("Event_passesSingleElTriggers", &passElTrig);
  chain.SetBranchAddress("Event_passesSingleMuTriggers", &passMuTrig);
  chain.SetBranchAddress("Event_passesHadronicTriggers", &passHadTrig);

  
  
  /********************************************************************************/
  /**************                    Histogram booking              ***************/
  /********************************************************************************/
  TH1F *h_cutFlow = new TH1F("h_cutFlow","cutflow",10,-0.5,9.5);
  TH1F *h_weight_sign = new TH1F("h_weight_sign","weight correction factor before all selections, gives the effective number of events",2,-2.0,2.0);
  chainNEvents.Project("h_weight_sign","Event_LHEWeight10/abs(Event_LHEWeight10)"); // I messed up and weight zero is not  saved in the full chain... but this one should do the trick as well! Just do h_weight_sign->GetMean() and multiply this to get the correct number :)
  double zeroWeight=h_weight_sign->GetMean(); // Will have to add a zeroWeight vector in the systematics for a proper treatment of the acceptance in systematics. For now we do it via hardcoding outside.
  cout << " mean  "<< zeroWeight<<endl;
    
  
  TH1F *h_weightZero[maxSysts]; systZero.initHistogramsSysts(h_weightZero,"h_weightZero","weight before all selections, normalized by the total number of events",2000,0,10.0);
  TH1F *h_nPV[maxSysts];        systZero.initHistogramsSysts(h_nPV,"h_nPV","nPV",102,0,51);
  TH1F *h_nGoodPV[maxSysts];    systZero.initHistogramsSysts(h_nGoodPV,"h_nGoodPV","nGoodPV",92,0,46);
  TH1F *h_nTruePV[maxSysts];    systZero.initHistogramsSysts(h_nTruePV,"h_nTruePV","nTruePV",80,0,40);
  
  TH1F *h_nJets[maxSysts];      systZero.initHistogramsSysts(h_nJets,"h_nJets","Number of tight jets",13,-0.5,12.5);
  TH1F *h_nbJets[maxSysts];     systZero.initHistogramsSysts(h_nbJets,"h_nbJets","Number of tight b-jets",11,-0.5,10.5);

  //2j0t 
  TH1F *h_2j0t_jetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j0t_jetpt40_1st,   "h_2j0t_jetpt40_leading","2j0t leading jet Pt ",250,0,500);
  TH1F *h_2j0t_jetpt40_2nd[maxSysts];    systZero.initHistogramsSysts(h_2j0t_jetpt40_2nd,   "h_2j0t_jetpt40_subleading","2j0t Sub. leading jet Pt ",250,0,500);
  TH1F *h_2j0t_mtw[maxSysts];            systZero.initHistogramsSysts(h_2j0t_mtw,           "h_2j0t_mtw",     "2j0t mtw ",250,0,500);
  
  //After mtw>50 cut
  TH1F *h_2j0t_mtwcut_jetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j0t_mtwcut_jetpt40_1st,   "h_2j0t_mtwcut_jetpt40_leading","2j0t leading jet Pt ",250,0,500);
  TH1F *h_2j0t_mtwcut_jetpt40_2nd[maxSysts];    systZero.initHistogramsSysts(h_2j0t_mtwcut_jetpt40_2nd,   "h_2j0t_mtwcut_jetpt40_subleading","2j0t Sub. leading jet Pt ",250,0,500);

  TH1F *h_2j0t_mtwcut_jeteta40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j0t_mtwcut_jeteta40_1st,   "h_2j0t_mtwcut_jeteta40_leading","2j0t leading jet eta ",200,-4.7,4.7);
  TH1F *h_2j0t_mtwcut_incljeteta[maxSysts];    systZero.initHistogramsSysts(h_2j0t_mtwcut_incljeteta,   "h_2j0t_mtwcut_incljeteta","2j0t leading inclusive jet eta ",200,-4.7,4.7);

  TH1F *h_2j0t_mtwcut_jetabseta40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j0t_mtwcut_jetabseta40_1st,   "h_2j0t_mtwcut_jetabseta40_leading","2j0t leading jet eta ",200,0,4.7);
  TH1F *h_2j0t_mtwcut_incljetabseta[maxSysts];    systZero.initHistogramsSysts(h_2j0t_mtwcut_incljetabseta,   "h_2j0t_mtwcut_incljetabseta","2j0t leading inclusive jet eta ",200,0,4.7);

  
  //2j1t Category
  TH1F *h_2j1t_bjetpt[maxSysts];    systZero.initHistogramsSysts(h_2j1t_bjetpt,  "h_2j1t_bjetpt",  "2j1t b jet pt",250,0,500);
  TH1F *h_2j1t_ljetpt[maxSysts];     systZero.initHistogramsSysts(h_2j1t_ljetpt,   "h_2j1t_ljetpt",   "2j1t light jet pt" ,250,0,500);
  TH1F *h_2j1t_ljeteta[maxSysts];   systZero.initHistogramsSysts(h_2j1t_ljeteta,  "h_2j1t_ljeteta",  "2j1t light jet eta",200,-4.7,4.7);
  TH1F *h_2j1t_nMu[maxSysts];       systZero.initHistogramsSysts(h_2j1t_nMu,     "h_2j1t_nMu",     "2j1t Number of tight Muons",13,-0.5,12.5);
  TH1F *h_2j1t_MuPt[maxSysts];      systZero.initHistogramsSysts(h_2j1t_MuPt,    "h_2j1t_MuPt",    "2j1t muon pt ",250,0,500);
  TH1F *h_2j1t_MuEta[maxSysts];     systZero.initHistogramsSysts(h_2j1t_MuEta,   "h_2j1t_MuEta",   "2j1t muon eta ",100,-2.1,2.1);
  TH1F *h_2j1t_MuPhi[maxSysts];     systZero.initHistogramsSysts(h_2j1t_MuPhi,   "h_2j1t_MuPhi",   "2j1t muon phi ",100, -3.2, 3.2);
  TH1F *h_2j1t_MuE[maxSysts];       systZero.initHistogramsSysts(h_2j1t_MuE,     "h_2j1t_MuE",     "2j1t muon e ",250,0,500);
  TH1F *h_2j1t_MuCharge[maxSysts];  systZero.initHistogramsSysts(h_2j1t_MuCharge,"h_2j1t_MuCharge","2j1t muon charge ",2,-1,1);
  TH1F *h_2j1t_mtw[maxSysts];       systZero.initHistogramsSysts(h_2j1t_mtw,     "h_2j1t_mtw",     "2j1t mtw ",250,0,500);
  TH1F *h_2j1t_topMass[maxSysts];   systZero.initHistogramsSysts(h_2j1t_topMass, "h_2j1t_topMass", "2j1t top mass",200,100,500);
  TH1F *h_2j1t_topPt[maxSysts];     systZero.initHistogramsSysts(h_2j1t_topPt,   "h_2j1t_topPt",   "2j1t top pt",250,0,500);
  
  TH1F *h_2j1t_etajprime[maxSysts];   systZero.initHistogramsSysts(h_2j1t_etajprime,  "h_2j1t_etajprime",  "2j1t light jet abs(eta)",100,0.,4.7);

  //Jets pt 20 spectra.
  TH1F *h_2j1t_nextrajets[maxSysts];   systZero.initHistogramsSysts(h_2j1t_nextrajets,  "h_2j1t_nextrajets",  "2j1t number of jets with 20<pt<40",5,-0.5,4.5);
  TH1F *h_2j1t_leadingextrajetpt[maxSysts];   systZero.initHistogramsSysts(h_2j1t_leadingextrajetpt,  "h_2j1t_leadingextrajetpt",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_leadingextrajetcsv[maxSysts];   systZero.initHistogramsSysts(h_2j1t_leadingextrajetcsv,  "h_2j1t_leadingextrajetcsv",  "2j1t light jet abs(eta)",100,0.,1.0);


  
  //After mtw>50 cut, label mtwcut
  TH1F *h_2j1t_mtwcut_topMass[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_topMass, "h_2j1t_mtwcut_topMass", "2j1t top mass",200,100,500);
  TH1F *h_2j1t_mtwcut_mtw[maxSysts];systZero.initHistogramsSysts(h_2j1t_mtwcut_mtw,     "h_2j1t_mtwcut_mtw",     "2j1t mtw ",250,0,500);
  TH1F *h_2j1t_mtwcut_bjetpt[maxSysts];    systZero.initHistogramsSysts(h_2j1t_mtwcut_bjetpt,  "h_2j1t_mtwcut_bjetpt",  "2j1t b jet pt",250,0,500);
  TH1F *h_2j1t_mtwcut_ljetpt[maxSysts];     systZero.initHistogramsSysts(h_2j1t_mtwcut_ljetpt,   "h_2j1t_mtwcut_ljetpt",   "2j1t light jet pt" ,250,0,500);
  TH1F *h_2j1t_mtwcut_ljeteta[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_ljeteta,  "h_2j1t_mtwcut_ljeteta",  "2j1t light jet eta",200,-4.7,4.7);
  TH1F *h_2j1t_mtwcut_topPt[maxSysts];     systZero.initHistogramsSysts(h_2j1t_mtwcut_topPt,   "h_2j1t_mtwcut_topPt",   "2j1t top pt",250,0,500);
  TH1F *h_2j1t_mtwcut_etajprime[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_etajprime,  "h_2j1t_mtwcut_etajprime",  "2j1t light jet abs(eta)",100,0.,4.7);


  TH1F *h_2j1t_mtwcut_nextrajets[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_nextrajets,  "h_2j1t_mtwcut_nextrajets",  "2j1t number of jets with 20<pt<40",5,-0.5,4.5);
  TH1F *h_2j1t_mtwcut_leadingextrajetpt[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt,  "h_2j1t_mtwcut_leadingextrajetpt",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_mtwcut_leadingextrajetcsv[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv,  "h_2j1t_mtwcut_leadingextrajetcsv",  "2j1t csv of the leading extra jet",100,0.,1.0);
  TH1F *h_2j1t_mtwcut_topMassExtra[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_topMassExtra, "h_2j1t_mtwcut_topMassExtra", "2j1t top mass with the first extra jet as b",200,100,500);

  
  //mtw > 50; |eta_j'|>2.5, label mtwcut_sr
  TH1F *h_2j1t_mtwcut_sr_topMass[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_topMass, "h_2j1t_mtwcut_sr_topMass", "2j1t top mass",200,100,500);
  TH1F *h_2j1t_mtwcut_sr_topPt[maxSysts];     systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_topPt,   "h_2j1t_mtwcut_sr_topPt",   "2j1t top pt",250,0,500);
  
  TH1F *h_2j1t_mtwcut_sr_nextrajets[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_nextrajets,  "h_2j1t_mtwcut_sr_nextrajets",  "2j1t number of jets with 20<pt<40",5,-0.5,4.5);
  TH1F *h_2j1t_mtwcut_sr_leadingextrajetpt[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt,  "h_2j1t_mtwcut_sr_leadingextrajetpt",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_mtwcut_sr_leadingextrajetcsv[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv,  "h_2j1t_mtwcut_sr_leadingextrajetcsv",  "2j1t csv of the leading extra jet",100,0.,1.0);
  TH1F *h_2j1t_mtwcut_sr_topMassExtra[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra, "h_2j1t_mtwcut_sr_topMassExtra", "2j1t top mass with the first extra jet as b",200,100,500);
  


  //Other possible cr:
  //Features of top quarks reconstructed with pt 20.
  
  ////
  //Additional t-channel enriched regions from TOP-12-038
  //Top mass SR 130<mtop<220, cut
  ////Top mass SR 130<mtop<220, mtw > 50 cuts
  //Top mass SB (130<mtop || mtop>220)
  ////

  //3j1t 
  TH1F *h_3j1t_bjetpt[maxSysts];   systZero.initHistogramsSysts(h_3j1t_bjetpt,  "h_3j1t_bjetpt",  "3j1t b jet pt ",250,0,500);
  TH1F *h_3j1t_MuPt[maxSysts];     systZero.initHistogramsSysts(h_3j1t_MuPt,    "h_3j1t_MuPt",    "muon pt ",250,0,500);
  TH1F *h_3j1t_MuEta[maxSysts];    systZero.initHistogramsSysts(h_3j1t_MuEta,   "h_3j1t_MuEta",   "muon eta ",100,-2.1,2.1);
  TH1F *h_3j1t_MuPhi[maxSysts];    systZero.initHistogramsSysts(h_3j1t_MuPhi,   "h_3j1t_MuPhi",   "muon phi ",100, -3.2, 3.2);
  TH1F *h_3j1t_MuE[maxSysts];      systZero.initHistogramsSysts(h_3j1t_MuE,     "h_3j1t_MuE",     "muon e ",250,0,500);
  TH1F *h_3j1t_mtw[maxSysts];      systZero.initHistogramsSysts(h_3j1t_mtw,     "h_3j1t_mtw",     "3j1t mtw ",250,0,500);
  
  //3j2t
  TH1F *h_3j2t_bjetpt[maxSysts];   systZero.initHistogramsSysts(h_3j2t_bjetpt,    "h_3j2t_bjetpt",     "3j2t b jet pt ",100,0,500);
  TH1F *h_3j2t_2ndbjetpt[maxSysts];systZero.initHistogramsSysts(h_3j2t_2ndbjetpt, "h_3j2t_2ndbjetpt",  "3j2t sub lead. b jet pt",100,0,500);
  TH1F *h_3j2t_ljetpt[maxSysts];   systZero.initHistogramsSysts(h_3j2t_ljetpt,    "h_3j2t_ljetpt",     "3j2t light jet pt ",100,0,500);
  TH1F *h_3j2t_ljeteta[maxSysts];   systZero.initHistogramsSysts(h_3j2t_ljeteta,    "h_3j2t_ljeteta",     "3j2t light jet eta ",100,0,500);
  TH1F *h_3j2t_etajprime[maxSysts];   systZero.initHistogramsSysts(h_3j2t_etajprime,    "h_3j2t_etajprime",     "3j2t eta j' ",100,0.,4.7);
  TH1F *h_3j2t_nextrajets[maxSysts];   systZero.initHistogramsSysts(h_3j2t_nextrajets,    "h_3j2t_nextrajets",     "3j2t number of extra jets ",5,-0.5,4.5);
  TH1F *h_3j2t_MuPt[maxSysts];     systZero.initHistogramsSysts(h_3j2t_MuPt,      "h_3j2t_MuPt",       "3j2t muon pt ",250,0,500);
  TH1F *h_3j2t_MuEta[maxSysts];    systZero.initHistogramsSysts(h_3j2t_MuEta,     "h_3j2t_MuEta",      "3j2t muon eta ",100,-2.1,2.1);
  TH1F *h_3j2t_MuPhi[maxSysts];    systZero.initHistogramsSysts(h_3j2t_MuPhi,     "h_3j2t_MuPhi",      "3j2t muon phi ",100, -3.2, 3.2);
  TH1F *h_3j2t_MuE[maxSysts];      systZero.initHistogramsSysts(h_3j2t_MuE,       "h_3j2t_MuE",        "3j2t muon e ",250,0,500);
  TH1F *h_3j2t_MuCharge[maxSysts]; systZero.initHistogramsSysts(h_3j2t_MuCharge,  "h_3j2t_MuCharge",   "3j2t muon charge ",2,-1,1);
  TH1F *h_3j2t_mtw[maxSysts];      systZero.initHistogramsSysts(h_3j2t_mtw,       "h_3j2t_mtw",        "3j2t mtw ",250,0,500);
  
  //After mtw>50 cut



  TopUtilities topUtils;
  
  //Pileup Reweighting
  edm::LumiReWeighting LumiWeights_, LumiWeightsUp_, LumiWeightsDown_;

  TFile* file_eff_trig = TFile::Open("data/SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root");
  TFile* file_sf_id = TFile::Open("data/MuonID_Z_RunBCD_prompt80X_7p65.root");
  TFile* file_sf_iso = TFile::Open("data/MuonIso_Z_RunBCD_prompt80X_7p65.root");
  
  //cout << "before getting weights "<<endl;

  Weights muTightTrigEff(file_eff_trig,"IsoMu22_OR_IsoTkMu22_PtEtaBins_Run273158_to_274093/efficienciesDATA/abseta_pt_DATA");
  Weights muTightIdSF(file_sf_id,"MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio");
  Weights muTightIsoSF(file_sf_iso,"MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio");

  //  cout << "got weights" <<endl;

  if(isData=="MC"){
    //LumiWeights_ = edm::LumiReWeighting("pu/MyDataPileupHistogram.root", "pu/MyDataPileupHistogram.root","pileup","pileup");
    //    LumiWeightsUp_ = edm::LumiReWeighting("pu/MyDataPileupHistogram.root", "pu/MyDataPileupHistogram.root","pileup","pileup");
    //    LumiWeightsDown_ = edm::LumiReWeighting("pu/MyDataPileupHistogram.root", "pu/MyDataPileupHistogram.root","pileup","pileup");

    LumiWeights_ = edm::LumiReWeighting("pu/MCPU.root", "pu/MyDataPileupHistogram.root","pileup","pileup");
    LumiWeightsUp_ = edm::LumiReWeighting("pu/MCPU.root", "pu/MyDataPileupHistogram_up.root","pileup","pileup");
    LumiWeightsDown_ = edm::LumiReWeighting("pu/MCPU.root", "pu/MyDataPileupHistogram_down.root","pileup","pileup");


    //LumiWeightsUp_ = edm::LumiReWeighting("pu/puMC.root", "pu/MyDataPileupHistogram_up.root","MC_pu","pileup");
    //LumiWeightsDown_ = edm::LumiReWeighting("pu/puMC.root", "pu/MyDataPileupHistogram_down.root","MC_pu","pileup");
  }
  
for(Int_t evt=0; evt<nEvents; evt++ ){
    if(evt%10000==1 ){
    cout<<"Info: Running on event: "<<evt<<endl; 
    }
  chain.GetEntry(evt);
  int maxJetLoop = min(15, jetSize);
  int maxJet20Loop = min(15, jet20Size);

  //  cout << " jet 20 loop size"<< maxJet20Loop <<endl;
  int maxMuLoop = min(6, muSize);
  if(channel == "muonantiiso") maxMuLoop = min(6, muAntiIsoSize); 
  int maxElLoop = min(6, elSize);
  
  //step 1 Trigger
  if(isData=="DATA"){
    TrigIsoMu20= false;
    //TrigIsoMu20= (runNumber>=272023 && runNumber<=274443 && (slTrigIsoMu20_v1 || slTrigIsoMu20_v2 || slTrigIsoMu20_v3));
    //TrigIsoMu22= (runNumber>=274954 && runNumber<=276811 && (slTrigIsoMu22_v1 || slTrigIsoMu22_v2 || slTrigIsoMu22_v3));  
    TrigIsoMu22= (runNumber>=0.0 && runNumber<=276811 && (slTrigIsoMu22_v1 || slTrigIsoMu22_v2 || slTrigIsoMu22_v3 ||
							  slTrigIsoTkMu22_v1 || slTrigIsoTkMu22_v2 || slTrigIsoTkMu22_v3));  
    TrigIsoMu24= (runNumber>=276824 && runNumber<=999999 && (slTrigIsoMu24_v1 || slTrigIsoMu24_v2 || slTrigIsoMu24_v3  ||
							     slTrigIsoTkMu24_v1 || slTrigIsoTkMu24_v2 || slTrigIsoTkMu24_v3
							     ));  
  }
  if(isData=="MC"){
    TrigIsoMu20=false;
    TrigIsoMu24=false;
    TrigIsoMu22= true;//Trigger efficiency:use the measured one atm --> always pass the trigger in MC// (slTrigIsoMu22_v1 || slTrigIsoMu22_v2 || slTrigIsoMu22_v3);  
  }

  // if (TrigIsoMu20) cout << "TrigIsoMu20 fired in Event = "<<evt<<" Run Number = "<<runNumber<<endl;
  //if (TrigIsoMu22) cout << "TrigIsoMu22 fired in Event = "<<evt<<" Run Number = "<<runNumber<<endl;
  //if (TrigIsoMu24) cout << "TrigIsoMu24 fired in Event = "<<evt<<" Run Number = "<<runNumber<<endl;
  
  bool muonTrigger = (TrigIsoMu20 || TrigIsoMu22 || TrigIsoMu24); 
  bool passesAnyTrigger = muonTrigger || muonTrigger;//REMINDER: ADD ALL TRIGGERS HERE 
  if(!passesAnyTrigger)continue;
  
  n_trig += w;
  nev_trig +=1;

  double puUpFact=1.0;//(LumiWeightsUp_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));
  double puDownFact=1.0;//(LumiWeightsDown_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));

  if(isData=="MC"){
    LHEWeightSign[0]=w_zero/fabs(w_zero);
    w = LHEWeightSign[0];
    //    w=1.;
    w_pu = LumiWeights_.weight(numTrueInt);
    //    cout << " ntrueint "<< numTrueInt <<" puweight "<<w_pu<<endl;
    w = w * w_pu;
    //    w = w  w_pu;
    
    if(sample=="TTNoTT"){
      w*=w_top/topWeight;
    }
    
    if(w_pu!=0){
      puUpFact=(LumiWeightsUp_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));
      puDownFact=(LumiWeightsDown_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));
    }
    //    cout << " ntrueint "<< numTrueInt <<" puweight "<<w_pu<<" wup "<< puUpFact<<" wdown "<<puDownFact<<endl;
    
  } 

    //if(numTrueInt>49){
    //  cout << " --> numTrueInt very high!!" << endl;
    //  puUpFact =0;
    //  puDownFact=0;
    //  }


  met = metPt[0];
  metpx = metPx[0];
  metpy = metPy[0];
  
  TLorentzVector lep1;
  TLorentzVector lep2;
  TLorentzVector lep, mu, el;

  vector<TLorentzVector> tightEl, tightMu;
  int nMu(0.), nEl(0.);//, nVetoEl(0.), nLooseMu(0.);
  for(int e= 0; e<maxElLoop;++e ){
      if(elPt[e]>20){
      el.SetPtEtaPhiE(elPt[e], elEta[e], elPhi[e],elE[e]);
      tightEl.push_back(el);
      }
  }//end loop on electrons     
  vector<float> selectedIso;
  if(channel!="muonantiiso"){
    for(int m= 0; m<maxMuLoop;++m ){
      if(muPt[m]>20){

	      mu.SetPtEtaPhiE(muPt[m], muEta[m], muPhi[m],muE[m]);
	      tightMu.push_back(mu);
      }
    }
  }
  else {
      for(int m= 0; m<maxMuLoop;++m ){
      if(muAntiIsoPt[m]>20 && muAntiIsoIso[m]>0.25){
	    selectedIso.push_back(muAntiIsoIso[m]);
	    mu.SetPtEtaPhiE(muAntiIsoPt[m], muAntiIsoEta[m], muAntiIsoPhi[m],muAntiIsoE[m]);
	    tightMu.push_back(mu);
      }
    }
  }
  //end loop on muons
  
  nMu = tightMu.size();
  nEl = tightEl.size();
  vector<lept> leptons;
  //  lepWeight1Mu(1.), lepWeight1MuLepUp(1.), lepWeight1MuLepDown(1.);
  //  lepWeight1Ele(1.), lepWeight1EleLepUp(1.), lepWeight1EleLepDown(1.);

  //Weights muTightTrigEff(file_eff_trig,"IsoMu22_OR_IsoTkMu22_PtEtaBins_Run273158_to_274093/efficienciesDATA/abseta_pt");
  //Weights muTightIdSF(file_eff_trig,"MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio");
  //Weights muTightIsoSF(file_eff_trig,"MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio");
  if(isData=="MC"){
    if(tightMu.size()==1){
      /*      cout << "event with 1 muon, "<< evt<< " getting trigg eff" <<endl;
      cout << muTightTrigEff.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())<<endl;
      cout << "getting id sf"<< endl;
      cout << muTightIdSF.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())<<endl;
      cout << "getting iso sf"<< endl;
      cout << muTightIsoSF.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())<<endl;
      cout << "lep sf retrieved successfully!"<<endl;*/
      lepWeight1Mu = muTightTrigEff.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())*
	muTightIdSF.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())*
	muTightIsoSF.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt());
      float errMu= muTightTrigEff.getErr(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())*
	muTightIdSF.getErr(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())*
	muTightIsoSF.getErr(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt());
      lepWeight1MuLepUp=lepWeight1Mu+errMu;
      lepWeight1MuLepDown=lepWeight1Mu-errMu;
    }
  }
  for (size_t e = 0; e < (size_t)tightEl.size(); ++e){
    lept lepton;
    lepton.p4 = tightEl[e];
    lepton.flavour = 11;
    leptons.push_back(lepton);
    }
  
  for (size_t m = 0; m < (size_t)tightMu.size(); ++m){
    lept lepton;
    lepton.p4 = tightMu[m];
    lepton.flavour = 13;
    leptons.push_back(lepton);
    }

  std::sort(leptons.begin(), leptons.end(), by_pt());
  
  float w_lep(1.0),sf_lepUp(1.0),sf_lepDown(1.0);
  if(isData=="MC"){
    if(channel=="muon"){
      w_lep=lepWeight1Mu;
      if(w_lep!=0){
	sf_lepUp=lepWeight1MuLepUp/w_lep;
	sf_lepDown=lepWeight1MuLepDown/w_lep;
      }
    }
    if(channel=="electron"){
      w_lep=lepWeight1Ele;
      if(w_lep!=0){
	sf_lepUp=lepWeight1EleLepUp/w_lep;
	sf_lepDown=lepWeight1EleLepDown/w_lep;
      }
    }
    
    w*=w_lep;
  }
  
  systZero.setWeight(0,1.);
  systZero.setWeight("btagUp",1.);
  systZero.setWeight("btagDown",1.);
  systZero.setWeight("mistagUp",1.);
  systZero.setWeight("mistagDown",1.);
  systZero.setWeight("puDown",1.);
  systZero.setWeight("puUp",1.);
  systZero.setWeight("lepDown",1.);
  systZero.setWeight("lepUp",1.);

  if(addPDF)systZero.setPDFWeights(w_pdfs,nPDF,w_zero,true);
  if(addQ2)systZero.setQ2Weights(w_q2up,w_q2down,w_zero,true);
  if(addTopPt)systZero.setTWeight(w_top,topWeight,true);

  

  syst0BM.copySysts(systZero);
  syst0BM.setWeight(0,bWeight0CSVM);
  syst0BM.setWeight("btagUp",bWeight0CSVMBTagUp);
  syst0BM.setWeight("btagDown",bWeight0CSVMBTagDown);
  syst0BM.setWeight("mistagUp",bWeight0CSVMMisTagUp);
  syst0BM.setWeight("mistagDown",bWeight0CSVMMisTagDown);
  syst0BM.setWeight("puUp",bWeight0CSVM*puUpFact);
  syst0BM.setWeight("puDown",bWeight0CSVM*puDownFact);
  syst0BM.setWeight("lepUp",bWeight0CSVM*sf_lepUp);
  syst0BM.setWeight("lepDown",bWeight0CSVM*sf_lepDown);

  syst1BM.copySysts(systZero);
  syst1BM.setWeight(0,bWeight1CSVM);
  syst1BM.setWeight("btagUp",bWeight1CSVMBTagUp);
  syst1BM.setWeight("btagDown",bWeight1CSVMBTagDown);
  syst1BM.setWeight("mistagUp",bWeight1CSVMMisTagUp);
  syst1BM.setWeight("mistagDown",bWeight1CSVMMisTagDown);
  syst1BM.setWeight("lepUp",bWeight1CSVM*sf_lepUp);
  syst1BM.setWeight("lepDown",bWeight1CSVM*sf_lepDown);
 
  syst2BM.copySysts(systZero);
  syst2BM.setWeight(0,bWeight2CSVM);
  syst2BM.setWeight("btagUp",bWeight2CSVMBTagUp);
  syst2BM.setWeight("btagDown",bWeight2CSVMBTagDown);
  syst2BM.setWeight("mistagUp",bWeight2CSVMMisTagUp);
  syst2BM.setWeight("mistagDown",bWeight2CSVMMisTagDown);
  syst2BM.setWeight("puUp",bWeight2CSVM*puUpFact);
  syst2BM.setWeight("puDown",bWeight2CSVM*puDownFact);
  syst2BM.setWeight("lepUp",bWeight2CSVM*sf_lepUp);
  syst2BM.setWeight("lepDown",bWeight2CSVM*sf_lepDown);
  
  vector<float> jetsPhi;
  vector<TLorentzVector> bjets, jets_nob, jets, jets20;    
  
  struct btag{
      TLorentzVector vect;
      float csv;
      };
  
  std::vector<btag> bvects,extrajets;
  struct by_csv{
      bool operator()(btag const &b1, btag const &b2){
      return b1.csv > b2.csv;
      }
      };

  struct by_pt_jet{
      bool operator()(btag const &jet1, btag const &jet2){
      return jet2.vect.Pt() < jet2.vect.Pt();
      }
  };
  
   for (int j = 0; j <maxJet20Loop;++j){	
     //     cout <<" in j20 loop"<<endl;
     TLorentzVector jet20;
     if(jet20Pt[j]<40){
       //       cout <<" # j= "<<j<< " pt "<< jet20Pt[j] <<endl;
       jet20.SetPtEtaPhiE(jet20Pt[j], jet20Eta[j], jet20Phi[j], jet20E[j]);
       jets20.push_back(jet20);
       btag bj;
       bj.vect = jet20;
       bj.csv = jet20ak4chs_csvv2[j];
       extrajets.push_back(bj);
     }
   }
   
   for (int j = 0; j <maxJetLoop;++j){	
      TLorentzVector jet, all_jets;
      
      all_jets.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);
      jets.push_back(all_jets);
      bool btagcond = jetIsCSVM[j]>0.;
      if(useCSVTSelection){
	btagcond = jetIsCSVT[j]>0.;
      }

  if(!btagcond){
      jet.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);
      jets_nob.push_back(jet);
      }
  
  jetsPhi.push_back(jetPhi[j]);
  TLorentzVector bjet;
  
  if( btagcond && abs(jetEta[j])<2.4){
     bjet.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);
     bjets.push_back(bjet);
     btag b;
     b.vect = bjet;
     b.csv = jetak4chs_csvv2[j];
     bvects.push_back(b);
    }  
  }// End of jet loop 
 
  nJets = jetsPhi.size();
  nCSVJets=bjets.size();
  //  cout <<"passes trig "<<endl; 
  
  //  std::sort(jets20.begin(), jets20.end(), by_pt()); 
  std::sort(bvects.begin(), bvects.end(), by_pt_jet()); 
  std::sort(extrajets.begin(), extrajets.end(), by_pt_jet()); 
  //  bool passmuon = muonTrigger && nMu == 1 && muLooseSize==1 && nEl==0 && elLooseSize == 0;
  bool passmuon = muonTrigger && nMu == 1 ;
  if(passmuon){n_lepton+=w;nev_lepton+=1;}

  passmuon = passmuon && muLooseSize==1;
  if(passmuon){ n_loose_veto+=w; nev_loose_veto+=1;}
  
  passmuon = passmuon && nEl==0 && elLooseSize == 0;
  if(passmuon){n_lepton_cross_veto+=w;nev_lepton_cross_veto+=1;}

  //bool passantiisomuon = muonTrigger && nMu == 1 && muLooseSize == 0 && nEl==0 && elLooseSize == 0;
  //bool passantiisomuon = muonTrigger && muAntiIsoSize==1 && nMu==1 && nEl==0 && elLooseSize == 0;
  bool passantiisomuon = muonTrigger && muAntiIsoSize==1 && nEl==0 && elLooseSize == 0;
  bool passelectron = false;
  bool passsinglelepton = passmuon || passelectron;


  if(channel=="muonantiiso" && !passantiisomuon) continue;
  if(channel=="muon" && !passmuon)continue; 
  if(channel=="electron" && !passelectron)continue;
  if(channel=="singlelepton" && !passsinglelepton)continue;
  
  if(jets.size() ==2){ n_2j+=w; nev_2j+=1;}
  if(jets.size() ==2 && bjets.size()==1){ n_2j1t+=w; nev_2j1t+=1;  
    if(doSynch){
      if((tightMu.size())==1){
	TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
	float phi_lmet = fabs(deltaPhi(tightMu[0].Phi(), metPhi[0]) );
	mt = sqrt(2* tightMu[0].Pt() * met* ( 1- cos(phi_lmet)));
      }
      if(mt> 50){
	n_2j1tmtw+=w;
	nev_2j1tmtw+=1.;
      }
    }
  }
  //2j0t 
  if((jets.size() ==2 && bjets.size() == 0)){
    
    for(size_t i = 0; i < (size_t)tightMu.size();++i ){
        if((tightMu.size())< 2){
          TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
          float phi_lmet = fabs(deltaPhi(tightMu[i].Phi(), metPhi[0]) );
          mt = sqrt(2* tightMu[i].Pt() * met* ( 1- cos(phi_lmet)));
	  //	  cout << "muon iso "<< selectedIso[i]<< " mtw "<< mt <<endl;
	  
          syst0BM.fillHistogramsSysts(h_2j0t_mtw,mt,w);
          }
      }
    for(size_t j= 0; j< (size_t)jets.size();++j ){
    //    cout <<runNumber<<"   "<<evtNumber<<"   "<<"jet["<<i<< "] "<<jets[i].Pt()<<"   bjetsize "<< bjets.size() <<"   Mu["<<i<< "] "<<tightMu[i].Pt()<< std::endl;
    //    cout << "weight 0 "  << bWeight0CSVM <<" w "<< w <<endl;
      if(j==0)syst0BM.fillHistogramsSysts(h_2j0t_jetpt40_1st,jets[j].Pt(),w);  
      if(j==1)syst0BM.fillHistogramsSysts(h_2j0t_jetpt40_2nd,jets[j].Pt(),w);  
      if(mt>50.0){
	if(j==0){
	  syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_jetpt40_1st,jets[j].Pt(),w);  
	  syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_jeteta40_1st,jets[j].Pt(),w);  
	  syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_jetabseta40_1st,fabs(jets[j].Eta()),w);  

	  syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_incljeteta,jets[j].Eta(),w);  
	  syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_incljetabseta,fabs(jets[j].Eta()),w);  
	}
	if(j==1){
	  syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_jetpt40_2nd,jets[j].Pt(),w);  
	  syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_incljeteta,jets[j].Eta(),w);  
	  syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_incljetabseta,fabs(jets[j].Eta()),w);  
	}
      }
    }
  }

  double etathreshold=2.5;
  //2j1t
  if((jets.size() == 2 && bjets.size() == 1)){
    for(size_t b= 0; b< (size_t)bjets.size();++b ){
      //      cout << " w "<< w <<" b weight 1 "<< bWeight1CSVM << endl;
      //      cout <<" pu weight = "<< w_pu<<endl; 
      //      cout <<" lep weight = "<< w_lep<<endl; 

      if(b==0){
	syst1BM.fillHistogramsSysts(h_2j1t_bjetpt,bjets[b].Pt(),w,NULL,false);
      //cout << bjets[i].Pt()<< endl;
      }
    }
    bool qcddepleted=false, signalenriched=false;
    //define QCD depletion condition:
    for(size_t i = 0; i < (size_t)tightMu.size();++i ){
      syst1BM.fillHistogramsSysts(h_2j1t_MuPt,tightMu[i].Pt(),w);
      syst1BM.fillHistogramsSysts(h_2j1t_MuEta,tightMu[i].Eta(),w);
      syst1BM.fillHistogramsSysts(h_2j1t_MuPhi,tightMu[i].Phi(),w);
      syst1BM.fillHistogramsSysts(h_2j1t_MuE,tightMu[i].E(),w);
      
      if((tightMu.size())<2 ){
        TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
        float phi_lmet = fabs(deltaPhi(tightMu[i].Phi(), metPhi[0]) );
        mt = sqrt(2* tightMu[i].Pt() * met* ( 1- cos(phi_lmet)));
        //cout <<" mt = "<<mt<<endl;
        //cout <<" calculate_mtw() = "<<calculate_mtw(metPt, metPhi,tightMu)<<endl;     
        syst1BM.fillHistogramsSysts(h_2j1t_mtw,mt,w); 
        if (calculate_mtw(metPt, metPhi,tightMu) > 50.0){
	  qcddepleted=true;
          syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_mtw,mt,w);

	} 
      }
    }
    //define signal enriching condition:
    for (size_t j= 0; j< (size_t)jets.size();++j ){
      bool btagcond = jetIsCSVM[j]>0.;
      if(useCSVTSelection){
	btagcond = jetIsCSVT[j]>0.;
      }
      if(!(btagcond)){
	syst1BM.fillHistogramsSysts(h_2j1t_ljetpt,jets[j].Pt(),w);
	syst1BM.fillHistogramsSysts(h_2j1t_ljeteta,jets[j].Eta(),w);
	syst1BM.fillHistogramsSysts(h_2j1t_etajprime,fabs(jets[j].Eta()),w);
	signalenriched=signalenriched||fabs(jets[j].Eta())>etathreshold;
	if(qcddepleted) {
	  syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_ljetpt,jets[j].Pt(),w);
	  syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_ljeteta,jets[j].Eta(),w);
	  syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_etajprime,fabs(jets[j].Eta()),w);
	}
      } 
    }
    //    cout << " filling extrajets "<<endl;
    syst1BM.fillHistogramsSysts(h_2j1t_nextrajets,extrajets.size(),w);
    //    cout << " filled nextrajets extrajets "<<endl;
    if(extrajets.size()){
      //      syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetpt,jets20.at(0).Pt(),w);
      //      cout << "method 1 ok "<<endl;
      syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetpt,extrajets.at(0).vect.Pt(),w);
      //      cout << "method 2 ok "<<endl;
      syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetcsv,extrajets.at(0).csv,w);
      //      cout << "csv ok"<<endl;
    }
    //    cout << "extrajets ok"<<endl;
    
    //Adding top related variables
    vector<TLorentzVector> tops= topUtils.top4Momenta(tightMu, bjets ,metpx, metpy);
    if(tops.size()!=0){
      double mass = tops.at(0).M();
      double pt = tops.at(0).Pt();
      syst1BM.fillHistogramsSysts(h_2j1t_topMass,mass,w);
      syst1BM.fillHistogramsSysts(h_2j1t_topPt,pt,w);
      //      if(calculate_mtw(metPt, metPhi,tightMu) > 50.0) {
      if(qcddepleted){ 
	syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_topMass,mass,w);
        syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_topPt,pt,w);
      }
      if(signalenriched && qcddepleted){
	syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_topMass,mass,w);
        syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_topPt,pt,w);
      }
    }
    //Adding extra-jet related variables: only for QCD and signal enriched regions.
    //    cout << "tops ok"<<endl;
    if(qcddepleted){
       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_nextrajets,extrajets.size(),w);
       if(signalenriched)syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_nextrajets,extrajets.size(),w);
       //6&       cout << "extrajet sizes ok"<<endl;

       if(extrajets.size()){
	 vector<TLorentzVector> topsextra= topUtils.top4Momenta(tightMu, jets20 ,metpx, metpy);
	 if(topsextra.size()!=0){
	   syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_topMassExtra,topsextra[0].M(),w);
	   syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt,extrajets.at(0).vect.Pt(),w);
	   syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv,extrajets.at(0).csv,w);
	   
	   //	   cout << "extrajet csvs ok"<<endl;

	 }
	 if(signalenriched){
	   //	   cout<< " filling histograms for the signal enriched region"<<endl;
	   syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra,topsextra[0].M(),w);
	   syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt,extrajets.at(0).vect.Pt(),w);
	   syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv,extrajets.at(0).csv,w);

	   //	   cout << "extrajet csvs sr ok"<<endl;
	 }
       }
    }
    //AOB for 2j1t?
  }
  //End 2j1t

  //3j1t   
  if((jets.size() == 3 && bjets.size() ==1)) {
  for(size_t b= 0; b< (size_t)bjets.size();++b ){
  if(b==0)systZero.fillHistogramsSysts(h_3j1t_bjetpt,bjets[b].Pt(),w);
   } 
  for(size_t i = 0; i < (size_t)tightMu.size();++i ){
        syst1BM.fillHistogramsSysts(h_3j1t_MuPt,tightMu[i].Pt(),w);
        syst1BM.fillHistogramsSysts(h_3j1t_MuEta,tightMu[i].Eta(),w);
        syst1BM.fillHistogramsSysts(h_3j1t_MuPhi,tightMu[i].Phi(),w);
        syst1BM.fillHistogramsSysts(h_3j1t_MuE,tightMu[i].E(),w);
        if((tightMu.size())<2 ){
          TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
          float phi_lmet = fabs(deltaPhi(tightMu[i].Phi(), metPhi[0]) );
          mt = sqrt(2* tightMu[i].Pt() * met* ( 1- cos(phi_lmet)));
          systZero.fillHistogramsSysts(h_3j1t_mtw,mt,w);
          }
      }
  //AOB for 3j1t?
  }
  //End 3j1t


  //3j2t
  if((jets.size() == 3 && bjets.size()==2)){
    //      cout << " w "<< w <<" b weight 2 "<< bWeight2CSVM << endl;
    //      cout <<" pu weight = "<< w_pu<<endl; 
    //      cout <<" lep weight = "<< w_lep<<endl; 

    //    cout << " event "<<evt<<endl;
    for(size_t b=0; b< (size_t)bjets.size();++b ){
      //      cout << " ,bWeight2CSV "<< bWeight2CSVM <<",bWeight2CSVMBTagUp "<<bWeight2CSVMBTagUp<<endl;
      // cout << " i is "<< i << endl;
      if(b==0)syst2BM.fillHistogramsSysts(h_3j2t_bjetpt,bjets[b].Pt(),w,NULL,false);
      if(b==1)syst2BM.fillHistogramsSysts(h_3j2t_2ndbjetpt,bjets[b].Pt(),w);
    }
    for(size_t i = 0; i < (size_t)tightMu.size();++i ){
      if(tightMu[i].Pt()>20){
        syst2BM.fillHistogramsSysts(h_3j2t_MuPt,tightMu[i].Pt(),w);
        syst2BM.fillHistogramsSysts(h_3j2t_MuEta,tightMu[i].Eta(),w);
        syst2BM.fillHistogramsSysts(h_3j2t_MuPhi,tightMu[i].Phi(),w);
        syst2BM.fillHistogramsSysts(h_3j2t_MuE,tightMu[i].E(),w);
        if((tightMu.size())<2 ){
          TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
          float phi_lmet = fabs(deltaPhi(tightMu[i].Phi(), metPhi[0]) );
          mt = sqrt(2* tightMu[i].Pt() * met* ( 1- cos(phi_lmet)));
          syst2BM.fillHistogramsSysts(h_3j2t_mtw,mt,w);
	}
      }
    }
    for(size_t j=0; j< (size_t)jets.size();++j ){
      //      bool btagcond=false;
      bool btagcond = jetIsCSVM[j]>0.;
      if(useCSVTSelection){
	btagcond = jetIsCSVT[j]>0.;
      }

      if(!btagcond){
	syst2BM.fillHistogramsSysts(h_3j2t_ljetpt,(jets.at(j).Pt()),w);
	syst2BM.fillHistogramsSysts(h_3j2t_ljeteta,(jets.at(j).Eta()),w);
	syst2BM.fillHistogramsSysts(h_3j2t_etajprime,(fabs(jets.at(j).Eta())),w);
      }
    }
    syst2BM.fillHistogramsSysts(h_3j2t_nextrajets,extrajets.size(),w);
    /* TH1F *h_3j2t_ljetpt[maxSysts];   systZero.initHistogramsSysts(h_3j2t_ljetpt,    "h_3j2t_ljetpt",     "3j2t b jet pt ",100,0,500);
    TH1F *h_3j2t_ljeteta[maxSysts];   systZero.initHistogramsSysts(h_3j2t_ljeteta,    "h_3j2t_ljeteta",     "3j2t b jet pt ",100,0,500);
    TH1F *h_3j2t_etajprime[maxSysts];   systZero.initHistogramsSysts(h_3j2t_etajprime,    "h_3j2t_etajprime",     "3j2t b jet pt ",100,-4.7,4.7);*/
					  ;

    //AOB for 3j2t?
  }
  //End 3j2t
  
  //To add here: 2l3j2t: 2 leptons, 3 jets, 2 tags for ttbar control

  //End signal and control regions.

  systZero.fillHistogramsSysts(h_nJets,nJets,w); 
  systZero.fillHistogramsSysts(h_nbJets,nCSVJets,w); 
  systZero.fillHistogramsSysts(h_nPV,nPV,w);
  systZero.fillHistogramsSysts(h_nGoodPV,nPV,w);
  systZero.fillHistogramsSysts(h_nTruePV,numTrueInt,w);
  
  //if(nTightElectrons != nEl)cout << "warning! problem with tight el"<<endl;
  //if(nTightMuons != nMu)cout << "warning! problem with tight mu"<<endl;
  
  nTightElectrons = nEl;
  nTightMuons = nMu;    


  if(doSynch)fileout<<std::fixed<<std::setprecision(0)<<runNumber<<"   "<<evtNumber<<"   "<<lumiSec<<"   "<<std::endl;
  
  }//end of loop over events 
  
 
 if(doSynch){
   //  int nev_trig(0),nev_lepton(0),nev_loose_veto(0),nev_lepton_cross_veto(0),nev_2j(0),nev_2j1t(0);
   fileout<<" ntrig "<< " nlep "<< " nlooseveto "<< " nloosevetoele "<< " n2jets "<< " n2jets1tag "<< " n2jets1tagmtw "<<endl;
   fileout<< " " << nev_trig << " " << nev_lepton << " " << nev_loose_veto << " " << nev_lepton_cross_veto << " " << nev_2j << " " << nev_2j1t <<" "<< nev_2j1tmtw <<endl;
   fileout<<" wtrig "<< " wlep "<< " wlooseveto "<<" wloosevetoele "<< " w2jets "<< " w2jets1tag "<<" w2jets1tagmtw "<<endl;
   fileout<< " " << n_trig << " " << n_lepton << " " << n_loose_veto <<n_lepton_cross_veto << " " << n_2j << " " << n_2j1t <<" "<< n_2j1tmtw<<endl;

  fileout.close();  
 }//return h

 //Step 0:get only the negative weights version, 
 if(isData=="MC"){
   //chainNEvents.Project("h_weight_sign","Event_LHEWeight10/abs(Event_LHEWeight10)"); // I messed up and weight zero is not  saved in the full chain... but this one should do the trick as well! Just do h_weight_sign->GetMean() and multiply this to get the correct number :)
   cout << " mean  "<< zeroWeight<<endl;
 }
 // TO FIX! ADD PROJECTION OF ALL PDFS TO THE EVENT!
 // systZero.projectAllPDF(h_weightZero,chainNEvents);

 bool hasLHEWeights=true; 
 cout << "haslheweights?" << hasLHEWeights<<endl;
 if(sample.find("QCDMu")!=std::string::npos || sample.find("_tW_")!=std::string::npos || isData=="DATA")hasLHEWeights=false;
 // Cut Flow
 
 if(hasLHEWeights) {
   h_cutFlow->SetBinContent(0,nEventsTot*zeroWeight);//Underflow: number of events pre-preselection.}
   cout << " nEvents tot "<< nEventsTot << " * weight"<< zeroWeight<<" = "<< h_cutFlow->GetBinContent(0)<<endl; 
 }
 else{ h_cutFlow->SetBinContent(0,nEventsTot);//Underflow: number of events pre-preselection.
   cout << " nEvents tot "<< nEventsTot <<" = "<< h_cutFlow->GetBinContent(0)<<endl; 
}

  h_cutFlow->SetBinContent(1,nEvents);
  h_cutFlow->GetXaxis()->SetBinLabel(1,"no selection");
  h_cutFlow->SetBinContent(2, n_trig);
  h_cutFlow->GetXaxis()->SetBinLabel(2, "trigger");
  h_cutFlow->SetBinContent(3,n_lepton);
  h_cutFlow->GetXaxis()->SetBinLabel(3, "=1 Iso. Muon");
  h_cutFlow->SetBinContent(4,n_lepton_cross_veto);
  h_cutFlow->GetXaxis()->SetBinLabel(4, "=0 Loose Ele.");
  h_cutFlow->SetBinContent(5,n_2j);
  h_cutFlow->GetXaxis()->SetBinLabel(5, "2-jet");
  h_cutFlow->SetBinContent(6,n_2j1t);
  h_cutFlow->GetXaxis()->SetBinLabel(6, "2-jet 1-tag");

  
  //Write the Histogramms here  
  systZero.writeSingleHistogramSysts(h_cutFlow, allMyFiles); 
  systZero.writeSingleHistogramSysts(h_weight_sign, allMyFiles); 
  systZero.writeHistogramsSysts(h_nPV, allMyFiles); 
  systZero.writeHistogramsSysts(h_nGoodPV, allMyFiles); 
  systZero.writeHistogramsSysts(h_nTruePV, allMyFiles); 
  
  systZero.writeHistogramsSysts(h_2j0t_jetpt40_1st, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_jetpt40_2nd, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtw,   allMyFiles); 

  systZero.writeHistogramsSysts(h_2j0t_mtwcut_jeteta40_1st,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_jetabseta40_1st,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_jetpt40_1st,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_jetpt40_2nd,   allMyFiles); 
  
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_incljeteta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_incljetabseta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_incljeteta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_incljetabseta,   allMyFiles); 
  
  //2j1t
  systZero.writeHistogramsSysts(h_2j1t_bjetpt,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_ljetpt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_ljeteta,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuPt,    allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuEta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuPhi,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuE,     allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuCharge,allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtw,     allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_topMass, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_topPt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_nextrajets,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_leadingextrajetpt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_leadingextrajetcsv,   allMyFiles); 

  systZero.writeHistogramsSysts(h_2j1t_mtwcut_mtw,     allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_topMass, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_topPt, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_bjetpt, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_ljetpt, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_ljeteta, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_etajprime, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_nextrajets,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_topMassExtra,   allMyFiles); 

  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_topMass,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_topPt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_nextrajets,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra,   allMyFiles); 

  
 


  //3j1t
  systZero.writeHistogramsSysts(h_3j1t_bjetpt,  allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_MuPt,    allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_MuEta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_MuPhi,   allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_MuE,     allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_mtw,     allMyFiles); 
  
  systZero.writeHistogramsSysts(h_3j2t_bjetpt,  allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_2ndbjetpt, allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuPt,    allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuEta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuPhi,   allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuE,     allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuCharge,allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_mtw,     allMyFiles); 

  systZero.writeHistogramsSysts(h_3j2t_ljetpt,  allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_ljeteta,allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_etajprime,allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_nextrajets,allMyFiles); 
  
  systZero.writeHistogramsSysts(h_nJets, allMyFiles); 
  systZero.writeHistogramsSysts(h_nbJets, allMyFiles); 


  
  
  systZero.closeFilesSysts(allMyFiles);
  
  std::cout<< "Info: Total No. of events for sample["<<sample<<"] : "<<nEvents<<std::endl;
}

void systWeights::closeFilesSysts(TFile *filesout[(int)MAXSYSTS]){
  int MAX= this->maxSystsNonPDF;
  //  int MAXTOT= this->maxSysts;
  
  bool useOnlyNominal = this->onlyNominal;
  for (int c = 0; c < this->nCategories; c++){
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      //cout << "c is now "<< c << " sy "<< sy << " location "<< sy+(MAXTOT+1)*c <<" is histo there? " << histo[sy+(MAXTOT+1)*c] << " file location "<<sy+(MAX+1)*c << " is file there "<< filesout[sy+(MAX+1)*c]<< endl;
      if(!(!useOnlyNominal || sy==0)) continue;
      filesout[(int)sy+(MAX+1)*(c)]->Close();
    }
    
    if(this->addPDF){
      if(!useOnlyNominal){
	filesout[MAX+(MAX+1)*(c)]->Close();
      }
    }
  }
}
void systWeights::setOnlyNominal(bool useOnlyNominal){
  this->onlyNominal=useOnlyNominal;
}
void systWeights::setSystValue(string name, double value, bool mult){
  float zerofact=1.0;
  if(mult)zerofact=this->weightedSysts[0];
  int MAX = this->maxSysts;
  for(int sy=0;sy<(int)MAX;++sy){
    if(this->weightedNames[(int)sy] ==name){
      this->weightedSysts[(int)sy] =value*zerofact;
    }
  }
}


void systWeights::setSystValue(int place, double value, bool mult){
  float zerofact=1.0;
  if(mult)zerofact=this->weightedSysts[0];
  this->weightedSysts[place] =value*zerofact;
}

void systWeights::setWeight(string name, double value, bool mult){
  this->setSystValue(name, value, mult);
}

void systWeights::setWeight(int place, double value, bool mult){
  this->setSystValue(place, value, mult);
}

void systWeights::initHistogramsSysts(TH1F** histo,TString name, TString title, int nbins, float min, float max){
  for (int c = 0; c < this->nCategories; c++){
    int MAX = this->maxSysts;
    bool useOnlyNominal = this->onlyNominal;
    TString cname= (this->categoriesNames[c]).c_str();
    for(int sy=0;sy<(int)MAX;++sy){
      TString ns= (this->weightedNames[sy]).c_str();
      if(sy==0){
          if(c==0) histo[sy+((MAX+1)*(c))]=new TH1F(name,title,nbins,min,max);
          else histo[sy+((MAX+1)*(c))]=new TH1F(name+"_"+cname,title,nbins,min,max);
          }
      
      if(sy!=0 && !useOnlyNominal) {
          if(c==0)histo[sy+(MAX+1)*c]=new TH1F(name+"_"+ns,title,nbins,min,max);
          else histo[sy+(MAX+1)*c]=new TH1F(name+"_"+ns+"_"+cname,title,nbins,min,max);
          }
          //cout << " initialized histogram "<< histo[sy+(MAX+1)*c]->GetName() <<" sy " << sy << " c  "<< c <<" location " << sy+(MAX+1)*c << endl;
    }
  }
}
void writeHistogramsSysts(TH1F* histo[(int)MAXSYSTS], TFile *filesout[(int)MAXSYSTS], bool useOnlyNominal=false){
  for(int sy=0;sy<(int)MAXSYSTS;++sy){
    //cout << " writing histo "<< histo[(int)sy]->GetName()<< " in file "<< filesout[(int)sy]->GetName()<<endl;;
    //TString ns= weightedSystsNames((weightedSysts)sy);
    filesout[(int)sy]->cd();
    histo[sy]->Write(histo[0]->GetName());
    //    histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
  }
}
void systWeights::writeHistogramsSysts(TH1F** histo, TFile **filesout){
  int MAX= this->maxSystsNonPDF;
  int MAXTOT= this->maxSysts;
  bool useOnlyNominal = this->onlyNominal;
  for (int c = 0; c < this->nCategories; c++){
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      //cout << "c is now "<< c << " sy "<< sy << " location "<< sy+(MAXTOT+1)*c <<" is histo there? " << histo[sy+(MAXTOT+1)*c] << " file location "<<sy+(MAX+1)*c << " is file there "<< filesout[sy+(MAX+1)*c]<< endl;
      
      //cout << " writing histo "<< histo[sy+(MAXTOT+1)*c]->GetName()<< " in file "<< filesout[sy+(MAX+1)*c]->GetName()<<endl;;
      
      //TString ns= weightedSystsNames((weightedSysts)sy);
      if(!(!useOnlyNominal || sy==0)) continue;
      filesout[(int)sy+(MAX+1)*(c)]->cd();
      if(this->addPDF){
      //if(this->weightedNames[sy]=="pdf_totalUp")calcPDFHisto(histo, histo[sy+(MAXTOT+1)*(c)],1.0,c);
      //if(this->weightedNames[sy]=="pdf_totalDown")calcPDFHisto(histo, histo[sy+(MAXTOT+1)*(c)],-1.0,c);
      ;      //this->
      }

      histo[sy+(MAXTOT+1)*c]->Write(histo[0]->GetName());
      //histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
      }
      
      if(this->addPDF){
      if(!useOnlyNominal){
        filesout[MAX+(MAX+1)*(c)]->cd();
        //cout << " file max is "<< filesout[MAX+(MAX+1)*c]->GetName()<<endl;
        //int npdf=this->maxSysts-this->maxSystsNonPdf;
        int MAXPDF=this->maxSysts;
        for(int sy=MAX;sy<MAXPDF;++sy){
        //    cout << " writing sy "<<sy+(MAXTOT+1)*c<<endl;
        //    cout << " histo is there? "<< histo[sy+(MAXTOT+1)*c]<<endl;
        histo[sy+(MAXTOT+1)*(c)]->Write();
        //    cout << " written sy "<< histo[sy+(MAXTOT+1)*c]->GetName()<<endl;
        }
      }
    }

  }

}

void systWeights::createFilesSysts(  TFile ** allFiles, TString basename, TString opt){
  for (int c = 0; c < this->nCategories; c++){
    int MAX = this->maxSystsNonPDF;
    int MAXTOT = this->maxSystsNonPDF;
    bool useOnlyNominal = this->onlyNominal;
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      TString ns= (this->weightedNames[(int)sy]);
      cout << " creating file for syst "<< ns<<endl;
      if (c!=0)     cout << " category is "<< c<<endl;
      cout << "onlynominal is "<<useOnlyNominal<<endl;

      if(sy==0){
	      cout<<" filename is "<< basename+ns+cname+".root"<<endl;
	      allFiles[sy+(MAX+1)*c]= TFile::Open((basename+ns+cname+".root"), opt);
        }
      
      else{
	      if(!useOnlyNominal){
	          //if((ns!="1lep") && (ns!="2lep")&& (ns!="0lep")){
	          cout<<" filename is "<< basename+ns+cname+".root"<<endl;
	          allFiles[sy+(MAX+1)*c]= TFile::Open((basename+"_"+ns+cname+".root"), opt);
	          }
          }
        //TFile *outTree = TFile::Open(("trees/tree_"+outFileName).c_str(), "RECREATE");
        //cout << " created file at c "<< c << " s "<< sy << " location "<< sy+(MAXTOT+1)*c<< " fname "<<allFiles[sy+(MAXTOT+1)*c]->GetName()<<endl;   
        }
      
      if(this->addPDF){
      if(!useOnlyNominal)allFiles[MAX+((MAX+1)*c)]= TFile::Open((basename+"_pdf"+cname+".root"), opt);
      //cout << " created file at c "<< c << " s "<< MAX+(MAX+1)*c << " location "<< MAX+(MAX+1)*c<<endl;
      cout<< " fname "<<allFiles[MAX+(MAXTOT+1)*c]->GetName()<<endl;
        }
      }
    //return allFiles;
}



 
void systWeights::fillHistogramsSysts(TH1F** histo, float v, float w, double * wcats, bool verbose ){
  if(wcats==NULL){
    wcats=this->wCats;
  }
  for (int c = 0; c < this->nCategories; c++){
    if(verbose){
      cout<< " in cat loop "<< c<<endl;
      cout<< " value "<< wcats[c] <<endl;
    }
    int MAX = this->maxSysts;
    bool useOnlyNominal = this->onlyNominal;
    for(int sy=0;sy<(int)MAX;++sy){
      if(verbose){
	cout<< " in syst loop "<< sy<< endl;
	cout<<" value "<< this->weightedSysts[(int)sy] <<endl ;
      }
      if(sy!=0&& useOnlyNominal)continue;
      float ws = (this->weightedSysts[(int)sy])*wcats[c];
      if(verbose)cout << " filling histogram "<< histo[(int)sy]->GetName() << " with value "<< v <<" and weight "<< w <<" ws "<< ws<<endl;
      histo[sy+(MAX+1)*(c)]->Fill(v, w*ws);
    }
  }
}

//void  initHistogramsSysts (TH1F* histo[(int)MAXSYSTS],TString name, TString title, int nbins, float min, float max, bool useOnlyNominal=false){
//    for(int sy=0;sy<(int)MAXSYSTS;++sy){
//          //TString ns= weightedSystsNames((weightedSysts)sy);
//          histo[sy]=new TH1F(name,title,nbins,min,max);
//          }
//        }

void systWeights::setWCats(double * wcats){
  for(int i =0;i<this->nCategories;++i){
    //    cout << "setting wcat #"<< i << " to be "<<wcats[i]<<endl;
    this->wCats[i]=wcats[i];
  }
 
}

void systWeights::copySysts(systWeights sys){
  for(int i =0; i < sys.maxSysts;++i){
    this->weightedNames[i]=sys.weightedNames[i];
    this->weightedSysts[i]=sys.weightedSysts[i];

  }
  this->setOnlyNominal(sys.onlyNominal);
  this->setMax(sys.maxSysts);
  this->setMaxNonPDF(sys.maxSystsNonPDF);
  this->nPDF=sys.nPDF;
  this->nCategories=sys.nCategories;  
  this->addQ2=sys.addQ2;
  this->addPDF=sys.addPDF;
  this->addTopPt=sys.addTopPt;
  this->addVHF=sys.addVHF;
  this->addTTSplit=sys.addTTSplit;
  this->setWCats(sys.wCats);

}

void systWeights::setMax(int max){
  this->maxSysts =  max;
  }
void systWeights::setMaxNonPDF(int max){
  this->maxSystsNonPDF =  max;
  }

void systWeights::prepareDefault(bool addDefault, bool addQ2, bool addPDF, bool addTopPt,bool addVHF, bool addTTSplit, int numPDF){ 
  this->addPDF=addPDF;
  this->addQ2=addQ2;
  this->addTopPt=addTopPt;
  this->addVHF=addVHF;
  this->addTTSplit=addTTSplit;
  this->nPDF=numPDF;
  this->nCategories=1;
  categoriesNames[0]="";
  this->wCats[0]=1.0;
  if(addDefault){
    this->weightedNames[0]="";
    this->weightedNames[1]="btagUp";
    this->weightedNames[2]="btagDown";
    this->weightedNames[3]="mistagUp";
    this->weightedNames[4]="mistagDown";
    this->weightedNames[5]="puUp";
    this->weightedNames[6]="puDown";
    this->weightedNames[7]="lepUp";
    this->weightedNames[8]="lepDown";
    this->setMax(9);
    this->setMaxNonPDF(9);
    this->weightedNames[this->maxSysts]="";
  }
  if(addQ2){
    this->weightedNames[this->maxSysts]= "q2Up";
    this->weightedNames[this->maxSysts+1]= "q2Down";
    this->setMax(this->maxSysts+2);
    this->setMaxNonPDF(this->maxSystsNonPDF+2);
    this->weightedNames[this->maxSysts]= "";
  }

  if(addTopPt){
    this->weightedNames[this->maxSysts]="topPtWeightUp";
    this->weightedNames[this->maxSysts+1]="topPtWeightDown";
    this->setMax(this->maxSysts+2);
    this->setMaxNonPDF(this->maxSystsNonPDF+2);
    this->weightedNames[this->maxSysts]= "";
  }

  if(addVHF){
    this->weightedNames[this->maxSysts]="VHFWeightUp";
    this->weightedNames[this->maxSysts+1]="VHFWeightDown";
    this->setMax(this->maxSysts+2);
    this->setMaxNonPDF(this->maxSystsNonPDF+2);
    this->weightedNames[this->maxSysts]= "";
  }

  if(addTTSplit){
    this->nCategories=4;
    categoriesNames[1]="TT0lep";
    categoriesNames[2]="TT1lep";
    categoriesNames[3]="TT2lep";
    this->wCats[1]=1.0;
    this->wCats[2]=1.0;
    this->wCats[3]=1.0;
    }

  if(addPDF){
    this->weightedNames[this->maxSysts]= "pdf_totalUp";
    this->weightedNames[this->maxSysts+1]= "pdf_totalDown";
    this->weightedNames[this->maxSysts+2]= "pdf_asUp";
    this->weightedNames[this->maxSysts+3]= "pdf_asDown";
    this->weightedNames[this->maxSysts+4]= "pdf_zmUp";
    this->weightedNames[this->maxSysts+5]= "pdf_zmDown";
    this->setMax(this->maxSysts+6);
    this->setMaxNonPDF(this->maxSystsNonPDF+6);
    int nPDF=this->nPDF;
    for(int i =0; i < nPDF;++i){
      stringstream ss;
      ss<< i+1;
      this->weightedNames[i+this->maxSysts]= "pdf"+ss.str();
      }
    this->setMax(maxSysts+nPDF);
    this->weightedNames[this->maxSysts]= "";
    }
}

TString  weightedSystsNames (weightedSysts sy){
  switch(sy){
  case NOSYST : return "";
  case BTAGUP : return "btagUp";
  case BTAGDOWN : return "btagDown";
  case MISTAGUP : return "mistagUp";
  case MISTAGDOWN : return "mistagDown";
  case PUUP : return "puUp";
  case PUDOWN : return "puDown";
  case LEPUP : return "lepUp";
  case LEPDOWN : return "lepDown";
  //case ISOUP : return "isoUp";
  //case ISODOWN : return "isoDown";
  //case TRIGUP : return "trigUp";
  //case TRIGDOWN : return "trigDown";
  case MAXSYSTS : return "";
  }
  return "noSyst";
}

void systWeights::writeSingleHistogramSysts(TH1F* histo, TFile **filesout){  
  int MAX= this->maxSystsNonPDF;
  bool useOnlyNominal = this->onlyNominal;
  for (int c = 0; c < this->nCategories; c++){
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      if(!(!useOnlyNominal || sy==0)) continue;
      //cout << " writing histo "<< histo->GetName()<< " in file "<< filesout[(int)sy]->GetName()<<endl;;
      filesout[(int)sy+(MAX+1)*c]->cd();
      histo->Write();
      //    histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
    }
    if(this->addPDF){
      if(!useOnlyNominal){
        filesout[MAX+(MAX+1)*c]->cd();
        int MAXPDF=this->maxSysts;
        for(int sy=MAX;sy<MAXPDF;++sy){
        //cout << " writing sy "<< histo[sy]->GetName()<<endl;
        histo->Write();
        //cout << " written sy "<< histo[sy]->GetName()<<endl;
        }
      }
    }
  }
}


void systWeights::setPDFWeights(float * wpdfs, int numPDFs, float wzero,bool mult){
  float zerofact=1.0;
  if(mult)zerofact=this->weightedSysts[0];
  for (int i = 0; i < numPDFs; ++i){
    this->setPDFValue(i,zerofact*wpdfs[i]/wzero);
  }
  this->setSystValue("pdf_asUp", this->getPDFValue(this->nPDF-2)/wzero);
  this->setSystValue("pdf_asDown", zerofact);
  this->setSystValue("pdf_zmUp", this->getPDFValue(this->nPDF-1)/wzero);
  this->setSystValue("pdf_zmDown", zerofact);
  this->setSystValue("pdf_totalUp", zerofact);
  this->setSystValue("pdf_totalDown", zerofact);
}

//void systWeights::setTWeight(float tweight, float totalweight){
void systWeights::setTWeight(float tweight, float wtotsample,bool mult){
  float zerofact=1.0;
  //  cout << " weighted syst 0 is "<< weightedSysts[0]<<endl;
  if(mult)zerofact=this->weightedSysts[0];
  this->setSystValue("topPtWeightUp", zerofact*tweight/wtotsample);
  this->setSystValue("topPtWeightDown", zerofact/tweight*wtotsample);
}

void systWeights::setVHFWeight(int vhf,bool mult,double shiftval){
  float zerofact=1.0;
  double w_shift=0.0;
  //  cout << "vhf is "<<vhf<<endl;
  if (vhf>1)w_shift=shiftval;
  //  cout << " weighted syst 0 is "<< weightedSysts[0]<<endl;
  if(mult)zerofact=this->weightedSysts[0];
  this->setSystValue("VHFWeightUp", zerofact*(1+w_shift));
  this->setSystValue("VHFWeightDown", zerofact*(1-w_shift));
}


void systWeights::setQ2Weights(float q2up, float q2down, float wzero, bool mult){
  float zerofact=1.0;
  if(mult){
    zerofact=this->weightedSysts[0];
    //    cout <<  "zerofact "<< zerofact << endl;
  }
  //  cout <<  "zerofact "<< zerofact << " q2up weight "<< q2up/wzero << " tot to fill "<< zerofact*q2up/wzero<<endl;
  //  cout <<  "zerofact "<< zerofact << " q2down weight "<< q2down/wzero << " tot to fill "<< zerofact*q2down/wzero<<endl;
  this->setSystValue("q2Up", zerofact*q2up/wzero);
  this->setSystValue("q2Down", zerofact*q2down/wzero);
}

double systWeights::getPDFValue(int numPDF){
  if(!addPDF){ cout << "error! No PDF used, this will do nothing."<<endl;return 0.;}
  int MIN = this->maxSystsNonPDF;
  return (double)this->weightedSysts[numPDF+MIN];

}
void systWeights::setPDFValue(int numPDF, double w){
  if(!addPDF){ cout << "error! No PDF used, this will do nothing."<<endl;return;}
  int MIN = this->maxSystsNonPDF;
  this->weightedSysts[numPDF+MIN]=w;

}

void systWeights::calcPDFHisto(TH1F** histo, TH1F* singleHisto, double scalefactor, int c){//EXPERIMENTAL
  if(!addPDF){ cout << "error! No PDF used, this will do nothing."<<endl;return;}
  int MAX = this->maxSysts;
  //  for (int c = 0; c < this->nCategories; c++){
    int MIN = this->maxSystsNonPDF+(MAX+1)*c;
    for(int b = 0; b< singleHisto->GetNbinsX();++b){
      float val = singleHisto->GetBinContent(b);
      //      cout << "bin # "<<b << " val "<<val<<endl;
      float mean = 0, devst=0;
      //      cout << "name is "<< singleHisto->GetName()<<endl;
      for(int i = 0; i<this->nPDF;++i ){
	//cout<< " now looking at pdf # "<<i<< " coordinate is "<< MIN+i<<endl;
	//	cout << "is histo there? "<< histo[i+MIN]<<endl;
	//	cout << " histo should be "<< (histo[i+MIN])->GetName()<<endl;
	mean = mean+ histo[i+MIN]->GetBinContent(b);
      }
      mean = mean/this->nPDF;
      //mean = val;//mean/this->nPDF;
      for(int i = 0; i<this->nPDF;++i ){
	devst+=(mean-histo[i+MIN]->GetBinContent(b))*(mean-histo[i+MIN]->GetBinContent(b));
      }
      devst= sqrt(devst/this->nPDF);
      singleHisto->SetBinContent(b,val+devst*scalefactor);
      //      singleHisto->SetBinContent(b,mean+devst*scalefactor);
    }
    //}
}

