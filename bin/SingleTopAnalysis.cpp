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
#include "SystematicsUtilities.h"
//Extra variables
#include "../src/DMTopVariables.h"
#include "../src/MT2Utility.h"
#include "../src/mt2w_bisect.h"
#include "../src/mt2bl_bisect.h"
#include "../src/Mt2Com_bisect.h"
//#include "./EquationSolver.h"
//#include "./DMTopVariables.h"
//#include "./Weights.h"


using namespace std;

typedef vector<double> vdouble;
typedef vector<float> vfloat;
typedef vector<int> vint;
typedef vector<bool> vbool;
typedef vector<string> vstring;



int main(int argc, char **argv) {

    std::cout<<"Info: Starting to run ..."<<endl;
    std::cout<<"--------------------"<<endl; 
    std::cout<<"Usage:\t "<<argv[0]<<" "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<" "<<argv[6]<<" "<<argv[7]<<" "<<endl;
    

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

    string treesDir(argv[8]);
    if(treesDir != "noTrees")std::cout<<"\t treesDir: "<<treesDir<<endl;

    string addMVA(argv[9]);
    if(addMVA != "noMVA")std::cout<<"\t adding mvas: "<<addMVA<<endl;
    
    string outputpath(argv[10]);
    std::cout<<"\t Writing output to : "<<outputpath<<endl;

    bool doMVA=false;
    if(addMVA != "noMVA")doMVA=true;
    
    TString path_ = path ; 
    std::cout<<"\t File to open: "<<path_<<endl;
    std::cout<<"--------------------"<<endl; 
    
    std::cout << "Info: Loading file collection from " << path << std::endl;
    TFileCollection fc(sample.c_str(),sample.c_str(),path.c_str());
    std::cout << "Info: Files found : " << fc.GetNFiles() << std::endl;

    bool useCMVATSelection = true;
    bool useCSVTSelection = false;
    bool useCSVMSelection = false;
    //  cout << "cat is "<<endl;
    if(cat=="CSVT"){
        useCSVMSelection = false;
        useCSVTSelection = true;
        useCMVATSelection = false;
        }

    if(cat=="CMVAT"){
        useCSVMSelection = false;
        useCSVTSelection = false;
        useCMVATSelection = true;
        }

    if(cat=="CSVM"){
        useCSVMSelection = true;
        useCSVTSelection = false;
        useCMVATSelection = false;
        }
    
    string btagalgo="CMVAT";
    string btagvarname= "CMVAv2";
    string btagreshapename= "CMVA";
    if(useCSVTSelection){
      btagvarname= "CSVv2";
      btagreshapename= "CSV";
    }

   
    if(useCSVMSelection) btagalgo= "CSVM";
    if(useCSVTSelection) btagalgo= "CSVT";
    if(useCMVATSelection) btagalgo= "CMVAT";

    string reportName = "SelectedEvents_"+channel+"_"+cat+"_"+sample+".txt";
    string reportName_step0 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step0.txt";
    string reportName_step1 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step1.txt";
    string reportName_step2 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step2.txt";
    string reportName_step3 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step3.txt";
    string reportName_step4 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step4.txt";
    string reportName_step5 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step5.txt";
    string reportName_step6 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step6.txt";
    //  string reportName = "SelectedEvents_"+channel+""cat"_"+sample+".txt";
    ofstream fileout;
    ofstream fileout_step0;
    ofstream fileout_step1;
    ofstream fileout_step2;
    ofstream fileout_step3;
    ofstream fileout_step4;
    ofstream fileout_step5;
    ofstream fileout_step6;

    bool doSynch=sync.find("sync")!=std::string::npos;
    if(doSynch){
        fileout.open(reportName.c_str(),ios::in | ios::out | ios::trunc);
        fileout<<"RunNumber EvtNumber Lumi "<<std::endl;
        fileout_step0.open(reportName_step0.c_str(),ios::in | ios::out | ios::trunc);
        fileout_step1.open(reportName_step1.c_str(),ios::in | ios::out | ios::trunc);
        //fileout_step1<<"RunNumber | EvtNumber | Lumi | #TightMu | #LooseMu | #LooseEl | #Jets | #BJets | TightMu Pt | TightMuRelIso | AddLoose.Mu Pt | Add.LooseEl Pt | 1stJetPt | 2ndJetPt | 1stB Dis.| 2ndB Dis.| MET | MTW | "<<std::endl;
        fileout_step2.open(reportName_step2.c_str(),ios::in | ios::out | ios::trunc);
        fileout_step3.open(reportName_step3.c_str(),ios::in | ios::out | ios::trunc);
        fileout_step4.open(reportName_step4.c_str(),ios::in | ios::out | ios::trunc);
        fileout_step5.open(reportName_step5.c_str(),ios::in | ios::out | ios::trunc);
        fileout_step6.open(reportName_step6.c_str(),ios::in | ios::out | ios::trunc);
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
    
    bool addTrees=false;
    if (treesDir!="noTrees") addTrees=true;
    
    
    TFile *outTreeFile = new TFile;
    if(addTrees){
      outTreeFile = TFile::Open((treesDir+"/trees_"+sample+"_"+channel+".root").c_str(), "RECREATE");
    }
    
    
    Int_t nEvents = (Int_t)chain.GetEntries();
    std::cout<<"Info: Number of Events: "<<nEvents<< endl;
    //nEvents = std::min(nEvents, 1000);
    
    TString weightedSystsNames (weightedSysts sy);
    
    systWeights systZero,syst0BM,syst1BM,syst2BM; 
    int maxSysts=0; 
    int sizeMax=50;
    int muSize, elSize, jetSize,jet20Size, muLooseSize,elLooseSize,muLooseIsoGE0p15Size,muTightIsoGE0p15Size;
    //float passTrigHT(0.), Ht(0.);
    float Ht(0.), mt(0.);
    float runNumber(0.), lumiSec(0.);
    double evtNumber(0.);
    float jetE[sizeMax], jetPt[sizeMax], jetPhi[sizeMax], jetEta[sizeMax];
    float jetReshapeFactorCSV[sizeMax],jetReshapeFactorCSV_SD[sizeMax];
    //    float jetReshapeFactorCMVA[sizeMax],jetReshapeFactorCMVA_SD[sizeMax];
    float jetHadronFlavour[sizeMax],jetPartonFlavour[sizeMax];
    
    float jetJesUpReshapeFactorCSV[sizeMax],jetJesUpReshapeFactorCSV_SD[sizeMax], jetJesUpReshapeFactorCMVA[sizeMax],jetJesUpReshapeFactorCMVA_SD[sizeMax], jetJesUpHadronFlavour[sizeMax],jetJesUpPartonFlavour[sizeMax];

    float jetJesDownReshapeFactorCSV[sizeMax],jetJesDownReshapeFactorCSV_SD[sizeMax], jetJesDownReshapeFactorCMVA[sizeMax],jetJesDownReshapeFactorCMVA_SD[sizeMax], jetJesDownHadronFlavour[sizeMax],jetJesDownPartonFlavour[sizeMax];

    float jet20E[sizeMax], jet20Pt[sizeMax], jet20Phi[sizeMax], jet20Eta[sizeMax];


    int jetJesUpSize, jetJesDownSize, jet20JesUpSize, jet20JesDownSize;
    
    float jetJesUpE[sizeMax], jetJesUpPt[sizeMax], jetJesUpPhi[sizeMax], jetJesUpEta[sizeMax];
    float jetJesDownE[sizeMax], jetJesDownPt[sizeMax], jetJesDownPhi[sizeMax], jetJesDownEta[sizeMax];

    float jet20JesUpE[sizeMax], jet20JesUpPt[sizeMax], jet20JesUpPhi[sizeMax], jet20JesUpEta[sizeMax];
    float jet20JesDownE[sizeMax], jet20JesDownPt[sizeMax], jet20JesDownPhi[sizeMax], jet20JesDownEta[sizeMax];



    float jetIsCSVL[sizeMax], jetIsCSVM[sizeMax], jetIsCSVT[sizeMax],jetIsLoose[sizeMax],jetIsTight[sizeMax],jetak4chs_csvv2[sizeMax],jetJesUpak4chs_csvv2[sizeMax],jetJesDownak4chs_csvv2[sizeMax];

    float jetJesUpIsCSVT[sizeMax],jet20JesUpIsCSVT[sizeMax],jetJesDownIsCSVT[sizeMax],jet20JesDownIsCSVT[sizeMax];
    
    float jetIsCMVAL[sizeMax], jetIsCMVAM[sizeMax], jetIsCMVAT[sizeMax],jetak4chs_cmvav2[sizeMax];
    float jetPassesB[sizeMax], jet20PassesB[sizeMax], jetJesUpPassesB[sizeMax], jet20JesUpPassesB[sizeMax], jetJesDownPassesB[sizeMax], jet20JesDownPassesB[sizeMax];
    float jet20IsCSVL[sizeMax], jet20IsCSVM[sizeMax], jet20IsCSVT[sizeMax],jet20IsLoose[sizeMax],jet20IsTight[sizeMax],jet20ak4chs_csvv2[sizeMax];
    float jet20ReshapeFactorCSV[sizeMax],jet20ReshapeFactorCSV_SD[sizeMax];
 
    float jet20JesUpReshapeFactorCSV[sizeMax],jet20JesUpReshapeFactorCSV_SD[sizeMax],jet20JesUpak4chs_csvv2[sizeMax];
    float jet20JesDownReshapeFactorCSV[sizeMax],jet20JesDownReshapeFactorCSV_SD[sizeMax],jet20JesDownak4chs_csvv2[sizeMax];
  

    //    float jet20ReshapeFactorCMVA[sizeMax],jet20ReshapeFactorCMVA_SD[sizeMax];
    float jet20HadronFlavour[sizeMax],jet20PartonFlavour[sizeMax];
    float jet20JesUpHadronFlavour[sizeMax],jet20JesUpPartonFlavour[sizeMax];
    float jet20JesDownHadronFlavour[sizeMax],jet20JesDownPartonFlavour[sizeMax];



    float nTightMuons, nTightElectrons, nVetoElectrons, nLooseMuons, nJets,nCSVJets;//, nCSVLJets;
    float nGoodPV, nPV, numTrueInt, w_pu;
    float metPt[1],metPhi[1],metPx[1],metPy[1];
    float metZeroPt[1],metZeroPhi[1],metZeroPx[1],metZeroPy[1];
    float metJesUpPt[1],metJesUpPhi[1],metJesUpPx[1],metJesUpPy[1];
    float metJesDownPt[1],metJesDownPhi[1],metJesDownPx[1],metJesDownPy[1];
    
    //double met,metpx,metpy;
    float met,metpx,metpy;
    float muE[sizeMax], muPt[sizeMax], muEta[sizeMax], muIso[sizeMax], muIsTight[sizeMax],muIsLoose[sizeMax], muPhi[sizeMax],elPt[sizeMax],muLoosePhi[sizeMax],muLooseEta[sizeMax], muLoosePt[sizeMax],muLooseE[sizeMax], muLooseIso[sizeMax];
    float muLooseIsoGE0p15Iso[sizeMax],muLooseIsoGE0p15Pt[sizeMax],muLooseIsoGE0p15Eta[sizeMax],muLooseIsoGE0p15Phi[sizeMax],muLooseIsoGE0p15E[sizeMax];
    float muTightIsoGE0p15Iso[sizeMax],muTightIsoGE0p15Pt[sizeMax],muTightIsoGE0p15Eta[sizeMax],muTightIsoGE0p15Phi[sizeMax],muTightIsoGE0p15E[sizeMax];

    float muAntiIsoE[sizeMax], muAntiIsoPt[sizeMax], muAntiIsoEta[sizeMax], muAntiIsoIso[sizeMax], muAntiIsoIsTight[sizeMax],muAntiIsoIsLoose[sizeMax], muAntiIsoPhi[sizeMax];
    int   muAntiIsoSize;
    float muCharge[sizeMax], elCharge[sizeMax];
    float muAntiIsoCharge[sizeMax];
    float topcharge,topsize,antitopsize;

    float elEta[sizeMax], elIso[sizeMax], elIsTight[sizeMax], elIsVeto[sizeMax], elE[sizeMax], elPhi[sizeMax], elPassesDRmu[sizeMax];    
    float lep1Pt(0), lep1Eta(0), lep1Phi(0), lep1E(0), lep2Pt(0), lep2Eta(0), lep2Phi(0),  lep2E(0), lep1Flavour(0), lep1Charge(0), lep2Flavour(0), lep2Charge(0);
    int nPDF=102;

    float slTrigIsoMu20_v1(0.), slTrigIsoMu20_v2(0.), slTrigIsoMu20_v3(0.);
    float slTrigIsoMu22_v1(0.), slTrigIsoMu22_v2(0.), slTrigIsoMu22_v3(0.);
    float slTrigIsoMu24_v1(0.), slTrigIsoMu24_v2(0.), slTrigIsoMu24_v3(0.),slTrigIsoMu24_v4(0.);

    float slTrigIsoTkMu20_v1(0.), slTrigIsoTkMu20_v2(0.), slTrigIsoTkMu20_v3(0.);
    float slTrigIsoTkMu22_v1(0.), slTrigIsoTkMu22_v2(0.), slTrigIsoTkMu22_v3(0.);
    float slTrigIsoTkMu24_v1(0.), slTrigIsoTkMu24_v2(0.), slTrigIsoTkMu24_v3(0.),slTrigIsoTkMu24_v4(0.);
    float slHLTEle32_eta2p1_WPTight_Gsf_v1(0.),slHLTEle32_eta2p1_WPTight_Gsf_v2(0.),slHLTEle32_eta2p1_WPTight_Gsf_v3(0.),slHLTEle32_eta2p1_WPTight_Gsf_v4(0.),slHLTEle32_eta2p1_WPTight_Gsf_v5(0.);
    float slTrigEle_v1(0.), slTrigEle_v2(0.);
    
    bool  TrigIsoMu20 = false;
    bool  TrigIsoMu22 = false;
    bool  TrigIsoMu24 = false;
    bool  TrigGSFEl32 = false;

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
   
    float lepWeightBCDEF1Mu(1.);
    float lepWeightGH1Mu(1.);

    Mt2Com_bisect *Mt2cal = new Mt2Com_bisect();
    //  float bWeight0CSVT, bWeight0CSVTBTagUp,bWeight0CSVTBTagDown,bWeight0CSVTMisTagUp,bWeight0CSVTMisTagDown;
    //  float bWeight1CSVT, bWeight1CSVTBTagUp,bWeight1CSVTBTagDown,bWeight1CSVTMisTagUp,bWeight1CSVTMisTagDown;
    //  float bWeight2CSVT, bWeight2CSVTBTagUp,bWeight2CSVTBTagDown,bWeight2CSVTMisTagUp,bWeight2CSVTMisTagDown;
    
    bool doTtoSDDecay=false;
    bool reweightall=false;

    bool onlyNominal=false;
    systZero.setOnlyNominal(onlyNominal);

    
    bool addJES=true,addJER=false;
    
    // addJES=false for sync exe
    if(doSynch) addJES = false;
    if(isData=="DATA") addJES=false;
    bool addPDF=false,addQ2=false,addTopPt=false,addVHF=false,addTTSplit=false;
    //addQ2=true;addPDF=true;
    //addJES=false;
    if(sample=="TT"){
        addTTSplit=false;
        }

    if(sample.find("_sd")!=std::string::npos){
        doTtoSDDecay=true;
        }

    systZero.prepareDefault(true, addQ2, addPDF, addTopPt,addJES,addJER,addVHF,addTTSplit);
    //Here must add all systematics we want to put there so that the size of the syst vector is set

    //this maxSysts is to set the number of systematics we do want IMPORTANT FOR ALL INITIALIZATIONS
    maxSysts= systZero.maxSysts;
    TFile * allMyFiles[maxSysts];
    
    systZero.createFilesSysts(allMyFiles,outputpath+"/res/"+sample + "_" +channel);
    
    //systZero.createFilesSysts(allMyFiles,"./res/"+sample + "_" +channel);
    //addWZNLO=false, 
    //addPDF=true;
    //addQ2=true;
    //addTopPt=true;
    //addVHF=false;
    //addWZNLO=false;
    //if(isData=="DATA"){addPDF=false, addQ2=false;addVHF=false;addTTSplit=false;} 

      cout << w<<endl;
      if(isData=="MC"){

        chain.SetBranchAddress("Event_LHEWeight0", &w_zero); 
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
            //cout << "Event_LHEWeight"+pstr<<w_pdfs<<endl;
            }
        chain.SetBranchAddress("Event_T_Weight",&w_top);
      //    chain.SetBranchAddress("Event_T_Weight",&w_top);
        }
    
    chain.SetBranchAddress("Event_T_size",&topsize);
    chain.SetBranchAddress("Event_Tbar_size",&antitopsize);

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
    chain.SetBranchAddress("muonsLoose_Phi", muLoosePhi);
    chain.SetBranchAddress("muonsLoose_Eta", muLooseEta);
    chain.SetBranchAddress("muonsLoose_Pt", muLoosePt);
    chain.SetBranchAddress("muonsLoose_E",  muLooseE);
    chain.SetBranchAddress("muonsLoose_Iso04", muLooseIso);
    chain.SetBranchAddress("muonsLoose_Iso04_0p15_GE_Iso04", muLooseIsoGE0p15Iso);
    chain.SetBranchAddress("muonsLoose_Iso04_0p15_GE_size", &muLooseIsoGE0p15Size);
    chain.SetBranchAddress("muonsLoose_Iso04_0p15_GE_Pt", muLooseIsoGE0p15Pt);
    chain.SetBranchAddress("muonsLoose_Iso04_0p15_GE_Eta", muLooseIsoGE0p15Eta);
    chain.SetBranchAddress("muonsLoose_Iso04_0p15_GE_Phi", muLooseIsoGE0p15Phi);
    chain.SetBranchAddress("muonsLoose_Iso04_0p15_GE_E", muLooseIsoGE0p15E);
    
    chain.SetBranchAddress("muonsTight_Iso04_0p15_GE_Iso04", muTightIsoGE0p15Iso);
    chain.SetBranchAddress("muonsTight_Iso04_0p15_GE_size", &muTightIsoGE0p15Size);
    chain.SetBranchAddress("muonsTight_Iso04_0p15_GE_Pt",    muTightIsoGE0p15Pt);
    chain.SetBranchAddress("muonsTight_Iso04_0p15_GE_Eta",   muTightIsoGE0p15Eta);
    chain.SetBranchAddress("muonsTight_Iso04_0p15_GE_Phi",   muTightIsoGE0p15Phi);
    chain.SetBranchAddress("muonsTight_Iso04_0p15_GE_E",     muTightIsoGE0p15E);
 
    chain.SetBranchAddress("muonsTightAntiIso_E", muAntiIsoE);
    chain.SetBranchAddress("muonsTightAntiIso_Phi", muAntiIsoPhi);
    chain.SetBranchAddress("muonsTightAntiIso_Eta", muAntiIsoEta);
    chain.SetBranchAddress("muonsTightAntiIso_Pt", muAntiIsoPt);
    chain.SetBranchAddress("muonsTightAntiIso_Charge", muAntiIsoCharge);
    chain.SetBranchAddress("muonsTightAntiIso_Iso04", muAntiIsoIso);
    chain.SetBranchAddress("muonsTightAntiIso_IsTightMuon", muAntiIsoIsTight);
    chain.SetBranchAddress("muonsTightAntiIso_IsLooseMuon", muAntiIsoIsLoose);
    chain.SetBranchAddress("muonsTightAntiIso_size", &muAntiIsoSize);

    if(!doSynch){
      chain.SetBranchAddress("jetsAK4CHSTight_CorrE",      &jetE);      chain.SetBranchAddress("jetsAK4CHSTight_CorrPt",     &jetPt);
      //chain.SetBranchAddress("jetsAK4CHSTight_E",      &jetE);      chain.SetBranchAddress("jetsAK4CHSTight_Pt",     &jetPt);

    }
    else{
      chain.SetBranchAddress("jetsAK4CHSTight_E",      &jetE);      chain.SetBranchAddress("jetsAK4CHSTight_Pt",     &jetPt);
    }
    //    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt",     &jetPt);

    if(!doSynch){
      chain.SetBranchAddress("metFull_CorrT1Pt",metZeroPt);        chain.SetBranchAddress("metFull_CorrT1Phi",metZeroPhi);        chain.SetBranchAddress("metFull_CorrT1Px",metZeroPx);	chain.SetBranchAddress("metFull_CorrT1Py",metZeroPy);            
      //chain.SetBranchAddress("metFull_Pt",metZeroPt);      chain.SetBranchAddress("metFull_Phi",metZeroPhi);      chain.SetBranchAddress("metFull_Px",metZeroPx);      chain.SetBranchAddress("metFull_Py",metZeroPy);

    }   
    else{
      chain.SetBranchAddress("metFull_Pt",metZeroPt);      chain.SetBranchAddress("metFull_Phi",metZeroPhi);      chain.SetBranchAddress("metFull_Px",metZeroPx);      chain.SetBranchAddress("metFull_Py",metZeroPy);
      
    }
//    chain.SetBranchAddress("metFull_Pt",metZeroPt);
    //    chain.SetBranchAddress("metFull_Phi",metZeroPhi);
    //    chain.SetBranchAddress("metFull_Px",metZeroPx);
    //    chain.SetBranchAddress("metFull_Py",metZeroPy);
    
    vector<string>scenarios;
    scenarios.push_back("nominal");
    
 

   //    addJES=false;
    if(addJES){
      scenarios.push_back("jesUp");
      scenarios.push_back("jesDown");

      chain.SetBranchAddress("jetsAK4CHSTightJESUp_size",      &jetJesUpSize);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_size",      &jetJesDownSize);
      chain.SetBranchAddress("jetsAK4CHSTightJESUp_CorrE",      &jetJesUpE);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_CorrE",      &jetJesDownE);
      chain.SetBranchAddress("jetsAK4CHSTightJESUp_CorrPt",      &jetJesUpPt);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_CorrPt",      &jetJesDownPt);
      chain.SetBranchAddress("jetsAK4CHSTightJESUp_Eta",      &jetJesUpEta);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_Eta",      &jetJesDownEta);
      chain.SetBranchAddress("jetsAK4CHSTightJESUp_Phi",      &jetJesUpPhi);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_Phi",      &jetJesDownPhi);
      chain.SetBranchAddress("jetsAK4CHSTightJESUp_IsCSVT",      &jetJesUpIsCSVT);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_IsCSVT",      &jetJesDownIsCSVT);

      chain.SetBranchAddress(("jetsAK4CHSTightJESUp_Is"+btagalgo).c_str(),  &jetJesUpPassesB);
      chain.SetBranchAddress(("jetsAK4CHSTightJESUp_CorrPt_20_Is"+btagalgo).c_str(),  &jet20JesUpPassesB);

     chain.SetBranchAddress(("jetsAK4CHSTightJESDown_Is"+btagalgo).c_str(),  &jetJesDownPassesB);
     chain.SetBranchAddress(("jetsAK4CHSTightJESDown_CorrPt_20_Is"+btagalgo).c_str(),  &jet20JesDownPassesB);
     // // //
      chain.SetBranchAddress("jetsAK4CHSTightJESUp_CorrPt_20_size",      &jet20JesUpSize);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_CorrPt_20_size",      &jet20JesDownSize);
      chain.SetBranchAddress("jetsAK4CHSTightJESUp_CorrPt_20_CorrE",      &jet20JesUpE);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_CorrPt_20_CorrE",      &jet20JesDownE);
      chain.SetBranchAddress("jetsAK4CHSTightJESUp_CorrPt_20_CorrPt",      &jet20JesUpPt);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_CorrPt_20_CorrPt",      &jet20JesDownPt);
      chain.SetBranchAddress("jetsAK4CHSTightJESUp_CorrPt_20_Eta",      &jet20JesUpEta);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_CorrPt_20_Eta",      &jet20JesDownEta);
      chain.SetBranchAddress("jetsAK4CHSTightJESUp_CorrPt_20_Phi",      &jet20JesUpPhi);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_CorrPt_20_Phi",      &jet20JesDownPhi);
      chain.SetBranchAddress("jetsAK4CHSTightJESUp_CorrPt_20_IsCSVT",      &jet20JesUpIsCSVT);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_CorrPt_20_IsCSVT",      &jet20JesDownIsCSVT);


      chain.SetBranchAddress(("jetsAK4CHSTightJESDown_CorrPt_20_"+btagvarname).c_str(),  &jet20JesDownak4chs_csvv2);
      chain.SetBranchAddress(("jetsAK4CHSTightJESDown_CorrPt_20_reshapeFactor"+btagreshapename).c_str(),  &jet20JesDownReshapeFactorCSV);
      chain.SetBranchAddress(("jetsAK4CHSTightJESDown_CorrPt_20_reshapeFactor"+btagreshapename+"_SD").c_str(),  &jet20JesDownReshapeFactorCSV_SD);

      chain.SetBranchAddress(("jetsAK4CHSTightJESUp_CorrPt_20_"+btagvarname).c_str(),  &jet20JesUpak4chs_csvv2);
      chain.SetBranchAddress(("jetsAK4CHSTightJESUp_CorrPt_20_reshapeFactor"+btagreshapename).c_str(),  &jet20JesUpReshapeFactorCSV);
      chain.SetBranchAddress(("jetsAK4CHSTightJESUp_CorrPt_20_reshapeFactor"+btagreshapename+"_SD").c_str(),  &jet20JesUpReshapeFactorCSV_SD);

      chain.SetBranchAddress(("jetsAK4CHSTightJESDown_"+btagvarname).c_str(),  &jetJesDownak4chs_csvv2);
      chain.SetBranchAddress(("jetsAK4CHSTightJESDown_reshapeFactor"+btagreshapename).c_str(),  &jetJesDownReshapeFactorCSV);
      chain.SetBranchAddress(("jetsAK4CHSTightJESDown_reshapeFactor"+btagreshapename+"_SD").c_str(),  &jetJesDownReshapeFactorCSV_SD);

      chain.SetBranchAddress(("jetsAK4CHSTightJESUp_"+btagvarname).c_str(),  &jetJesUpak4chs_csvv2);
      chain.SetBranchAddress(("jetsAK4CHSTightJESUp_reshapeFactor"+btagreshapename).c_str(),  &jetJesUpReshapeFactorCSV);
      chain.SetBranchAddress(("jetsAK4CHSTightJESUp_reshapeFactor"+btagreshapename+"_SD").c_str(),  &jetJesUpReshapeFactorCSV_SD);


      chain.SetBranchAddress("jetsAK4CHSTightJESUp_HadronFlavour", &jetJesUpHadronFlavour);
      chain.SetBranchAddress("jetsAK4CHSTightJESUp_PartonFlavour",  &jetJesUpPartonFlavour);

      chain.SetBranchAddress("jetsAK4CHSTightJESDown_HadronFlavour", &jetJesDownHadronFlavour);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_PartonFlavour",  &jetJesDownPartonFlavour);

      chain.SetBranchAddress("jetsAK4CHSTightJESUp_CorrPt_20_HadronFlavour", &jetJesUpHadronFlavour);
      chain.SetBranchAddress("jetsAK4CHSTightJESUp_CorrPt_20_PartonFlavour",  &jetJesUpPartonFlavour);

      chain.SetBranchAddress("jetsAK4CHSTightJESDown_CorrPt_20_HadronFlavour", &jetJesDownHadronFlavour);
      chain.SetBranchAddress("jetsAK4CHSTightJESDown_CorrPt_20_PartonFlavour",  &jetJesDownPartonFlavour);


      chain.SetBranchAddress("metFull_CorrT1PtJESUp",metJesUpPt);
      chain.SetBranchAddress("metFull_CorrT1PhiJESUp",metJesUpPhi);
      chain.SetBranchAddress("metFull_CorrT1PxJESUp",metJesUpPx);
      chain.SetBranchAddress("metFull_CorrT1PyJESUp",metJesUpPy);

      chain.SetBranchAddress("metFull_CorrT1PtJESDown",metJesDownPt);
      chain.SetBranchAddress("metFull_CorrT1PhiJESDown",metJesDownPhi);
      chain.SetBranchAddress("metFull_CorrT1PxJESDown",metJesDownPx);
      chain.SetBranchAddress("metFull_CorrT1PyJESDown",metJesDownPy);

    }
    chain.SetBranchAddress("jetsAK4CHSTight_Phi",    &jetPhi);
    chain.SetBranchAddress("jetsAK4CHSTight_Eta",    &jetEta);
    chain.SetBranchAddress("jetsAK4CHSTight_size",   &jetSize);
    chain.SetBranchAddress("jetsAK4CHSTight_IsCSVL", &jetIsCSVL);
    chain.SetBranchAddress("jetsAK4CHSTight_IsCSVM", &jetIsCSVM);
    chain.SetBranchAddress("jetsAK4CHSTight_IsCSVT", &jetIsCSVT);
    chain.SetBranchAddress("jetsAK4CHSTight_IsLoose",&jetIsLoose);
    chain.SetBranchAddress("jetsAK4CHSTight_IsTight",&jetIsTight);
    //btagvarname= "CSVv2";
    //      btagreshapename= "CSV";

    chain.SetBranchAddress(("jetsAK4CHSTight_"+btagvarname).c_str(),  &jetak4chs_csvv2);
    chain.SetBranchAddress(("jetsAK4CHSTight_reshapeFactor"+btagreshapename).c_str(),  &jetReshapeFactorCSV);
    chain.SetBranchAddress(("jetsAK4CHSTight_reshapeFactor"+btagreshapename+"_SD").c_str(),  &jetReshapeFactorCSV_SD);
    //    chain.SetBranchAddress("jetsAK4CHSTight_reshapeFactorCMVA",  &jetReshapeFactorCMVA);
    //    chain.SetBranchAddress("jetsAK4CHSTight_reshapeFactorCMVA_SD",  &jetReshapeFactorCMVA_SD);
    
    chain.SetBranchAddress("jetsAK4CHSTight_HadronFlavour", &jetHadronFlavour);
    chain.SetBranchAddress("jetsAK4CHSTight_PartonFlavour",  &jetPartonFlavour);

    //  float jetReshapeFactorCSV[sizeMax],jetReshapeFactorCSV_SD[sizeMax];
     chain.SetBranchAddress("jetsAK4CHSTight_IsCMVAL", &jetIsCMVAL);
     chain.SetBranchAddress("jetsAK4CHSTight_IsCMVAM", &jetIsCMVAM);
     chain.SetBranchAddress("jetsAK4CHSTight_IsCMVAT", &jetIsCMVAT);
     //     chain.SetBranchAddress("jetsAK4CHSTight_CMVAv2",  &jetak4chs_cmvav2);

     chain.SetBranchAddress(("jetsAK4CHSTight_Is"+btagalgo).c_str(),  &jetPassesB);
     chain.SetBranchAddress(("jetsAK4CHSTight_CorrPt_20_Is"+btagalgo).c_str(),  &jet20PassesB);
        

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
    chain.SetBranchAddress(("jetsAK4CHSTight_CorrPt_20_"+btagvarname).c_str(),  &jet20ak4chs_csvv2);
    chain.SetBranchAddress(("jetsAK4CHSTight_CorrPt_20_reshapeFactor"+btagreshapename).c_str(),  &jet20ReshapeFactorCSV);
    chain.SetBranchAddress(("jetsAK4CHSTight_CorrPt_20_reshapeFactor"+btagreshapename+"_SD").c_str(),  &jet20ReshapeFactorCSV_SD);
    //    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_reshapeFactorCSV",  &jet20ReshapeFactorCSV);
    //    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_reshapeFactorCSV_SD",  &jet20ReshapeFactorCSV_SD);
    //    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_reshapeFactorCSV",  &jet20ReshapeFactorCMVA);
    //    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_reshapeFactorCSV_SD",  &jet20ReshapeFactorCMVA_SD);

    
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_HadronFlavour", &jet20HadronFlavour);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_PartonFlavour",  &jet20PartonFlavour);  

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
      
      
     
      chain.SetBranchAddress(("Event_bWeight2"+btagalgo+"Tight").c_str(), &bWeight2CSVM);
      chain.SetBranchAddress(("Event_bWeightMisTagUp2"+btagalgo+"Tight").c_str(), &bWeight2CSVMMisTagUp);
      chain.SetBranchAddress(("Event_bWeightMisTagDown2"+btagalgo+"Tight").c_str(), &bWeight2CSVMMisTagDown);
      chain.SetBranchAddress(("Event_bWeightBTagUp2"+btagalgo+"Tight").c_str(), &bWeight2CSVMBTagUp);
      chain.SetBranchAddress(("Event_bWeightBTagDown2"+btagalgo+"Tight").c_str(), &bWeight2CSVMBTagDown);
      
      chain.SetBranchAddress(("Event_bWeight0"+btagalgo+"Tight").c_str(), &bWeight0CSVM);
      
      chain.SetBranchAddress(("Event_bWeightMisTagUp0"+btagalgo+"Tight").c_str(), &bWeight0CSVMMisTagUp);
      chain.SetBranchAddress(("Event_bWeightMisTagDown0"+btagalgo+"Tight").c_str(), &bWeight0CSVMMisTagDown);
      chain.SetBranchAddress(("Event_bWeightBTagUp0"+btagalgo+"Tight").c_str(), &bWeight0CSVMBTagUp);
      chain.SetBranchAddress(("Event_bWeightBTagDown0"+btagalgo+"Tight").c_str(), &bWeight0CSVMBTagDown);
      
      chain.SetBranchAddress(("Event_bWeight1"+btagalgo+"Tight").c_str(), &bWeight1CSVM);
      chain.SetBranchAddress(("Event_bWeightMisTagUp1"+btagalgo+"Tight").c_str(), &bWeight1CSVMMisTagUp);
      chain.SetBranchAddress(("Event_bWeightMisTagDown1"+btagalgo+"Tight").c_str(), &bWeight1CSVMMisTagDown);
      chain.SetBranchAddress(("Event_bWeightBTagUp1"+btagalgo+"Tight").c_str(), &bWeight1CSVMBTagUp);
      chain.SetBranchAddress(("Event_bWeightBTagDown1"+btagalgo+"Tight").c_str(), &bWeight1CSVMBTagDown);
      
    }


    
    //Muon trigger
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

    chain.SetBranchAddress("Event_passesHLT_Ele32_eta2p1_WPTight_Gsf_v1OB", &slHLTEle32_eta2p1_WPTight_Gsf_v1);
    chain.SetBranchAddress("Event_passesHLT_Ele32_eta2p1_WPTight_Gsf_v2", &slHLTEle32_eta2p1_WPTight_Gsf_v2);
    chain.SetBranchAddress("Event_passesHLT_Ele32_eta2p1_WPTight_Gsf_v3", &slHLTEle32_eta2p1_WPTight_Gsf_v3);
    chain.SetBranchAddress("Event_passesHLT_Ele32_eta2p1_WPTight_Gsf_v4", &slHLTEle32_eta2p1_WPTight_Gsf_v4);
    chain.SetBranchAddress("Event_passesHLT_Ele32_eta2p1_WPTight_Gsf_v8", &slHLTEle32_eta2p1_WPTight_Gsf_v5);
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
    
    chain.SetBranchAddress("Event_passesHLT_IsoMu24_v4", &slTrigIsoMu24_v4);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu24_v4", &slTrigIsoTkMu24_v4);
        

  
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
    /**************                    MiniTree branching              ***************/
    /********************************************************************************/
    if(addTrees){
    outTreeFile->cd();
    }
    TTree * trees1T[9];
    TTree * trees2T[6];
    //    float * flo = new float(1.0); 
    //    float * mtw_2j1t = new float(-999.0); 
    
    //
    //Standard selection: 2j1t + mtw cut
    map< string, float * >namesToVars;
    
    //Standard selection: 2j1t
    float * w_2j1t = new float(-999.0); 
    float * etajprime_2j1t  = new float(-999.0); namesToVars["etajprime"] = etajprime_2j1t;
    float * jprimeflavour_2j1t  = new float(-999.0); namesToVars["jprimeflavour"] = jprimeflavour_2j1t;
    float * nextrajets_2j1t = new float(-999.0); namesToVars["nextrajets"] = nextrajets_2j1t;
    
    //b-jet variables
    float * bjetflavour_2j1t  = new float(-999.0); namesToVars["bjetflavour"] = bjetflavour_2j1t;
    float * bjeteta_2j1t  = new float(-999.0); namesToVars["bjeteta"] = bjeteta_2j1t;
    float * bjetpt_2j1t  = new float(-999.0); namesToVars["bjetpt"] = bjetpt_2j1t;
    
    //Top variables
    float * topMass_2j1t    = new float(-999.0); namesToVars["topMass"] = topMass_2j1t;
    float * topY_2j1t = new float(-999.0); namesToVars["topY"] = topY_2j1t;
    float * topPt_2j1t = new float(-999.0); namesToVars["topPt"] = topPt_2j1t;
    float * topMt_2j1t = new float(-999.0); namesToVars["topMt"] = topMt_2j1t;
    float * topEta_2j1t = new float(-999.0); namesToVars["topEta"] = topEta_2j1t;
    float * costhetapol_2j1t = new float(-999.0); namesToVars["costhetapol"] = costhetapol_2j1t;
    float * costhetael_2j1t = new float(-999.0); namesToVars["costhetael"] = costhetael_2j1t;
    
    //Mass of lepton +jet
    float * mlb_2j1t = new float(-999.0); namesToVars["mlb"] = mlb_2j1t;
    float * mljprime_2j1t = new float(-999.0); namesToVars["mljprime"] = mljprime_2j1t;
    float * mljextra_2j1t = new float(-999.0); namesToVars["mljextra"] = mljextra_2j1t;
     
    //global and met vvariables variables
    float * mindeltaphi_2j1t = new float(-999.0); namesToVars["mindeltaphi"] = mindeltaphi_2j1t;
    float * mindeltaphi20_2j1t = new float(-999.0); namesToVars["mindeltaphi20"] = mindeltaphi20_2j1t;
    float * mt2w_2j1t = new float(-999.0); namesToVars["mt2w"] = mt2w_2j1t;
    float * mtw_2j1t = new float(-999.0); namesToVars["mtw"] = mtw_2j1t;
    float * met_2j1t = new float(-999.0); namesToVars["MET"] = met_2j1t;

    
    //extra jet 
    float * leadingextrajeteta_2j1t = new float(-999.0); namesToVars["leadingextrajeteta"] = leadingextrajeteta_2j1t;
    float * leadingextrajetflavour_2j1t = new float(-999.0);namesToVars["leadingextrajetflavour"] = leadingextrajetflavour_2j1t;
    float * leadingextrajetpt_2j1t = new float(-999.0); namesToVars["leadingextrajetpt"] = leadingextrajetpt_2j1t;
    float * leadingextrajetcsv_2j1t = new float(-999.0); namesToVars["leadingextrajetcsv"] = leadingextrajetcsv_2j1t;
    float * leadingextrajetcsvweight_2j1t = new float(-999.0); namesToVars["leadingextrajetcsvweight"] = leadingextrajetcsvweight_2j1t;
    float * leadingextrajetcsvweight_sd_2j1t = new float(-999.0); namesToVars["leadingextrajetcsvweight"] = leadingextrajetcsvweight_2j1t;

    //secondary top quark
    float * topYExtra_2j1t = new float(-999.0); namesToVars["topYExtra"] = topYExtra_2j1t;
    float * topMassExtra_2j1t = new float(-999.0); namesToVars["topMassExtra"] = topMassExtra_2j1t;
    float * topPtExtra_2j1t = new float(-999.0); namesToVars["topPtExtra"] = topPtExtra_2j1t;
    float * topEtaExtra_2j1t = new float(-999.0); namesToVars["topEtaExtra"] = topEtaExtra_2j1t;
    float * topMtExtra_2j1t = new float(-999.0); namesToVars["topMtExtra"] = topMtExtra_2j1t;
    float * costhetapolExtra_2j1t = new float(-999.0); namesToVars["costhetapolExtra"] = costhetapolExtra_2j1t;
    float * costhetaelExtra_2j1t = new float(-999.0); namesToVars["costhetaelExtra"] = costhetaelExtra_2j1t; 

    //mvas
    float * bdt_st_vs_qcd_2j1t = new float(-999.0); 
    
    //MVAS:
    //ST vs other:
    float * bdt_st_vs_vj_2j1t_mtwcut = new float(-999.0); 
    float * bdt_st_vs_tt_2j1t_mtwcut = new float(-999.0); 
    
    //ST sd vs other:
    float * bdt_stsd_vs_st_2j1t_mtwcut = new float(-999.0); 
    float * bdt_stsd_vs_vj_2j1t_mtwcut = new float(-999.0); 
    
    //3j2t selection
    float * w_3j2t = new float(-999.0); 
    float * mtw_3j2t = new float(-999.0); 
    float * etajprime_3j2t = new float(-999.0); 
    float * topMassLeading_3j2t = new float(-999.0); 
    float * nextrajets_3j2t = new float(-999.0); 
    float * mt2w_3j2t = new float(-999.0); 
    
    float * mljprime_3j2t = new float(-999.0);
    float * mlbleading_3j2t = new float(-999.0);
    float * mlbsecond_3j2t = new float(-999.0);
    
    float * deltaEtabb_3j2t = new float(-999.0);
    float * topMtLeading_3j2t = new float(-999.0);
    float * topPtLeading_3j2t = new float(-999.0);
    float * topYLeading_3j2t = new float(-999.0);
    float * topEtaLeading_3j2t = new float(-999.0);
    float * costhetapolLeading_3j2t = new float(-999.0);
    float * costhetaelLeading_3j2t = new float(-999.0);
    
    float * topMassSecond_3j2t = new float(-999.0);
    float * topMtSecond_3j2t = new float(-999.0);
    float * topPtSecond_3j2t = new float(-999.0);
    float * topYSecond_3j2t = new float(-999.0);
    float * topEtaSecond_3j2t = new float(-999.0);
    float * costhetapolSecond_3j2t = new float(-999.0);
    float * costhetaelSecond_3j2t = new float(-999.0); 
    
    
    //3j1t selection
    float * w_3j1t = new float(-999.0);
    float * mtw_3j1t = new float(-999.0);
    float * mt2w_3j1t = new float(-999.0); 
    
    float * etajprime_3j1t = new float(-999.0);
    float * jprimeflavour_3j1t = new float(-999.0);
    
    float * mlb_3j1t = new float(-999.0);
    float * mljprime_3j1t = new float(-999.0);
    float * mljextra_3j1t = new float(-999.0);
    
    float * bjetflavour_3j1t = new float(-999.0);
    float * bjeteta_3j1t = new float(-999.0);
    float * bjetpt_3j1t = new float(-999.0); 
     
    float * topMass_3j1t = new float(-999.0);
    float * topPt_3j1t = new float(-999.0);
    float * topMt_3j1t = new float(-999.0);
    float * topY_3j1t = new float(-999.0);
    float * topEta_3j1t = new float(-999.0);
    float * costhetael_3j1t = new float(-999.0);
    float * costhetapol_3j1t = new float(-999.0);
    
    float * topMassExtra_3j1t = new float(-999.0);
    float * topMtExtra_3j1t = new float(-999.0);
    float * topPtExtra_3j1t = new float(-999.0);
    float * topYExtra_3j1t = new float(-999.0);
    float * topEtaExtra_3j1t = new float(-999.0);
    float * costhetapolExtra_3j1t = new float(-999.0);
    float * costhetaelExtra_3j1t = new float(-999.0); 
    
     //extra jet 
     float * leadingextrajeteta_3j1t = new float(-999.0);
     float * leadingextrajetflavour_3j1t = new float(-999.0);
     float * leadingextrajetpt_3j1t = new float(-999.0); 
     float * leadingextrajetcsv_3j1t = new float(-999.0);
     float * leadingextrajetcsvweight_3j1t = new float(-999.0);
     float * leadingextrajetcsvweight_sd_3j1t = new float(-999.0);
     
    syst0BM.copySysts(systZero);
    syst1BM.copySysts(systZero);
    syst2BM.copySysts(systZero);
    if(addTrees){
      cout << "adding trees "<<endl;
      syst1BM.addSelection("2j1t");
      syst1BM.addSelection("3j1t");
      syst1BM.initTreesSysts(trees1T, outTreeFile);
      
      syst2BM.addSelection("3j2t");
      syst2BM.initTreesSysts(trees2T, outTreeFile);
      cout << " bef branch trees" <<endl;
      
      //without mtw cut
      
      syst1BM.branchTreesSysts(trees1T,"2j1t","w", outTreeFile, w_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","bdt_st_vs_qcd", outTreeFile, bdt_st_vs_qcd_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","etajprime", outTreeFile, etajprime_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","jprimeflavour", outTreeFile, jprimeflavour_2j1t);
      
      syst1BM.branchTreesSysts(trees1T,"2j1t","bjetflavour", outTreeFile, bjetflavour_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","bjeteta", outTreeFile, bjeteta_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","bjetpt", outTreeFile, bjetpt_2j1t);
      
      syst1BM.branchTreesSysts(trees1T,"2j1t","topMass",outTreeFile, topMass_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","topY",   outTreeFile, topY_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","topPt",  outTreeFile, topPt_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","topMt",  outTreeFile, topMt_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","topEta", outTreeFile, topEta_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","costhetapol", outTreeFile, costhetapol_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","costhetael", outTreeFile, costhetael_2j1t);
      
      syst1BM.branchTreesSysts(trees1T,"2j1t","mlb", outTreeFile, mlb_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","mljprime", outTreeFile, mljprime_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","mljextra", outTreeFile, mljextra_2j1t);
      
      syst1BM.branchTreesSysts(trees1T,"2j1t","mt2w", outTreeFile, mt2w_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","mindeltaphi", outTreeFile, mindeltaphi_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","mindeltaphi20", outTreeFile, mindeltaphi20_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","mtw", outTreeFile, mtw_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","MET", outTreeFile, met_2j1t);
      
      syst1BM.branchTreesSysts(trees1T,"2j1t","topMassExtra", outTreeFile, topMassExtra_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","topMtExtra", outTreeFile, topMtExtra_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","topYExtra", outTreeFile, topYExtra_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","topPtExtra", outTreeFile, topPtExtra_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","topEtaExtra", outTreeFile, topEtaExtra_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","costhetapolExtra", outTreeFile, costhetapolExtra_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","costhetaelExtra", outTreeFile, costhetaelExtra_2j1t);
      
      syst1BM.branchTreesSysts(trees1T,"2j1t","nextrajets", outTreeFile, nextrajets_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","leadingextrajetcsv", outTreeFile, leadingextrajetcsv_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","leadingextrajetflavour", outTreeFile, leadingextrajetflavour_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","leadingextrajetpt", outTreeFile, leadingextrajetpt_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","leadingextrajeteta", outTreeFile, leadingextrajeteta_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","leadingextrajetcsvweight", outTreeFile, leadingextrajetcsvweight_2j1t);
      syst1BM.branchTreesSysts(trees1T,"2j1t","leadingextrajetcsvweight_sd", outTreeFile, leadingextrajetcsvweight_sd_2j1t);
      
      
      //3j1t
      syst1BM.branchTreesSysts(trees1T,"3j1t","w", outTreeFile, w_3j1t);
      
      syst1BM.branchTreesSysts(trees1T,"3j1t","mtw", outTreeFile, mtw_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","mt2w", outTreeFile, mt2w_3j1t);
      
      syst1BM.branchTreesSysts(trees1T,"3j1t","etajprime", outTreeFile, etajprime_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","jprimeflavour", outTreeFile, jprimeflavour_3j1t);
      
      syst1BM.branchTreesSysts(trees1T,"3j1t","mlb", outTreeFile, mlb_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","mljprime", outTreeFile, mljprime_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","mljextra", outTreeFile, mljextra_3j1t);
      
      syst1BM.branchTreesSysts(trees1T,"3j1t","bjetflavour", outTreeFile, bjetflavour_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","bjeteta", outTreeFile, bjeteta_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","bjetpt", outTreeFile, bjetpt_3j1t);
      
      syst1BM.branchTreesSysts(trees1T,"3j1t","topMass",outTreeFile, topMass_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","topY",   outTreeFile, topY_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","topPt",  outTreeFile, topPt_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","topMt",  outTreeFile, topMt_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","topEta", outTreeFile, topEta_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","costhetael", outTreeFile, costhetael_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","costhetapol", outTreeFile, costhetapol_3j1t);
      
      syst1BM.branchTreesSysts(trees1T,"3j1t","topMassExtra",outTreeFile, topMassExtra_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","topYExtra",   outTreeFile, topYExtra_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","topPtExtra",  outTreeFile, topPtExtra_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","topMtExtra",  outTreeFile, topMtExtra_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","topEtaExtra", outTreeFile, topEtaExtra_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","costhetaelExtra", outTreeFile, costhetaelExtra_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","costhetapolExtra", outTreeFile, costhetapolExtra_3j1t);
      
      syst1BM.branchTreesSysts(trees1T,"3j1t","leadingextrajetcsv", outTreeFile, leadingextrajetcsv_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","leadingextrajetpt", outTreeFile, leadingextrajetpt_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","leadingextrajeteta", outTreeFile, leadingextrajeteta_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","leadingextrajetcsvweight", outTreeFile, leadingextrajetcsvweight_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","leadingextrajetcsvweight_sd", outTreeFile, leadingextrajetcsvweight_sd_3j1t);
      syst1BM.branchTreesSysts(trees1T,"3j1t","leadingextrajetflavour", outTreeFile, leadingextrajetflavour_3j1t);

      
      syst2BM.branchTreesSysts(trees2T,"3j2t","w", outTreeFile, w_3j2t);

      syst2BM.branchTreesSysts(trees2T,"3j2t","mtw", outTreeFile, mtw_3j2t);
      syst2BM.branchTreesSysts(trees2T,"3j2t","etajprime", outTreeFile, etajprime_3j2t);
      syst2BM.branchTreesSysts(trees2T,"3j2t","topMassLeading", outTreeFile, topMassLeading_3j2t);
      syst2BM.branchTreesSysts(trees2T,"3j2t","topMassSecond", outTreeFile, topMassSecond_3j2t);
      syst2BM.branchTreesSysts(trees2T,"3j2t","mt2w", outTreeFile, mt2w_3j2t);


       syst2BM.branchTreesSysts(trees2T,"3j2t","deltaEtabb", outTreeFile, deltaEtabb_3j2t);
 
       syst2BM.branchTreesSysts(trees2T,"3j2t","mlbleading", outTreeFile, mlbleading_3j2t);
       syst2BM.branchTreesSysts(trees2T,"3j2t","mlbsecond", outTreeFile, mlbsecond_3j2t);
       syst2BM.branchTreesSysts(trees2T,"3j2t","mljprime", outTreeFile, mljprime_3j2t);
 
       syst2BM.branchTreesSysts(trees2T,"3j2t","topYLeading",   outTreeFile, topYLeading_3j2t);
       syst2BM.branchTreesSysts(trees2T,"3j2t","topPtLeading",  outTreeFile, topPtLeading_3j2t);
       syst2BM.branchTreesSysts(trees2T,"3j2t","topMtLeading",  outTreeFile, topMtLeading_3j2t);
       syst2BM.branchTreesSysts(trees2T,"3j2t","topEtaLeading", outTreeFile, topEtaLeading_3j2t);
       syst2BM.branchTreesSysts(trees2T,"3j2t","costhetaelLeading", outTreeFile, costhetaelLeading_3j2t);
       syst2BM.branchTreesSysts(trees2T,"3j2t","costhetapolLeading", outTreeFile, costhetapolLeading_3j2t);
 
       syst2BM.branchTreesSysts(trees2T,"3j2t","topYSecond",   outTreeFile, topYSecond_3j2t);
       syst2BM.branchTreesSysts(trees2T,"3j2t","topPtSecond",  outTreeFile, topPtSecond_3j2t);
       syst2BM.branchTreesSysts(trees2T,"3j2t","topMtSecond",  outTreeFile, topMtSecond_3j2t);
       syst2BM.branchTreesSysts(trees2T,"3j2t","topEtaSecond", outTreeFile, topEtaSecond_3j2t);
       syst2BM.branchTreesSysts(trees2T,"3j2t","costhetaelSecond", outTreeFile, costhetaelSecond_3j2t);
       syst2BM.branchTreesSysts(trees2T,"3j2t","costhetapolSecond", outTreeFile, costhetapolSecond_3j2t);

    }

    //doMVA=true;
    // MVA input variables
    string STsd_ST_name="STsd_vs_ST", STsd_VJ_name="STsd_vs_VJ",STsd_TT_name= "STsd_vs_TT";
    string ST_VJ_name="ST_vs_VJ", ST_TT_name="ST_vs_TT";
    string ST_QCDMu_name="ST_vs_QCDMu";

    TMVA::Reader * STsd_VJ = new TMVA::Reader(STsd_VJ_name);//STsd_VJ(STsd_VJ_name),STsd_TT(STsd_TT_name);//Single top Vtsd discriminators 
    TMVA::Reader * STsd_ST = new TMVA::Reader(STsd_ST_name);//STsd_VJ(STsd_VJ_name),STsd_TT(STsd_TT_name);//Single top Vtsd discriminators 
    TMVA::Reader * ST_TT = new TMVA::Reader(ST_TT_name);
    TMVA::Reader * ST_VJ = new TMVA::Reader(ST_VJ_name);//Single top Vtb
    TMVA::Reader * ST_QCDMu = new TMVA::Reader(ST_QCDMu_name);//QCD
    vector<TMVA::Reader *> AllReaders;
    vector<string> AllReaderNames;
    if(doMVA){
      AllReaders.push_back(ST_VJ); AllReaderNames.push_back(ST_VJ_name);
      AllReaders.push_back(ST_TT); AllReaderNames.push_back(ST_TT_name);
      AllReaders.push_back(STsd_ST); AllReaderNames.push_back(STsd_ST_name);
      AllReaders.push_back(STsd_VJ); AllReaderNames.push_back(STsd_VJ_name);
      AllReaders.push_back(ST_QCDMu); AllReaderNames.push_back(ST_QCDMu_name);
    }
    if(doMVA){
      //      STsd_ST.AddSpectator("etajprime",&mvavars1);
      for(size_t rd = 0; rd<AllReaders.size();++rd){
	
	ifstream varCfgFile;
	varCfgFile.open("MVA/cfg"+AllReaderNames.at(rd)+".txt",ifstream::in);
	string bas,option;
	while(std::getline(varCfgFile, bas)) {
	  if( bas.find("#") != std::string::npos) {
	    if(bas.find("#") == 0) continue;
	    bas = bas.substr(0, bas.find("#"));
	  }
	  if( bas.empty() ) continue;
	  if( bas.substr(0,1) == "[" && bas.substr(bas.length() -1 , bas.length()) == "]" ) {
	    option = bas; 
	    cout << "Option " << option << " is initialized" << endl;
	    continue; 
	  }
	  if(option == "[variables]" ) {
	    AllReaders.at(rd)->AddVariable(bas,namesToVars[bas]);
	    cout << "variable: " << bas << endl;
	    continue;
	  }
	}
	//TString weightsig="",weightbkg=""; 
	//STsd_ST.AddVariable("etajprime", &etajprime_2j1t_mtwcut); 
	//STsd_ST.AddVariable("etajprime", &etajprime_2j1t_mtwcut); 
	const std::string weightsfile = std::string("MVA/weights/TMVAMu"+(AllReaderNames.at(rd))+"_BDT.weights.xml");
	AllReaders.at(rd)->BookMVA("BDTG", weightsfile.c_str()); 
      }
      
      //      const std::string cmssw_base = getenv("CMSSW_BASE");
      //      const std::string weightsfile = cmssw_base + std::string("/src/ttDM/TopTagResolved/data/toptrainingbits_prob.root_BDTG.weights.xml");
    }

    
    /********************************************************************************/
    /**************                    Histogram booking              ***************/
    /********************************************************************************/
    TH1F *h_cutFlow = new TH1F("h_cutFlow","cutflow",10,-0.5,9.5);
    TH1F *h_weight_sign = new TH1F("h_weight_sign","weight correction factor before all selections, gives the effective number of events",2,-2.0,2.0);
    chainNEvents.Project("h_weight_sign","Event_LHEWeight10/abs(Event_LHEWeight10)"); // I messed up and weight zero is not  saved in the full chain... but this one should do the trick as well! Just do h_weight_sign->GetMean() and multiply this to get the correct number :)
    double zeroWeight=h_weight_sign->GetMean(); // Will have to add a zeroWeight vector in the systematics for a proper treatment of the acceptance in systematics. For now we do it via hardcoding outside.
    cout << " mean  "<< zeroWeight<<endl;
    
  
    TH1F *h_weightZero[maxSysts]; systZero.initHistogramsSysts(h_weightZero,"h_weightZero","weight before all selections, normalized by the total number of events",2000,0,10.0);
    TH1F *h_nPV[maxSysts];        systZero.initHistogramsSysts(h_nPV,"h_nPV","nPV",34,0,74);
    TH1F *h_nGoodPV[maxSysts];    systZero.initHistogramsSysts(h_nGoodPV,"h_nGoodPV","GoodPV",34,0,74);
    TH1F *h_nTruePV[maxSysts];    systZero.initHistogramsSysts(h_nTruePV,"h_nTruePV","TruePV",34,0,74);
    TH1F *h_nJets[maxSysts];      systZero.initHistogramsSysts(h_nJets,"h_nJets","Number of tight jets",13,-0.5,12.5);
    TH1F *h_nbJets[maxSysts];     systZero.initHistogramsSysts(h_nbJets,"h_nbJets","Number of tight b-jets",11,-0.5,10.5);

    //2j0t 
    TH1F *h_2j0t_dR_lepjetpt40_1st[maxSysts];      systZero.initHistogramsSysts(h_2j0t_dR_lepjetpt40_1st,   "h_2j0t_dR_lepjetpt40_1st",  "dR lep-jet ",50,0,10);
    TH1F *h_2j0t_dPhi_lepjetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j0t_dPhi_lepjetpt40_1st,   "h_2j0t_dPhi_lepjetpt40_1st",  "dPhi lep-jet ",100,0.0,3.2);
    TH1F *h_2j0t_dEta_lepjetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j0t_dEta_lepjetpt40_1st,   "h_2j0t_dEta_lepjetpt40_1st",  "dPhi lep-jet ",50,-3.2,3.2);
    TH1F *h_2j0t_MuEta[maxSysts];          systZero.initHistogramsSysts(h_2j0t_MuEta,       "h_2j0t_MuEta",      "muon eta ",100,-2.4,2.4); 
    TH1F *h_2j0t_MuPhi[maxSysts];          systZero.initHistogramsSysts(h_2j0t_MuPhi,       "h_2j0t_MuPhi",      "muon phi ",100,-3.2,3.2); 
    TH1F *h_2j0t_jetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j0t_jetpt40_1st,   "h_2j0t_jetpt40_leading","2j0t leading jet Pt ",250,0,500);
    TH1F *h_2j0t_jetpt40_2nd[maxSysts];    systZero.initHistogramsSysts(h_2j0t_jetpt40_2nd,   "h_2j0t_jetpt40_subleading","2j0t Sub. leading jet Pt ",250,0,500);
    TH1F *h_2j0t_mtw[maxSysts];            systZero.initHistogramsSysts(h_2j0t_mtw,           "h_2j0t_mtw",     "2j0t mtw ",250,0,500);
    TH1F *h_2j0t_muIso[maxSysts];          systZero.initHistogramsSysts(h_2j0t_muIso,"h_2j0t_muIso","Muon Isolation",40,0,1);
    TH1F *h_2j0t_MuPt[maxSysts];           systZero.initHistogramsSysts(h_2j0t_MuPt,    "h_2j0t_MuPt",    "2j0t muon pt",500,0,500);

    //After mtw>50 cut
    TH1F *h_2j0t_mtwcut_jetpt40_1st[maxSysts];   systZero.initHistogramsSysts(h_2j0t_mtwcut_jetpt40_1st,   "h_2j0t_mtwcut_jetpt40_leading","2j0t leading jet Pt ",250,0,500);
    TH1F *h_2j0t_mtwcut_jetpt40_2nd[maxSysts];   systZero.initHistogramsSysts(h_2j0t_mtwcut_jetpt40_2nd,   "h_2j0t_mtwcut_jetpt40_subleading","2j0t Sub. leading jet Pt ",250,0,500);
    TH1F *h_2j0t_mtwcut_jeteta40_1st[maxSysts];  systZero.initHistogramsSysts(h_2j0t_mtwcut_jeteta40_1st,   "h_2j0t_mtwcut_jeteta40_leading","2j0t leading jet eta ",200,-4.7,4.7);
    TH1F *h_2j0t_mtwcut_incljeteta[maxSysts];    systZero.initHistogramsSysts(h_2j0t_mtwcut_incljeteta,   "h_2j0t_mtwcut_incljeteta","2j0t leading inclusive jet eta ",200,-4.7,4.7);

    TH1F *h_2j0t_mtwcut_jetabseta40_1st[maxSysts]; systZero.initHistogramsSysts(h_2j0t_mtwcut_jetabseta40_1st,   "h_2j0t_mtwcut_jetabseta40_leading","2j0t leading jet eta ",200,0,4.7);
    TH1F *h_2j0t_mtwcut_incljetabseta[maxSysts];  systZero.initHistogramsSysts(h_2j0t_mtwcut_incljetabseta,   "h_2j0t_mtwcut_incljetabseta","2j0t leading inclusive jet eta ",200,0,4.7);
    TH1F *h_2j0t_incljetabseta[maxSysts];          systZero.initHistogramsSysts(h_2j0t_incljetabseta,   "h_2j0t_incljetabseta","2j0t leading inclusive jet eta ",200,0,4.7);
    TH1F *h_2j0t_mtwcut_MuPt[maxSysts];           systZero.initHistogramsSysts(h_2j0t_mtwcut_MuPt,    "h_2j0t_mtwcut_MuPt",    "2j0t muon pt with mtwcut",500,0,500);
    TH1F *h_2j0t_mtwcut_muIso[maxSysts];          systZero.initHistogramsSysts(h_2j0t_mtwcut_muIso,"h_2j0t_mtwcut_muIso","Muon Isolation with mtwcut",100,0,1);

  
  //2j1t Category
  TH1F *h_2j1t_bjetpt[maxSysts];    systZero.initHistogramsSysts(h_2j1t_bjetpt,  "h_2j1t_bjetpt",  "2j1t b jet pt",250,0,500);
  TH1F *h_2j1t_ljetpt[maxSysts];    systZero.initHistogramsSysts(h_2j1t_ljetpt,  "h_2j1t_ljetpt",   "2j1t light jet pt" ,250,0,500);
  TH1F *h_2j1t_ljeteta[maxSysts];   systZero.initHistogramsSysts(h_2j1t_ljeteta, "h_2j1t_ljeteta",  "2j1t light jet eta",200,-4.7,4.7);
  TH1F *h_2j1t_nMu[maxSysts];       systZero.initHistogramsSysts(h_2j1t_nMu,     "h_2j1t_nMu",     "2j1t Number of tight Muons",13,-0.5,12.5);
  TH1F *h_2j1t_MuPt[maxSysts];      systZero.initHistogramsSysts(h_2j1t_MuPt,    "h_2j1t_MuPt",    "2j1t muon pt ",500,0,500);
  TH1F *h_2j1t_MuEta[maxSysts];     systZero.initHistogramsSysts(h_2j1t_MuEta,   "h_2j1t_MuEta",   "2j1t muon eta ",100,-2.4,2.4);
  TH1F *h_2j1t_MuPhi[maxSysts];     systZero.initHistogramsSysts(h_2j1t_MuPhi,   "h_2j1t_MuPhi",   "2j1t muon phi ",100, -3.2, 3.2);
  TH1F *h_2j1t_MuE[maxSysts];       systZero.initHistogramsSysts(h_2j1t_MuE,     "h_2j1t_MuE",     "2j1t muon e ",250,0,500);
  TH1F *h_2j1t_MuCharge[maxSysts];  systZero.initHistogramsSysts(h_2j1t_MuCharge,"h_2j1t_MuCharge","2j1t muon charge ",2,-1,1);
  TH1F *h_2j1t_mtw[maxSysts];       systZero.initHistogramsSysts(h_2j1t_mtw,     "h_2j1t_mtw",     "2j1t mtw ",250,0,500);
  TH1F *h_2j1t_topMass[maxSysts];   systZero.initHistogramsSysts(h_2j1t_topMass, "h_2j1t_topMass", "2j1t top mass",200,100,500);
  TH1F *h_2j1t_topPt[maxSysts];     systZero.initHistogramsSysts(h_2j1t_topPt,   "h_2j1t_topPt",   "2j1t top pt",250,0,500);
  TH1F *h_2j1t_etajprime[maxSysts];   systZero.initHistogramsSysts(h_2j1t_etajprime,  "h_2j1t_etajprime",  "2j1t light jet abs(eta)",100,0.,4.7);

   TH1F *h_2j1t_dR_lepjetpt40_1st[maxSysts];      systZero.initHistogramsSysts(h_2j1t_dR_lepjetpt40_1st,   "h_2j1t_dR_lepjetpt40_1st",  "dR lep-jet ",50,0,10);
   TH1F *h_2j1t_dPhi_lepjetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j1t_dPhi_lepjetpt40_1st, "h_2j1t_dPhi_lepjetpt40_1st",  "dPhi lep-jet ",100,0.0,3.2);
   TH1F *h_2j1t_dEta_lepjetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j1t_dEta_lepjetpt40_1st, "h_2j1t_dEta_lepjetpt40_1st",  "dPhi lep-jet ",50,-3.2,3.2);
   TH1F *h_2j1t_muIso[maxSysts];          systZero.initHistogramsSysts(h_2j1t_muIso,"h_2j1t_muIso","Muon Isolation",100,0,1);

  TH1F *h_2j1t_jetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j1t_jetpt40_1st,   "h_2j1t_jetpt40_1st","2j1t leading jet Pt ",250,0,500);
  TH1F *h_2j1t_jetpt40_2nd[maxSysts];    systZero.initHistogramsSysts(h_2j1t_jetpt40_2nd,   "h_2j1t_jetpt40_2nd","2j1t 2nd leading jet Pt ",250,0,500);
  
  //Jets pt 20 spectra.
  TH1F *h_2j1t_nextrajets[maxSysts];   systZero.initHistogramsSysts(h_2j1t_nextrajets,  "h_2j1t_nextrajets",  "2j1t number of jets with 20<pt<40",5,-0.5,4.5);
  TH1F *h_2j1t_leadingextrajetpt[maxSysts];   systZero.initHistogramsSysts(h_2j1t_leadingextrajetpt,  "h_2j1t_leadingextrajetpt",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_leadingextrajetcsv[maxSysts];   systZero.initHistogramsSysts(h_2j1t_leadingextrajetcsv,  "h_2j1t_leadingextrajetcsv",  "2j1t leading extra jet csv value",100,0.,1.0);
  TH1F *h_2j1t_leadingextrajetpt_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_leadingextrajetpt_sd_b,  "h_2j1t_leadingextrajetpt_sd_b",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_leadingextrajetcsv_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_leadingextrajetcsv_sd_b,  "h_2j1t_leadingextrajetcsv_sd_b",  "2j1t leading extra jet csv value for b from top decays",100,0.,1.0);
  
  TH1F *h_2j1t_leadingextrajetcsv_reshape_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_leadingextrajetcsv_reshape_sd_b,  "h_2j1t_leadingextrajetcsv_reshape_sd_b",  "2j1t leading extra jet csv reshaping from b to sd for b from top decays",1000,0.,10.0);


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

  TH1F *h_2j1t_mtwcut_leadingextrajetpt_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt_sd_b,  "h_2j1t_mtwcut_leadingextrajetpt_sd_b",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_mtwcut_leadingextrajetcsv_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv_sd_b,  "h_2j1t_mtwcut_leadingextrajetcsv_sd_b",  "2j1t csv of the leading extra jet  for b from top decays",100,0.,1.0);
  TH1F *h_2j1t_mtwcut_topMassExtra_sd_b[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_topMassExtra_sd_b, "h_2j1t_mtwcut_topMassExtra_sd_b", "2j1t top mass with the first extra jet as b, for b-quarks from top",200,100,500);

  TH1F *h_2j1t_mtwcut_leadingextrajetcsv_reshape_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv_reshape_sd_b,  "h_2j1t_mtwcut_leadingextrajetcsv_reshape_sd_b",  "2j1t leading extra jet csv reshaping from b to sd for b from top decays",1000,0.,10.0);
   TH1F *h_2j1t_mtwcut_MuPt[maxSysts];     systZero.initHistogramsSysts(h_2j1t_mtwcut_MuPt,    "h_2j1t_mtwcut_MuPt",    "2j1t muon pt with mtwcut",500,0,500);
   TH1F *h_2j1t_mtwcut_muIso[maxSysts];          systZero.initHistogramsSysts(h_2j1t_mtwcut_muIso,"h_2j1t_mtwcut_muIso","Muon Isolation with mtwcut",100,0,1);

   //mtwcut + eta < 2.5
   TH1F *h_2j1t_mtwcut_le_sr_BDT_STvsVJ[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_le_sr_BDT_STvsVJ, "h_2j1t_mtwcut_le_sr_BDT_STvsVJ", "BDT ST-vs-VJ 2j1t with mtwcut_le_sr", 50, -1, 1);
   TH1F *h_2j1t_mtwcut_le_sr_BDT_STvsTT[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_le_sr_BDT_STvsTT, "h_2j1t_mtwcut_le_sr_BDT_STvsTT", "BDT ST-vs-VJ 2j1t with mtwcut_le_sr", 50, -1, 1); 
   TH1F *h_2j1t_mtwcut_le_sr_BDT_STsdvsVJ[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_le_sr_BDT_STsdvsVJ, "h_2j1t_mtwcut_le_sr_BDT_STsdvsVJ","BDT STsd-vs-VJ 2j1t with mtwcut_le_sr", 50, -1, 1);
   TH1F *h_2j1t_mtwcut_le_sr_BDT_STsdvsST[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_le_sr_BDT_STsdvsST, "h_2j1t_mtwcut_le_sr_BDT_STsdvsST","BDT STsd-vs-VJ 2j1t with mtwcut_le_sr", 50, -1, 1);
 
   //mtwcut + eta > 2.5
   TH1F *h_2j1t_mtwcut_sr_BDT_STvsVJ[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_BDT_STvsVJ, "h_2j1t_mtwcut_sr_BDT_STvsVJ", "BDT ST-vs-VJ 2j1t with mtwcut_ge_sr", 50, -1, 1);
   TH1F *h_2j1t_mtwcut_sr_BDT_STvsTT[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_BDT_STvsTT, "h_2j1t_mtwcut_sr_BDT_STvsTT", "BDT ST-vs-VJ 2j1t with mtwcut_ge_sr", 50, -1, 1); 
   TH1F *h_2j1t_mtwcut_sr_BDT_STsdvsVJ[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_BDT_STsdvsVJ, "h_2j1t_mtwcut_sr_BDT_STsdvsVJ","BDT STsd-vs-VJ 2j1t with mtwcut_ge_sr", 50, -1, 1);
   TH1F *h_2j1t_mtwcut_sr_BDT_STsdvsST[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_BDT_STsdvsST, "h_2j1t_mtwcut_sr_BDT_STsdvsST","BDT STsd-vs-VJ 2j1t with mtwcut_ge_sr", 50, -1, 1);

  
  //mtw > 50; |eta_j'|>2.5, label mtwcut_sr
  TH1F *h_2j1t_mtwcut_sr_topMass[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_topMass, "h_2j1t_mtwcut_sr_topMass", "2j1t top mass",200,100,500);
  TH1F *h_2j1t_mtwcut_sr_topPt[maxSysts];     systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_topPt,   "h_2j1t_mtwcut_sr_topPt",   "2j1t top pt",250,0,500);
  TH1F *h_2j1t_mtwcut_sr_nextrajets[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_nextrajets,  "h_2j1t_mtwcut_sr_nextrajets",  "2j1t number of jets with 20<pt<40",5,-0.5,4.5);

  TH1F *h_2j1t_mtwcut_sr_leadingextrajetpt[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt,  "h_2j1t_mtwcut_sr_leadingextrajetpt",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_mtwcut_sr_leadingextrajetcsv[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv,  "h_2j1t_mtwcut_sr_leadingextrajetcsv",  "2j1t csv of the leading extra jet",100,0.,1.0);
  TH1F *h_2j1t_mtwcut_sr_topMassExtra[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra, "h_2j1t_mtwcut_sr_topMassExtra", "2j1t top mass with the first extra jet as b",200,100,500);
  TH1F *h_2j1t_mtwcut_sr_leadingextrajetcsv_reshape_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv_reshape_sd_b,  "h_2j1t_mtwcut_sr_leadingextrajetcsv_reshape_sd_b",  "2j1t leading extra jet csv reshaping from b to sd for b from top decays",1000,0.,10.0);
  TH1F *h_2j1t_mtwcut_sr_MuPt[maxSysts];      systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_MuPt,    "h_2j1t_mtwcut_sr_MuPt",    "2j1t muon pt with mtwcut sr",500,0,500);

  
  TH1F *h_2j1t_mtwcut_sr_leadingextrajetpt_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt_sd_b,  "h_2j1t_mtwcut_sr_leadingextrajetpt_sd_b",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b,  "h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b",  "2j1t csv of the leading extra jet",100,0.,1.0);
  TH1F *h_2j1t_mtwcut_sr_topMassExtra_sd_b[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra_sd_b, "h_2j1t_mtwcut_sr_topMassExtra_sd_b", "2j1t top mass with the first extra jet as b",200,100,500);
  TH1F *h_2j1t_mtwcut_sr_etajprime[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_etajprime,"h_2j1t_mtwcut_sr_etajprime","2j1t ljet eta>2.5",100,2.5,5.0);  
  //BDT plots 
  TH1F *h_2j1t_BDT_STvsQCDMu[maxSysts]; systZero.initHistogramsSysts(h_2j1t_BDT_STvsQCDMu, "h_2j1t_BDT_STvsQCDMu", "BDT ST-vs-QCD 2j1t", 50, -1, 1);
  TH1F *h_2j1t_mtwcut_BDT_STvsVJ[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_BDT_STvsVJ, "h_2j1t_mtwcut_BDT_STvsVJ", "BDT ST-vs-VJ 2j1t with mtwcut", 50, -1, 1);
  TH1F *h_2j1t_mtwcut_BDT_STvsTT[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_BDT_STvsTT, "h_2j1t_mtwcut_BDT_STvsTT", "BDT ST-vs-VJ 2j1t with mtwcut", 50, -1, 1); 
  
  TH1F *h_2j1t_mtwcut_BDT_STsdvsVJ[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_BDT_STsdvsVJ, "h_2j1t_mtwcut_BDT_STsdvsVJ", "BDT STsd-vs-VJ 2j1t with mtwcut", 50, -1, 1);
  TH1F *h_2j1t_mtwcut_BDT_STsdvsST[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_BDT_STsdvsST, "h_2j1t_mtwcut_BDT_STsdvsST", "BDT STsd-vs-VJ 2j1t with mtwcut", 50, -1, 1);
  
  //Other possible cr:
  TH1F *h_2j1t_mtw_le50[maxSysts];       systZero.initHistogramsSysts(h_2j1t_mtw_le50,     "h_2j1t_mtw_le50",     "2j1t mtw less than 50",50,0,50);
  TH1F *h_2j1t_mtwcut_le_sr_etajprime[maxSysts]; systZero.initHistogramsSysts(h_2j1t_mtwcut_le_sr_etajprime,"h_2j1t_mtwcut_le_sr_etajprime","2j1t ljet abs(eta)<2p5",100,0.,2.5);

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
  TH1F *h_3j1t_MuEta[maxSysts];    systZero.initHistogramsSysts(h_3j1t_MuEta,   "h_3j1t_MuEta",   "muon eta ",100,-2.4,2.4);
  TH1F *h_3j1t_MuPhi[maxSysts];    systZero.initHistogramsSysts(h_3j1t_MuPhi,   "h_3j1t_MuPhi",   "muon phi ",100, -3.2, 3.2);
  TH1F *h_3j1t_MuE[maxSysts];      systZero.initHistogramsSysts(h_3j1t_MuE,     "h_3j1t_MuE",     "muon e ",250,0,500);
  TH1F *h_3j1t_mtw[maxSysts];      systZero.initHistogramsSysts(h_3j1t_mtw,     "h_3j1t_mtw",     "3j1t mtw ",250,0,500);
  TH1F *h_3j1t_dR_lepjetpt40_1st[maxSysts];      systZero.initHistogramsSysts(h_3j1t_dR_lepjetpt40_1st,   "h_3j1t_dR_lepjetpt40_1st",  "dR lep-jet ",50,0,10);
  TH1F *h_3j1t_dPhi_lepjetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_3j1t_dPhi_lepjetpt40_1st, "h_3j1t_dPhi_lepjetpt40_1st",  "dPhi lep-jet ",100,0.0,3.2);
  TH1F *h_3j1t_dEta_lepjetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_3j1t_dEta_lepjetpt40_1st, "h_3j1t_dEta_lepjetpt40_1st",  "dPhi lep-jet ",50,-3.2,3.2);
  TH1F *h_3j1t_muIso[maxSysts];          systZero.initHistogramsSysts(h_3j1t_muIso,"h_3j1t_muIso","Muon Isolation",100,0,1);
  TH1F *h_3j1t_mtwcut_MuPt[maxSysts];     systZero.initHistogramsSysts(h_3j1t_mtwcut_MuPt,    "h_3j1t_mtwcut_MuPt",    "muon pt ",500,0,500);
  TH1F *h_3j1t_mtwcut_muIso[maxSysts];          systZero.initHistogramsSysts(h_3j1t_mtwcut_muIso,"h_3j1t_mtwcut_muIso","Muon Isolation with mtwcut",100,0,1);
  
  TH1F *h_3j1t_jetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_3j1t_jetpt40_1st,   "h_3j1t_jetpt40_1st","3j1t leading jet Pt", 250,0,500);
  TH1F *h_3j1t_jetpt40_2nd[maxSysts];    systZero.initHistogramsSysts(h_3j1t_jetpt40_2nd,   "h_3j1t_jetpt40_2nd","3j1t 2nd leading jet Pt", 250,0,500);
  TH1F *h_3j1t_jetpt40_3rd[maxSysts];    systZero.initHistogramsSysts(h_3j1t_jetpt40_3rd,   "h_3j1t_jetpt40_3rd","3j1t 3rd leading jet Pt", 250,0,500);
  
  //3j2t
  TH1F *h_3j2t_bjetpt[maxSysts];   systZero.initHistogramsSysts(h_3j2t_bjetpt,    "h_3j2t_bjetpt",     "3j2t b jet pt ",100,0,500);
  TH1F *h_3j2t_2ndbjetpt[maxSysts];systZero.initHistogramsSysts(h_3j2t_2ndbjetpt, "h_3j2t_2ndbjetpt",  "3j2t sub lead. b jet pt",100,0,500);
  TH1F *h_3j2t_ljetpt[maxSysts];   systZero.initHistogramsSysts(h_3j2t_ljetpt,    "h_3j2t_ljetpt",     "3j2t light jet pt ",100,0,500);
  TH1F *h_3j2t_ljeteta[maxSysts];   systZero.initHistogramsSysts(h_3j2t_ljeteta,    "h_3j2t_ljeteta",     "3j2t light jet eta ",100,0,500);
  TH1F *h_3j2t_etajprime[maxSysts];   systZero.initHistogramsSysts(h_3j2t_etajprime,    "h_3j2t_etajprime",     "3j2t eta j' ",100,0.,4.7);
  TH1F *h_3j2t_nextrajets[maxSysts];   systZero.initHistogramsSysts(h_3j2t_nextrajets,    "h_3j2t_nextrajets",     "3j2t number of extra jets ",5,-0.5,4.5);
  TH1F *h_3j2t_MuPt[maxSysts];     systZero.initHistogramsSysts(h_3j2t_MuPt,      "h_3j2t_MuPt",       "3j2t muon pt ",500,0,500);
  TH1F *h_3j2t_MuEta[maxSysts];    systZero.initHistogramsSysts(h_3j2t_MuEta,     "h_3j2t_MuEta",      "3j2t muon eta ",100,-2.4,2.4);
  TH1F *h_3j2t_MuPhi[maxSysts];    systZero.initHistogramsSysts(h_3j2t_MuPhi,     "h_3j2t_MuPhi",      "3j2t muon phi ",100, -3.2, 3.2);
  TH1F *h_3j2t_MuE[maxSysts];      systZero.initHistogramsSysts(h_3j2t_MuE,       "h_3j2t_MuE",        "3j2t muon e ",250,0,500);
  TH1F *h_3j2t_MuCharge[maxSysts]; systZero.initHistogramsSysts(h_3j2t_MuCharge,  "h_3j2t_MuCharge",   "3j2t muon charge ",2,-1,1);
  TH1F *h_3j2t_mtw[maxSysts];      systZero.initHistogramsSysts(h_3j2t_mtw,       "h_3j2t_mtw",        "3j2t mtw ",250,0,500);

  TH1F *h_3j2t_muIso[maxSysts];                  systZero.initHistogramsSysts(h_3j2t_muIso,"h_3j2t_muIso","Muon Isolation",100,0,1);
  TH1F *h_3j2t_dR_lepjetpt40_1st[maxSysts];      systZero.initHistogramsSysts(h_3j2t_dR_lepjetpt40_1st,   "h_3j2t_dR_lepjetpt40_1st",  "dR lep-jet ",50,0,10);
  TH1F *h_3j2t_dPhi_lepjetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_3j2t_dPhi_lepjetpt40_1st, "h_3j2t_dPhi_lepjetpt40_1st",  "dPhi lep-jet ",100,0.0,3.2);
  TH1F *h_3j2t_dEta_lepjetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_3j2t_dEta_lepjetpt40_1st, "h_3j2t_dEta_lepjetpt40_1st",  "dPhi lep-jet ",50,-3.2,3.2); 
  TH1F *h_3j2t_mtwcut_MuPt[maxSysts];            systZero.initHistogramsSysts(h_3j2t_mtwcut_MuPt,      "h_3j2t_mtwcut_MuPt",       "3j2t muon pt with mtwcut",500,0,500);
  TH1F *h_3j2t_mtwcut_muIso[maxSysts];           systZero.initHistogramsSysts(h_3j2t_mtwcut_muIso,"h_3j2t_mtwcut_muIso","Muon Isolation with mtwcut",100,0,1);
  
  TH1F *h_3j2t_jetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_3j2t_jetpt40_1st,   "h_3j2t_jetpt40_1st","3j2t leading jet Pt", 250,0,500);
  TH1F *h_3j2t_jetpt40_2nd[maxSysts];    systZero.initHistogramsSysts(h_3j2t_jetpt40_2nd,   "h_3j2t_jetpt40_2nd","3j2t 2nd leading jet Pt", 250,0,500);
  TH1F *h_3j2t_jetpt40_3rd[maxSysts];    systZero.initHistogramsSysts(h_3j2t_jetpt40_3rd,   "h_3j2t_jetpt40_3rd","3j2t 3rd leading jet Pt", 250,0,500);
  //After mtw>50 cut

  TopUtilities topUtils;
  
  //Pileup Reweighting
  edm::LumiReWeighting LumiWeights_, LumiWeightsUp_, LumiWeightsDown_;
  //New-SF files
  TFile* file_eff_trig_bcdef = TFile::Open("data/EfficienciesAndSF_BCDEF_Trigger_SF.root");
  TFile* file_eff_trig_gh    = TFile::Open("data/EfficienciesAndSF_GH_Trigger_SF.root");
  
  TFile* file_sf_id_bcdef = TFile::Open("data/EfficienciesAndSF_BCDEF_ID_SF.root");
  TFile* file_sf_id_gh    = TFile::Open("data/EfficienciesAndSF_GH_ID_SF.root");
  
  TFile* file_sf_iso_bcdef = TFile::Open("data/EfficienciesAndSF_BCDEF_ISO_SF.root");
  TFile* file_sf_iso_gh = TFile::Open("data/EfficienciesAndSF_GH_ISO_SF.root"); 

  Weights muTightTrigEff_BCDEF(file_eff_trig_bcdef,"IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA");
  Weights muTightTrigEff_GH(file_eff_trig_gh,"IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA");
  
  Weights muTightIdSF_BCDEF(file_sf_id_bcdef,"MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
  Weights muTightIdSF_GH(file_sf_id_gh,"MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
  
  Weights muTightIsoSF_BCDEF(file_sf_iso_bcdef,"TightISO_TightID_pt_eta/abseta_pt_ratio");
  Weights muTightIsoSF_GH(file_sf_iso_gh,"TightISO_TightID_pt_eta/abseta_pt_ratio"); 

  //TFile* file_eff_trig = TFile::Open("data/SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root");
  //TFile* file_sf_id = TFile::Open("data/MuonID_Z_RunBCD_prompt80X_7p65.root");
  //TFile* file_sf_iso = TFile::Open("data/MuonIso_Z_RunBCD_prompt80X_7p65.root");
  //cout << "before getting weights "<<endl;
  //Weights muTightTrigEff(file_eff_trig,"IsoMu22_OR_IsoTkMu22_PtEtaBins_Run273158_to_274093/efficienciesDATA/abseta_pt_DATA");
  //Weights muTightIdSF(file_sf_id,"MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio");
  //Weights muTightIsoSF(file_sf_iso,"MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio");
  //  cout << "got weights" <<endl;

    if(isData=="MC"){
        LumiWeights_ = edm::LumiReWeighting("pu/MCPU.root", "pu/MyDataPileupHistogram.root","pileup","pileup");
        LumiWeightsUp_ = edm::LumiReWeighting("pu/MCPU.root", "pu/MyDataPileupHistogram_up.root","pileup","pileup");
        LumiWeightsDown_ = edm::LumiReWeighting("pu/MCPU.root", "pu/MyDataPileupHistogram_down.root","pileup","pileup");
        }
  
    for(Int_t evt=0; evt<nEvents; evt++ ){
        if(evt%10000==1 ){
        cout<<"Info: Running on event: "<<evt<<endl; 
        }
    chain.GetEntry(evt);
    
    //    if (evtNumber!=1751573)continue;
    //    cout << "found it!"<<endl;
    int maxJetLoop = min(15, jetSize);
    int maxJet20Loop = min(15, jet20Size);
    
    //  cout << " jet 20 loop size"<< maxJet20Loop <<endl;
    
    int maxMuLoop = min(6, muSize);
    if(channel == "muonantiiso") maxMuLoop = min(20, muTightIsoGE0p15Size);
    int maxElLoop = min(6, elSize);
  
    //step 1 Trigger
    if(isData=="DATA"){
    TrigIsoMu20= false;
    TrigIsoMu24 = slTrigIsoMu24_v1 || slTrigIsoMu24_v2|| slTrigIsoMu24_v3 || slTrigIsoMu24_v4 || slTrigIsoTkMu24_v1 || slTrigIsoTkMu24_v2|| slTrigIsoTkMu24_v3|| slTrigIsoTkMu24_v4;
    }

    if(isData=="MC"){
        TrigIsoMu20=false;
        TrigIsoMu22=false;
        TrigIsoMu24=true;//Trigger efficiency:use the measured one atm --> always pass the trigger in MC// (slTrigIsoMu22_v1 || slTrigIsoMu22_v2 || slTrigIsoMu22_v3);  
	TrigGSFEl32=true;
        if (doSynch  || channel == "muonantiiso"){
	  TrigIsoMu24=false;
	  TrigGSFEl32=false;
	  //TrigIsoMu24 = true; 
	  TrigIsoMu24 = slTrigIsoMu24_v1 || slTrigIsoMu24_v2|| slTrigIsoMu24_v3 || slTrigIsoMu24_v4 || slTrigIsoTkMu24_v1 || slTrigIsoTkMu24_v2|| slTrigIsoTkMu24_v3|| slTrigIsoTkMu24_v4;
	  TrigGSFEl32= slHLTEle32_eta2p1_WPTight_Gsf_v1||slHLTEle32_eta2p1_WPTight_Gsf_v2||slHLTEle32_eta2p1_WPTight_Gsf_v3||slHLTEle32_eta2p1_WPTight_Gsf_v4||slHLTEle32_eta2p1_WPTight_Gsf_v5;
	}
    }
    bool muonTrigger = (TrigIsoMu20 || TrigIsoMu22 || TrigIsoMu24); 
    bool electronTrigger = (TrigGSFEl32); 
    //bool passesAnyTrigger = ( muonTrigger && channel== "muon") || (electronTrigger&& channel == "electron");//REMINDER: ADD ALL TRIGGERS HERE 
    bool passesAnyTrigger = ( muonTrigger && (channel== "muon" || channel== "muonantiiso")) || (electronTrigger&& channel == "electron");//REMINDER: ADD ALL TRIGGERS HERE 
   
    if(!passesAnyTrigger)continue;

    if(doSynch && passesAnyTrigger){
      fileout_step0<<std::fixed<<std::setprecision(0)<<evtNumber<<std::endl;
    }
    
    n_trig += w;
    nev_trig +=1;

    double puUpFact=1.0;//(LumiWeightsUp_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));
    double puDownFact=1.0;//(LumiWeightsDown_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));


    if(isData=="MC"){
      //cout <<" step 1 "<< w<<" wzero "<< w_zero<< " LHEWeightSign "<< LHEWeightSign[0] << endl;
      LHEWeightSign[0]=w_zero;
      if(w_zero!=0) LHEWeightSign[0]=w_zero/fabs(w_zero);
      w = LHEWeightSign[0];
      
	//cout <<" step 2 "<< w<<endl;
	w_pu = LumiWeights_.weight(numTrueInt);
	
        w = w * w_pu;
    if(sample=="TTNoTT"){
        w*=w_top/topWeight;
        }
    
    if(w_pu!=0){
        puUpFact=(LumiWeightsUp_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));
        puDownFact=(LumiWeightsDown_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));
        }
        } 

    TLorentzVector lep1;
    TLorentzVector lep2;
    TLorentzVector lep, mu, el, muloose;

    vector<TLorentzVector> tightEl, tightMu,looseMu, tightLep;
    int nMu(0.), nEl(0.);
    
    for(int e= 0; e<maxElLoop;++e ){
      if(elPt[e]>35 && fabs(elEta[e])<2.1){
	  el.SetPtEtaPhiE(elPt[e], elEta[e], elPhi[e],elE[e]);
	  tightEl.push_back(el);
      }
    }//end loop on electrons     
  
    vector<float> selectedIso;

    vector<float> mindeltaRmuonloop;
    if(channel!="muonantiiso"){
      for(int m= 0; m<maxMuLoop;++m ){
	if(muPt[m]>26 && abs(muEta[m])<2.4 && muIso[m]<0.06){
	  selectedIso.push_back(muIso[m]);
	  mu.SetPtEtaPhiE(muPt[m], muEta[m], muPhi[m],muE[m]);
	  tightMu.push_back(mu);
	  for (int j = 0; j <maxJetLoop;++j){
	    TLorentzVector all_jets;
	    all_jets.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);
	    float dRtemp = deltaR(all_jets.Eta(),all_jets.Phi(),mu.Eta(),mu.Phi());
	    mindeltaRmuonloop.push_back(dRtemp);
	  } 
	}
      }
    }
    else{
      for(int m= 0; m<maxMuLoop;++m ){
	if(muTightIsoGE0p15Pt[m]>26 && abs(muTightIsoGE0p15Eta[m])<2.4 && muTightIsoGE0p15Iso[m]>0.20){
	  selectedIso.push_back(muTightIsoGE0p15Iso[m]);
	  mu.SetPtEtaPhiE(muTightIsoGE0p15Pt[m], muTightIsoGE0p15Eta[m], muTightIsoGE0p15Phi[m],muTightIsoGE0p15E[m]);
	  tightMu.push_back(mu);
	}
      }
    }
    

    nMu = tightMu.size();
    nEl = tightEl.size();
    vector <lept> leptons;

    if(isData=="MC"){
      if(tightMu.size()==1){
        //--------new efficiencies
        lepWeightBCDEF1Mu = muTightTrigEff_BCDEF.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())*
        muTightIdSF_BCDEF.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())*
        muTightIsoSF_BCDEF.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt());

        lepWeightGH1Mu = muTightTrigEff_GH.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())*
        muTightIdSF_GH.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())*
        muTightIsoSF_GH.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt());
       
        float lumi_bcdef = 18.041; //lumi in fb 
        float lumi_gh = 16.159;
        
        lepWeight1Mu = (lepWeightBCDEF1Mu*lumi_bcdef + lepWeightGH1Mu*lumi_gh)/(lumi_bcdef + lumi_gh);

        float errMu_BCDEF = sqrt(pow(muTightTrigEff_BCDEF.getErr(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt()),2)+
				 pow(muTightIdSF_BCDEF.getErr(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt()),2)+
				 pow(muTightIsoSF_BCDEF.getErr(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt()),2));
       
        float errMu_GH = sqrt(pow(muTightTrigEff_GH.getErr(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt()),2)+
			      pow(muTightIdSF_GH.getErr(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt()),2)+
			      pow(muTightIsoSF_GH.getErr(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt()),2));
        
        float errMu = sqrt(pow(errMu_BCDEF,2)+pow(errMu_GH,2));
	
        lepWeight1MuLepUp=lepWeight1Mu+errMu;
        lepWeight1MuLepDown=lepWeight1Mu-errMu;
        if (channel == "muonantiiso"){
         lepWeightBCDEF1Mu = 1.;
         lepWeightGH1Mu = 1.; 
         lepWeight1MuLepUp = 0.;
         lepWeight1MuLepDown = 0.;
        }
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

    if(channel=="muon" || channel == "muonantiiso")for (size_t m =0;m<tightMu.size();++m){tightLep.push_back(tightMu.at(m));}
    if(channel=="electron")for (size_t e =0;e<tightEl.size();++e){tightLep.push_back(tightEl.at(e));}
    

        
    for (size_t scen=0; scen<scenarios.size();++scen){
      string scenario=scenarios.at(scen);
      
      if(scenario=="jesUp"){ maxJetLoop = min(15, jetJesUpSize);maxJet20Loop = min(15, jet20JesUpSize); 
	metPt[0]=metJesUpPt[0];	metPhi[0]=metJesUpPhi[0];metPx[0]=metJesUpPx[0];metPy[0]=metJesUpPy[0];
      }
      else if(scenario=="jesDown"){ maxJetLoop = min(15, jetJesDownSize);maxJet20Loop = min(15, jet20JesDownSize); 
	metPt[0]=metJesDownPt[0];metPhi[0]=metJesDownPhi[0];metPx[0]=metJesDownPx[0];metPy[0]=metJesDownPy[0];
      }
      else {
	metPt[0]=metZeroPt[0];	metPhi[0]=metZeroPhi[0];metPx[0]=metZeroPx[0];metPy[0]=metZeroPy[0];
	maxJetLoop = min(15, jetSize);
	maxJet20Loop = min(15, jet20Size);
      }

      met = metPt[0];
      metpx = metPx[0];
      metpy = metPy[0];

      //cout << "uffa8 "<<endl;
      systZero.setWeight(0,1.);
      systZero.setWeight("btagUp",1.);
      systZero.setWeight("btagDown",1.);
      systZero.setWeight("mistagUp",1.);
      systZero.setWeight("mistagDown",1.);
      systZero.setWeight("puDown",1.);
      systZero.setWeight("puUp",1.);
      systZero.setWeight("lepDown",1.);
      systZero.setWeight("lepUp",1.);
      
      systZero.setScenario(scenario);
      systZero.setEventBasedDefault();
      
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
      syst1BM.setWeight("puUp",bWeight1CSVM*puUpFact);
      syst1BM.setWeight("puDown",bWeight1CSVM*puDownFact);
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
      
      if(addPDF){
	syst0BM.setPDFWeights(w_pdfs,nPDF,w_zero,true);
	syst1BM.setPDFWeights(w_pdfs,nPDF,w_zero,true);
	syst2BM.setPDFWeights(w_pdfs,nPDF,w_zero,true);
	
      }
      if(addQ2){
	syst0BM.setQ2Weights(w_q2up,w_q2down,w_zero,true);
	syst1BM.setQ2Weights(w_q2up,w_q2down,w_zero,true);
	syst2BM.setQ2Weights(w_q2up,w_q2down,w_zero,true);
	
      }
      if(addTopPt){
	syst0BM.setTWeight(w_top,topWeight,true);
	syst1BM.setTWeight(w_top,topWeight,true);
	syst2BM.setTWeight(w_top,topWeight,true);
      }
      
      systZero.setEventBasedDefault();
      syst0BM.setEventBasedDefault();
      syst1BM.setEventBasedDefault();
      syst2BM.setEventBasedDefault();
      
    vector<float> jetsPhi;
    
    vector<TLorentzVector> bjets, jets_nob, jets, jets20, jets_nob20, bjets20;    

    struct btag{
        TLorentzVector vect;
        int passesB;
        float csv;
        float csvreweight;
        float csvreweightsd;
        float partonFlavour;
        float hadronFlavour;
        };
  
    std::vector<btag> bvects,standardjets,extrajets;
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

    struct by_eta_jet{
      bool operator()(TLorentzVector const &jet1, TLorentzVector const &jet2){
	return fabs(jet1.Eta()) < fabs(jet1.Pt());
      }
    };
    
    for (int j = 0; j <maxJet20Loop;++j){
      //     cout <<" in j20 loop"<<endl;
      TLorentzVector jet20;

      float jetLocEta=jet20Eta[j], jetLocPt=jet20Pt[j], jetLocPhi=jet20Phi[j], jetLocE=jet20E[j], jetLocPassesB=jet20PassesB[j],jetLocak4chs_csvv2=jet20ak4chs_csvv2[j],jetLocReshapeFactorCSV=jet20ReshapeFactorCSV[j],jetLocReshapeFactorCSV_SD=jet20ReshapeFactorCSV_SD[j],jetLocPartonFlavour=jet20PartonFlavour[j];
      
      if(scenario=="jesUp"){jetLocEta=jet20JesUpEta[j];jetLocPt=jet20JesUpPt[j]; jetLocPhi=jet20JesUpPhi[j]; jetLocE=jet20JesUpE[j]; jetLocPassesB=jet20JesUpPassesB[j]; jetLocak4chs_csvv2=jet20JesUpak4chs_csvv2[j];jetLocReshapeFactorCSV=jet20JesUpReshapeFactorCSV[j];jetLocReshapeFactorCSV_SD=jet20JesUpReshapeFactorCSV_SD[j];jetLocPartonFlavour=jet20JesUpPartonFlavour[j];}
      if(scenario=="jesDown"){jetLocEta=jet20JesDownEta[j];jetLocPt=jet20JesDownPt[j]; jetLocPhi=jet20JesDownPhi[j]; jetLocE=jet20JesDownE[j]; jetLocPassesB=jet20JesDownPassesB[j]; jetLocak4chs_csvv2=jet20JesDownak4chs_csvv2[j];jetLocReshapeFactorCSV=jet20JesDownReshapeFactorCSV[j];jetLocReshapeFactorCSV_SD=jet20JesDownReshapeFactorCSV_SD[j];jetLocPartonFlavour=jet20JesDownPartonFlavour[j];}
      
      jet20.SetPtEtaPhiE(jetLocPt, jetLocEta, jetLocPhi, jetLocE);
	bool btagcond=jetLocPassesB>0.&&fabs(jetLocEta)<2.4;
	if(btagcond)bjets20.push_back(jet20);
	else jets_nob20.push_back(jet20);
	
	if(jetLocPt<40 && fabs(jetLocEta)<2.4){
            //cout <<" # j= "<<j<< " pt "<< jet20Pt <<endl;
            jets20.push_back(jet20);
            btag bj;
            bj.vect = jet20;
            bj.csv = jetLocak4chs_csvv2;
            bj.csvreweight = jetLocReshapeFactorCSV;
            bj.csvreweightsd = 1.;//jet20:ak4chsReweightFactorCSV_SD[j];
            if( doTtoSDDecay){
                bj.csvreweightsd = jetLocReshapeFactorCSV_SD;
		bj.partonFlavour=jetLocPartonFlavour;
	    }
            extrajets.push_back(bj);
            }
        }
   
    for (int j = 0; j <maxJetLoop;++j){
      TLorentzVector jet, all_jets;
      float jetLocEta=jetEta[j], jetLocPt=jetPt[j], jetLocPhi=jetPhi[j], jetLocE=jetE[j], jetLocPassesB=jetPassesB[j], jetLocak4chs_csvv2=jetak4chs_csvv2[j], jetLocReshapeFactorCSV=jetReshapeFactorCSV[j],jetLocReshapeFactorCSV_SD=jetReshapeFactorCSV_SD[j],jetLocPartonFlavour=jetPartonFlavour[j];
      
      if(scenario=="jesUp"){jetLocEta=jetJesUpEta[j]; jetLocPt=jetJesUpPt[j]; jetLocPhi=jetJesUpPhi[j]; jetLocE=jetJesUpE[j]; jetLocPassesB=jetJesUpPassesB[j]; jetLocak4chs_csvv2=jetJesUpak4chs_csvv2[j]; jetLocReshapeFactorCSV=jetJesUpReshapeFactorCSV[j]; jetLocReshapeFactorCSV_SD=jetJesUpReshapeFactorCSV_SD[j]; jetLocPartonFlavour=jetJesUpPartonFlavour[j];}
      
      if(scenario=="jesDown"){jetLocEta=jetJesDownEta[j]; jetLocPt=jetJesDownPt[j]; jetLocPhi=jetJesDownPhi[j]; jetLocE=jetJesDownE[j]; jetLocPassesB=jetJesDownPassesB[j]; jetLocak4chs_csvv2=jetJesDownak4chs_csvv2[j]; jetLocReshapeFactorCSV=jetJesDownReshapeFactorCSV[j]; jetLocReshapeFactorCSV_SD=jetJesDownReshapeFactorCSV_SD[j]; jetLocPartonFlavour=jetJesDownPartonFlavour[j];}
      
      all_jets.SetPtEtaPhiE(jetLocPt, jetLocEta, jetLocPhi, jetLocE);
      
      
      bool passesDR=true;
      double dR=0; 
      //--------------------------------------------------------------------------------
      
      
      if(channel=="muonantiiso"){
	for (size_t i = 0; i<(size_t)tightMu.size();++i){
	  dR= deltaR(all_jets.Eta(),all_jets.Phi(),tightMu[i].Eta(),tightMu[i].Phi());
	  if (dR<0.4) passesDR=false; //throw away the jet if deltaR btw jet and tightmu in less than 0.4
	}
      } 
      if(!passesDR)continue;
      
      
      jets.push_back(all_jets);
      bool btagcond = jetLocPassesB>0. && fabs(jetLocEta)<2.4;
      
      if(!btagcond){
  	jet.SetPtEtaPhiE(jetLocPt, jetLocEta, jetLocPhi, jetLocE);
  	jets_nob.push_back(jet);
      }
      
      jetsPhi.push_back(jetLocPhi);
      TLorentzVector bjet;
      
      bjet.SetPtEtaPhiE(jetLocPt, jetLocEta, jetLocPhi, jetLocE);
      btag b;
      b.vect = bjet;
      //      cout << "csv before loop "<< jetLocak4chs_csvv2<< endl;
      b.csv = jetLocak4chs_csvv2;
      b.passesB = jetLocPassesB;
      b.partonFlavour=jetLocPartonFlavour;
      b.csvreweightsd = jetLocReshapeFactorCSV;

      if( doTtoSDDecay){
  	b.csvreweightsd = jetLocReshapeFactorCSV_SD;
      }
      
      if( btagcond && fabs(jetLocEta)<2.4){
  	bjets.push_back(bjet);
  	bvects.push_back(b);
      }
      standardjets.push_back(b);
    }//End of jet loop
    
    nJets = jets.size();
    nCSVJets=bjets.size();
    //  cout <<"passes trig "<<endl; 
    std::sort(bvects.begin(), bvects.end(), by_pt_jet()); 
    std::sort(extrajets.begin(), extrajets.end(), by_pt_jet()); 
    bool passmuon = muonTrigger && nMu == 1 ;
    bool passelectron = electronTrigger && nEl == 1 ;

    bool passlepton = (passmuon && ((channel=="muon" || channel=="muonantiiso"))) || ((passelectron && channel=="electron"));
    if(passlepton  ==  true ){n_lepton+=w;nev_lepton+=1;}
    //    cout << "passlepton" <<passlepton;
    if(doSynch && passlepton && scenario=="nominal"){
      fileout_step1<<std::fixed<<std::setprecision(0)<<evtNumber<<endl;
    }
    if(channel=="muonantiiso") passlepton = passlepton && ( ( (channel=="muon" || channel=="muonantiiso")  && muLooseSize==0) || (channel =="electron" && elLooseSize==1) );
    else passlepton = passlepton && ( ( (channel=="muon" || channel=="muonantiiso")  && muLooseSize==1) || (channel =="electron" && elLooseSize==1) );
    //    cout << "passlepton after lep step" <<passlepton;
    if(passlepton){ n_loose_veto+=w; nev_loose_veto+=1;}
    if(doSynch && passlepton){
      fileout_step2<<std::fixed<<std::setprecision(0)<<evtNumber<<std::endl;
    }
    passlepton = passlepton && ( ( (channel=="muon" || channel=="muonantiiso") && elLooseSize == 0 ) || (channel =="electron" && muLooseSize==0) );
    //    cout << "passlepton after lep 2 step" <<passlepton;
    if(passlepton){n_lepton_cross_veto+=w;nev_lepton_cross_veto+=1;}
    if(doSynch && passlepton && scenario=="nominal"){
      fileout_step3<<std::fixed<<std::setprecision(0)<<evtNumber<<std::endl;
    }
    
    if(topsize==1)topcharge=1.;
    else if(antitopsize==1)topcharge=-1.;
    bool passsinglelepton = passlepton;

    if(channel=="muonantiiso" && !passlepton) continue;
    if((channel=="muon" || channel=="electron") && !passlepton)continue; 
    //cout <<"Evt No. after "<<evt<<" tightLepons "<<nMu<<" Loose Muons "<<muLooseSize<<endl; 
    if(channel=="singlelepton" && !passsinglelepton)continue;
    
    
    if(jets.size() ==2){
      n_2j+=w; nev_2j+=1;
      if(doSynch && scenario=="nominal"){
	fileout_step4<<std::fixed<<std::setprecision(0)<<evtNumber<<std::endl;
      }
    }
    mt = -999;
    if(jets.size() ==2 && bjets.size()==1){ n_2j1t+=w; nev_2j1t+=1;  

      if(doSynch && scenario=="nominal"){
	fileout_step5<<std::fixed<<std::setprecision(0)<<evtNumber<<std::endl;
      }
      if(doSynch&& scenario=="nominal"){
	if((tightLep.size())==1){
	  TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
	  float phi_lmet = fabs(deltaPhi(tightLep[0].Phi(), metPhi[0]) );
	  mt = sqrt(2* tightLep[0].Pt() * met* ( 1- cos(phi_lmet)));
	}
	if(((channel == "muon" || channel=="muonantiiso") && mt> 50) || (channel == "electron" && metPt[0] >50)){
          n_2j1tmtw+=w;
	  nev_2j1tmtw+=1.;
	  
	  if(doSynch){
	    fileout_step6<<std::fixed<<std::setprecision(0)<<evtNumber<<std::endl;
	  }
	}
      }
    }
    

    systZero.fillHistogramsSysts(h_nJets,nJets,w); 
    systZero.fillHistogramsSysts(h_nbJets,nCSVJets,w); 
    systZero.fillHistogramsSysts(h_nPV,nPV,w);
    systZero.fillHistogramsSysts(h_nGoodPV,nPV,w);
    systZero.fillHistogramsSysts(h_nTruePV,numTrueInt,1);
    

    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------    
    
    //2j0t 
    if((jets.size() == 2 && bjets.size() == 0)){
      for(size_t i = 0; i < (size_t)tightLep.size();++i ){
	if((tightLep.size()) == 1){
            TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
            float phi_lmet = fabs(deltaPhi(tightLep[i].Phi(), metPhi[0]) );
            mt = sqrt(2* tightLep[i].Pt() * met* ( 1- cos(phi_lmet)));
	    //cout << "muon iso "<< selectedIso[i]<< " mtw "<< mt <<endl;
            syst0BM.fillHistogramsSysts(h_2j0t_mtw,mt,w);
            }
        }
        for(size_t j= 0; j< (size_t)jets.size();++j ){
        //cout <<"Evt No. "<<evt<<"   "<<"   "<<"jet["<<j<< "] "<<jets[j].Pt()<<"   bjetsize "<< bjets.size() <<"   Mu[0] "<<tightLep[0].Pt()<< std::endl;
	  if(j==0){
	    syst0BM.fillHistogramsSysts(h_2j0t_jetpt40_1st,jets[j].Pt(),w);
	    syst0BM.fillHistogramsSysts(h_2j0t_MuEta,tightLep[0].Eta(),w);
	    syst0BM.fillHistogramsSysts(h_2j0t_MuPhi,tightLep[0].Phi(),w);     
	    syst0BM.fillHistogramsSysts(h_2j0t_dR_lepjetpt40_1st,deltaR(jets[j].Eta(),jets[j].Phi(),tightLep[0].Eta(),tightLep[0].Phi()),w);
	    syst0BM.fillHistogramsSysts(h_2j0t_dPhi_lepjetpt40_1st,fabs(deltaPhi(jets[j].Phi(),tightLep[0].Phi())),w);
	    syst0BM.fillHistogramsSysts(h_2j0t_dEta_lepjetpt40_1st,jets[j].Eta()-tightLep[0].Eta(),w);
	    syst0BM.fillHistogramsSysts(h_2j0t_MuPt,tightLep[0].Pt(),w);
	  }
	  
	  if(j==1)syst0BM.fillHistogramsSysts(h_2j0t_jetpt40_2nd,jets[j].Pt(),w);
	  syst0BM.fillHistogramsSysts(h_2j0t_incljetabseta,fabs(jets[j].Eta()),w);  
	  
	  if(mt>50.0){
	    if(j==0){
	      syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_jetpt40_1st,jets[j].Pt(),w);  
	      syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_jeteta40_1st,jets[j].Pt(),w);  
	      syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_jetabseta40_1st,fabs(jets[j].Eta()),w);  
	      syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_incljeteta,jets[j].Eta(),w);  
	      syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_incljetabseta,fabs(jets[j].Eta()),w);
	      syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_MuPt,tightLep[0].Phi(),w);
	      if(channel=="muon" || channel == "muonantiiso")syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_muIso,selectedIso[0],w);
	    }
	    
            if(j==1){
                syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_jetpt40_2nd,jets[j].Pt(),w);  
                syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_incljeteta,jets[j].Eta(),w);  
                syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_incljetabseta,fabs(jets[j].Eta()),w);  
	    }
	  }
        }
    }
    
    //end 2j0t
    
    double etathreshold=2.5;
    
    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------    
    
    //2j1t
    if((jets.size() == 2 && bjets.size() == 1)){
      //      cout << "syst1BM "<< " scenario "<< scenario<< " weight zero " << syst1BM.weightedSysts[0]<< " alt "<< syst1BM.getSystValue(0)<<" weight jesUp "<< syst1BM.getSystValue("jesUp")<<endl;
      
      for(size_t j= 0; j< (size_t)jets.size();++j ){
	if(j==0){
	  syst1BM.fillHistogramsSysts(h_2j1t_dR_lepjetpt40_1st,deltaR(jets[j].Eta(),jets[j].Phi(),tightLep[0].Eta(),tightLep[0].Phi()),w);
	  syst1BM.fillHistogramsSysts(h_2j1t_dPhi_lepjetpt40_1st,fabs(deltaPhi(jets[j].Phi(),tightLep[0].Phi())),w);
	  syst1BM.fillHistogramsSysts(h_2j1t_dEta_lepjetpt40_1st,jets[j].Eta()-tightLep[0].Eta(),w);
	  syst1BM.fillHistogramsSysts(h_2j1t_jetpt40_1st,jets[j].Pt(),w);
	}
	if(j==1){
	  syst1BM.fillHistogramsSysts(h_2j1t_jetpt40_2nd,jets[j].Pt(),w);
	}
      }
      for(size_t b= 0; b< (size_t)bjets.size();++b ){
	if(b==0){
          syst1BM.fillHistogramsSysts(h_2j1t_bjetpt,bjets[b].Pt(),w,NULL,false);
	}
      }
      
      bool qcddepleted=false, signalenriched=false;
      //define QCD depletion condition:
      for(size_t i = 0; i < (size_t)tightLep.size();++i ){
        syst1BM.fillHistogramsSysts(h_2j1t_MuPt,tightLep[i].Pt(),w);
        syst1BM.fillHistogramsSysts(h_2j1t_MuEta,tightLep[i].Eta(),w);
        syst1BM.fillHistogramsSysts(h_2j1t_MuPhi,tightLep[i].Phi(),w);
        syst1BM.fillHistogramsSysts(h_2j1t_MuE,tightLep[i].E(),w);
	if((tightLep.size()) == 1 ){
            TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
            float phi_lmet = fabs(deltaPhi(tightLep[i].Phi(), metPhi[0]) );
            mt = sqrt(2* tightLep[i].Pt() * met* ( 1 - cos(phi_lmet)));
            syst1BM.fillHistogramsSysts(h_2j1t_mtw,mt,w);
            if(channel=="muon" || channel == "muonantiiso")syst1BM.fillHistogramsSysts(h_2j1t_muIso,selectedIso[i],w);
	    *mtw_2j1t = mt;
	    *met_2j1t = met;
            if (calculate_mtw(metPt, metPhi,tightLep) > 50.0){
	      qcddepleted=true;
                syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_mtw,mt,w);
                syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_MuPt,tightLep[i].Pt(),w);
                if(channel=="muon" || channel == "muonantiiso")syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_muIso,selectedIso[i],w);
            } 
            if (!qcddepleted) syst1BM.fillHistogramsSysts(h_2j1t_mtw_le50,mt,w);
        }
      }
      
      
      TLorentzVector jetprime_2j1t;
      //define signal enriching condition:
      *bjetflavour_2j1t=-9.5;
      for (size_t j= 0; j< (size_t)jets.size();++j ){
	bool btagcond = standardjets[j].passesB>0. && fabs(jets[j].Eta())<2.4;
	if(!(btagcond)){
	  syst1BM.fillHistogramsSysts(h_2j1t_ljetpt,jets[j].Pt(),w);
	  syst1BM.fillHistogramsSysts(h_2j1t_ljeteta,jets[j].Eta(),w);
	  syst1BM.fillHistogramsSysts(h_2j1t_etajprime,fabs(jets[j].Eta()),w);
	  signalenriched=signalenriched||fabs(jets[j].Eta())>etathreshold;
	  *etajprime_2j1t = fabs(jets[j].Eta());
	  *jprimeflavour_2j1t=standardjets[j].partonFlavour;
	  jetprime_2j1t = jets[j];
	  *mljprime_2j1t=topUtils.mlj(tightLep.at(0), jets.at(j));
	  if(qcddepleted) {
	    syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_ljetpt,jets[j].Pt(),w);
	    syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_ljeteta,jets[j].Eta(),w);
	    syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_etajprime,fabs(jets[j].Eta()),w);
	  }
	  if(qcddepleted && !signalenriched) syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_le_sr_etajprime,fabs(jets[j].Eta()),w);
	  if(qcddepleted && signalenriched)  syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_etajprime,fabs(jets[j].Eta()),w);
	}
	else{
	  *bjetflavour_2j1t=standardjets[j].partonFlavour;
	  *bjeteta_2j1t=jets[j].Eta();
	  *bjetpt_2j1t=jets[j].Pt();
	}
      }
      syst1BM.fillHistogramsSysts(h_2j1t_nextrajets,extrajets.size(),w);
      reweightall=true;
      if(extrajets.size()){
	
	float weightreshape=extrajets.at(0).csvreweight;
	float weight_sd_b=1.0;
	bool fillcond=true;
	
	if( doTtoSDDecay){
	  fillcond = !(extrajets.at(0).partonFlavour ==5 && extrajets.at(0).partonFlavour*topcharge>0&&fabs(extrajets.at(0).vect.Eta())<2.4);
	}
	if(fillcond){
	  syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetpt,extrajets.at(0).vect.Pt(),w*weightreshape);
	  syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetcsv,extrajets.at(0).csv,w*weightreshape);
	}
	else{
	  weight_sd_b=extrajets.at(0).csvreweightsd;
	  syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetcsv_sd_b,extrajets.at(0).csv,w*weightreshape*weight_sd_b);
	  syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetcsv_reshape_sd_b, weight_sd_b,w);
	  if(reweightall){syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetpt_sd_b,extrajets.at(0).vect.Pt(),w*weightreshape*weight_sd_b);}
	  else {syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetpt,extrajets.at(0).vect.Pt(),w);}
	}
      }
      //    cout << "extrajets ok"<<endl;
      
      //Adding top related variables
      vector<TLorentzVector> tops= topUtils.top4Momenta(tightLep, bjets ,metpx, metpy);
      
      *mindeltaphi_2j1t = topUtils.mindeltaphi(metPhi[0],jets);
      *mindeltaphi20_2j1t = topUtils.mindeltaphi(metPhi[0],jets20);
      
      TVector2 metv( metpx, metpy);
      double Mt2w = Mt2cal->calculateMT2w(jets_nob20,bjets20,tightLep.at(0), metv,"MT2w"); 
      *mt2w_2j1t = Mt2w;
      
      if(tops.size()!=0){
	double mass = tops.at(0).M();
	double pt = tops.at(0).Pt();
	
	
	*topMass_2j1t=mass;
	*topY_2j1t=tops.at(0).Y();
	*topPt_2j1t=pt;
	*topEta_2j1t=tops.at(0).Eta();
	*costhetael_2j1t=topUtils.costhetael(tops.at(0), tightLep.at(0), metpx, metpy);
	*costhetapol_2j1t=topUtils.costhetapol(tightLep.at(0),jetprime_2j1t, tops.at(0));
	*mlb_2j1t=topUtils.mlj(tightLep.at(0), bjets.at(0));
	
	*topMt_2j1t=topUtils.topMtw(tightLep.at(0), bjets.at(0) ,metpx, metpy);
      
	*topPt_2j1t=pt;
	*topEta_2j1t=tops.at(0).Eta();
	*mlb_2j1t=topUtils.mlj(tightLep.at(0), bjets.at(0));
	syst1BM.fillHistogramsSysts(h_2j1t_topMass,mass,w);
	syst1BM.fillHistogramsSysts(h_2j1t_topPt,pt,w);
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
      *nextrajets_2j1t=extrajets.size();
      *leadingextrajetpt_2j1t=-5;
      *leadingextrajetcsv_2j1t=-5;
      *topMassExtra_2j1t=-5;
      *topYExtra_2j1t=-5;
      *topPtExtra_2j1t=-5;
      *topEtaExtra_2j1t=-5;
      *topMtExtra_2j1t=-5;
      *leadingextrajeteta_2j1t=-5;
      *leadingextrajetcsvweight_2j1t=0;
      *leadingextrajetcsvweight_sd_2j1t=0;
      *leadingextrajetflavour_2j1t=-9.5;
      *costhetaelExtra_2j1t=-5;
      *costhetapolExtra_2j1t=-5;
      *mljextra_2j1t=-5;
      
      if(extrajets.size()){
	vector<TLorentzVector> topsextra= topUtils.top4Momenta(tightLep, jets20 ,metpx, metpy);
	
	float weightreshape=extrajets.at(0).csvreweight;
	bool fillcond=true;
	if( doTtoSDDecay){
	  fillcond = (!(extrajets.at(0).partonFlavour ==5
			&& (extrajets.at(0).partonFlavour*topcharge>0)
			&& fabs(extrajets.at(0).vect.Eta())<2.4));		 
	}
	if(topsextra.size()!=0){
	  *leadingextrajetcsvweight_2j1t=extrajets.at(0).csvreweight;
	  if(fillcond){
	  *leadingextrajetcsvweight_sd_2j1t=0;
	  }
	  else{
	    float weight_sd_b=extrajets.at(0).csvreweightsd;
	    *leadingextrajetcsvweight_sd_2j1t=weight_sd_b;
	  }
	  //      cout << "extrajet csvs ok"<<endl;
	  *leadingextrajetpt_2j1t=extrajets.at(0).vect.Pt();
	  *leadingextrajetflavour_2j1t=extrajets.at(0).partonFlavour;
	  *leadingextrajeteta_2j1t=extrajets.at(0).vect.Eta();
	  *leadingextrajetcsv_2j1t=extrajets.at(0).csv;
	  *topMassExtra_2j1t=topsextra[0].M();
	  *topYExtra_2j1t=topsextra[0].Y();
	  *topEtaExtra_2j1t=topsextra[0].Eta();
	  *topPtExtra_2j1t=topsextra[0].Pt();
	  *topMtExtra_2j1t = topUtils.topMtw(tightLep.at(0), jets20.at(0) ,metpx, metpy);
	  *costhetaelExtra_2j1t=topUtils.costhetael(topsextra.at(0), tightLep.at(0), metpx, metpy);
	  *costhetapolExtra_2j1t=topUtils.costhetapol(tightLep.at(0),jetprime_2j1t, topsextra.at(0));
	  *mljextra_2j1t=topUtils.mlj(tightLep.at(0), extrajets.at(0).vect);
	}
      }
      *w_2j1t=w;
      
      if(doMVA){
	for(size_t rd =0; rd < AllReaders.size();++rd){
	  if(AllReaderNames.at(rd)==("ST_vs_QCDMu")) {
	  *bdt_st_vs_qcd_2j1t = AllReaders.at(rd)->EvaluateMVA("BDTG"); 
	  syst1BM.fillHistogramsSysts(h_2j1t_BDT_STvsQCDMu,*bdt_st_vs_qcd_2j1t,w);
	  }
	}
      }
      if(addTrees)syst1BM.fillTreesSysts(trees1T,"2j1t");
      
    
      //End part for generic 2j1t, follow 2j1t specific CRs
      
      if(qcddepleted){
	syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_nextrajets,extrajets.size(),w);

       if(signalenriched){
	 syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_nextrajets,extrajets.size(),w);
       }
       if(extrajets.size()){
	 vector<TLorentzVector> topsextra= topUtils.top4Momenta(tightLep, jets20 ,metpx, metpy);
	 float weightreshape=extrajets.at(0).csvreweight;
	 bool fillcond=true;
	 if( doTtoSDDecay){
	   fillcond = (!(extrajets.at(0).partonFlavour ==5
			 && (extrajets.at(0).partonFlavour*topcharge>0)
			 && fabs(extrajets.at(0).vect.Eta())<2.4));		 
	 }
	 if(topsextra.size()!=0){
	   if(fillcond){
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_topMassExtra,topsextra[0].M(),w);
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt,extrajets.at(0).vect.Pt(),w);
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv,extrajets.at(0).csv,w*weightreshape);
	   }
	   else{
	     float weight_sd_b=extrajets.at(0).csvreweightsd;
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv_sd_b,extrajets.at(0).csv,w*weightreshape*weight_sd_b);
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv_reshape_sd_b, weight_sd_b,w);
	     if(reweightall){
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_topMassExtra_sd_b,topsextra[0].M(),w*weight_sd_b);
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt_sd_b,extrajets.at(0).vect.Pt(),w*weight_sd_b);
	     }
	     else{
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_topMassExtra,topsextra[0].M(),w);
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt,extrajets.at(0).vect.Pt(),w);
	     }
	   }
	   //	   cout << "extrajet csvs ok"<<endl;
	   
	 }
	 if(signalenriched){
	   if(fillcond){
	   //	   cout<< " filling histograms for the signal enriched region"<<endl;
	   syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra,topsextra[0].M(),w);
	   syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt,extrajets.at(0).vect.Pt(),w);
	   syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv,extrajets.at(0).csv,w*weightreshape);
	   }
	   else{
	     float weight_sd_b=extrajets.at(0).csvreweightsd;
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b,extrajets.at(0).csv,w*weightreshape*weight_sd_b);
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv_reshape_sd_b, weight_sd_b,w);
	     if(reweightall){
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra_sd_b,topsextra[0].M(),w*weight_sd_b);
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt_sd_b,extrajets.at(0).vect.Pt(),w*weight_sd_b);
	     }
	     else{
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra_sd_b,topsextra[0].M(),w);
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt_sd_b,extrajets.at(0).vect.Pt(),w);
	     }
	   }
	   //	   cout << "extrajet csvs sr ok"<<endl;
	 }
       }
      }
      if(qcddepleted){
	if(doMVA){
	  for(size_t rd =0; rd < AllReaders.size();++rd){
	    //cout << AllReaders.at(rd)->EvaluateMVA("BDTG")<<endl;
	    if(AllReaderNames.at(rd)==("ST_vs_VJ")) {
	      *bdt_st_vs_vj_2j1t_mtwcut = AllReaders.at(rd)->EvaluateMVA("BDTG");
	      syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_BDT_STvsVJ,*bdt_st_vs_vj_2j1t_mtwcut,w);
	    }
	    if(AllReaderNames.at(rd)==("ST_vs_TT")) {
	    *bdt_st_vs_tt_2j1t_mtwcut = AllReaders.at(rd)->EvaluateMVA("BDTG");
	    syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_BDT_STvsTT,*bdt_st_vs_tt_2j1t_mtwcut,w);
	    }
	    if(AllReaderNames.at(rd)==("STsd_vs_VJ")) {
	      *bdt_stsd_vs_vj_2j1t_mtwcut = AllReaders.at(rd)->EvaluateMVA("BDTG");
	      syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_BDT_STsdvsVJ,*bdt_stsd_vs_vj_2j1t_mtwcut,w);
	    }
	    if(AllReaderNames.at(rd)==("STsd_vs_ST")) {
	    *bdt_stsd_vs_st_2j1t_mtwcut = AllReaders.at(rd)->EvaluateMVA("BDTG");
	    syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_BDT_STsdvsST,*bdt_stsd_vs_st_2j1t_mtwcut,w);
	    }
	  }
	  
	  if(!signalenriched){
	    syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_le_sr_BDT_STvsVJ,*bdt_st_vs_vj_2j1t_mtwcut,w);
	    syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_le_sr_BDT_STvsTT,*bdt_st_vs_tt_2j1t_mtwcut,w);
	    syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_le_sr_BDT_STsdvsVJ,*bdt_stsd_vs_vj_2j1t_mtwcut,w);
	    syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_le_sr_BDT_STsdvsST,*bdt_stsd_vs_st_2j1t_mtwcut,w);
	  }
	  if(signalenriched){
	    syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_BDT_STvsVJ,*bdt_st_vs_vj_2j1t_mtwcut,w);
	    syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_BDT_STvsTT,*bdt_st_vs_tt_2j1t_mtwcut,w);
	    syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_BDT_STsdvsVJ,*bdt_stsd_vs_vj_2j1t_mtwcut,w);
	    syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_BDT_STsdvsST,*bdt_stsd_vs_st_2j1t_mtwcut,w);
	  }
	}
      }
      //AOB for 2j1t?
    }
    //End 2j1t

    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------    
    
    //3j1t   
    if((jets.size() == 3 && bjets.size() ==1)) {
      *w_3j1t=w;
      for(size_t j= 0; j< (size_t)jets.size();++j ){
	if(j==0){
	  syst1BM.fillHistogramsSysts(h_3j1t_dR_lepjetpt40_1st,deltaR(jets[j].Eta(),jets[j].Phi(),tightLep[0].Eta(),tightLep[0].Phi()),w);
	  syst1BM.fillHistogramsSysts(h_3j1t_dPhi_lepjetpt40_1st,fabs(deltaPhi(jets[j].Phi(),tightLep[0].Phi())),w);
	  syst1BM.fillHistogramsSysts(h_3j1t_dEta_lepjetpt40_1st,jets[j].Eta()-tightLep[0].Eta(),w);
	  syst1BM.fillHistogramsSysts(h_3j1t_jetpt40_1st,jets[j].Pt(),w);
	}
	if(j==1){syst1BM.fillHistogramsSysts(h_3j1t_jetpt40_2nd,jets[j].Pt(),w); }
	if(j==2){syst1BM.fillHistogramsSysts(h_3j1t_jetpt40_3rd,jets[j].Pt(),w); }
      }
      for(size_t b= 0; b< (size_t)bjets.size();++b ){
	if(b==0)syst1BM.fillHistogramsSysts(h_3j1t_bjetpt,bjets[b].Pt(),w);
      } 
      for(size_t i = 0; i < (size_t)tightLep.size();++i ){
        syst1BM.fillHistogramsSysts(h_3j1t_MuPt,tightLep[i].Pt(),w);
        syst1BM.fillHistogramsSysts(h_3j1t_MuEta,tightLep[i].Eta(),w);
        syst1BM.fillHistogramsSysts(h_3j1t_MuPhi,tightLep[i].Phi(),w);
        syst1BM.fillHistogramsSysts(h_3j1t_MuE,tightLep[i].E(),w);
        if((tightLep.size()) == 1 ){
          TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
          float phi_lmet = fabs(deltaPhi(tightLep[i].Phi(), metPhi[0]) );
          mt = sqrt(2* tightLep[i].Pt() * met* ( 1- cos(phi_lmet)));
	  *mtw_3j1t=mt;
	  syst1BM.fillHistogramsSysts(h_3j1t_mtw,mt,w);
	  if(channel=="muon" || channel == "muonantiiso")syst1BM.fillHistogramsSysts(h_3j1t_muIso,selectedIso[i],w);
	  if(mt>50){
	    syst1BM.fillHistogramsSysts(h_3j1t_mtwcut_MuPt,tightLep[i].Pt(),w);
	    if(channel=="muon" || channel == "muonantiiso")syst1BM.fillHistogramsSysts(h_3j1t_mtwcut_muIso,selectedIso[i],w);
	  }
	}
      }
      float maxEta=0., minEta=6.;
      TLorentzVector jetprime_3j1t,bjet;
      btag  jetnobnoprime;
      vector<TLorentzVector> jetstopMass;
      if(bjets.size()!=0){
	for(size_t i=0; i<(size_t)standardjets.size(); ++i){
	  bool btagcond = standardjets[i].passesB>0. && fabs(jets[i].Eta())<2.4;
	  if(!btagcond){
	    if(fabs(standardjets[i].vect.Eta())>maxEta){
	      maxEta=fabs(standardjets[i].vect.Eta());
	      jetprime_3j1t=standardjets[i].vect;
	      *etajprime_3j1t=maxEta;
	      *jprimeflavour_3j1t=standardjets[i].partonFlavour;
	    }
	    if(fabs(standardjets[i].vect.Eta())<minEta){
	      minEta=fabs(standardjets[i].vect.Eta());
	      jetnobnoprime=standardjets[i];
	      *leadingextrajetflavour_3j1t=standardjets[i].partonFlavour;
	    }
	  }
	  else{
	    bjet=standardjets[i].vect;
	    *bjetflavour_3j1t = standardjets[i].partonFlavour; 
	    *bjetpt_3j1t = standardjets[i].vect.Pt(); 
	    *bjeteta_3j1t = standardjets[i].vect.Eta(); 
	  }
	}}
      else{
	for(size_t i=0; i<(size_t)standardjets.size(); ++i){
	  cout <<"passeB"<< standardjets[i].passesB << "eta jet" << fabs(jets[i].Eta()) << endl;
	}}
      jetstopMass.push_back(bjet);
      jetstopMass.push_back(jetnobnoprime.vect);
      
      *mlb_3j1t=topUtils.mlj(tightLep.at(0), bjet);
      *mljprime_3j1t=topUtils.mlj(tightLep.at(0), jetprime_3j1t);
      *mljextra_3j1t=topUtils.mlj(tightLep.at(0), jetnobnoprime.vect);
      
      vector<TLorentzVector> tops = topUtils.top4Momenta(tightLep, jetstopMass, metpx, metpy);
      *topMass_3j1t=tops.at(0).M();
      *topMt_3j1t=topUtils.topMtw(tops.at(0), bjet, metpx, metpy);
      *topPt_3j1t=tops.at(0).Pt();
      *topEta_3j1t=tops.at(0).Eta();
      *topY_3j1t=tops.at(0).Y();
      *costhetapol_3j1t = topUtils.costhetapol(tightLep.at(0), jetprime_3j1t, tops.at(0));
      *costhetael_3j1t = topUtils.costhetael(tops.at(0), tightLep.at(0), metpx, metpy);
      
      TVector2 metv( metpx, metpy);
      double Mt2w = Mt2cal->calculateMT2w(jets_nob20,bjets20,tightLep.at(0), metv,"MT2w"); 
      *mt2w_3j1t = Mt2w;

      float weightreshape=jetnobnoprime.csvreweight;
      bool fillcond=true;
      if(doTtoSDDecay){
	fillcond = (!(jetnobnoprime.partonFlavour==5
		      && (jetnobnoprime.partonFlavour*topcharge>0)
		      && fabs(jetnobnoprime.vect.Eta())<2.4));		 
      }
      *leadingextrajetcsvweight_3j1t=jetnobnoprime.csvreweight;
      *leadingextrajetcsv_3j1t=jetnobnoprime.csv;
      if(fillcond){
	*leadingextrajetcsvweight_sd_3j1t=0;
      }
      else{
	float weight_sd_b=jetnobnoprime.csvreweightsd;
	*leadingextrajetcsvweight_sd_3j1t=weight_sd_b;
      }
      *leadingextrajetpt_3j1t=jetnobnoprime.vect.Pt();
      *leadingextrajeteta_3j1t=jetnobnoprime.vect.Eta();
      
      *topMassExtra_3j1t=tops.at(1).M();
      *topMtExtra_3j1t=topUtils.topMtw(tops.at(1), jetnobnoprime.vect, metpx, metpy);
      *topPtExtra_3j1t=tops.at(1).Pt();
      *topEtaExtra_3j1t=tops.at(1).Eta();
      *topYExtra_3j1t=tops.at(1).Y();
      *costhetapolExtra_3j1t = topUtils.costhetapol(tightLep.at(0), jetprime_3j1t, tops.at(1));
      *costhetaelExtra_3j1t = topUtils.costhetael(tops.at(1), tightLep.at(0), metpx, metpy);
      
      
      if(addTrees)syst1BM.fillTreesSysts(trees1T,"3j1t");
      //AOB for 3j1t?
    }    
    //End 3j1t
    
    //----------------------------------------------------------------------------------    
    //----------------------------------------------------------------------------------    
    
    //3j2t
    if((jets.size() == 3 && bjets.size()==2)){
      for(size_t j= 0; j< (size_t)jets.size();++j ){
	if(j==0){
	  syst2BM.fillHistogramsSysts(h_3j2t_dR_lepjetpt40_1st,deltaR(jets[j].Eta(),jets[j].Phi(),tightLep[0].Eta(),tightLep[0].Phi()),w);
	  syst2BM.fillHistogramsSysts(h_3j2t_dPhi_lepjetpt40_1st,fabs(deltaPhi(jets[j].Phi(),tightLep[0].Phi())),w);
	  syst2BM.fillHistogramsSysts(h_3j2t_dEta_lepjetpt40_1st,jets[j].Eta()-tightLep[0].Eta(),w);
	  syst2BM.fillHistogramsSysts(h_3j2t_jetpt40_1st,jets[j].Pt(),w);
	}
	if(j==1){syst2BM.fillHistogramsSysts(h_3j2t_jetpt40_2nd,jets[j].Pt(),w); }
	if(j==2){syst2BM.fillHistogramsSysts(h_3j2t_jetpt40_3rd,jets[j].Pt(),w); }
      }
      for(size_t b=0; b< (size_t)bjets.size();++b ){
	if(b==0)syst2BM.fillHistogramsSysts(h_3j2t_bjetpt,bjets[b].Pt(),w,NULL,false);
	if(b==1)syst2BM.fillHistogramsSysts(h_3j2t_2ndbjetpt,bjets[b].Pt(),w);
      }
      for(size_t i = 0; i < (size_t)tightLep.size();++i ){
	if(tightLep[i].Pt()>20){
	  syst2BM.fillHistogramsSysts(h_3j2t_MuPt,tightLep[i].Pt(),w);
	  syst2BM.fillHistogramsSysts(h_3j2t_MuEta,tightLep[i].Eta(),w);
	  syst2BM.fillHistogramsSysts(h_3j2t_MuPhi,tightLep[i].Phi(),w);
	  syst2BM.fillHistogramsSysts(h_3j2t_MuE,tightLep[i].E(),w);
	  if((tightLep.size()) == 1 ){
	    TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
	    float phi_lmet = fabs(deltaPhi(tightLep[i].Phi(), metPhi[0]) );
	    mt = sqrt(2* tightLep[i].Pt() * met* ( 1- cos(phi_lmet)));
	    syst2BM.fillHistogramsSysts(h_3j2t_mtw,mt,w);
	    *mtw_3j2t=mt;
	    if(channel=="muon" || channel == "muonantiiso")syst2BM.fillHistogramsSysts(h_3j2t_muIso,selectedIso[i],w);
	    if(mt>50){
	      syst2BM.fillHistogramsSysts(h_3j2t_mtwcut_MuPt,tightLep[i].Pt(),w);
	      if(channel=="muon" || channel == "muonantiiso")syst2BM.fillHistogramsSysts(h_3j2t_mtwcut_muIso,selectedIso[i],w);
	    }
	  }
	}
      }
      TLorentzVector jetprime_3j2t;
      for(size_t j=0; j< (size_t)jets.size();++j ){
	bool btagcond = standardjets[j].passesB>0. && fabs(jets[j].Eta())<2.4;
	
	if(!btagcond){
	  syst2BM.fillHistogramsSysts(h_3j2t_ljetpt,(jets.at(j).Pt()),w);
	  syst2BM.fillHistogramsSysts(h_3j2t_ljeteta,(jets.at(j).Eta()),w);
	  syst2BM.fillHistogramsSysts(h_3j2t_etajprime,(fabs(jets.at(j).Eta())),w);
	  *etajprime_3j2t=fabs(jets.at(j).Eta());
	  jetprime_3j2t=jets.at(j);
	  
	}
      }
      
      *mljprime_3j2t=topUtils.mlj(tightLep.at(0), jetprime_3j2t);
      *mlbleading_3j2t=topUtils.mlj(tightLep.at(0), bjets.at(0));
      *mlbsecond_3j2t=topUtils.mlj(tightLep.at(0), bjets.at(1));
      *deltaEtabb_3j2t=fabs(bjets.at(0).Eta()-bjets.at(1).Eta());
      
      vector<TLorentzVector> tops= topUtils.top4Momenta(tightLep, bjets, metpx, metpy);
      if(tops.size()!=0){
	double massLeading = tops.at(0).M();
	double massSecond = tops.at(1).M();
	*topMassLeading_3j2t=massLeading;
	*topMtLeading_3j2t=topUtils.topMtw(tops.at(0), bjets.at(0), metpx, metpy);
	*topPtLeading_3j2t=tops.at(0).Pt();
	*topEtaLeading_3j2t= tops.at(0).Eta();
	*topYLeading_3j2t= tops.at(0).Y();
	*costhetapolLeading_3j2t = topUtils.costhetapol(tightLep.at(0), jetprime_3j2t, tops.at(0));	
	*costhetaelLeading_3j2t = topUtils.costhetael(tops.at(0), tightLep.at(0), metpx, metpy);	
	
	*topMassSecond_3j2t=massSecond;
	*topMtSecond_3j2t=topUtils.topMtw(tops.at(1), bjets.at(1), metpx, metpy);
	*topPtSecond_3j2t=tops.at(1).Pt();
	*topEtaSecond_3j2t= tops.at(1).Eta();
	*topYSecond_3j2t= tops.at(1).Y();
	*costhetapolSecond_3j2t = topUtils.costhetapol(tightLep.at(0), jetprime_3j2t, tops.at(1));
	*costhetaelSecond_3j2t = topUtils.costhetael(tops.at(1), tightLep.at(0), metpx, metpy);	
      }
      syst2BM.fillHistogramsSysts(h_3j2t_nextrajets,extrajets.size(),w);
      *nextrajets_3j2t= extrajets.size();
      TVector2 met( metpx, metpy);
      double Mt2w = Mt2cal->calculateMT2w(jets_nob20,bjets20,tightLep.at(0), met,"MT2w"); 
     *mt2w_3j2t = Mt2w;
      
      *w_3j2t=w;
      if(addTrees)syst2BM.fillTreesSysts(trees2T,"3j2t");
      //AOB for 3j2t?
    }
    //End 3j2t

    //----------------------------------------------------------------------------------    
    //----------------------------------------------------------------------------------    
    
    //To add here: 2l3j2t: 2 leptons, 3 jets, 2 tags for ttbar control
    
    }
    //End signal and control regions.
    
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
   //   cout << " mean  "<< zeroWeight<<endl;
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

  //Rescale the extra jet csv by the average of the sf
  systZero.rescaleHistograms(h_2j1t_leadingextrajetcsv_sd_b,h_2j1t_leadingextrajetcsv_reshape_sd_b);
  systZero.addHistograms(h_2j1t_leadingextrajetcsv,h_2j1t_leadingextrajetcsv_sd_b);

  systZero.rescaleHistograms(h_2j1t_mtwcut_leadingextrajetcsv_sd_b,h_2j1t_mtwcut_leadingextrajetcsv_reshape_sd_b);
  systZero.addHistograms(h_2j1t_mtwcut_leadingextrajetcsv,h_2j1t_mtwcut_leadingextrajetcsv_sd_b);

  systZero.rescaleHistograms(h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b,h_2j1t_mtwcut_sr_leadingextrajetcsv_reshape_sd_b);
  systZero.addHistograms(h_2j1t_mtwcut_sr_leadingextrajetcsv,h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b);


  //Write the Histogramms here  
  systZero.writeSingleHistogramSysts(h_cutFlow, allMyFiles); 
  systZero.writeSingleHistogramSysts(h_weight_sign, allMyFiles); 
  systZero.writeHistogramsSysts(h_nPV, allMyFiles); 
  systZero.writeHistogramsSysts(h_nGoodPV, allMyFiles); 
  systZero.writeHistogramsSysts(h_nTruePV, allMyFiles); 
  
  systZero.writeHistogramsSysts(h_2j0t_dR_lepjetpt40_1st,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_dPhi_lepjetpt40_1st,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_dEta_lepjetpt40_1st,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_MuEta,      allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_MuPhi,      allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_MuPt,      allMyFiles);
  
  systZero.writeHistogramsSysts(h_2j0t_jetpt40_1st, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_jetpt40_2nd, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtw,   allMyFiles); 

  systZero.writeHistogramsSysts(h_2j0t_mtwcut_MuPt,      allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_jeteta40_1st,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_jetabseta40_1st,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_jetpt40_1st,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_jetpt40_2nd,   allMyFiles); 
  
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_incljeteta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_incljetabseta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_incljetabseta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_muIso,   allMyFiles);
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_muIso,   allMyFiles);
  
  //2j1t
  systZero.writeHistogramsSysts(h_2j1t_bjetpt,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_ljetpt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_ljeteta,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_etajprime,  allMyFiles); 
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
  systZero.writeHistogramsSysts(h_2j1t_dR_lepjetpt40_1st,    allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_dPhi_lepjetpt40_1st,  allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_dEta_lepjetpt40_1st,  allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_muIso,   allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_muIso,   allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_jetpt40_1st,allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_jetpt40_2nd,allMyFiles);

  if(doTtoSDDecay){    systZero.writeHistogramsSysts(h_2j1t_leadingextrajetpt_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_leadingextrajetcsv_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_leadingextrajetcsv_reshape_sd_b,   allMyFiles);   }


  //BDT Histos  
  systZero.writeHistogramsSysts(h_2j1t_BDT_STvsQCDMu, allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_BDT_STvsVJ, allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_BDT_STvsTT, allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_BDT_STsdvsVJ, allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_BDT_STsdvsST, allMyFiles);
  
   //mtw + eta < 2.5
   systZero.writeHistogramsSysts(h_2j1t_mtwcut_le_sr_BDT_STvsVJ, allMyFiles);
   systZero.writeHistogramsSysts(h_2j1t_mtwcut_le_sr_BDT_STvsTT, allMyFiles);
   systZero.writeHistogramsSysts(h_2j1t_mtwcut_le_sr_BDT_STsdvsVJ, allMyFiles);
   systZero.writeHistogramsSysts(h_2j1t_mtwcut_le_sr_BDT_STsdvsST, allMyFiles);
 
   //mtw + eta > 2.5
   systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_BDT_STvsVJ, allMyFiles);
   systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_BDT_STvsTT, allMyFiles);
   systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_BDT_STsdvsVJ, allMyFiles);
   systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_BDT_STsdvsST, allMyFiles);


  systZero.writeHistogramsSysts(h_2j1t_dR_lepjetpt40_1st,    allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_dPhi_lepjetpt40_1st,  allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_dEta_lepjetpt40_1st,  allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_muIso,   allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_muIso,   allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_MuPt,    allMyFiles);
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_MuPt,   allMyFiles);
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
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_topMassExtra,   allMyFiles); 
  if(doTtoSDDecay){    systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_topMassExtra_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv_reshape_sd_b,   allMyFiles);   }
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_topMass,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_topPt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_nextrajets,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra,   allMyFiles); 
    if(doTtoSDDecay){
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt_sd_b,   allMyFiles); //Part for the reweighted csv
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv_reshape_sd_b,   allMyFiles);   }

  // Other control region plots  
  systZero.writeHistogramsSysts(h_2j1t_mtw_le50,     allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_le_sr_etajprime, allMyFiles); 
  

  //3j1t
  systZero.writeHistogramsSysts(h_3j1t_bjetpt,  allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_MuPt,    allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_MuEta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_MuPhi,   allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_MuE,     allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_mtw,     allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_dR_lepjetpt40_1st,    allMyFiles);
  systZero.writeHistogramsSysts(h_3j1t_dPhi_lepjetpt40_1st,  allMyFiles);
  systZero.writeHistogramsSysts(h_3j1t_dEta_lepjetpt40_1st,  allMyFiles);
  systZero.writeHistogramsSysts(h_3j1t_muIso,   allMyFiles);  
  systZero.writeHistogramsSysts(h_3j1t_mtwcut_MuPt, allMyFiles);
  systZero.writeHistogramsSysts(h_3j1t_mtwcut_muIso,   allMyFiles);  
  systZero.writeHistogramsSysts(h_3j1t_jetpt40_1st, allMyFiles);
  systZero.writeHistogramsSysts(h_3j1t_jetpt40_2nd, allMyFiles);
  systZero.writeHistogramsSysts(h_3j1t_jetpt40_3rd, allMyFiles);

  //3j2t
  systZero.writeHistogramsSysts(h_3j2t_bjetpt,  allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_2ndbjetpt, allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuPt,    allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuEta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuPhi,   allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuE,     allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuCharge,allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_mtw,     allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_muIso,     allMyFiles);
  systZero.writeHistogramsSysts(h_3j2t_dR_lepjetpt40_1st,    allMyFiles);
  systZero.writeHistogramsSysts(h_3j2t_dPhi_lepjetpt40_1st,  allMyFiles);
  systZero.writeHistogramsSysts(h_3j2t_dEta_lepjetpt40_1st,  allMyFiles);
  systZero.writeHistogramsSysts(h_3j2t_mtwcut_MuPt,    allMyFiles);
  systZero.writeHistogramsSysts(h_3j2t_mtwcut_muIso,     allMyFiles);
  systZero.writeHistogramsSysts(h_3j2t_ljetpt,  allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_ljeteta,allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_etajprime,allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_nextrajets,allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_jetpt40_1st, allMyFiles);
  systZero.writeHistogramsSysts(h_3j2t_jetpt40_2nd, allMyFiles);
  systZero.writeHistogramsSysts(h_3j2t_jetpt40_3rd, allMyFiles);
  
  systZero.writeHistogramsSysts(h_nJets, allMyFiles); 
  systZero.writeHistogramsSysts(h_nbJets, allMyFiles);
  
  if(addTrees){
    syst1BM.writeTreesSysts(trees1T,outTreeFile);
    syst2BM.writeTreesSysts(trees2T,outTreeFile);
  }
  
  systZero.closeFilesSysts(allMyFiles);
  
  std::cout<< "Info: Total No. of events for sample["<<sample<<"] : "<<nEvents<<std::endl;
}

