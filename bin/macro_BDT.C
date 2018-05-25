#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>


using namespace std;
void macro_BDT(TString channel, TString njmt, TString training_name, TString syst, TString lep)
{
  if((lep=="electron"||lep=="electronantiiso")&& channel=="QCDMuPt20toInf") return;
  if((lep=="muonantiiso"||lep=="electronantiiso")&& channel=="TT_sd") return;
  if((channel.Contains("hdamp")||channel.Contains("psq2"))&&((lep=="muonantiiso"||lep=="electronantiiso")||(syst.EqualTo("jesUp")||syst.EqualTo("jesDown")||(syst.EqualTo("jerUp")||syst.EqualTo("jerDown"))))) return;
  TFile * f = TFile::Open("trees_lumi/"+lep+"/trees_"+channel+"_"+lep+".root","UPDATE");
  TString treename;
  if(syst.EqualTo("")) treename = "events_"+njmt;
  if(syst.EqualTo("jesUp")||syst.EqualTo("jesDown")||syst.EqualTo("jerUp")||syst.EqualTo("jerDown")) treename = "events_"+njmt+"_"+syst;
  TTree * t =(TTree*)f->Get(treename);
  Float_t BDT;
  TString branch_name = "BDT_"+training_name;
  TBranch * newbranch = t->Branch(branch_name, &BDT, branch_name+"/F");
  TMVA::Reader * Reader = new TMVA::Reader(training_name);
  map<string, float> namesToVars;
  vector<string> variables_names;
  vector<float> variables; 
  ifstream varCfgFile;
  varCfgFile.open("MVA/cfg"+training_name+".txt",ifstream::in);
  string bas,option;
  while(std::getline(varCfgFile, bas))       {
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
      variables_names.push_back(bas);
      namesToVars[bas]=-99999.;
      variables.push_back(namesToVars[bas]);
      cout << "variable: " << bas << endl;
      continue;
    }
  }
  varCfgFile.close();
  //cout << "variables size is " << variables.size() << endl;
  for(Int_t i=0; i<variables.size(); i++) {
    t->SetBranchAddress(variables_names[i].c_str(), &(variables[i]));
    Reader->AddVariable(variables_names[i], &(variables[i]));
  }
  const std::string weightsfile = std::string("MVA/weights/TMVAMu"+training_name+"_BDT.weights.xml");
  Reader->BookMVA("BDTG", weightsfile.c_str());
  
  Int_t nentries = t->GetEntries();
  for(Int_t i=0; i<nentries; i++){
  //  for(Int_t i=0; i<30000; i++){
    if(i%100000==1) cout << " event # " << i << " completion: " << (float)i/nentries*100. << endl;
    t->GetEntry(i);
    BDT = Reader->EvaluateMVA("BDTG");
    newbranch->Fill();
  }  
  f->cd();
  t->Write();
  delete f;
}
