void PrepareBDiff(TString dir="./",TString treename="DMTreesDumper/ttDM__noSyst"){
//PrepareBDiff(TString dir,TString treename){

  bool ext_file = true;
  ifstream samples_txt;
  string infile = "/afs/cern.ch/work/a/adeiorio/CMSSW_8_0_26_patch1/src/Analysis/NAAnaFW/bin/files/renamed/ST_T_tch_sd.txt";
  //  string infile = "/afs/cern.ch/work/a/adeiorio/CMSSW_8_0_26_patch1/src/Analysis/NAAnaFW/bin/files/renamed/ST_T_tch_sd_5.txt";
  string bas;
  vector<string> signal_files;

  if(ext_file){
    samples_txt.open(infile,ifstream::in);
    while(std::getline(samples_txt, bas))       {
      if( bas.empty() ) continue;
      if( bas.substr(0,4) == "root" ) {
	signal_files.push_back(bas);
	//      cout << "signal is: " << bas << endl;
	continue;
      }
    }
  }

  TFile * file0 = TFile::Open("testErrors.root","UPDATE");
  TChain *chain = new TChain(treename);

  if(!ext_file){
    dir = "/afs/cern.ch/work/w/wajid/public/xWajid/sdsamples";
    chain->Add(dir+"/*");
  }

  if(ext_file){
    for(int i=0; i < signal_files.size(); i++){
      chain->Add(signal_files[i].c_str());
      cout << signal_files[i].c_str() << endl;
    }
  }
  int nbins=100;
  float lowend=-1.,highend=1.;
  TH1F *nom = new TH1F("nom","nom",nbins,lowend,highend);
  TH1F *den = new TH1F("den","den",nbins,lowend,highend);

  TString cut ="",cutname="";
  //cut = "&&(jetsAK4CHSTight_CorrPt<60 )"; cutname= "PtLE60";
  TH1F *res = new TH1F("res"+cutname,"res"+cutname,nbins,lowend,highend);
  
  nom->Sumw2();
  den->Sumw2();
  res->Sumw2();
  chain->Project("nom","jetsAK4CHSTight_CMVAv2","(abs(jetsAK4CHSTight_PartonFlavour)<4 && abs(jetsAK4CHSTight_Eta)<2.4"+cut+")");
  chain->Project("den","jetsAK4CHSTight_CMVAv2","jetsAK4CHSTight_reshapeFactorCMVA_SD*(jetsAK4CHSTight_PartonFlavour==5 && abs(jetsAK4CHSTight_Eta)<2.4"+cut+")");
  //nom->DrawNormalized();
  //den->DrawNormalized("samee");
  nom->Scale(1/nom->Integral());
  den->Scale(1/den->Integral());
  res = (TH1F*)nom->Clone("res"+cutname);
  res->Divide(den);
  nom->Write();
  den->Write();
  res->Draw();
  res->Write();
}
