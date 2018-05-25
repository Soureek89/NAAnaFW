void trainTMVA(string inputname, string cfgfile){

  TString postfix_bdt = inputname;
  TString filename = "tmva"+postfix_bdt+".root";
  cout <<  postfix_bdt<<" filename " << filename <<endl;
    
  TFile outFile(filename, "RECREATE");
  TMVA::Factory factory("TMVAMu"+postfix_bdt, &outFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
  int nTrainSignal=0, nTestSignal=0, nTrainBackground=0, nTestBackground=0;
  ifstream varCfgFile;

  varCfgFile.open(cfgfile, ifstream::in);
  
  
  string bas, option="";
  TString weightexpr = "";
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
       factory.AddVariable(bas);
       cout << "variable: " << bas << endl;
       continue;
    }
    TString weightsig="",weightbkg=""; 
    if(option == "[signalFiles]" ) {
       const TString file = bas.substr(0, bas.find(" "));
       bas = bas.substr(bas.find(" ")+1, bas.length());
       const TString tree = bas.substr(0, bas.find(" "));
       bas = bas.substr(bas.find(" ")+1, bas.length());
       const TString type = bas.substr(0, bas.find(" "));
       bas = bas.substr(bas.find(" ")+1, bas.length());
       TString weight = bas.substr(0, bas.length());
       if(weight == "") weight="1";
       cout << " sig weight "<< weight <<endl;
       TFile *theFile = new TFile(file);
       TTree *theTree = (TTree*) theFile->Get(tree);

       cout << " treee "<< tree << " file "<< theTree << " file "<<  file <<" theFile "<< theFile << endl;
       factory.AddSignalTree ( theTree, weight.Atof(), type);
       if(type == "Training") nTrainSignal += theTree->GetEntries();
       if(type == "Test") nTestSignal += theTree->GetEntries();
       continue;
    }
    if(option == "[bgFiles]" ) {
       const TString file = bas.substr(0, bas.find(" "));
       bas = bas.substr(bas.find(" ")+1, bas.length());
       const TString tree = bas.substr(0, bas.find(" "));
       bas = bas.substr(bas.find(" ")+1, bas.length());
       const TString type = bas.substr(0, bas.find(" "));
       bas = bas.substr(bas.find(" ")+1, bas.length());
       TString weight = bas.substr(0, bas.length());
       if(weight == "") weight="1";
       TFile *theFile = new TFile(file);
       TTree *theTree = (TTree*) theFile->Get(tree);
       cout << " bkg weight "<< weight <<endl;
       factory.AddBackgroundTree ( theTree, weight.Atof(), type);
       if(type == "Training") nTrainBackground += theTree->GetEntries();
       if(type == "Test") nTestBackground += theTree->GetEntries();
       continue;
    }
    if(option == "[weightExpr]" ) {
      weightexpr = bas.substr(0, bas.find(" "));
    }
  }
  varCfgFile.close();
  cout << ">> The config file has been read and loaded." << endl;
  //  TCut sigCut = "isSignal > 0.5";
  //  TCut bkgCut = "isSignal < 0.5";
  //  factory.SetInputTrees(sampleTree, sigCut, bkgCut);
  
  
  
  
  //  factory.SetWeightExpression("1.0");  
  factory.SetWeightExpression("w_nominal");  
  cout << "weight expression is "<<weightexpr<<endl;
  if(weightexpr!="")factory.SetWeightExpression(weightexpr);  


  TCut preselection_sig = "";//weightcut;//"nextrajets<0.10"; 
  TCut preselection_bkg = "";//weightcut;//"nextrajets<0.10";

  //  TCut preselection_sig_nojets = "nextrajets<0.10"; 
  //  TCut preselection_bkg_nojets = "nextrajets<0.10";


  factory.PrepareTrainingAndTestTree(preselection_sig, preselection_bkg, TString("nTrain_Signal=nTrainSignal:nTrain_Background=nTrainBackground:nTest_Signal=nTestSignal:nTest_Background=nTestBackground:SplitMode=Block:VerboseLevel=Info"));

  factory.BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=800:MinNodeSize=0.05:BoostType=AdaBoost:SeparationType=GiniIndex:PruneMethod=CostComplexity:MaxDepth=3:PruningValFraction=0.3:PruneStrength=-1:UseRandomisedTrees=True");

  //Normal settings
  //factory.BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=800:MinNodeSize=0.05:BoostType=AdaBoost:SeparationType=GiniIndex:PruneMethod=CostComplexity:MaxDepth=3:PruningValFraction=0.3:PruneStrength=-1");
  
  //Adding UseRandomisedTrees=True
  //factory.BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=800:MinNodeSize=0.05:BoostType=AdaBoost:SeparationType=GiniIndex:PruneMethod=CostComplexity:MaxDepth=3:PruningValFraction=0.3:PruneStrength=-1:UseRandomisedTrees=True");
  
  //Adding AdaBoostBeta=0.5
  //factory.BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=800:MinNodeSize=0.05:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:PruneMethod=CostComplexity:MaxDepth=3:PruningValFraction=0.3:PruneStrength=-1:UseRandomisedTrees=True");
  
  //Adding UseNvars=5 default value 2
  //factory.BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=800:MinNodeSize=0.05:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:PruneMethod=CostComplexity:MaxDepth=3:PruningValFraction=0.3:PruneStrength=-1:UseRandomisedTrees=True:UseNvars=5");
 
  factory.TrainAllMethods();
  factory.TestAllMethods();
  factory.EvaluateAllMethods();

  /*
  outFile.Close();
  if (!gROOT->IsBatch()) {
    gROOT ->Execute(".L TMVAGui.C"); 
    TMVAGui( outFile );
  }
  */
}

