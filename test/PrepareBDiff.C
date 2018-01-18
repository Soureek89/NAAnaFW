void PrepareBDiff(TString dir="./",TString treename="DMTreesDumper/ttDM__noSyst"){
//PrepareBDiff(TString dir,TString treename){


  TFile * file0 = TFile::Open("testErrors.root","UPDATE");
  int nbins=100;
  float lowend=-1.,highend=1.;
  TH1D nom("nom","nom",nbins,lowend,highend);
  TH1D den("den","den",nbins,lowend,highend);

  TString cut ="",cutname="";
  //cut = "&&(jetsAK4CHSTight_CorrPt<60 )"; cutname= "PtLE60";
  //cut = "&&(jetsAK4CHSTight_CorrPt>60 )"; cutname= "PtGE60";
  //  cut = "&&(jetsAK4CHSTight_CorrPt>80 )"; cutname= "PtGE80";
  TH1D res("res"+cutname,"res"+cutname,nbins,lowend,highend);
  
  TChain chain(treename);
  chain.Add(dir+"/*");
  chain.Project("nom","jetsAK4CHSTight_CMVAv2","(abs(jetsAK4CHSTight_PartonFlavour)<4 && abs(jetsAK4CHSTight_Eta)<2.4"+cut+")");
  chain.Project("den","jetsAK4CHSTight_CMVAv2","jetsAK4CHSTight_reshapeFactorCMVA_SD*(jetsAK4CHSTight_PartonFlavour==5 && abs(jetsAK4CHSTight_Eta)<2.4"+cut+")");
  //nom.DrawNormalized();
  //den.DrawNormalized("samee");
  nom.Scale(1/nom.Integral());
  den.Scale(1/den.Integral());
  for(int B=1;B<=nbins;++B){
    if(den.GetBinContent(B)!=0){
      float re=(nom.GetBinContent(B)/den.GetBinContent(B));
      cout << "bin " << B<< " re " <<re<<endl;
      res.SetBinContent(B,re);
    }
    else res.SetBinContent(B,0.);
  }
  nom.Write();
  den.Write();
  res.Draw();
  res.Write();
}
