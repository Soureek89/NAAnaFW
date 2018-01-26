#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"

using namespace std;
void macro_plot(TString channel, TString variabile, TString title, TString peso, int nbins, int xmin, int xmax, TString njmt, TString syst, TString cut_tag, TString lep) {
  if((lep=="electron"||lep=="electronantiiso")&& channel=="QCDMuPt20toInf") return;
  if((lep=="muonantiiso"||lep=="electronantiiso")&& channel=="TT_sd") return;
  if((channel.Contains("hdamp")||channel.Contains("psq2"))&&((lep=="muonantiiso"||lep=="electronantiiso")||(syst!=""))) return;
  TH1F *h1_uncut= new TH1F;
  TFile *f1= TFile::Open("trees_lumi/"+lep+"/trees_"+channel+"_"+lep+".root");
  TFile *f1_uncut;
  TString histoname;
  TString treename = "events_"+njmt;
  TString taglio, foutput;
  if(cut_tag.EqualTo(""))    histoname = "h_"+njmt+"_"+variabile;
  else    histoname = "h_"+njmt+"_"+variabile+"_"+cut_tag;
  TH1F *h1= new TH1F(histoname, variabile+"_"+channel, nbins, xmin, xmax);
  if(syst.EqualTo("jesUp") || syst.EqualTo("jesDown") || syst.EqualTo("jerUp") || syst.EqualTo("jerDown")){
    treename = "events_"+njmt+"_"+syst;
    taglio ="w_nominal*"+peso;
    foutput = "Plot/"+lep+"/"+channel+"_"+lep+"_"+syst+".root"; 
  }
  else if(syst.EqualTo("")){  
    taglio ="w_nominal*"+peso;
    foutput = "Plot/"+lep+"/"+channel+"_"+lep+".root";
  }
  else{
    taglio = syst+"*"+peso;
    foutput = "Plot/"+lep+"/"+channel+"_"+lep+"_"+syst+".root"; 
  }
  //cout << taglio << endl;
  ((TTree*)f1->Get(treename))->Project(histoname,variabile,taglio);
  //cout << h1->Integral() << endl;
  TFile *fout = TFile::Open(foutput, "UPDATE"); 
  
  fout->cd();
  /*  TLegend *leg= new TLegend(0.75,0.57,0.90,0.90);
  gStyle->SetOptStat(0);
  leg->AddEntry(h1,(variabile+channel),"l");
  h1->GetXaxis()->SetTitle(title);
  h1->GetYaxis()->SetTitle("Events/bin");
  h1->Draw();
  leg->Draw("SAME");*/
  h1->Write();
  fout->Close();
  f1->Close();
}
