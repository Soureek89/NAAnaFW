#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include <sstream>
//#include ""

using namespace std;

std::string Convert (float number){
  std::ostringstream buff;
  buff<<number;
  return buff.str();   
}

void macro_2Dto1D(TString channel, TString njmt, TString syst, TString region, Float_t cut1x, Float_t cut2x, Float_t cut1y, Float_t cut2y, TString lep){
  TString infilename, outfilename, histo1Dname, histo2Dname;
  if((lep=="muonantiiso"||lep=="electronantiiso")&& (channel=="QCDMuPt20toInf"|| channel=="TT_sd")) return;
  if(lep=="electron" && channel=="QCDMuPt20toInf") return;
  if((channel.Contains("hdamp")||channel.Contains("psq2"))&&((lep=="muonantiiso"||lep=="electronantiiso")||(syst!=""))) return;
  infilename = "Plot/histo2D_"+lep+".root";  
  histo1Dname = "h1D_"+njmt+"_"+region;
  if(syst.EqualTo("")){
    histo2Dname = "h2D_"+njmt+"_"+channel+"_"+region;
    cout << histo2Dname << endl;
    outfilename = "Plot/"+lep+"/"+channel+"_"+lep+".root";
  }
  else if(syst.EqualTo("jesUp")||syst.EqualTo("jesDown")||syst.EqualTo("jerUp")||syst.EqualTo("jerDown")){
    histo2Dname = "h2D_"+njmt+"_"+channel+"_"+region+"_"+syst;
    outfilename = "Plot/"+lep+"/"+channel+"_"+lep+"_"+syst+".root";
  }
  else{
    histo2Dname = "h2D_"+njmt+"_"+channel+"_"+region+"_"+syst;
    outfilename = "Plot/"+lep+"/"+channel+"_"+lep+"_"+syst+".root";
  }
  TFile * fin = TFile::Open(infilename);
  TH2F * h2 = (TH2F*)fin->Get(histo2Dname);
  cout << " opening file "<< infilename<< endl;
  Int_t nbins=10;
  Int_t xmin=-1;
  Int_t xmax=1;
  Int_t ymin=-1; 
  Int_t ymax=1;
  Int_t bin1x, bin2x, bin1y, bin2y, global_bin;
  Float_t h1content=0; 
  bin1x = (1+cut1x)*nbins/(xmax-xmin)+1;
  //  cout << bin1x << endl;
  bin2x = (1+cut2x)*nbins/(xmax-xmin);
  //  cout << bin2x << endl;
  bin1y = (1+cut1y)*nbins/(ymax-ymin)+1;
  //  cout << bin1y << endl;
  bin2y = (1+cut2y)*nbins/(ymax-ymin);
  //  cout << bin2y << endl;
  Int_t n = (bin2x-bin1x+1)*(bin2y-bin1y+1);
  //  cout << n <<endl;
  //h1->SetBins(n,1,n+1);
  //  h1->SetName(histo1Dname+code.c_str());

  string code = Convert(bin1x)+Convert(bin2x)+Convert(bin1y)+Convert(bin2y);
  //  cout << code << endl;
  TH1F * h1 = new TH1F(histo1Dname+"_"+code.c_str(),histo1Dname+"_"+code.c_str(),n,1,n+1);
  //  cout << h2->GetName() << endl;
  Int_t a=0;
  for(Int_t i=bin1y; i<=bin2y; i++){
    for(Int_t j=bin1x; j<=bin2x; j++){
      a++;
      global_bin = h2->GetBin(j,i);
      //      cout << " gbin is "<< global_bin<<endl;
      
      h1content = h2->GetBinContent(global_bin);
      //      cout << i << "_" << j << "_" << h1content <<endl;
      h1->SetBinContent(a, h1content);
    }
  }
  TFile * fout = TFile::Open(outfilename,"UPDATE");
  h1->Write();
  cout << h1->GetName() << " integral is:" << h1->Integral() << endl;
  fin->Close();
  fout->Close();
}
