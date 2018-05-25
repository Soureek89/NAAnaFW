#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
//#include "TTree.h"

using namespace std;
void macro_2Dhisto(TString channel, TString njmt, TString syst, TString region, TString BDT1_name, TString BDT2_name, TString lep){
  if((lep=="electron"||lep=="electronantiiso")&& channel=="QCDMuPt20toInf") return;
  if((lep=="muonantiiso"||lep=="electronantiiso") && channel=="TT_sd") return;
  if((channel.Contains("hdamp")||channel.Contains("psq2"))&&((lep=="muonantiiso"||lep=="electronantiiso")||(syst!=""))) return;
  TFile * f=TFile::Open("trees_lumi/"+lep+"/trees_"+channel+"_"+lep+".root");
  TString treename, histoname;
  TH2F * h2 = new TH2F("h2", "h2",10,-1,1,10,-1,1);
  //  Float_t BDT1, BDT2, etajprime, mtw, w, w_nominal, var_syst; 
  gStyle->SetOptStat(0);
  TString outfilename, taglio, selection;
  selection = "BDT_"+BDT1_name+":BDT_"+BDT2_name;
  if(syst.EqualTo("jesUp")||syst.EqualTo("jesDown")||syst.EqualTo("jerUp")||syst.EqualTo("jerDown")) treename = "events_"+njmt+"_"+syst;
  else treename = "events_"+njmt;
  TTree * t =(TTree*)f->Get(treename);
  if(region=="cr") taglio = "(mtw>50 && etajprime<2.4";
  else taglio = "(mtw>50 && etajprime>2.4";
  if(lep=="electronantiiso") taglio += " &&mlb>30";
  if(syst.EqualTo("")){
    histoname = "h2D_"+njmt+"_"+channel+"_"+region;
    t->Project("h2",selection,taglio+")*w*w_nominal");
  }
  else{
    if(syst.EqualTo("jesUp")||syst.EqualTo("jesDown")||syst.EqualTo("jerUp")||syst.EqualTo("jerDown")){
      histoname = "h2D_"+njmt+"_"+channel+"_"+region+"_"+syst;
      t->Project("h2",selection,taglio+")*w*w_nominal");
    }
    else{
      histoname = "h2D_"+njmt+"_"+channel+"_"+region+"_"+syst;
      t->Project("h2",selection,taglio+")*w*"+syst);
    }
  }
  h2->SetName(histoname);
  TCanvas * c = new TCanvas( "c1"," ");
  Float_t newMargin1 = 0.13;
  Float_t newMargin2 = 0.15;
  c->SetGrid();
  c->SetTicks();
  c->SetLeftMargin  ( newMargin2 );
  c->SetBottomMargin( newMargin2 );
  c->SetRightMargin ( newMargin1 );
  c->SetTopMargin   ( newMargin1 );
  gStyle->SetPalette( 1, 0 );
  gStyle->SetPaintTextFormat( "3g" );

  h2->SetMarkerSize( 1.5 );
  h2->SetMarkerColor( 0 );
  Float_t labelSize = 0.040;
  h2->GetXaxis()->SetLabelSize( labelSize );
  h2->GetYaxis()->SetLabelSize( labelSize );
  h2->GetXaxis()->SetTitle( "BDT_"+BDT1_name );
  h2->GetYaxis()->SetTitle( "BDT_"+BDT2_name );
  //  h2->LabelsOption( "d" );
  h2->SetLabelOffset( 0.011 );// label offset on x axis    
  h2->Draw("colz"); // color pads   
  c->Update();
  // modify properties of paletteAxis
  TPaletteAxis * paletteAxis = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject( "palette" );
  paletteAxis->SetLabelSize( 0.03 );
  paletteAxis->SetX1NDC( paletteAxis->GetX1NDC() + 0.02 );

  h2->Draw("textsame");  // add text
  // TMVAGlob::plot_logo( );
  c->Update();
  outfilename = "Plot/histo2D_"+lep+".root";
  TFile * fout = TFile::Open(outfilename, "UPDATE");
  h2->Write();
  c->SaveAs("Plot2D/"+histoname+"_"+lep+".pdf");
}
