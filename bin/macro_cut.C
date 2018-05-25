#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"

using namespace std;
void macro_cut(TString variabile, TString peso, int nbins, int xmin, int xmax, TString njmt, TString taglio1, TString tag1, TString taglio2, TString tag2)
{
  TH1F *h_sig1 = new TH1F("hsig1", variabile, nbins, xmin, xmax); 
  TH1F *h_bkg1 = new TH1F("hbkg1", variabile, nbins, xmin, xmax); 
  TH1F *h_sig2 = new TH1F("hsig2", variabile, nbins, xmin, xmax); 
  TH1F *h_bkg2 = new TH1F("hbkg2", variabile, nbins, xmin, xmax); 
  TString treename = "events"+njmt;
  TString channel[17]={"ST_T_tch","ST_T_tch_sd","ST_Tbar_tch","ST_Tbar_tch_sd","DYJetsToLL","TT","WToLNu0J","WToLNu1J","WToLNu2J", "ST_T_sch","ST_T_tW", "ST_Tbar_tW","WWTo1L1Nu2Q","WWTo2L2Nu","WZTo1L1Nu2Q","WZTo2L2Q","ZZTo2L2Q"};
  TFile *file[17];
  TH1F* histo[17];
  TH1F* histoST[4];
  TH1F *tmp= new TH1F("tmp", "tmp", nbins, xmin, xmax);
  TString taglio="w_nominal*"+peso+"*("+taglio1+")";
  //cout << taglio << endl;
  for(Int_t j=0; j<17; j++)
    { 
      file[j]=TFile::Open("treesSmall_lumi/trees_"+channel[j]+"_muon.root");
    }
  for(Int_t i=0; i<17; i++)
    {
      file[i]->cd();
      TH1F *h = (TH1F*)tmp->Clone("h");
      /*
      //---------------------------------------------------------------------
       //da commentare una volta avuto i tree funzionanti del t-channel
      if(i<4) 
	{
	  taglio ="w_nominal*1.8*("+taglio1+")";
	  ((TTree*)file[i]->Get(treename))->Project("h",variabile,taglio); //
	}
      else
	{
	  taglio ="w_nominal*"+peso+"*("+taglio1+")";
	  ((TTree*)file[i]->Get(treename))->Project("h",variabile,taglio); //
	}
      //---------------------------------------------------------------------
*/	 
      ((TTree*)file[i]->Get(treename))->Project("h",variabile,taglio); //
      //cout << "Integrali tmp-h " << tmp->Integral() << "-" << h->Integral()  << endl;
      tmp->Reset("ICES");
      histo[i]=(TH1F*)h->Clone(channel[i]);
      if(i<4) h_sig1->Add(histo[i]);
      else h_bkg1->Add(histo[i]);
      //cout << "Integrali histo-hsig-hbkg " << histo[i]->Integral() << "-" << h_sig1->Integral() << "-" << h_bkg1->Integral() << endl;
      delete h;
    }
  TCanvas *c1 = new TCanvas("Confronto_bkg_vs_ST_tch",variabile);
  h_sig1->SetLineColor(kRed);
  h_bkg1->SetLineColor(kBlue);
  Double_t fom1 = h_sig1->Integral()/(pow(h_sig1->Integral()+h_bkg1->Integral(),0.5));
  cout << "la fom per il taglio1 è: " << fom1 << endl;
  h_sig1->Scale(1/h_sig1->Integral());
  h_bkg1->Scale(1/h_bkg1->Integral());
  h_sig1->Draw();
  h_bkg1->Draw("SAME");

  TLegend *leg1 = new TLegend(0.80,0.52,0.95,0.85);
  leg1->SetBorderSize(0);
  leg1->SetLineStyle(0);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.03);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->AddEntry(h_bkg1,"Fondo","l");
  leg1->AddEntry(h_sig1,"Segnale","l");
  leg1->Draw();


  taglio ="w_nominal*"+peso+"*("+taglio1+" && "+taglio2+")";
  //taglio ="w_nominal*1.8*("+taglio1+" && "+taglio2+")";
  for(Int_t i=0; i<4; i++)
    {
      file[i]->cd();
      TH1F *h = (TH1F*)tmp->Clone("h");
      ((TTree*)file[i]->Get(treename))->Project("h",variabile,taglio);//
      tmp->Reset("ICES");
      histoST[i]=(TH1F*)h->Clone(channel[i]);
      if(i==1 || i==3) h_sig2->Add(histoST[i]);
      else h_bkg2->Add(histoST[i]);
      delete h;
    }
  TCanvas *c2 = new TCanvas("Confronto_ST_tch_vs_ST_tch_sd",variabile);
  h_sig2->SetLineColor(kRed);
  h_bkg2->SetLineColor(kBlue);
  Double_t fom2 = h_sig2->Integral()/(pow(h_sig2->Integral()+h_bkg2->Integral(),0.5));
  cout << "la fom per il taglio2  è: " << fom2 << endl;
  h_sig2->Scale(1/h_sig2->Integral());
  h_bkg2->Scale(1/h_bkg2->Integral());
  h_sig2->Draw();
  h_bkg2->Draw("SAME");

  TLegend *leg2 = new TLegend(0.80,0.52,0.95,0.85);
  leg2->SetBorderSize(0);
  leg2->SetLineStyle(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.03);
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->AddEntry(h_bkg2,"Fondo","l");
  leg2->AddEntry(h_sig2,"Segnale","l");
  leg2->Draw(); 

  TString fileout = "Var_discr/"+variabile+"_"+tag1+"_"+tag2+".root";
  TFile *fout = TFile::Open(fileout, "RECREATE");
  c1->Write();
  c2->Write();
}
