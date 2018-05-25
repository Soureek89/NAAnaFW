#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"

using namespace std;
void macro_lumi(TString channel, double lumi, double sigma, TString lep)
{
  if((lep=="electron"||lep=="electronantiiso")&& channel=="QCDMuPt20toInf") return;
  if((lep=="muonantiiso"||lep=="electronantiiso")&& channel=="TT_sd") return;
  if((lep=="muonantiiso"||lep=="electronantiiso")&& (channel.Contains("hdamp")||channel.Contains("psq2"))) return;
  Float_t w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15;
  TFile * f1_uncut= TFile::Open("res/"+lep+"/"+channel+"_"+lep+".root");
  h1_uncut =(TH1F*)f1_uncut->Get("h_cutFlow");
  Float_t n_uncut = h1_uncut->GetBinContent(0);
  TFile * f1 = TFile::Open("trees/"+lep+"/trees_"+channel+"_"+lep+".root");

  TTree * t1 =(TTree*)f1->Get("events_2j1t");
  TTree * t2 =(TTree*)f1->Get("events_2j1t_jesUp");
  TTree * t3 =(TTree*)f1->Get("events_2j1t_jesDown");
  TTree * t4 =(TTree*)f1->Get("events_2j1t_jerUp");
  TTree * t5 =(TTree*)f1->Get("events_2j1t_jerDown");
  TTree * t6 =(TTree*)f1->Get("events_3j1t");
  TTree * t7 =(TTree*)f1->Get("events_3j1t_jesUp");
  TTree * t8 =(TTree*)f1->Get("events_3j1t_jesDown");
  TTree * t9 =(TTree*)f1->Get("events_3j1t_jerUp");
  TTree * t10 =(TTree*)f1->Get("events_3j1t_jerDown");
  TTree * t11 =(TTree*)f1->Get("events_3j2t");
  TTree * t12 =(TTree*)f1->Get("events_3j2t_jesUp");
  TTree * t13 =(TTree*)f1->Get("events_3j2t_jesDown");
  TTree * t14 =(TTree*)f1->Get("events_3j2t_jerUp");
  TTree * t15 =(TTree*)f1->Get("events_3j2t_jerDown");

  t1->SetBranchAddress("w", &w1); 
  t2->SetBranchAddress("w", &w2); 
  t3->SetBranchAddress("w", &w3); 
  t4->SetBranchAddress("w", &w4); 
  t5->SetBranchAddress("w", &w5); 
  t6->SetBranchAddress("w", &w6); 
  t7->SetBranchAddress("w", &w7); 
  t8->SetBranchAddress("w", &w8); 
  t9->SetBranchAddress("w", &w9); 
  t10->SetBranchAddress("w", &w10); 
  t11->SetBranchAddress("w", &w11); 
  t12->SetBranchAddress("w", &w12); 
  t13->SetBranchAddress("w", &w13); 
  t14->SetBranchAddress("w", &w14); 
  t15->SetBranchAddress("w", &w15); 

  TFile *f1_new = TFile::Open("trees_lumi/"+lep+"/trees_"+channel+"_"+lep+".root","RECREATE");  

  TTree *t1_new =t1->CloneTree(0); 
  TTree *t2_new =t2->CloneTree(0); 
  TTree *t3_new =t3->CloneTree(0); 
  TTree *t4_new =t4->CloneTree(0); 
  TTree *t5_new =t5->CloneTree(0); 
  TTree *t6_new =t6->CloneTree(0); 
  TTree *t7_new =t7->CloneTree(0); 
  TTree *t8_new =t8->CloneTree(0); 
  TTree *t9_new =t9->CloneTree(0); 
  TTree *t10_new =t10->CloneTree(0); 
  TTree *t11_new =t11->CloneTree(0); 
  TTree *t12_new =t12->CloneTree(0); 
  TTree *t13_new =t13->CloneTree(0); 
  TTree *t14_new =t14->CloneTree(0); 
  TTree *t15_new =t15->CloneTree(0); 

  Int_t nentries1 = (Int_t) t1->GetEntries();
  for (Int_t i = 0; i<nentries1; i++)
    {
      t1->GetEntry(i);
      w1*=(lumi*sigma*1000/n_uncut);
      t1_new->Fill();
    }
  Int_t nentries6 = (Int_t) t6->GetEntries();
  for (Int_t i = 0; i<nentries6; i++)
    {
      t6->GetEntry(i);
      w6*=(lumi*sigma*1000/n_uncut);
      t6_new->Fill();
    }
  Int_t nentries11 = (Int_t) t11->GetEntries();
  for (Int_t i = 0; i<nentries11; i++)
    {
      t11->GetEntry(i);
      w11*=(lumi*sigma*1000/n_uncut);
      t11_new->Fill();
    }
  if(!(channel.Contains("hdamp")||channel.Contains("psq2"))){
    Int_t nentries2 = (Int_t) t2->GetEntries();
    for (Int_t i = 0; i<nentries2; i++)
      {
	t2->GetEntry(i);
	w2*=(lumi*sigma*1000/n_uncut);
	t2_new->Fill();
      }
    Int_t nentries3 = (Int_t) t3->GetEntries();
    for (Int_t i = 0; i<nentries3; i++)
      {
	t3->GetEntry(i);
	w3*=(lumi*sigma*1000/n_uncut);
	t3_new->Fill();
      }
    Int_t nentries4 = (Int_t) t4->GetEntries();
    for (Int_t i = 0; i<nentries4; i++)
      {
	t4->GetEntry(i);
	w4*=(lumi*sigma*1000/n_uncut);
	t4_new->Fill();
      }
    Int_t nentries5 = (Int_t) t5->GetEntries();
    for (Int_t i = 0; i<nentries5; i++)
      {
	t5->GetEntry(i);
	w5*=(lumi*sigma*1000/n_uncut);
	t5_new->Fill();
      }
    Int_t nentries7 = (Int_t) t7->GetEntries();
    for (Int_t i = 0; i<nentries7; i++)
      {
	t7->GetEntry(i);
	w7*=(lumi*sigma*1000/n_uncut);
	t7_new->Fill();
      }
    Int_t nentries8 = (Int_t) t8->GetEntries();
    for (Int_t i = 0; i<nentries8; i++)
      {
	t8->GetEntry(i);
	w8*=(lumi*sigma*1000/n_uncut);
	t8_new->Fill();
      }
    Int_t nentries9 = (Int_t) t9->GetEntries();
    for (Int_t i = 0; i<nentries9; i++)
      {
	t9->GetEntry(i);
	w9*=(lumi*sigma*1000/n_uncut);
	t9_new->Fill();
      }
    Int_t nentries10 = (Int_t) t10->GetEntries();
    for (Int_t i = 0; i<nentries10; i++)
      {
	t10->GetEntry(i);
	w10*=(lumi*sigma*1000/n_uncut);
	t10_new->Fill();
      }
    Int_t nentries12 = (Int_t) t12->GetEntries();
    for (Int_t i = 0; i<nentries12; i++)
      {
	t12->GetEntry(i);
	w12*=(lumi*sigma*1000/n_uncut);
	t12_new->Fill();
      }
    Int_t nentries13 = (Int_t) t13->GetEntries();
    for (Int_t i = 0; i<nentries13; i++)
      {
	t13->GetEntry(i);
	w13*=(lumi*sigma*1000/n_uncut);
	t13_new->Fill();
      }
    Int_t nentries14 = (Int_t) t14->GetEntries();
    for (Int_t i = 0; i<nentries14; i++)
      {
	t14->GetEntry(i);
	w14*=(lumi*sigma*1000/n_uncut);
	t14_new->Fill();
      }
    Int_t nentries15 = (Int_t) t15->GetEntries();
    for (Int_t i = 0; i<nentries15; i++)
      {
	t15->GetEntry(i);
	w15*=(lumi*sigma*1000/n_uncut);
	t15_new->Fill();
      }
  }
  f1_new->Write();
}
