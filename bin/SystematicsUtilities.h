#ifndef _Syst_Weights_
#define _Syst_Weights_

enum weightedSysts { NOSYST=0, BTAGUP = 1,BTAGDOWN=2,MISTAGUP=3,MISTAGDOWN=4, PUUP=5, PUDOWN=6, LEPUP=7, LEPDOWN=8, MAXSYSTS=9};

struct systWeights{
  void initHistogramsSysts(TH1F** histo, TString name, TString, int, float, float);
  void createFilesSysts(TFile ** allFiles, TString basename, TString opt="RECREATE");
  void fillHistogramsSysts(TH1F** histo, float v, float W, double *wcats= NULL,bool verbose=false);
  void fillHistogramsSysts(TH1F** histo, float v, float W,  float *systWeights, int nFirstSysts=(int)MAXSYSTS, double *wcats=NULL, bool verbose=false);
  void rescaleHistograms(TH1F** histo, TH1F ** histoweights);
  void addHistograms(TH1F** histo, TH1F ** histoweights);

  //Trees parts
  void initTreesSysts(TTree ** trees, TFile * file);
  void initTreeSysts(TTree * tree,bool isEventBasedSyst=false);
  
  void branchTreesSysts(TTree ** trees, TString selection, TString name, TFile * file, float * f);
  void fillTreesSysts(TTree ** trees, TString selection);
  void writeTreesSysts(TTree ** trees, TFile * file);
  

  void closeFilesSysts(TFile ** allFiles);
  void writeHistogramsSysts(TH1F** histo, TFile ** allFiles );
  void writeSingleHistogramSysts(TH1F* histo,TFile ** allMyFiles);
  void setMax(int max);
  void setMaxNonPDF(int max);
  //  TFile** initFilesSysts();
  void setSystValue(string name, double value, bool mult=false);
  void setSystValue(int systPlace, double value, bool mult=false);

  float getSystValue(string name);
  float getSystValue(int systPlace);

  void setPDFWeights(float * wpdfs, int numPDFs, float wzero=1.0,bool mult=true);
  void setQ2Weights(float q2up, float q2down, float wzero=1.0,bool mult=true);
  void setTWeight(float tweight, float wtotsample=1.0,bool mult=true);
  void setVHFWeight(int vhf,bool mult=true, double shiftval=0.65);
  void setPDFValue(int numPDF, double value);
  double getPDFValue(int numPDF);
  void setWeight(string name, double value, bool mult=false);
  void setWeight(int systPlace, double value, bool mult=false);
  void prepareDefault(bool addDefault, bool addPDF, bool addQ2, bool addTopPt, bool addJES, bool addJER, bool addVHF, bool addTTSplit, int numPDF=102);
  void addEventBasedSyst(string name);
  void addSyst(string name);
  void addSystNonPDF(string name);
  void setWCats(double *wcats);
  void setSelectionsNames(string * selections);
  
  //Selections
  void addSelection(string name);

  void addkFact(string name);
  void setkFact(string name,float kfact_nom, float kfact_up,float kfact_down,  bool mult=true);

  void copySysts(systWeights sys,bool copySelections=false);
  void calcPDFHisto(TH1F** histos, TH1F* singleHisto, double scalefactor=1.0, int c = 0);
  void setOnlyNominal(bool useOnlyNominal=false);
  //Event based systs:
  bool isEventBasedSyst(int sy);
  bool isEventBasedSelection(int sy);
  //  int getBaseSelection(int sy);
  void setScenario(string scenario);
  void setEventBasedDefault();

  void setData(bool isData);

  bool onlyNominal;
  bool addPDF, addQ2, addTopPt, addVHF, addTTSplit,addJER,addJES;
  bool shortPDFFiles;
  bool isData;
  int maxSysts, maxSystsNonPDF;
  int nPDF;
  int nCategories;
  int nSelections;
  int nEventBasedSysts;
  float weightedSysts[150];
  string eventBasedScenario;
  double wCats[10];
  string eventBasedNames[10];
  int baseSelections[20];
  string weightedNames[150];
  string selectionsNames[20];
  string categoriesNames[10];
  };

struct lept{
  TLorentzVector p4;
  int flavour;
  };

struct by_pt{
    bool operator()(lept const &l1, lept const &l2){
    return l1.p4.Pt() > l2.p4.Pt();
    }
  };

bool pt_ordered (double i,double j){ 
  return (j<i);
  }

float calculate_mtw(float met[0], float metphi[0], vector<TLorentzVector> & mu) {
  TVector2 met_(met[0]*cos(metphi[0]),met[0]*sin(metphi[0]));
  float phi_lmet = fabs(deltaPhi(mu.at(0).Phi(),metphi[0])); 
  float mtw=sqrt(2*mu.at(0).Pt()*met[0]*(1-cos(phi_lmet)));
  return mtw;
  }

void callme(){
  std::cout<<" NaN value"<<std::endl;
  }



void systWeights::closeFilesSysts(TFile *filesout[(int)MAXSYSTS]){
  int MAX= this->maxSystsNonPDF;
  //  int MAXTOT= this->maxSysts;
  
  bool useOnlyNominal = this->onlyNominal;
  for (int c = 0; c < this->nCategories; c++){
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      //cout << "c is now "<< c << " sy "<< sy << " location "<< sy+(MAXTOT+1)*c <<" is histo there? " << histo[sy+(MAXTOT+1)*c] << " file location "<<sy+(MAX+1)*c << " is file there "<< filesout[sy+(MAX+1)*c]<< endl;
      if(!(!useOnlyNominal || sy==0)) continue;
      filesout[(int)sy+(MAX+1)*(c)]->Close();
    }
    
    if(this->addPDF && !this->shortPDFFiles){
      if(!useOnlyNominal){
	filesout[MAX+(MAX+1)*(c)]->Close();
      }
    }
  }
}
void systWeights::setOnlyNominal(bool useOnlyNominal){
  this->onlyNominal=useOnlyNominal;
}
void systWeights::setSystValue(string name, double value, bool mult){
  float zerofact=1.0;
  if(mult){
    zerofact=this->weightedSysts[0];
    //    cout << " mult is"  << endl;
  }
  int MAX = this->maxSysts;
  for(int sy=0;sy<(int)MAX;++sy){
    if(this->weightedNames[(int)sy] ==name){
      this->weightedSysts[(int)sy] =value*zerofact;
    }
  }
}


void systWeights::setSystValue(int place, double value, bool mult){
  float zerofact=1.0;
  if(mult)zerofact=this->weightedSysts[0];
  this->weightedSysts[place] =value*zerofact;
}


float systWeights::getSystValue(string name){
  float val=0.0;
  int MAX = this->maxSysts;
  for(int sy=0;sy<(int)MAX;++sy){
    if(this->weightedNames[(int)sy] ==name){
      val = this->weightedSysts[(int)sy];
    }
  }
  return val;
}

float systWeights::getSystValue(int place){
  //  float val=0.0;
  return this->weightedSysts[place] ;
}


void systWeights::setWeight(string name, double value, bool mult){
  this->setSystValue(name, value, mult);
}

void systWeights::setWeight(int place, double value, bool mult){
  this->setSystValue(place, value, mult);
}


void systWeights::initTreesSysts(TTree ** trees, TFile * file){
  file->cd();
  for(int s =0;s<this->nSelections;++s){
    trees[s]= new TTree(("events_"+this->selectionsNames[s]).c_str(),"");
    bool isEventBasedSelection=this->isEventBasedSelection(s);
    this->initTreeSysts(trees[s],isEventBasedSelection);
  }
}

void systWeights::initTreeSysts(TTree * tree, bool isEventBasedSyst){
  int MAX = this->maxSysts;
  if(this->shortPDFFiles) MAX= this->maxSystsNonPDF;
  for(int sy=0;sy<(int)MAX;++sy){
    if(isEventBasedSyst&&sy>0)continue;
    TString ns= (this->weightedNames[sy]).c_str();
    if(sy==0)ns = "w_nominal";
    tree->Branch(ns,&(this->weightedSysts[(int)sy]));
  }
  for (int c = 0; c < this->nCategories; c++){
    TString cname=  (this->categoriesNames[c]).c_str();
    tree->Branch(cname, &(this->wCats[c]));
  }  

}

void systWeights::branchTreesSysts(TTree ** trees, TString selection, TString name, TFile * file, float * f){
  file->cd();
  for(int s =0;s<this->nSelections;++s){
    if(selection == this->selectionsNames[s]) trees[s]->Branch(name,f);
    //    if((this->selectionsNames[s].find(selection)!=std::string::npos)&&this->isEventBasedSelection(s)){
    if(this->isEventBasedSelection(s)){
      //      cout<<" selection is "<< selection <<" s is "<<  this->selectionsNames[s] << " base selection " << this->selectionsNames[this->baseSelections[s]]<<endl;
      if(selection == this->selectionsNames[this->baseSelections[s]])trees[s]->Branch(name,f);
    }//need to improve this
  }
}
void systWeights::fillTreesSysts(TTree ** trees, TString selection){
  for(int s =0;s<this->nSelections;++s){
    if(selection == this->selectionsNames[s] && !this->isEventBasedSelection(s) && this->eventBasedScenario=="nominal")trees[s]->Fill();
    //    cout << " selection is "<< this->selectionsNames[s]<< " is event based? "<< this->isEventBasedSelection(s)<<endl;
    //    cout << " contains the scenario? " << (selection.Contains(this->eventBasedScenario))<<endl;
    //    if(selection == this->selectionsNames[this->baseSelections[s]])trees[s]->Branch(name,f);
    if(this->isEventBasedSelection(s)){	
      if(this->selectionsNames[s].find(this->eventBasedScenario)!=std::string::npos &&
	 selection == this->selectionsNames[this->baseSelections[s]])trees[s]->Fill();
    }
  }
}
void systWeights::writeTreesSysts(TTree ** trees, TFile * file){
  file->cd();
  for(int s =0;s<this->nSelections;++s){
    trees[s]->Write();
  }

}

void systWeights::initHistogramsSysts(TH1F** histo,TString name, TString title, int nbins, float min, float max){
  for (int c = 0; c < this->nCategories; c++){
    int MAX = this->maxSysts;
    bool useOnlyNominal = this->onlyNominal;
    TString cname= (this->categoriesNames[c]).c_str();
    for(int sy=0;sy<(int)MAX;++sy){
      TString ns= (this->weightedNames[sy]).c_str();
      if(sy==0){
          if(c==0) histo[sy+((MAX+1)*(c))]=new TH1F(name,title,nbins,min,max);
          else histo[sy+((MAX+1)*(c))]=new TH1F(name+"_"+cname,title,nbins,min,max);
          }
      
      if(sy!=0 && !useOnlyNominal) {
          if(c==0)histo[sy+(MAX+1)*c]=new TH1F(name+"_"+ns,title,nbins,min,max);
          else histo[sy+(MAX+1)*c]=new TH1F(name+"_"+ns+"_"+cname,title,nbins,min,max);
          }
      //      cout << " initialized histogram "<< histo[sy+(MAX+1)*c]->GetName() <<" sy " << sy << " c  "<< c <<" location " << sy+(MAX+1)*c << endl;
    }
  }
}
void writeHistogramsSysts(TH1F* histo[(int)MAXSYSTS], TFile *filesout[(int)MAXSYSTS], bool useOnlyNominal=false){
  for(int sy=0;sy<(int)MAXSYSTS;++sy){
    //cout << " writing histo "<< histo[(int)sy]->GetName()<< " in file "<< filesout[(int)sy]->GetName()<<endl;;
    //TString ns= weightedSystsNames((weightedSysts)sy);
    filesout[(int)sy]->cd();
    histo[sy]->Write(histo[0]->GetName());
    //    histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
  }
}


void systWeights:: addHistograms(TH1F** histo, TH1F ** histoweights){
int MAX= this->maxSystsNonPDF;
  int MAXTOT= this->maxSysts;
  bool useOnlyNominal = this->onlyNominal;
  for (int c = 0; c < this->nCategories; c++){
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      if(!(!useOnlyNominal || sy==0)) continue;
      histo[sy+(MAXTOT+1)*c]->Add(histoweights[sy+(MAXTOT+1)*c]);
      //histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
      }
      
      if(this->addPDF){
      if(!useOnlyNominal){
        //cout << " file max is "<< filesout[MAX+(MAX+1)*c]->GetName()<<endl;
        //int npdf=this->maxSysts-this->maxSystsNonPdf;
        int MAXPDF=this->maxSysts;
        for(int sy=MAX;sy<MAXPDF;++sy){
        //    cout << " writing sy "<<sy+(MAXTOT+1)*c<<endl;
        //    cout << " histo is there? "<< histo[sy+(MAXTOT+1)*c]<<endl;
	  histo[sy+(MAXTOT+1)*(c)]->Add(histoweights[sy+(MAXTOT+1)*c]);
	  //    cout << " written sy "<< histo[sy+(MAXTOT+1)*c]->GetName()<<endl;
	}
      }
    }

  }
}

void systWeights:: rescaleHistograms(TH1F** histo, TH1F ** histoweights){
int MAX= this->maxSystsNonPDF;
  int MAXTOT= this->maxSysts;
  bool useOnlyNominal = this->onlyNominal;
  for (int c = 0; c < this->nCategories; c++){
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      if(!(!useOnlyNominal || sy==0)) continue;
      if(histoweights[sy+(MAXTOT+1)*c]->GetMean()!=0) histo[sy+(MAXTOT+1)*c]->Scale(1./histoweights[sy+(MAXTOT+1)*c]->GetMean());
      //histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
      }
      
      if(this->addPDF){
      if(!useOnlyNominal){
        //cout << " file max is "<< filesout[MAX+(MAX+1)*c]->GetName()<<endl;
        //int npdf=this->maxSysts-this->maxSystsNonPdf;
        int MAXPDF=this->maxSysts;
        for(int sy=MAX;sy<MAXPDF;++sy){
        //    cout << " writing sy "<<sy+(MAXTOT+1)*c<<endl;
        //    cout << " histo is there? "<< histo[sy+(MAXTOT+1)*c]<<endl;
	  if(histoweights[sy+(MAXTOT+1)*c]->GetMean()!=0) histo[sy+(MAXTOT+1)*(c)]->Scale(1./histoweights[sy+(MAXTOT+1)*c]->GetMean());
	  //    cout << " written sy "<< histo[sy+(MAXTOT+1)*c]->GetName()<<endl;
	}
      }
    }

  }
}
void systWeights::writeHistogramsSysts(TH1F** histo, TFile **filesout){
  int MAX= this->maxSystsNonPDF;
  int MAXTOT= this->maxSysts;
  bool useOnlyNominal = this->onlyNominal;
  for (int c = 0; c < this->nCategories; c++){
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      //cout << "c is now "<< c << " sy "<< sy << " location "<< sy+(MAXTOT+1)*c <<" is histo there? " << histo[sy+(MAXTOT+1)*c] << " file location "<<sy+(MAX+1)*c << " is file there "<< filesout[sy+(MAX+1)*c]<< endl;
      
      //      cout << " writing histo "<< histo[sy+(MAXTOT+1)*c]->GetName()<< " in file "<< filesout[sy+(MAX+1)*c]->GetName()<<endl;;
      
      //      TString ns= weightedSystsNames((weightedSysts)sy);
      if(!(!useOnlyNominal || sy==0)) continue;
      filesout[(int)sy+(MAX+1)*(c)]->cd();
      if(this->addPDF){
	if(this->weightedNames[sy]=="pdf_totalUp")calcPDFHisto(histo, histo[sy+(MAXTOT+1)*(c)],1.0,c);
	if(this->weightedNames[sy]=="pdf_totalDown")calcPDFHisto(histo, histo[sy+(MAXTOT+1)*(c)],-1.0,c);
	
	;      //this->
      }
      
      histo[sy+(MAXTOT+1)*c]->Write(histo[0]->GetName());
      //histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
    }
    
    //    cout << " histo pdf now "<<endl;
    if(this->addPDF && !this->shortPDFFiles){
      if(!useOnlyNominal){
        filesout[MAX+(MAX+1)*(c)]->cd();
        //cout << " file max is "<< filesout[MAX+(MAX+1)*c]->GetName()<<endl;
        //int npdf=this->maxSysts-this->maxSystsNonPdf;
        int MAXPDF=this->maxSysts;
        for(int sy=MAX;sy<MAXPDF;++sy){
	  //    cout << " writing sy "<<sy+(MAXTOT+1)*c<<endl;
        //    cout << " histo is there? "<< histo[sy+(MAXTOT+1)*c]<<endl;
	  histo[sy+(MAXTOT+1)*(c)]->Write();
	  //    cout << " written sy "<< histo[sy+(MAXTOT+1)*c]->GetName()<<endl;
        }
      }
    }
  }
  //  cout << " end writing histo "<<endl;
}

void systWeights::createFilesSysts(  TFile ** allFiles, TString basename, TString opt){
  for (int c = 0; c < this->nCategories; c++){
    int MAX = this->maxSystsNonPDF;
    int MAXTOT = this->maxSystsNonPDF;
    bool useOnlyNominal = this->onlyNominal;
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      TString ns= (this->weightedNames[(int)sy]);
      //cout << " creating file for syst "<< ns<<endl;
      //      if (c!=0)     cout << " category is "<< c<<endl;
      //      cout << "onlynominal is "<<useOnlyNominal<<endl;

      if(sy==0){
	//	      cout<<" filename is "<< basename+ns+cname+".root"<<endl;
	      allFiles[sy+(MAX+1)*c]= TFile::Open((basename+ns+cname+".root"), opt);
        }
      
      else{
	      if(!useOnlyNominal){
	          //if((ns!="1lep") && (ns!="2lep")&& (ns!="0lep")){
		//        cout<<" filename is "<< basename+ns+cname+".root"<<endl;
	          allFiles[sy+(MAX+1)*c]= TFile::Open((basename+"_"+ns+cname+".root"), opt);
	          }
          }
        //TFile *outTree = TFile::Open(("trees/tree_"+outFileName).c_str(), "RECREATE");
        //cout << " created file at c "<< c << " s "<< sy << " location "<< sy+(MAXTOT+1)*c<< " fname "<<allFiles[sy+(MAXTOT+1)*c]->GetName()<<endl;   
        }
      
      if(this->addPDF){
      if(!useOnlyNominal && !this->shortPDFFiles)allFiles[MAX+((MAX+1)*c)]= TFile::Open((basename+"_pdf"+cname+".root"), opt);
      //cout << " created file at c "<< c << " s "<< MAX+(MAX+1)*c << " location "<< MAX+(MAX+1)*c<<endl;
      //      cout<< " fname "<<allFiles[MAX+(MAXTOT+1)*c]->GetName()<<endl;
        }
      }
    //return allFiles;
}



 
void systWeights::fillHistogramsSysts(TH1F** histo, float v, float w, double * wcats, bool verbose ){
  if(wcats==NULL){
    wcats=this->wCats;
  }
  for (int c = 0; c < this->nCategories; c++){
    if(verbose){
      //      cout<< " in cat loop "<< c<<endl;
      //      cout<< " value "<< wcats[c] <<endl;
    }
    int MAX = this->maxSysts;
    bool useOnlyNominal = this->onlyNominal;
    for(int sy=0;sy<(int)MAX;++sy){
      if(verbose){
	cout<< " in syst loop "<< sy<< endl;
	cout<<" value "<< this->weightedSysts[(int)sy] <<endl ;
      }
      if(this->eventBasedScenario!="nominal"){
       	if(!this->isEventBasedSyst(sy))continue;
	if(this->weightedNames[sy]!=this->eventBasedScenario)continue;
      }
      else if(this->isEventBasedSyst(sy))continue;
      if(sy!=0&& useOnlyNominal)continue;
      float ws = (this->weightedSysts[(int)sy])*wcats[c];
      if(verbose)cout << " filling histogram "<< histo[(int)sy]->GetName() << " with value "<< v <<" and weight "<< w <<" ws "<< ws<<endl;
      histo[sy+(MAX+1)*(c)]->Fill(v, w*ws);
    }
  }
}

//void  initHistogramsSysts (TH1F* histo[(int)MAXSYSTS],TString name, TString title, int nbins, float min, float max, bool useOnlyNominal=false){
//    for(int sy=0;sy<(int)MAXSYSTS;++sy){
//          //TString ns= weightedSystsNames((weightedSysts)sy);
//          histo[sy]=new TH1F(name,title,nbins,min,max);
//          }
//        }

bool systWeights::isEventBasedSelection(int sy){
  bool isEventBased=false;
  for(int e = 0; e < this->nEventBasedSysts;++e){
    if(this->selectionsNames[sy].find(this->eventBasedNames[e])!=std::string::npos){
      isEventBased=true;
      return true;
    }
  }
  return isEventBased;
}
bool systWeights::isEventBasedSyst(int sy){
  bool isEventBased=false;
  for(int e = 0; e < this->nEventBasedSysts;++e){
    if(this->weightedNames[sy]==this->eventBasedNames[e]){
      isEventBased=true;
      return true;
    }
  }
    return isEventBased;
}
  



void systWeights::setWCats(double * wcats){
  for(int i =0;i<this->nCategories;++i){
    //   cout << "setting wcat #"<< i << " to be "<<wcats[i]<<endl;
    this->wCats[i]=wcats[i];
  }
 
}

void systWeights::addSyst(string name){
  this->weightedNames[this->maxSysts]= name;
  this->setMax(maxSysts+1);
  if(name.find("pdf")!=std::string::npos)this->setMaxNonPDF(maxSysts+1);
  this->weightedNames[this->maxSysts]= "";
}

void systWeights::addSystNonPDF(string name){
  this->weightedNames[this->maxSystsNonPDF]= name;
  this->setMaxNonPDF(maxSystsNonPDF+1);
  int nPDF=this->nPDF;
  for(int i =0; i < nPDF;++i){
    stringstream ss;
    ss<< i+1;
    this->weightedNames[i+this->maxSystsNonPDF]= "pdf"+ss.str();
  }
  this->setMax(maxSystsNonPDF+nPDF);
  this->weightedNames[this->maxSysts]= "";
}

void systWeights::addEventBasedSyst(string name){
  this->addSystNonPDF(name);
  this->eventBasedNames[this->nEventBasedSysts]=name;
  this->nEventBasedSysts=this->nEventBasedSysts+1;
}

void systWeights::setScenario(string scenario){
  this->eventBasedScenario=scenario;
  //  cout << " event based scenario is "<<this->eventBasedScenario<<endl;
}

void systWeights::addSelection(string selection){
  this->selectionsNames[this->nSelections]=selection;
  int initSelection=this->nSelections;
  this->nSelections = this->nSelections+1;
  for(int sc=0; sc<this->nEventBasedSysts;sc++){
    this->selectionsNames[this->nSelections]=selection+"_"+this->eventBasedNames[sc];
    this->baseSelections[this->nSelections]=initSelection;
    this->nSelections = this->nSelections+1;
  }
}

void systWeights::setSelectionsNames(string * selections){
  for(int s =0;s<this->nSelections;++s){
    this->selectionsNames[s]=selections[s];
    //   cout << "setting wcat #"<< i << " to be "<<wcats[i]<<endl;
  }
}

void systWeights::copySysts(systWeights sys, bool copySelections){
  for(int i =0; i < sys.maxSysts;++i){
    this->weightedNames[i]=sys.weightedNames[i];
    this->weightedSysts[i]=sys.weightedSysts[i];
  }
  this->nEventBasedSysts=sys.nEventBasedSysts;
  for(int i =0; i < sys.nEventBasedSysts;++i){
    this->eventBasedNames[i]=sys.eventBasedNames[i];
  }
  this->eventBasedScenario=sys.eventBasedScenario;
  this->setOnlyNominal(sys.onlyNominal);
  this->setMax(sys.maxSysts);
  this->setMaxNonPDF(sys.maxSystsNonPDF);
  this->nPDF=sys.nPDF;
  this->nCategories=sys.nCategories;  
  if(copySelections)this->nSelections=sys.nSelections;  
  this->addQ2=sys.addQ2;
  this->addJES=sys.addJES;
  this->addJER=sys.addJER;
  this->addPDF=sys.addPDF;
  this->addTopPt=sys.addTopPt;
  this->addVHF=sys.addVHF;
  this->addTTSplit=sys.addTTSplit;
  this->setWCats(sys.wCats);
  if(copySelections)this->setSelectionsNames(sys.selectionsNames);
  this->shortPDFFiles=sys.shortPDFFiles;

}

void systWeights::setMax(int max){
  this->maxSysts =  max;
  }
void systWeights::setMaxNonPDF(int max){
  this->maxSystsNonPDF =  max;
  }
void systWeights::setData(bool isData){
  this->isData=isData;
  if(this->isData==true){
    this->addQ2=false; this->addPDF=false;
    this->addJES=false; this->addVHF=false; this->addJER=false;
    this->addTTSplit=false;
    this->nPDF=0;
    this->nCategories=1;
  }
}
void systWeights::prepareDefault(bool addDefault, bool addQ2, bool addPDF, bool addTopPt, bool addJES, bool addJER, bool addVHF, bool addTTSplit, int numPDF){ 
  if(this->isData==true){
    addDefault=false;  addQ2=false; addPDF=false;
    addJES=false; addVHF=false; addJER=false;
    addTTSplit=false;
    numPDF=0;
    this->weightedNames[0]="";
    this->setMax(1);
    this->setMaxNonPDF(1);
  }
  this->addPDF=addPDF;
  this->shortPDFFiles=true;
  this->addQ2=addQ2;
  this->addJES=addJES;
  this->addJER=addJER;
  this->addTopPt=addTopPt;
  this->addVHF=addVHF;
  this->addTTSplit=addTTSplit;
  this->nPDF=numPDF;
  this->nCategories=1;
  categoriesNames[0]="";
  this->wCats[0]=1.0;
  this->nSelections=0;
  this->eventBasedScenario="nominal";
  nEventBasedSysts=0;
  if(addDefault){
    this->weightedNames[0]="";
    this->weightedNames[1]="btagUp";
    this->weightedNames[2]="btagDown";
    this->weightedNames[3]="mistagUp";
    this->weightedNames[4]="mistagDown";
    this->weightedNames[5]="puUp";
    this->weightedNames[6]="puDown";
    this->weightedNames[7]="lepUp";
    this->weightedNames[8]="lepDown";
    this->setMax(9);
    this->setMaxNonPDF(9);
    this->weightedNames[this->maxSysts]="";
  }
  if(addQ2){
    this->weightedNames[this->maxSysts]= "q2Up";
    this->weightedNames[this->maxSysts+1]= "q2Down";
    this->setMax(this->maxSysts+2);
    this->setMaxNonPDF(this->maxSystsNonPDF+2);
    this->weightedNames[this->maxSysts]= "";
    
  }
  if(addJES){
    this->weightedNames[this->maxSysts]= "jesUp";
    this->weightedNames[this->maxSysts+1]= "jesDown";
    this->setMax(this->maxSysts+2);
    this->setMaxNonPDF(this->maxSystsNonPDF+2);
    this->weightedNames[this->maxSysts]= "";
    //Consider that it is event based:
    this->eventBasedNames[this->nEventBasedSysts]="jesUp";
    this->eventBasedNames[this->nEventBasedSysts+1]="jesDown";
    nEventBasedSysts=this->nEventBasedSysts+2;

  }
  if(addJER){
    this->weightedNames[this->maxSysts]= "jerUp";
    this->weightedNames[this->maxSysts+1]= "jerDown";
    this->setMax(this->maxSysts+2);
    this->setMaxNonPDF(this->maxSystsNonPDF+2);
    this->weightedNames[this->maxSysts]= "";
    //Consider that it is event based:
    this->eventBasedNames[this->nEventBasedSysts]="jerUp";
    this->eventBasedNames[this->nEventBasedSysts+1]="jerDown";
    nEventBasedSysts=this->nEventBasedSysts+2;
  }
  if(addTopPt){
    this->weightedNames[this->maxSysts]="topPtWeightUp";
    this->weightedNames[this->maxSysts+1]="topPtWeightDown";
    this->setMax(this->maxSysts+2);
    this->setMaxNonPDF(this->maxSystsNonPDF+2);
    this->weightedNames[this->maxSysts]= "";
  }

  if(addVHF){
    this->weightedNames[this->maxSysts]="VHFWeightUp";
    this->weightedNames[this->maxSysts+1]="VHFWeightDown";
    this->setMax(this->maxSysts+2);
    this->setMaxNonPDF(this->maxSystsNonPDF+2);
    this->weightedNames[this->maxSysts]= "";
  }

  if(addTTSplit){
    this->nCategories=4;
    categoriesNames[1]="TT0lep";
    categoriesNames[2]="TT1lep";
    categoriesNames[3]="TT2lep";
    this->wCats[1]=1.0;
    this->wCats[2]=1.0;
    this->wCats[3]=1.0;
    }
  if(addPDF){
    this->weightedNames[this->maxSysts]= "pdf_totalUp";
    this->weightedNames[this->maxSysts+1]= "pdf_totalDown";
    this->setMax(this->maxSysts+2);
    this->setMaxNonPDF(this->maxSystsNonPDF+2);
    if(!this->shortPDFFiles){
      this->weightedNames[this->maxSysts]= "pdf_asUp";
      this->weightedNames[this->maxSysts+1]= "pdf_asDown";
      this->weightedNames[this->maxSysts+2]= "pdf_zmUp";
      this->weightedNames[this->maxSysts+3]= "pdf_zmDown";
      this->setMax(this->maxSysts+4);
      this->setMaxNonPDF(this->maxSystsNonPDF+4);
    }
    int nPDF=this->nPDF;
    for(int i =0; i < nPDF;++i){
      stringstream ss;
      ss<< i+1;
      this->weightedNames[i+this->maxSysts]= "pdf"+ss.str();
      }
    this->setMax(maxSysts+nPDF);
    this->weightedNames[this->maxSysts]= "";
    }
  this->setEventBasedDefault();
}

TString  weightedSystsNames (weightedSysts sy){
  switch(sy){
  case NOSYST : return "";
  case BTAGUP : return "btagUp";
  case BTAGDOWN : return "btagDown";
  case MISTAGUP : return "mistagUp";
  case MISTAGDOWN : return "mistagDown";
  case PUUP : return "puUp";
  case PUDOWN : return "puDown";
  case LEPUP : return "lepUp";
  case LEPDOWN : return "lepDown";
  //case ISOUP : return "isoUp";
  //case ISODOWN : return "isoDown";
  //case TRIGUP : return "trigUp";
  //case TRIGDOWN : return "trigDown";
  case MAXSYSTS : return "";
  }
  return "noSyst";
}

void systWeights::writeSingleHistogramSysts(TH1F* histo, TFile **filesout){  
  int MAX= this->maxSystsNonPDF;
  bool useOnlyNominal = this->onlyNominal;
  for (int c = 0; c < this->nCategories; c++){
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      if(!(!useOnlyNominal || sy==0)) continue;
      //cout << " writing histo "<< histo->GetName()<< " in file "<< filesout[(int)sy]->GetName()<<endl;;
      filesout[(int)sy+(MAX+1)*c]->cd();
      histo->Write();
      //    histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
    }
    if(this->addPDF && !this->shortPDFFiles){
      if(!useOnlyNominal){
        filesout[MAX+(MAX+1)*c]->cd();
        int MAXPDF=this->maxSysts;
        for(int sy=MAX;sy<MAXPDF;++sy){
        //cout << " writing sy "<< histo[sy]->GetName()<<endl;
        histo->Write();
        //cout << " written sy "<< histo[sy]->GetName()<<endl;
        }
      }
    }
  }
}

void systWeights::setEventBasedDefault(){
  int MAX= this->maxSysts;
  for(int sy=0;sy<(int)MAX;++sy){
    if(this->isEventBasedSyst(sy)){
      //      bool multi=true;
      //      cout << "value 0 is "<< this->weightedSysts[0];
      this->setSystValue(sy,1,true);//this->weightedSysts[0]);
      //      cout << " value gotten for syst " << this->weightedNames[sy]<<" is "<< this->weightedSysts[sy]<<endl;
    }
  }
      
      
}


void systWeights::setPDFWeights(float * wpdfs, int numPDFs, float wzero,bool mult){
  float zerofact=1.0;
  if(mult)zerofact=this->weightedSysts[0];
  float rms=0,mean=0;
  for (int i = 0; i < numPDFs; ++i){
    this->setPDFValue(i,zerofact*wpdfs[i]/wzero);
    mean += wpdfs[i]/wzero;
  }
  mean = mean/numPDFs;
  for (int i = 0; i < numPDFs; ++i){
    rms += (wpdfs[i]/wzero-mean)*(wpdfs[i]/wzero-mean);
  }  
  rms= sqrt(rms/numPDFs);
  if(!this->shortPDFFiles){
    this->setSystValue("pdf_asUp", this->getPDFValue(this->nPDF-2)/wzero);
    this->setSystValue("pdf_asDown", zerofact);
    this->setSystValue("pdf_zmUp", this->getPDFValue(this->nPDF-1)/wzero);
    this->setSystValue("pdf_zmDown", zerofact);
  }
  this->setSystValue("pdf_totalUp", zerofact*(1+rms));
  this->setSystValue("pdf_totalDown", zerofact*(1-rms));
}

//void systWeights::setTWeight(float tweight, float totalweight){
void systWeights::setTWeight(float tweight, float wtotsample,bool mult){
  float zerofact=1.0;
  //  cout << " weighted syst 0 is "<< weightedSysts[0]<<endl;
  if(mult)zerofact=this->weightedSysts[0];
  this->setSystValue("topPtWeightUp", zerofact*tweight/wtotsample);
  this->setSystValue("topPtWeightDown", zerofact/tweight*wtotsample);
}

void systWeights::setVHFWeight(int vhf,bool mult,double shiftval){
  float zerofact=1.0;
  double w_shift=0.0;
  //  cout << "vhf is "<<vhf<<endl;
  if (vhf>1)w_shift=shiftval;
  //  cout << " weighted syst 0 is "<< weightedSysts[0]<<endl;
  if(mult)zerofact=this->weightedSysts[0];
  this->setSystValue("VHFWeightUp", zerofact*(1+w_shift));
  this->setSystValue("VHFWeightDown", zerofact*(1-w_shift));
}

void systWeights::setQ2Weights(float q2up, float q2down, float wzero, bool mult){
  float zerofact=1.0;
  if(mult){
    zerofact=this->weightedSysts[0];
    //    cout <<  "zerofact "<< zerofact << endl;
  }
  //  cout <<  "zerofact "<< zerofact << " q2up weight "<< q2up/wzero << " tot to fill "<< zerofact*q2up/wzero<<endl;
  //  cout <<  "zerofact "<< zerofact << " q2down weight "<< q2down/wzero << " tot to fill "<< zerofact*q2down/wzero<<endl;
  this->setSystValue("q2Up", zerofact*q2up/wzero);
  this->setSystValue("q2Down", zerofact*q2down/wzero);
}

double systWeights::getPDFValue(int numPDF){
  if(!addPDF){ cout << "error! No PDF used, this will do nothing."<<endl;return 0.;}
  int MIN = this->maxSystsNonPDF;
  return (double)this->weightedSysts[numPDF+MIN];

}
void systWeights::setPDFValue(int numPDF, double w){
  if(!addPDF){ cout << "error! No PDF used, this will do nothing."<<endl;return;}
  int MIN = this->maxSystsNonPDF;
  this->weightedSysts[numPDF+MIN]=w;

}

//if(this->weightedNames[sy]=="pdf_totalUp")calcPDFHisto(histo, histo[sy+(MAXTOT+1)*(c)],1.0,c);
void systWeights::calcPDFHisto(TH1F** histo, TH1F* singleHisto, double scalefactor, int c){//EXPERIMENTAL
  if(!addPDF){ cout << "error! No PDF used, this will do nothing."<<endl;return;}
  int MAX = this->maxSysts;
  //  for (int c = 0; c < this->nCategories; c++){
    int MIN = this->maxSystsNonPDF+(MAX+1)*c;
    for(int b = 0; b< singleHisto->GetNbinsX();++b){
      float val = singleHisto->GetBinContent(b);
      //      cout << "bin # "<<b << " val "<<val<<endl;
      float mean = 0, devst=0;
      //      cout << "name is "<< singleHisto->GetName()<<endl;
      for(int i = 0; i<this->nPDF;++i ){
	//cout<< " now looking at pdf # "<<i<< " coordinate is "<< MIN+i<<endl;
	//	cout << "is histo there? "<< histo[i+MIN]<<endl;
	//	cout << " histo should be "<< (histo[i+MIN])->GetName()<<endl;
	mean = mean+ histo[i+MIN]->GetBinContent(b);
      }
      mean = mean/this->nPDF;
      //mean = val;//mean/this->nPDF;
      for(int i = 0; i<this->nPDF;++i ){
	devst+=(mean-histo[i+MIN]->GetBinContent(b))*(mean-histo[i+MIN]->GetBinContent(b));
      }
      devst= sqrt(devst/this->nPDF);
      singleHisto->SetBinContent(b,val+devst*scalefactor);
      //      singleHisto->SetBinContent(b,mean+devst*scalefactor);
    }
    //}
}



#endif
