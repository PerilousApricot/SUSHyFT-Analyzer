//#ifndef ROOT_TDISTRIBUTION
//#define ROOT_TDISTRIBUTION

#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

class TDistribution {

 public: 

  TDistribution();
  virtual ~TDistribution();
  
  TString sample;
  TString legentry; 
  TString var;
  int line_color;
  int line_width;
  int fill_color;
  int marker_style;
  int marker_color;
  TH1F *h;
  TString filename;
  TFile *f;
  TTree *tree;
  TString error;
  TString condition;
};

//#endif

