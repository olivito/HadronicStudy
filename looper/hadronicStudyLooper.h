#ifndef hadronicStudyLooper_h
#define hadronicStudyLooper_h

#include <vector>
#include <list>
#include <string>
#include <map>
#include <set>
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"

#include "TChain.h"
#include "TChainElement.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */


class hadronicStudyLooper
{
public: 
  hadronicStudyLooper();
  ~hadronicStudyLooper() {}

  int  ScanChain(TChain *chain, const TString& prefix = "" );
  void BookHistos (const TString& prefix);
  void InitBaby();
  float dRbetweenVectors(LorentzVector vec1, LorentzVector vec2 );

  // Set globals
  void set_createTree   (bool  b)    { g_createTree   = b; }
  void set_version      (const char* v)    { g_version      = v; }
  void set_json         (const char* v)    { g_json         = v; }        

  // Baby ntuple methods
  void makeOutput (const TString& prefix);
  void closeOutput ();

private:

  // Globals
  bool  g_createTree;
  const char* g_version;
  const char* g_json;      
  bool initialized;
  bool isdata_;
  TFile* outFile;

  // histograms
  TH1F* h_met;

};

#endif
