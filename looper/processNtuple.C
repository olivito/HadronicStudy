
#ifndef __CINT__
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "hadronicStudyLooper.h"

#include <iostream>
#endif

void pickSkimIfExists( TChain *ch, const std::string& base, const std::string& skimPrefix = "" )
{
  TChain *dummy = new TChain("Events");

  int nFiles = 0;
  if (dummy->Add(base.c_str())) {
    nFiles = ch->Add(base.c_str());
    std::cout << "Main " <<base.c_str() << " exists: use it. Loaded " 
              << nFiles << " files" << std::endl;
  } else
    std::cout << "Couldn't find " << base << std::endl;

  // be paranoid
  if (nFiles == 0) {
    std::cout << "ERROR: expected to read files " 
              << base.c_str() << "  but found none" << std::endl;
    //assert(0);
    exit(0);
  }

  return;
}

void processNtuple( TString outfileid = "tt_test", TString infile = "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/merged_ntuple_11.root" )
{

  //---------------------------------------------------------------
  // choose version, output will be written to output/[version]
  //---------------------------------------------------------------
  
  const char* version    = "V00-00-02";
  const char* jsonfile   = "jsons/Cert_198050-207279_8TeV_19p47ifb_Collisions12_JSON_goodruns.txt";
  const bool  useMCSkims = true;

  cout << "Version : " << version     << endl;
  cout << "json    : " << jsonfile    << endl;

  // Load Everything
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");

  gSystem->Load("libMiniFWLite.so");
  gSystem->Load("libhadronicStudyCORE.so");
  gSystem->Load("libhadronicStudyLooper.so");

  hadronicStudyLooper* looper = new hadronicStudyLooper();

  //make baby ntuple
  //  looper->set_createTree(1);
  //set version
  looper->set_version(version);
  //set json
  looper->set_json( jsonfile );
  //batch mode
  //  looper->SetBatchMode(true);

  TChain *chain = new TChain("Events");
  pickSkimIfExists(chain, infile.Data());

  //-------------------------------------
  //set name to get comprehensible output
  //-------------------------------------
  char* sample;
  //MC
  if (infile.Contains("TTJets_MassiveBinDECAY_TuneZ2star_8TeV"))     sample = Form("ttall_%s",  	 outfileid.Data());
  else if (infile.Contains("GluGluToHToGG_M-125")) sample = Form("hgg_m125_%s",  	 outfileid.Data());
  //otherwise
  else sample = Form("boiade_%s", outfileid.Data());

  cout<<"sample is "<<sample<<endl;

  //--------------------------------
  //set luminosity to scale to
  //--------------------------------

  looper->ScanChain(chain, TString(sample));

  //  gSystem->Exit(0);
 
}
