#include "hadronicStudyLooper.h"
#include "TTreeCache.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include "../Tools/goodrun.h"
#include "../Tools/vtxreweight.h"
//#include "../Tools/pfjetMVAtools.h"

#include "../CORE/eventSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/susySelections.h"
#include "../CORE/muonSelections.h"
#include "../CORE/jetSelections.h"
#include "../CORE/metSelections.h"
#include "../CORE/jetSmearingTools.h"
#include "../CORE/jetcorr/JetCorrectionUncertainty.h"

bool verbose              = false;
bool doTenPercent         = false;

using namespace std;
using namespace tas;

//--------------------------------------------------------------------

float hadronicStudyLooper::dRbetweenVectors(LorentzVector vec1, LorentzVector vec2 ) { 

  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();

  return sqrt(dphi*dphi + deta*deta);

}

//--------------------------------------------------------------------

hadronicStudyLooper::hadronicStudyLooper()
{

  std::cout << " construct " << std::endl;
  g_createTree   = false;
  initialized = false;
}

//--------------------------------------------------------------------

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return run < other.run;
  if (event != other.event)
    return event < other.event;
  if(lumi != other.lumi)
    return lumi < other.lumi;
  return false;
}

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return false;
  if (event != other.event)
    return false;
  return true;
}

//--------------------------------------------------------------------

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

//--------------------------------------------------------------------

// void hadronicStudyLooper::InitBaby(){
// }

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void hadronicStudyLooper::closeOutput()
{
  outFile->cd();
  //  outTree->Write();
  outFile->Write();
  outFile->Close();
  delete outFile;
}

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */

int hadronicStudyLooper::ScanChain(TChain* chain, const TString& prefix)

{

  //  cout << "ciao " << isData << endl;

  bool isData = false;
  if( prefix.Contains("data") || prefix.Contains("2012") 
      || prefix.Contains("dimu") || prefix.Contains("diel")
      || prefix.Contains("mueg") ){
    cout << "DATA!!!" << endl;
    isData       = true;
    doTenPercent = false;
  }

  cout << "IS DATA: " << isData << endl;

  if( doTenPercent ) cout << "Processing 10% of MC" << endl;

  //------------------------------------------------------------------------------------------------------
  // set json\, vertex reweighting function and msugra cross section files
  //------------------------------------------------------------------------------------------------------
  
  if( !initialized ){

    //set json
    cout << "setting json " << g_json << endl;
    set_goodrun_file( g_json );

    //    if( prefix.Contains("ttall_massivebin") ) 
    set_vtxreweight_rootfile("vtxreweight/vtxreweight_Summer12MC_PUS10_19fb_Zselection.root",true);

    initialized = true;
  }

  //------------------------------------------------------------------------------------------------------
  // latest-and-greatest JEC
  //------------------------------------------------------------------------------------------------------

  std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3;
  FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3;

  jetcorr_filenames_pfL1FastJetL2L3.clear();
  
  string pfUncertaintyFile;
  //string caloUncertaintyFile;

  if ( isData ) {
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_P_V39_AN3_L1FastJet_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_P_V39_AN3_L2Relative_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_P_V39_AN3_L3Absolute_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_P_V39_AN3_L2L3Residual_AK5PF.txt");

    pfUncertaintyFile = "jetCorrections/GR_P_V39_AN3_Uncertainty_AK5PF.txt";
  } 
  else {
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/START53_V21_L1FastJet_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/START53_V21_L2Relative_AK5PF.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/START53_V21_L3Absolute_AK5PF.txt");
    
    pfUncertaintyFile = "jetCorrections/START53_V21_Uncertainty_AK5PF.txt";
  }

  jet_corrector_pfL1FastJetL2L3  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);

  JetCorrectionUncertainty *pfUncertainty   = new JetCorrectionUncertainty( pfUncertaintyFile   );

  MetCorrector *met_corrector_pfL1FastJetL2L3 = new MetCorrector(jetcorr_filenames_pfL1FastJetL2L3);

  // /*
  //  *  Jet Smearer Object to obtain the jet pt uncertainty.
  //  */

  // std::vector<std::string> list_of_file_names;
  // list_of_file_names.push_back("jetSmearData/Spring10_PtResolution_AK5PF.txt");
  // list_of_file_names.push_back("jetSmearData/Spring10_PhiResolution_AK5PF.txt");
  // list_of_file_names.push_back("jetSmearData/jet_resolutions.txt");
  // JetSmearer *jetSmearer = makeJetSmearer(list_of_file_names);
 
  // -----------------------------------------------------------

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  makeOutput(prefix);
  //  if(g_createTree) makeTree(prefix, doFakeApp, frmode);

  //  BookHistos(prefix);
  BookHistos("h");

  cout << " done with initialization "  << endl;
  
  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;
  unsigned int nEventsPass = 0;
  // map isn't needed for this purpose, vector is sufficient
  // better would be to use a struct with run, lb, event
  map<int,int> m_events;

  // loop over files
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TChainElement* currentFile = 0;

  // test chain
  if (!chain)
    {
      throw std::invalid_argument("at::ScanChain: chain is NULL!");
    }
  if (chain->GetListOfFiles()->GetEntries()<1)
    {
      throw std::invalid_argument("at::ScanChain: chain has no files!");
    }
  if (not chain->GetFile())
    {
      throw std::invalid_argument("at::ScanChain: chain has no files or file path is invalid!");
    }
  int nSkip_els_conv_dist = 0;

  while((currentFile = (TChainElement*)fileIter.Next())) {
    TFile* f = new TFile(currentFile->GetTitle());

    cout << currentFile->GetTitle() << endl;

    if (!f || f->IsZombie()) {
      throw std::runtime_error(Form("ERROR::File from TChain is invalid or corrupt: %s", currentFile->GetTitle()));
    }
    
    // get the trees in each file
    // TTree *tree = (TTree*)f->Get("Events");
    TTree *tree = dynamic_cast<TTree*>(f->Get("Events"));
    if (!tree || tree->IsZombie()) {
      throw std::runtime_error(Form("ERROR::File from TChain has an invalid TTree or is corrupt: %s", currentFile->GetTitle()));
    }

    //Matevz
    TTreeCache::SetLearnEntries(100);
    tree->SetCacheSize(128*1024*1024);

    cms2.Init(tree);
      
    unsigned int nEntries = tree->GetEntries();
    for(unsigned int z = 0; z < nEntries; ++z) {
      ++nEventsTotal;

      /////////      cout << nEventsTotal << endl;

      if( doTenPercent ){
	if( !(nEventsTotal%10==0) ) continue;
      }

      // progress feedback to user
      if (nEventsTotal % 1000 == 0){
        
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)){
                
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }

      //Matevz
      tree->LoadTree(z);

      cms2.GetEntry(z);

      if( evt_ww_rho_vor() != evt_ww_rho_vor() ){
	cout << "Skipping event with rho = nan!!!" << endl;
	continue;
      }

      //      InitBaby();

      isdata_ = isData ? 1 : 0;

      if( verbose ){
	cout << "-------------------------------------------------------"   << endl;
	cout << "Event " << z                                               << endl;
	cout << "File  " << currentFile->GetTitle()                         << endl;
	cout << evt_dataset().at(0) << " " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
	cout << "-------------------------------------------------------"   << endl;
      }

      TString datasetname(evt_dataset().at(0));

      //---------------------------------------------
      // event cleaning and good run list
      //---------------------------------------------

      //      if( !cleaning_goodVertexApril2011() )                          continue;
      if( isData && !goodrun(evt_run(), evt_lumiBlock()) ) continue;

      //---------------------
      // skip duplicates
      //---------------------

      if( isData ) {
        DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
        if (is_duplicate(id) ){
          continue;
        }
      }

      //-------------------------------------
      // skip events with bad els_conv_dist
      //-------------------------------------

      bool skipEvent = false;
      for( unsigned int iEl = 0 ; iEl < els_conv_dist().size() ; ++iEl ){
        if( els_conv_dist().at(iEl) != els_conv_dist().at(iEl) ){
          skipEvent = true;
        }
        if( els_sigmaIEtaIEta().at(iEl) != els_sigmaIEtaIEta().at(iEl) ){
          skipEvent = true;
        }
        if( els_sigmaIEtaIEtaSC().at(iEl) != els_sigmaIEtaIEtaSC().at(iEl) ){
          skipEvent = true;
        }
      }
             
      if( skipEvent ){
        nSkip_els_conv_dist++;
        continue;
      }
   
      //-------------------------------------
      // loop over jets
      //-------------------------------------



      //      outTree->Fill();
    
    } // entries

    delete f;
  } // currentFile

  if( nSkip_els_conv_dist > 0 )
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist" << endl;

  cout << endl;
  cout << "Sample: " << prefix << endl;
  cout << endl;
  cout << "Processed events: " << nEventsTotal << endl;
  cout << "Passed events: " << nEventsPass << endl;
  cout << endl;

  closeOutput();
  //  if(g_createTree) closeTree();
  
  already_seen.clear();

  if (nEventsChain != nEventsTotal) 
    std::cout << "ERROR: number of events from files (" << nEventsChain 
	      << ") is not equal to total number of processed events (" << nEventsTotal << ")" << std::endl;
  
  return 0;

}


//--------------------------------------------------------------------
 
void hadronicStudyLooper::BookHistos(const TString& prefix)
{
  // Prefix comes from the sample and it is passed to the scanning function
  // Suffix is "ee" "em" "em" "all" which depends on the final state
  // For example: histogram named tt_hnJet_ee would be the Njet distribution
  // for the ee final state in the ttbar sample.
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  cout << "Begin book histos..." << endl;

  // TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  // rootdir->cd();
  if (outFile) outFile->cd();

  h_met = new TH1F(Form("%s_met",prefix.Data()),";E_{T}^{miss} [GeV]",1000,0,1000.);

  cout << "End book histos..." << endl;
}// CMS2::BookHistos()


void hadronicStudyLooper::makeOutput(const TString& prefix){
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  TString tpsuffix = "";
  if( doTenPercent ) tpsuffix = "_tenPercent";

  outFile   = new TFile(Form("output/%s_smallTree%s.root",prefix.Data(),tpsuffix.Data()), "RECREATE");
  //  outFile   = new TFile(Form("output/%s/%s_smallTree%s%s.root",g_version,prefix.Data(),frsuffix,tpsuffix), "RECREATE");
  //  outFile   = new TFile("baby.root","RECREATE");
  outFile->cd();
  //  outTree = new TTree("t","Tree");

  //Set branch addresses
  //variables must be declared in hadronicStudyLooper.h
}

