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
#include "../CORE/jetcorr/FactorizedJetCorrector.h"

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

  //  MetCorrector *met_corrector_pfL1FastJetL2L3 = new MetCorrector(jetcorr_filenames_pfL1FastJetL2L3);

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
  unsigned int nEventsPassDup = 0;
  unsigned int nEventsPassFilters = 0;
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

      ++nEventsPassDup;

      //-------------------------------------
      // MET filters
      //-------------------------------------

      float rhovor = evt_ww_rho_vor();

      // rho requirement
      if ( rhovor<0. || rhovor>=40. ) continue;

      // met filters
      if ( evt_cscTightHaloId()       != 0 ) continue;
      if ( evt_hbheFilter()           != 1 ) continue;
      if ( filt_hcalLaser()           != 1 ) continue;
      if ( filt_ecalTP()              != 1 ) continue;
      if ( filt_trackingFailure()     != 1 ) continue;
      if ( isData && ( filt_eeBadSc() != 1)) continue;
      if ( passHBHEFilter()           != 1 ) continue;

      //      ecallasernew_   = passECALLaserFilter();
      ++nEventsPassFilters;
   
      //----------------------------------------
      // Vertex counting
      //----------------------------------------

      int nvtx = 0;
      for (unsigned int ivtx = 0; ivtx < vtxs_position().size(); ++ivtx){
	if(isGoodVertex(ivtx)) ++nvtx;
      }

      //-------------------------------------
      // loop over muons and apply veto
      //-------------------------------------

      VofP4 muons_p4;

      for (unsigned int imu = 0; imu < mus_p4().size(); ++imu) {
	// cut on pt and eta
	if (mus_p4().at(imu).pt() < 10.) continue;
	if (fabs(mus_p4().at(imu).eta()) > 2.4) continue;
	// ID and iso
	if( !muonId( imu , ZMet2012_v1 ))         continue;
	muons_p4.push_back(mus_p4().at(imu));
      }

      if (muons_p4.size()) continue;

      //-------------------------------------
      // loop over electrons and apply veto
      //-------------------------------------

      VofP4 electrons_p4;

      for (unsigned int iel = 0; iel < els_p4().size(); ++iel) {
	// cut on pt and eta
	if (els_p4().at(iel).pt() < 10.) continue;
	if (fabs(els_p4().at(iel).eta()) > 2.4) continue;
	// ID and iso
	//  currently not vetoing transition (first false argument)
        if( !passElectronSelection_Stop2012_v3( iel , false,true,false) )  continue;
	electrons_p4.push_back(els_p4().at(iel));
      }

      if (electrons_p4.size()) continue;

      ++nEventsPass;

      //-------------------------------------
      // loop over jets
      //-------------------------------------

      VofP4 vpfrawjets_p4;
      vpfrawjets_p4.clear();

      vector<float> fullcors;
      vector<float> l2l3cors;
      vector<float> rescors;
      vector<float> l1cors;
      fullcors.clear();
      l2l3cors.clear();
      rescors.clear();
      l1cors.clear();

      float dmetx  = 0.0;
      float dmety  = 0.0;
      float jetptx = 0.0;
      float jetpty = 0.0;

      LorentzVector mht30_p4;
      LorentzVector mht15_p4;
      float ht = 0.;
      int njets = 0;
      int nbtags_med = 0;
      
      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ++ijet) {
	
	// skip jets with |eta| > 5.0
	if( fabs( pfjets_p4().at(ijet).eta() ) > 5.0 ) {
	  continue;
	}

	bool jetlep_overlap = false;
	// skip jets within dR < 0.4 of selected lepton
	for (unsigned int imu = 0; imu < muons_p4.size(); ++imu) {
	  float dr = dRbetweenVectors(muons_p4.at(imu),pfjets_p4().at(ijet));
	  if (dr < 0.4) {
	    jetlep_overlap = true;
	    break;
	  }
	}
	if (jetlep_overlap) continue;
	for (unsigned int iel = 0; iel < electrons_p4.size(); ++iel) {
	  float dr = dRbetweenVectors(electrons_p4.at(iel),pfjets_p4().at(ijet));
	  if (dr < 0.4) {
	    jetlep_overlap = true;
	    break;
	  }
	}
	if (jetlep_overlap) continue;

	// get L1FastL2L3Residual total correction
	jet_corrector_pfL1FastJetL2L3->setRho   ( evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( pfjets_area().at(ijet)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( pfjets_p4().at(ijet).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( pfjets_p4().at(ijet).eta() );
	double corr = jet_corrector_pfL1FastJetL2L3->getCorrection();

	// get L1Fast, L2, L3, Residual individual corrections
	jet_corrector_pfL1FastJetL2L3->setRho   ( evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( pfjets_area().at(ijet)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( pfjets_p4().at(ijet).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( pfjets_p4().at(ijet).eta() );
	vector<float> factors = jet_corrector_pfL1FastJetL2L3->getSubCorrections();

	// get residual correction only
	float rescorr = 1;
	if( isData ){
	  if( factors.size() == 4 ) rescorr = factors.at(3) / factors.at(2);
	  else                      cout << "ERROR! " << factors.size() << " jetSubCorrections" << endl;
	}

	LorentzVector vjet      = corr    * pfjets_p4().at(ijet);

	//---------------------------------------------------------------------------
	// get JES uncertainty
	//---------------------------------------------------------------------------
	
	pfUncertainty->setJetEta(vjet.eta());
	pfUncertainty->setJetPt(vjet.pt());   // here you must use the CORRECTED jet pt
	double unc = pfUncertainty->getUncertainty(true);

	LorentzVector vjetUp   = corr * pfjets_p4().at(ijet) * ( 1 + unc );
	//	LorentzVector vjetDown = corr * pfjets_p4().at(ijet) * ( 1 - unc );

	// PFJetID
	if( !passesPFJetID(ijet) ) continue;

	// pileup beta
	float beta = pfjet_beta( ijet, 2, 0.5 );

        // store raw pfjet p4's and corrections for type1 pfmet
        // using corr jet pT > 10 GeV and adding |eta|<4.7 as in AN2011/459 to avoid problems with
        // erroneously large corrections
        if( vjet.pt() > 10 && fabs(vjet.eta()) < 4.7 ){
	  float l1cor = factors.at(0);
          vpfrawjets_p4.push_back( pfjets_p4().at(ijet) );
          fullcors.push_back( corr );
	  l2l3cors.push_back( corr / l1cor );
          rescors.push_back( rescorr );
          l1cors.push_back( l1cor );
	  // cout<<"l1corr: "<<l1cor<<" (rho: "<<rhovor_
	  //     <<", nvts: "<<ndavtx_
	  //     <<") l2l3resi: "<<corr/l1cor
	  //     <<" all: "<<corr
	  //     <<" pt: "<<vjet.pt()<<" eta: "<<vjet.eta()
	  //     <<endl;
        }

	//------------------------------------------------------------------------------------------------------------
	// MET correction quantities
	// here we store 2 quantities:
	// the delta(METx/y) you get by varying the jet with pT > 10 GeV by their uncertainties (dmetx,dmety
	// the vector sum of pT > 10 GeV selected jets (jetptx,jetpty) --> use this to calculate unclustered energy
	//------------------------------------------------------------------------------------------------------------

	if( vjet.pt() > 10 ){
	  dmetx  += vjetUp.px() - vjet.px();
	  dmety  += vjetUp.py() - vjet.py();
	  jetptx += vjet.px();
	  jetpty += vjet.py();
	}

	// calculate MHT quantities: vector sum of jets above some threshold
	//  require beta > 0.2 for pileup rejection
	if ( (fabs( vjet.eta() ) < 2.4) && (beta > 0.2) ) {
	  if ( vjet.pt() > 15. ) {
	    mht15_p4 -= vjet;
	    if (vjet.pt() > 30. ) {
	      mht30_p4 -= vjet;
	    }
	  }
	}

	// njets: L1FastL2L3Residual, pt > 30 GeV
	if(       vjet.pt()    < 30. )           continue;
	if( fabs( vjet.eta() ) > 2.4 )           continue;
	++njets;
	ht += vjet.pt();

	//-------------------------------------
	// b-tag counting
	//-------------------------------------
	
	// btag variables: CSVM
	float discrimcsv = pfjets_combinedSecondaryVertexBJetTag().at(ijet);
	bool isbtagcsvm = ( discrimcsv > 0.679 ) ? true: false;
	if (isbtagcsvm) {  
	  ++nbtags_med;
	}

	// // btag variables: CSVL
	// bool isbtagcsvl = ( discrimcsv > 0.244 ) ? true: false;
	// if (isbtagcsvl)     nbtagscsvl_++;

	// // btag variables: CSVT
	// bool isbtagcsvt = ( discrimcsv > 0.898 ) ? true: false;
	// if (isbtagcsvt)     nbtagscsvt_++;

      } // loop over jets


      //-------------------------------------
      // corrected MET calculations
      //-------------------------------------
	
      // type1 met's
      pair<float, float> p_t1met10     = Type1PFMET( vpfrawjets_p4 , fullcors , l1cors , 10.0 );
      float t1met10        = p_t1met10.first;
      float t1met10phi     = p_t1met10.second;

      // pair<float, float> p_t1met20     = Type1PFMET( vpfrawjets_p4 , fullcors , l1cors , 20.0 );
      // float t1met20        = p_t1met20.first;
      // float t1met20phi     = p_t1met20.second;	  
      // pair<float, float> p_t1met30     = Type1PFMET( vpfrawjets_p4 , fullcors , l1cors , 30.0 );
      // float t1met30        = p_t1met30.first;
      // float t1met30phi     = p_t1met30.second;	  

      //phi-corrected type1 met
      pair<float, float> p_t1metphicorr = 
	getPhiCorrMET( t1met10, t1met10phi, nvtx, !isData);
      float t1metphicorr    = p_t1metphicorr.first;
      float t1metphicorrphi = p_t1metphicorr.second;

      //---------------------------------------    
      // now calculate METup and METdown
      //---------------------------------------

      // float pfmetx = t1metphicorr * cos( t1metphicorrphi );
      // float pfmety = t1metphicorr * sin( t1metphicorrphi );

      //--------------------------------------------------------
      // calculate unclustered energy x and y components
      // unclustered energy = -1 X ( MET + jets + leptons )
      //--------------------------------------------------------

      // float unclustered_x = -1 * ( pfmetx + jetptx );
      // float unclustered_y = -1 * ( pfmety + jetpty );

      // for( unsigned int imu = 0 ; imu < muons_p4.size() ; ++imu ){
      // 	unclustered_x -= muons_p4.at(imu).px();
      // 	unclustered_y -= muons_p4.at(imu).py();
      // }
      // for( unsigned int iel = 0 ; iel < electrons_p4.size() ; ++iel ){
      // 	unclustered_x -= electrons_p4.at(iel).px();
      // 	unclustered_y -= electrons_p4.at(iel).py();
      // }
      
      //------------------------------------------------------------------------------
      // now vary jets according to JEC uncertainty, vary unclustered energy by 10%
      //------------------------------------------------------------------------------

      // float pfmetx_up = pfmetx - dmetx - 0.1 * unclustered_x; 
      // float pfmety_up = pfmety - dmety - 0.1 * unclustered_y; 

      // // pfmet DOWN
      // float t1metphicorrup    = sqrt( pfmetx_up * pfmetx_up + pfmety_up * pfmety_up );
      // float t1metphicorrphiup = atan2( pfmety_up , pfmetx_up );

      // float pfmetx_dn = pfmetx + dmetx + 0.1 * unclustered_x; 
      // float pfmety_dn = pfmety + dmety + 0.1 * unclustered_y; 

      // // pfmet UP
      // float t1metphicorrdn    = sqrt( pfmetx_dn * pfmetx_dn + pfmety_dn * pfmety_dn );
      // float t1metphicorrphidn = atan2( pfmety_dn , pfmetx_dn );

      //------------------------------------------------------------------------------
      // Fill histograms
      //------------------------------------------------------------------------------

      float st = ht + t1metphicorr;
      float meff = mht30_p4.M();

      h_met->Fill(t1metphicorr);
      h_njets->Fill(njets);
      h_nbtags->Fill(nbtags_med);
      h_ht->Fill(ht);
      h_st->Fill(st);
      h_mht30->Fill(mht30_p4.pt());
      h_mht15->Fill(mht15_p4.pt());
      h_meff->Fill(meff);

      float mht30rel = mht30_p4.pt()/ht;
      float metrel = t1metphicorr/ht;

      h_mht30rel_vs_ht->Fill(ht,mht30rel);
      h_metrel_vs_ht->Fill(ht,metrel);
      h_mht30rel_vs_st->Fill(st,mht30rel);
      h_metrel_vs_st->Fill(st,metrel);
      h_mht30rel_vs_meff->Fill(meff,mht30rel);
      h_metrel_vs_meff->Fill(meff,metrel);

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
  cout << "Passed events before filters: " << nEventsPassDup << endl;
  cout << "Passed events after filters: " << nEventsPassFilters << endl;
  cout << "Passed events after lep vetoes: " << nEventsPass << endl;
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
  h_njets = new TH1F(Form("%s_njets",prefix.Data()),";N(jets)",10,0,10.);
  h_nbtags = new TH1F(Form("%s_nbtags",prefix.Data()),";N(CSVM btags)",10,0,10.);
  h_ht = new TH1F(Form("%s_ht",prefix.Data()),";H_{T} [GeV]",2000,0,2000.);
  h_st = new TH1F(Form("%s_st",prefix.Data()),";S_{T} [GeV]",2000,0,2000.);
  h_mht30 = new TH1F(Form("%s_mht30",prefix.Data()),";Missing H_{T} [GeV]",1000,0,1000.);
  h_mht15 = new TH1F(Form("%s_mht15",prefix.Data()),";Missing H_{T} [GeV]",1000,0,1000.);
  h_meff = new TH1F(Form("%s_meff",prefix.Data()),";M_{eff} [GeV]",2000,0,2000.);

  h_mht30rel_vs_ht = new TH2F(Form("%s_mht30rel_vs_ht",prefix.Data()),";H_{T} [GeV];Missing H_{T} / H_{T} ",2000,0.,2000.,100.,0,1.);
  h_metrel_vs_ht = new TH2F(Form("%s_metrel_vs_ht",prefix.Data()),";H_{T} [GeV];E_{T}^{miss} / H_{T} ",2000,0.,2000.,100.,0,1.);
  h_mht30rel_vs_st = new TH2F(Form("%s_mht30rel_vs_st",prefix.Data()),";S_{T} [GeV];Missing H_{T} / H_{T} ",2000,0.,2000.,100.,0,1.);
  h_metrel_vs_st = new TH2F(Form("%s_metrel_vs_st",prefix.Data()),";S_{T} [GeV];E_{T}^{miss} / H_{T} ",2000,0.,2000.,100.,0,1.);
  h_mht30rel_vs_meff = new TH2F(Form("%s_mht30rel_vs_meff",prefix.Data()),";M_{eff} [GeV];Missing H_{T} / H_{T} ",2000,0.,2000.,100.,0,1.);
  h_metrel_vs_meff = new TH2F(Form("%s_metrel_vs_meff",prefix.Data()),";M_{eff} [GeV];E_{T}^{miss} / H_{T} ",2000,0.,2000.,100.,0,1.);


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

//--------------------------------------------------------------------                                                                                                                                               

std::pair<float,float> hadronicStudyLooper::Type1PFMET( VofP4 jets_p4 , std::vector<float> cors , std::vector<float> l1cors , float minpt ){

  float metx = evt_pfmet() * cos( evt_pfmetPhi() );
  float mety = evt_pfmet() * sin( evt_pfmetPhi() );

  assert( jets_p4.size() == cors.size() );

  for( unsigned int i = 0 ; i < jets_p4.size() ; ++i ){
    float corrpt = jets_p4.at(i).pt() * cors.at(i);
    if( corrpt < minpt ) continue;
    float l1corr = (l1cors.size()==0) ? 1. : l1cors.at(i);
    metx += jets_p4.at(i).px() * l1corr - jets_p4.at(i).px() * cors.at(i);
    mety += jets_p4.at(i).py() * l1corr - jets_p4.at(i).py() * cors.at(i);
  }

  std::pair<float, float> type1met = make_pair( sqrt( metx*metx + mety*mety ), atan2( mety , metx ) );
  return type1met;
}

// --------------------------------------------------------

std::pair<float,float>  hadronicStudyLooper::getPhiCorrMET( float met, float metphi, int nvtx, bool ismc ) {

  //using met phi corrections from C. Veelken (revision 1.6)
  //functions are available here:                                                                                                            
  //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py                               

  float metx = met * cos( metphi );
  float mety = met * sin( metphi );

  float shiftx = 0.;
  float shifty = 0.;

  //use correction for data vs. mc                                                                                                                  
  shiftx = ismc ? (0.1166 + 0.0200*nvtx)
    : (0.2661 + 0.3217*nvtx);
  shifty = ismc ? (0.2764 - 0.1280*nvtx)
    : (-0.2251 - 0.1747*nvtx);

  metx -= shiftx;
  mety -= shifty;

  std::pair<float, float> phicorrmet = make_pair( sqrt( metx*metx + mety*mety ), atan2( mety , metx ) );
  return phicorrmet;
}

//--------------------------------------------------------------------                                                                                
