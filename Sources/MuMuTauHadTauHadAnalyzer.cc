#define MuMuTauHadTauHadAnalyzer_cxx
#include "MuMuTauHadTauHadAnalyzer.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <iomanip>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>
using namespace std;

void MuMuTauHadTauHadAnalyzer::Loop()
{
   TString outputfileName = createOutputFileName();
   TFile* outputFile = new TFile(outputfileName, "RECREATE");
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   if (nMaxEvents >= 0 && nMaxEvents  < nentries) nentries = nMaxEvents;
   cout << "We will run on " << nentries << " events" << endl;


   int DDTV_matched_count = 0;
   int DDTV_01_count = 0;
   int DDTV_015_count = 0;
   int DDTV_02_count = 0;
   int DDTV_025_count = 0;
   int DDTV_03_count = 0;
   int DDTV_035_count = 0;
   int DDTV_04_count = 0;
   int DDTV_045_count = 0;
   int DDTV_05_count = 0;
   int DDTV_055_count = 0;
   int DDTV_06_count = 0;
   int DDTV_065_count = 0;
   int DDTV_07_count = 0;
   int DDTV_075_count = 0;
   int DDTV_08_count = 0;
   int DDTV_085_count = 0;
   int DDTV_09_count = 0;
   int DDTV_095_count = 0;

   int DDTV_MD_matched_count = 0;
   int DDTV_MD_01_count = 0;
   int DDTV_MD_015_count = 0;
   int DDTV_MD_02_count = 0;
   int DDTV_MD_025_count = 0;
   int DDTV_MD_03_count = 0;
   int DDTV_MD_035_count = 0;
   int DDTV_MD_04_count = 0;
   int DDTV_MD_045_count = 0;
   int DDTV_MD_05_count = 0;
   int DDTV_MD_055_count = 0;
   int DDTV_MD_06_count = 0;
   int DDTV_MD_065_count = 0;
   int DDTV_MD_07_count = 0;
   int DDTV_MD_075_count = 0;
   int DDTV_MD_08_count = 0;
   int DDTV_MD_085_count = 0;
   int DDTV_MD_09_count = 0;
   int DDTV_MD_095_count = 0;



   int DeepVsjetraw_count = 0;
   int DeepVsjetvvvloose_count = 0;
   int DeepVsjetvvloose_count = 0;
   int DeepVsjetvloose_count = 0;
   int DeepVsjetloose_count = 0;
   int DeepVsjetmedium_count = 0;
   int DeepVsjettight_count = 0;
   int DeepVsjetvtight_count = 0;
   int DeepVsjetvvtight_count = 0;


   int DeepVSjet_vvvloose_count=0;

   Long64_t nbytes = 0, nb = 0;



   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      if (jentry % 1000 == 0 && jentry > 0) cout << "*** Processing #Events: " << jentry << endl;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

      // ---- define varibles that will be used to be filled into histograms ---
      TLorentzVector Mu1;
      TLorentzVector Mu2;
      //TLorentzVector Tau1;
      //TLorentzVector Tau2;
      //      TLorentzVector DeepDiTaujet;

      float Mu1Iso;
      float Mu2Iso;
      float Tau1Iso;
      float Tau1DM;
      float Tau2Iso;
      float Tau2DM;

      unsigned int indexGenTau = -1;
      unsigned int indexMu1 = -1;
      // =============================================================================
      //      std::cout << "******************" << std::endl;
      // ---- start loop on muon candidates for mu1 ----
      bool findMu1 = false;
      for (unsigned int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
      {
          if (recoMuonTriggerFlag->at(iMuon) == 1 && recoMuonIsolation->at(iMuon) < 0.25) 
          {
              Mu1.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
              Mu1Iso = recoMuonIsolation->at(iMuon);
              indexMu1 = iMuon;
              findMu1 = true;
              break;
          } // end if there is any matched Mu1 candidiate
      } // end loop for mu1
      // std::cout << "event line 65" << std::endl;
      if (!findMu1) continue;
      float smallestDR = 1.0; // dR cut between Mu1 and Mu2
      bool findMu2 = false;
      //std::cout << "event line 69" << std::endl;
      // ---- start loop on muon candidates for mu2 ----
      for (unsigned int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
      {
          if (iMuon == indexMu1) continue;
          if ((invertedMu2Iso == false && recoMuonIsolation->at(iMuon) > Mu2IsoThreshold) || (invertedMu2Iso == true && recoMuonIsolation->at(iMuon) < Mu2IsoThreshold)) continue;
	  //	  std::cout << "event line 75"<< std::endl;

          TLorentzVector Mu2Cand; // prepare this variable for dR(Mu1,Mu2) implementation
          Mu2Cand.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));

          if((Mu1.DeltaR(Mu2Cand) < smallestDR) && (recoMuonPDGId->at(indexMu1) == (-1) * recoMuonPDGId->at(iMuon)) && ((Mu1+Mu2Cand).M() > diMuonMassLowThreshold) && ((Mu1+Mu2Cand).M() < diMuonMassHighThreshold))
          {
              Mu2.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
              Mu2Iso = recoMuonIsolation->at(iMuon);
              smallestDR = Mu1.DeltaR(Mu2);
              findMu2 = true;
          } // end if pair candidates
      } // end loop for mu2
      // std::cout << "event here line 87" << std::endl;
          
      if (!findMu2) continue;


      // find gen matches for muons 


      TLorentzVector GenMu1;
      TLorentzVector GenMu2;

      bool findMatchedRecGenMu1 = false;
      bool findMatchedRecGenMu2 = false;
         
      unsigned int indexGenMu1 = -1;
      unsigned int indexGenMu2 = -1;

      smallestDR = 0.15;
      for (unsigned int iGenMu=0; iGenMu<genMuonPt->size(); iGenMu++)
	{
	  TLorentzVector GenMuCand;
	  GenMuCand.SetPtEtaPhiM(genMuonPt->at(iGenMu), genMuonEta->at(iGenMu), genMuonPhi->at(iGenMu), genMuonMass->at(iGenMu));
	  if (Mu1.DeltaR(GenMuCand) <= smallestDR)
	    {
	      smallestDR = Mu1.DeltaR(GenMuCand);
	      findMatchedRecGenMu1 = true;
	      GenMu1 = GenMuCand;
	      indexGenMu1 = iGenMu;
	    } // end if Mu1.DeltaR(GenMuCand) <= smallestDR                                                                                                                    
	} // end for loop on GenMu1                                                                                                                                            

      // --------- search for matched genMu2 for Mu2 --------------                                                                                                            
      smallestDR = 0.15;
      for (unsigned int iGenMu=0; iGenMu<genMuonPt->size(); iGenMu++)
	{
	  TLorentzVector GenMuCand;
	  GenMuCand.SetPtEtaPhiM(genMuonPt->at(iGenMu), genMuonEta->at(iGenMu), genMuonPhi->at(iGenMu), genMuonMass->at(iGenMu));
	  if (Mu2.DeltaR(GenMuCand) <= smallestDR && iGenMu != indexGenMu1)
	    {
	      smallestDR = Mu2.DeltaR(GenMuCand);
	      findMatchedRecGenMu2 = true;
	      GenMu2 = GenMuCand;
	      indexGenMu2 = iGenMu;
	    } // end if Mu2.DeltaR(GenMuCand) <= smallestDR && iGenMu != indexGenMu1                                                                                           
	}// end for loop on GenMu2                                                                                                                                             
      // } // gen Muonpt->Size>0 

      
      // ------- start loop on tau candidates -------
      TLorentzVector Tau1;
      TLorentzVector Tau2;    

      unsigned int indexTau1 = -1;
      unsigned int indexTau2 = -1;

      bool findTau1 = false;

      bool condTauDeepVsjetVVVLoose_Tau1 = false;
      bool condTauDeepVsjetVVLoose_Tau1 = false;
      bool condTauDeepVsjetVLoose_Tau1 = false;
      bool condTauDeepVsjetLoose_Tau1 = false;
      bool condTauDeepVsjetMedium_Tau1 = false;
      bool condTauDeepVsjetTight_Tau1 = false;
      bool condTauDeepVsjetVTight_Tau1 = false;

      bool condTauDeepVsjetVVVLoose_Tau2 = false;
      bool condTauDeepVsjetVVLoose_Tau2 = false;
      bool condTauDeepVsjetVLoose_Tau2 = false;
      bool condTauDeepVsjetLoose_Tau2 = false;
      bool condTauDeepVsjetMedium_Tau2 = false;
      bool condTauDeepVsjetTight_Tau2 = false;
      bool condTauDeepVsjetVTight_Tau2 = false;


      for (unsigned int iTau=0; iTau<recoTauPt->size(); iTau++)
      {

	Tau1.SetPtEtaPhiE(recoTauPt->at(iTau), recoTauEta->at(iTau), recoTauPhi->at(iTau), recoTauEnergy->at(iTau));
	indexTau1 = iTau;
	findTau1 = true;
	if (recoTauDeepVSjetVVVLoose->at(iTau)>0){condTauDeepVsjetVVVLoose_Tau1 = true;}
	if (recoTauDeepVSjetVVLoose->at(iTau)>0){condTauDeepVsjetVVLoose_Tau1 = true;}
	if (recoTauDeepVSjetVLoose->at(iTau)>0){condTauDeepVsjetVLoose_Tau1 = true;}
	if (recoTauDeepVSjetLoose->at(iTau)>0){condTauDeepVsjetLoose_Tau1 = true;}
	if (recoTauDeepVSjetMedium->at(iTau)>0){condTauDeepVsjetMedium_Tau1 = true;}
	if (recoTauDeepVSjetTight->at(iTau)>0){condTauDeepVsjetTight_Tau1 = true;}
	if (recoTauDeepVSjetVTight->at(iTau)>0){condTauDeepVsjetVTight_Tau1 = true;}

	break;
      }

      float smallestdRrecotaus = 1.0;
      bool findTau2 = false;
      if(findTau1){
	for (unsigned int iTau=0; iTau<recoTauPt->size(); iTau++){
	  if(iTau == indexTau1) continue;

	  TLorentzVector Tau2Cand; // prepare this variable for dR(Mu1,Mu2) implementation                                                                               
          Tau2Cand.SetPtEtaPhiE(recoTauPt->at(iTau), recoTauEta->at(iTau), recoTauPhi->at(iTau), recoTauEnergy->at(iTau));

          if((Tau1.DeltaR(Tau2Cand) < smallestdRrecotaus) && (recoTauPDGId->at(indexTau1) == (-1) * recoTauPDGId->at(iTau))){
              Tau2.SetPtEtaPhiE(recoTauPt->at(iTau), recoTauEta->at(iTau), recoTauPhi->at(iTau), recoTauEnergy->at(iTau));
              //Mu2Iso = recoMuonIsolation->at(iTau);
              smallestdRrecotaus = Tau1.DeltaR(Tau2);
              findTau2 = true;

	      if (recoTauDeepVSjetVVVLoose->at(iTau)>0){condTauDeepVsjetVVVLoose_Tau2 = true;}
	      if (recoTauDeepVSjetVVLoose->at(iTau)>0){condTauDeepVsjetVVLoose_Tau2 = true;}
	      if (recoTauDeepVSjetVLoose->at(iTau)>0){condTauDeepVsjetVLoose_Tau2 = true;}
	      if (recoTauDeepVSjetLoose->at(iTau)>0){condTauDeepVsjetLoose_Tau2 = true;}
	      if (recoTauDeepVSjetMedium->at(iTau)>0){condTauDeepVsjetMedium_Tau2 = true;}
	      if (recoTauDeepVSjetTight->at(iTau)>0){condTauDeepVsjetTight_Tau2 = true;}
	      if (recoTauDeepVSjetVTight->at(iTau)>0){condTauDeepVsjetVTight_Tau2 = true;}

	  } // if matching the reco tau pair
	} // loop over reco taus for tau 2
      }// if findTau1



      
      //	if (deepTauID && recoTauDeepVSjetraw->size() > 0){

	  //change the counts here to be checks for tau 1 and tau 2

	  

	  if (condTauDeepVsjetVVVLoose_Tau1 && condTauDeepVsjetVVVLoose_Tau2){
	    DeepVsjetvvvloose_count++;
	  }

	  if (condTauDeepVsjetVVLoose_Tau1 && condTauDeepVsjetVVLoose_Tau2){
            DeepVsjetvvloose_count++;
          }

	 
	  if (condTauDeepVsjetVLoose_Tau1 && condTauDeepVsjetVLoose_Tau2){
            DeepVsjetvloose_count++;
          }


	  if (condTauDeepVsjetLoose_Tau1 && condTauDeepVsjetLoose_Tau2){
            DeepVsjetloose_count++;
          }


	  if (condTauDeepVsjetMedium_Tau1 && condTauDeepVsjetMedium_Tau2){
            DeepVsjetmedium_count++;
          }


	  if (condTauDeepVsjetTight_Tau1 && condTauDeepVsjetTight_Tau2){
            DeepVsjettight_count++;
          }

	  if (condTauDeepVsjetVTight_Tau1 && condTauDeepVsjetVTight_Tau2){
            DeepVsjetvtight_count++;
          }

   
	/*
    if (deepTauID && recoTauDeepVSjetraw->size() > 0)
          {
              // -------------------------------------------------------------------------------
              bool condTauDeepVSeLoose = deepTauVSele == "LOOSE" && recoTauDeepVSeLoose->at(iTau)>0;
              bool condTauDeepVSjetLoose = deepTauVSjet == "LOOSE" && recoTauDeepVSjetLoose->at(iTau)>0;
              bool condTauDeepVSmuLoose = deepTauVSmu == "LOOSE" && recoTauDeepVSmuLoose->at(iTau)>0;
              bool condTauDeepVSeMedium = deepTauVSele == "MEDIUM" && recoTauDeepVSeMedium->at(iTau)>0;
              bool condTauDeepVSjetMedium = deepTauVSjet == "MEDIUM" && recoTauDeepVSjetMedium->at(iTau)>0;
              bool condTauDeepVSmuMedium = deepTauVSmu == "MEDIUM" && recoTauDeepVSmuMedium->at(iTau)>0;
              bool condTauDeepVSeTight = deepTauVSele == "TIGHT" && recoTauDeepVSeTight->at(iTau)>0;
              bool condTauDeepVSjetTight = deepTauVSjet == "TIGHT" && recoTauDeepVSjetTight->at(iTau)>0;
              bool condTauDeepVSmuTight = deepTauVSmu == "TIGHT" && recoTauDeepVSmuTight->at(iTau)>0;
              bool condTauDeepVSeVLoose = deepTauVSele == "VLOOSE" && recoTauDeepVSeVLoose->at(iTau)>0;
              bool condTauDeepVSjetVLoose = deepTauVSjet == "VLOOSE" && recoTauDeepVSjetVLoose->at(iTau)>0;
              bool condTauDeepVSmuVLoose = deepTauVSmu == "VLOOSE" && recoTauDeepVSmuVLoose->at(iTau)>0;
              bool condTauDeepVSeVTight = deepTauVSele == "VTIGHT" && recoTauDeepVSeVTight->at(iTau)>0;
              bool condTauDeepVSjetVTight = deepTauVSjet == "VTIGHT" && recoTauDeepVSjetVTight->at(iTau)>0;


 bool condTauDeepVSeVVLoose = deepTauVSele == "VVLOOSE" && recoTauDeepVSeVVLoose->at(iTau)>0;
              bool condTauDeepVSjetVVLoose = deepTauVSjet == "VVLOOSE" && recoTauDeepVSjetVVLoose->at(iTau)>0;
              
              bool condTauDeepVSeVVTight = deepTauVSele == "VVTIGHT" && recoTauDeepVSeVVTight->at(iTau)>0;
              bool condTauDeepVSjetVVTight = deepTauVSjet == "VVTIGHT" && recoTauDeepVSjetVVTight->at(iTau)>0;
              bool condTauDeepVSeVVVLoose = deepTauVSele == "VVVLOOSE" && recoTauDeepVSeVVVLoose->at(iTau)>0;
              bool condTauDeepVSjetVVVLoose = deepTauVSjet == "VVVLOOSE" && recoTauDeepVSjetVVVLoose->at(iTau)>0;
              bool condTauDeepVSeNull = deepTauVSele != "LOOSE" && deepTauVSele != "MEDIUM" && deepTauVSele != "TIGHT" && deepTauVSele != "VLOOSE" && deepTauVSele != "VTIGHT" && deepTauVSele != "VVLOOSE" && deepTauVSele != "VVTIGHT" && deepTauVSele != "VVVLOOSE";
              bool condTauDeepVSmuNull = deepTauVSmu != "LOOSE" && deepTauVSmu != "MEDIUM" && deepTauVSmu != "TIGHT" && deepTauVSmu != "VLOOSE";
              // -------------------------------------------------------------------------------
              bool passCondTauDeepVSele = (condTauDeepVSeLoose || condTauDeepVSeMedium || condTauDeepVSeTight || condTauDeepVSeVLoose || condTauDeepVSeVTight || condTauDeepVSeVVLoose || condTauDeepVSeVVTight || condTauDeepVSeVVVLoose || condTauDeepVSeNull);
              bool passCondTauDeepVSjet = (condTauDeepVSjetLoose || condTauDeepVSjetMedium || condTauDeepVSjetTight || condTauDeepVSjetVLoose || condTauDeepVSjetVTight || condTauDeepVSjetVVLoose || condTauDeepVSjetVVTight || condTauDeepVSjetVVVLoose);
              bool passCondTauDeepVSmu = (condTauDeepVSmuLoose || condTauDeepVSmuMedium || condTauDeepVSmuTight || condTauDeepVSmuVLoose || condTauDeepVSmuNull);
              bool passCondTauDeep = passCondTauDeepVSele && passCondTauDeepVSjet && passCondTauDeepVSmu;
              // -------------------- inverted deep Tau ID -----------------------------
              bool condInvertTauDeepVSjetLoose = deepTauVSjet == "LOOSE" && recoTauDeepVSjetLoose->at(iTau)<=0;
              bool condInvertTauDeepVSjetMedium = deepTauVSjet == "MEDIUM" && recoTauDeepVSjetMedium->at(iTau)<=0;
              bool condInvertTauDeepVSjetTight = deepTauVSjet == "TIGHT" && recoTauDeepVSjetTight->at(iTau)<=0;
              bool condInvertTauDeepVSjetVLoose = deepTauVSjet == "VLOOSE" && recoTauDeepVSjetVLoose->at(iTau)<=0;
              bool condInvertTauDeepVSjetVTight = deepTauVSjet == "VTIGHT" && recoTauDeepVSjetVTight->at(iTau)<=0;
              bool condInvertTauDeepVSjetVVLoose = deepTauVSjet == "VVLOOSE" && recoTauDeepVSjetVVLoose->at(iTau)<=0;
              bool condInvertTauDeepVSjetVVTight = deepTauVSjet == "VVTIGHT" && recoTauDeepVSjetVVTight->at(iTau)<=0;
              bool condInvertTauDeepVSjetVVVLoose = recoTauDeepVSjetVVVLoose->at(iTau)>0;
              bool condInvertTauDeepVSeVVVLoose = recoTauDeepVSeVVVLoose->at(iTau)>0;
              bool condInvertTauDeepVSmuVLoose = recoTauDeepVSmuVLoose->at(iTau)>0;
              // -------------------------------------------------------------------------------
              bool passCondInvertTauDeepVSjet = ((condInvertTauDeepVSjetLoose || condInvertTauDeepVSjetMedium || condInvertTauDeepVSjetTight || condInvertTauDeepVSjetVLoose || condInvertTauDeepVSjetVTight || condInvertTauDeepVSjetVVLoose || condInvertTauDeepVSjetVVTight) && condInvertTauDeepVSjetVVVLoose && condInvertTauDeepVSeVVVLoose && condInvertTauDeepVSmuVLoose);
              // -------------------------------------------------------------------------------
              if ((!invertedTauIso && !passCondTauDeep) || (invertedTauIso && !passCondInvertTauDeepVSjet)) continue;
          } // end if deepTauID && recoTauDeepVSjetraw->size() > 0

  else{
              bool condTauMVARaw = tauMVAIsoRawORWP == true && recoTauIsoMVArawValue->at(iTau) > tauMVAIsoRawThreshold;
              bool condTauMVAWPVVLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VVLOOSE" && recoTauIsoMVAVVLoose->at(iTau)>0;
              bool condTauMVAWPVLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VLOOSE" && recoTauIsoMVAVLoose->at(iTau)>0;
              bool condTauMVAWPLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "LOOSE" && recoTauIsoMVALoose->at(iTau)>0;
              bool condTauMVAWPMedium = tauMVAIsoRawORWP == false && tauMVAIsoWP == "MEDIUM" && recoTauIsoMVAMedium->at(iTau)>0;
              bool condTauMVAWPTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "TIGHT" && recoTauIsoMVATight->at(iTau)>0;
              bool condTauMVAWPVTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VTIGHT" && recoTauIsoMVAVTight->at(iTau)>0;
              bool condTauMVAWPVVTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VVTIGHT" && recoTauIsoMVAVVTight->at(iTau)>0;
              bool passCondTauMVA = (condTauMVARaw || condTauMVAWPVVLoose || condTauMVAWPVLoose || condTauMVAWPLoose || condTauMVAWPMedium || condTauMVAWPTight || condTauMVAWPVTight || condTauMVAWPVVTight);
              // -------------------------------------------------------------------------------------------------
              bool condInvertTauMVARaw = recoTauIsoMVArawValue->at(iTau) > tauMVAIsoRawThreshold;
              bool condInvertTauMVAWPVVLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VVLOOSE" && recoTauIsoMVAVVLoose->at(iTau)<=0;
              bool condInvertTauMVAWPVLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VLOOSE" && recoTauIsoMVAVLoose->at(iTau)<=0;
              bool condInvertTauMVAWPLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "LOOSE" && recoTauIsoMVALoose->at(iTau)<=0;
              bool condInvertTauMVAWPMedium = tauMVAIsoRawORWP == false && tauMVAIsoWP == "MEDIUM" && recoTauIsoMVAMedium->at(iTau)<=0;
              bool condInvertTauMVAWPTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "TIGHT" && recoTauIsoMVATight->at(iTau)<=0;
              bool condInvertTauMVAWPVTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VTIGHT" && recoTauIsoMVAVTight->at(iTau)<=0;
              bool condInvertTauMVAWPVVTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VVTIGHT" && recoTauIsoMVAVVTight->at(iTau)<=0;
              // ------ always require tau candidates pass vvvloose MVA id in order to have similar dynamic shape as real tau
              bool passCondInvertTauMVA = (condInvertTauMVARaw && (condInvertTauMVAWPVVLoose || condInvertTauMVAWPVLoose || condInvertTauMVAWPLoose || condInvertTauMVAWPMedium || condInvertTauMVAWPTight || condInvertTauMVAWPVTight || condInvertTauMVAWPVVTight));
              // -------------------------------------------------------------------------------------------------
              bool condTauAntiMuMVALoose = tauAntiMuDisc == "LOOSE" && recoTauAntiMuMVALoose->at(iTau)>0;
              bool condTauAntiMuMVATight = tauAntiMuDisc == "TIGHT" && recoTauAntiMuMVATight->at(iTau)>0; 
              bool condTauAntiMuMVANull = tauAntiMuDisc != "LOOSE" && tauAntiMuDisc != "TIGHT";
              bool passCondTauAntiMuMVA = (condTauAntiMuMVALoose || condTauAntiMuMVATight || condTauAntiMuMVANull);
              // -------------------------------------------------------------------------------------------------
              bool condTauAntiEleMVALoose = tauAntiEleDisc == "LOOSE" && recoTauAntiEleMVALoose->at(iTau)>0;
              bool condTauAntiEleMVAMedium = tauAntiEleDisc == "MEDIUM" && recoTauAntiEleMVAMedium->at(iTau)>0;
              bool condTauAntiEleMVATight = tauAntiEleDisc == "TIGHT" && recoTauAntiEleMVATight->at(iTau)>0; 
              bool condTauAntiEleMVANull = tauAntiEleDisc != "LOOSE" && tauAntiEleDisc != "MEDIUM" && tauAntiEleDisc != "TIGHT";
              bool passCondTauAntiEleMVA = (condTauAntiEleMVALoose || condTauAntiEleMVAMedium || condTauAntiEleMVATight || condTauAntiEleMVANull);
              // -------------------------------------------------------------------------------------------------

 if ((!invertedTauIso && !passCondTauMVA) || (invertedTauIso && !passCondInvertTauMVA) || !passCondTauAntiMuMVA || !passCondTauAntiEleMVA) continue;
          } // end if !deepTauID (tauMVAID)
          TLorentzVector TauCand;
          TauCand.SetPtEtaPhiE(recoTauPt->at(iTau), recoTauEta->at(iTau), recoTauPhi->at(iTau), recoTauEnergy->at(iTau));
          if (TauCand.DeltaR(Mu1) < 0.8 || TauCand.DeltaR(Mu2) < 0.8) continue;
          if ((recoTauDecayMode->at(iTau) != tauDecayModeThreshold) && (tauDecayModeThreshold == 0 || tauDecayModeThreshold == 1 || tauDecayModeThreshold == 10)) continue;
          Tau1.SetPtEtaPhiE(recoTauPt->at(iTau), recoTauEta->at(iTau), recoTauPhi->at(iTau), recoTauEnergy->at(iTau));
          Tau1Iso = deepTauID ? recoTauDeepVSjetraw->at(iTau) : recoTauIsoMVArawValue->at(iTau);
          Tau1DM = recoTauDecayMode->at(iTau);
          float smallestDR = 0.8; // dR cut between tau1 and tau2
          bool findTau2 = false;
          for (unsigned int iTau2=0; iTau2<recoTauPt->size(); iTau2++)
          {
              if (iTau2 == iTau) continue;
              if (deepTauID && recoTauDeepVSjetraw->size() > 0)
              {
  // -------------------------------------------------------------------------------
                  bool condTauDeepVSeLoose = deepTauVSele == "LOOSE" && recoTauDeepVSeLoose->at(iTau2)>0;
                  bool condTauDeepVSjetLoose = deepTauVSjet == "LOOSE" && recoTauDeepVSjetLoose->at(iTau2)>0;
                  bool condTauDeepVSmuLoose = deepTauVSmu == "LOOSE" && recoTauDeepVSmuLoose->at(iTau2)>0;
                  bool condTauDeepVSeMedium = deepTauVSele == "MEDIUM" && recoTauDeepVSeMedium->at(iTau2)>0;
                  bool condTauDeepVSjetMedium = deepTauVSjet == "MEDIUM" && recoTauDeepVSjetMedium->at(iTau2)>0;
                  bool condTauDeepVSmuMedium = deepTauVSmu == "MEDIUM" && recoTauDeepVSmuMedium->at(iTau2)>0;
                  bool condTauDeepVSeTight = deepTauVSele == "TIGHT" && recoTauDeepVSeTight->at(iTau2)>0;
                  bool condTauDeepVSjetTight = deepTauVSjet == "TIGHT" && recoTauDeepVSjetTight->at(iTau2)>0;
                  bool condTauDeepVSmuTight = deepTauVSmu == "TIGHT" && recoTauDeepVSmuTight->at(iTau2)>0;
                  bool condTauDeepVSeVLoose = deepTauVSele == "VLOOSE" && recoTauDeepVSeVLoose->at(iTau2)>0;
                  bool condTauDeepVSjetVLoose = deepTauVSjet == "VLOOSE" && recoTauDeepVSjetVLoose->at(iTau2)>0;
                  bool condTauDeepVSmuVLoose = deepTauVSmu == "VLOOSE" && recoTauDeepVSmuVLoose->at(iTau2)>0;
                  bool condTauDeepVSeVTight = deepTauVSele == "VTIGHT" && recoTauDeepVSeVTight->at(iTau2)>0;
                  bool condTauDeepVSjetVTight = deepTauVSjet == "VTIGHT" && recoTauDeepVSjetVTight->at(iTau2)>0;
                  bool condTauDeepVSeVVLoose = deepTauVSele == "VVLOOSE" && recoTauDeepVSeVVLoose->at(iTau2)>0;
                  bool condTauDeepVSjetVVLoose = deepTauVSjet == "VVLOOSE" && recoTauDeepVSjetVVLoose->at(iTau2)>0;
                  
                  bool condTauDeepVSeVVTight = deepTauVSele == "VVTIGHT" && recoTauDeepVSeVVTight->at(iTau2)>0;
                  bool condTauDeepVSjetVVTight = deepTauVSjet == "VVTIGHT" && recoTauDeepVSjetVVTight->at(iTau2)>0;
                  bool condTauDeepVSeVVVLoose = deepTauVSele == "VVVLOOSE" && recoTauDeepVSeVVVLoose->at(iTau2)>0;
                  bool condTauDeepVSjetVVVLoose = deepTauVSjet == "VVVLOOSE" && recoTauDeepVSjetVVVLoose->at(iTau2)>0;
                  bool condTauDeepVSeNull = deepTauVSele != "LOOSE" && deepTauVSele != "MEDIUM" && deepTauVSele != "TIGHT" && deepTauVSele != "VLOOSE" && deepTauVSele != "VTIGHT" && deepTauVSele != "VVLOOSE" && deepTauVSele != "VVTIGHT" && deepTauVSele != "VVVLOOSE";
                  bool condTauDeepVSmuNull = deepTauVSmu != "LOOSE" && deepTauVSmu != "MEDIUM" && deepTauVSmu != "TIGHT" && deepTauVSmu != "VLOOSE";
                  // -------------------------------------------------------------------------------
                  bool passCondTauDeepVSele = (condTauDeepVSeLoose || condTauDeepVSeMedium || condTauDeepVSeTight || condTauDeepVSeVLoose || condTauDeepVSeVTight || condTauDeepVSeVVLoose || condTauDeepVSeVVTight || condTauDeepVSeVVVLoose || condTauDeepVSeNull);
                  bool passCondTauDeepVSjet = (condTauDeepVSjetLoose || condTauDeepVSjetMedium || condTauDeepVSjetTight || condTauDeepVSjetVLoose || condTauDeepVSjetVTight || condTauDeepVSjetVVLoose || condTauDeepVSjetVVTight || condTauDeepVSjetVVVLoose);
                  bool passCondTauDeepVSmu = (condTauDeepVSmuLoose || condTauDeepVSmuMedium || condTauDeepVSmuTight || condTauDeepVSmuVLoose || condTauDeepVSmuNull);
                  bool passCondTauDeep = passCondTauDeepVSele && passCondTauDeepVSjet && passCondTauDeepVSmu;

  // -------------------- inverted deep Tau ID -----------------------------
                  bool condInvertTauDeepVSjetLoose = deepTauVSjet == "LOOSE" && recoTauDeepVSjetLoose->at(iTau2)<=0;
                  bool condInvertTauDeepVSjetMedium = deepTauVSjet == "MEDIUM" && recoTauDeepVSjetMedium->at(iTau2)<=0;
                  bool condInvertTauDeepVSjetTight = deepTauVSjet == "TIGHT" && recoTauDeepVSjetTight->at(iTau2)<=0;
                  bool condInvertTauDeepVSjetVLoose = deepTauVSjet == "VLOOSE" && recoTauDeepVSjetVLoose->at(iTau2)<=0;
                  bool condInvertTauDeepVSjetVTight = deepTauVSjet == "VTIGHT" && recoTauDeepVSjetVTight->at(iTau2)<=0;
                  bool condInvertTauDeepVSjetVVLoose = deepTauVSjet == "VVLOOSE" && recoTauDeepVSjetVVLoose->at(iTau2)<=0;
                  bool condInvertTauDeepVSjetVVTight = deepTauVSjet == "VVTIGHT" && recoTauDeepVSjetVVTight->at(iTau2)<=0;
                  bool condInvertTauDeepVSjetVVVLoose = recoTauDeepVSjetVVVLoose->at(iTau2)>0;
                  bool condInvertTauDeepVSeVVVLoose = recoTauDeepVSeVVVLoose->at(iTau2)>0;
                  bool condInvertTauDeepVSmuVLoose = recoTauDeepVSmuVLoose->at(iTau2)>0;
                  // -------------------------------------------------------------------------------
                  bool passCondInvertTauDeepVSjet = ((condInvertTauDeepVSjetLoose || condInvertTauDeepVSjetMedium || condInvertTauDeepVSjetTight || condInvertTauDeepVSjetVLoose || condInvertTauDeepVSjetVTight || condInvertTauDeepVSjetVVLoose || condInvertTauDeepVSjetVVTight) && condInvertTauDeepVSjetVVVLoose && condInvertTauDeepVSeVVVLoose && condInvertTauDeepVSmuVLoose);
                  // -------------------------------------------------------------------------------
                  if ((!invertedTauIso && !passCondTauDeep) || (invertedTauIso && !passCondInvertTauDeepVSjet)) continue;
              } // end if deepTauID && recoTauDeepVSjetraw->size() > 0
              else{
                  bool condTauMVARaw = tauMVAIsoRawORWP == true && recoTauIsoMVArawValue->at(iTau2) > tauMVAIsoRawThreshold;
                  bool condTauMVAWPVVLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VVLOOSE" && recoTauIsoMVAVVLoose->at(iTau2)>0;
                  bool condTauMVAWPVLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VLOOSE" && recoTauIsoMVAVLoose->at(iTau2)>0;
                  bool condTauMVAWPLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "LOOSE" && recoTauIsoMVALoose->at(iTau2)>0;
                  bool condTauMVAWPMedium = tauMVAIsoRawORWP == false && tauMVAIsoWP == "MEDIUM" && recoTauIsoMVAMedium->at(iTau2)>0;
                  bool condTauMVAWPTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "TIGHT" && recoTauIsoMVATight->at(iTau2)>0;
                  bool condTauMVAWPVTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VTIGHT" && recoTauIsoMVAVTight->at(iTau2)>0;
                  bool condTauMVAWPVVTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VVTIGHT" && recoTauIsoMVAVVTight->at(iTau2)>0;
                  bool passCondTauMVA = (condTauMVARaw || condTauMVAWPVVLoose || condTauMVAWPVLoose || condTauMVAWPLoose || condTauMVAWPMedium || condTauMVAWPTight || condTauMVAWPVTight || condTauMVAWPVVTight);
                  // -------------------------------------------------------------------------------------------------

  bool condInvertTauMVARaw = recoTauIsoMVArawValue->at(iTau2) > tauMVAIsoRawThreshold;
                  bool condInvertTauMVAWPVVLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VVLOOSE" && recoTauIsoMVAVVLoose->at(iTau2)<=0;
                  bool condInvertTauMVAWPVLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VLOOSE" && recoTauIsoMVAVLoose->at(iTau2)<=0;
                  bool condInvertTauMVAWPLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "LOOSE" && recoTauIsoMVALoose->at(iTau2)<=0;
                  bool condInvertTauMVAWPMedium = tauMVAIsoRawORWP == false && tauMVAIsoWP == "MEDIUM" && recoTauIsoMVAMedium->at(iTau2)<=0;
                  bool condInvertTauMVAWPTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "TIGHT" && recoTauIsoMVATight->at(iTau2)<=0;
                  bool condInvertTauMVAWPVTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VTIGHT" && recoTauIsoMVAVTight->at(iTau2)<=0;
                  bool condInvertTauMVAWPVVTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VVTIGHT" && recoTauIsoMVAVVTight->at(iTau2)<=0;
                  bool passCondInvertTauMVA = (condInvertTauMVARaw && (condInvertTauMVAWPVVLoose || condInvertTauMVAWPVLoose || condInvertTauMVAWPLoose || condInvertTauMVAWPMedium || condInvertTauMVAWPTight || condInvertTauMVAWPVTight || condInvertTauMVAWPVVTight));
                  // -------------------------------------------------------------------------------------------------
                  bool condTauAntiMuMVALoose = tauAntiMuDisc == "LOOSE" && recoTauAntiMuMVALoose->at(iTau2)>0;
                  bool condTauAntiMuMVATight = tauAntiMuDisc == "TIGHT" && recoTauAntiMuMVATight->at(iTau2)>0; 
                  bool condTauAntiMuMVANull = tauAntiMuDisc != "LOOSE" && tauAntiMuDisc != "TIGHT";
                  bool passCondTauAntiMuMVA = (condTauAntiMuMVALoose || condTauAntiMuMVATight || condTauAntiMuMVANull);
                  // -------------------------------------------------------------------------------------------------
                  bool condTauAntiEleMVALoose = tauAntiEleDisc == "LOOSE" && recoTauAntiEleMVALoose->at(iTau2)>0;
                  bool condTauAntiEleMVAMedium = tauAntiEleDisc == "MEDIUM" && recoTauAntiEleMVAMedium->at(iTau2)>0;
                  bool condTauAntiEleMVATight = tauAntiEleDisc == "TIGHT" && recoTauAntiEleMVATight->at(iTau2)>0; 
                  bool condTauAntiEleMVANull = tauAntiEleDisc != "LOOSE" && tauAntiEleDisc != "MEDIUM" && tauAntiEleDisc != "TIGHT";
                  bool passCondTauAntiEleMVA = (condTauAntiEleMVALoose || condTauAntiEleMVAMedium || condTauAntiEleMVATight || condTauAntiEleMVANull);
                  // -------------------------------------------------------------------------------------------------
                  if ((!invertedTauIso && !passCondTauMVA) || (invertedTauIso && !passCondInvertTauMVA) || !passCondTauAntiMuMVA || !passCondTauAntiEleMVA) continue;
              } // end if !deepTauID (tauMVAID)

                                                                                                                                                                                          if ((recoTauDecayMode->at(iTau2) != tauDecayModeThreshold) && (tauDecayModeThreshold == 0 || tauDecayModeThreshold == 1 || tauDecayModeThreshold == 10)) continue;
              TLorentzVector Tau2Cand; // prepare this variable for dR(tau1, tau2) implementation
              Tau2Cand.SetPtEtaPhiE(recoTauPt->at(iTau2), recoTauEta->at(iTau2), recoTauPhi->at(iTau2), recoTauEnergy->at(iTau2));
              if ((Tau1.DeltaR(Tau2Cand) < smallestDR) && (recoTauPDGId->at(iTau) == (-1) * recoTauPDGId->at(iTau2)) && ((Tau1+Tau2Cand).M() < 60.0) && (Tau2Cand.DeltaR(Mu1) > 0.8) && (Tau2Cand.DeltaR(Mu2) > 0.8))
              {
                  Tau2.SetPtEtaPhiE(recoTauPt->at(iTau2), recoTauEta->at(iTau2), recoTauPhi->at(iTau2), recoTauEnergy->at(iTau2));
                  Tau2Iso = deepTauID ? recoTauDeepVSjetraw->at(iTau2) : recoTauIsoMVArawValue->at(iTau2);
                  Tau2DM = recoTauDecayMode->at(iTau2);
                  smallestDR = Tau1.DeltaR(Tau2);
                  findTau2 = true;
              } // end if find tau2 with tau1 matched
          } // end loop for tau2
          if (!findTau2) continue;
          else{
              findTauTauPair = true;
              break;
          } // end if findTau2
	*/



      bool findTauTauPair = false;
      bool findgenTaugenTauPair = false;

      bool matchrecojetditau = false;
      TLorentzVector GenTauHad;
      TLorentzVector GenTauHad2;
      TLorentzVector matchedjet;     
      TLorentzVector combinedTaus;

      double smallestDRjettau = 0.2;
      float deepvalue;
      float deepvalueMD;



      // ---- prepare event weight info ----                                                                                                                                               
      double weight = 1;
      if (isMC == true)
	{
          weight *= genEventWeight;
	} // end if isMC == true   


      for (unsigned int iGenTau=0; iGenTau<genTauHadPt->size(); iGenTau++){                                                                                         
	indexGenTau = iGenTau; 
	GenTauHad.SetPtEtaPhiM(genTauHadPt->at(iGenTau),genTauHadEta->at(iGenTau),genTauHadPhi->at(iGenTau),genTauHadMass->at(iGenTau));
	//float smallestDRtwotaus = 0.8;
	bool findGenTauHad2 = false;
	for (unsigned int iGenTau2=0; iGenTau2<genTauHadPt->size(); iGenTau2++){
	  TLorentzVector GenTauHad2Cand;
	  if (iGenTau2 != indexGenTau){
	    if (genTauHadPDGId->at(iGenTau) != genTauHadPDGId->at(iGenTau2)){
	      GenTauHad2Cand.SetPtEtaPhiM(genTauHadPt->at(iGenTau2),genTauHadEta->at(iGenTau2),genTauHadPhi->at(iGenTau2),genTauHadMass->at(iGenTau2));
	      //if (GenTauHad2Cand.DeltaR(GenTauHad) < smallestDRtwotaus){
		GenTauHad2.SetPtEtaPhiM(genTauHadPt->at(iGenTau2),genTauHadEta->at(iGenTau2),genTauHadPhi->at(iGenTau2),genTauHadMass->at(iGenTau2));
		//	smallestDRtwotaus = GenTauHad2.DeltaR(GenTauHad);                                                                                                
		dRgenTaugenTau->Fill(GenTauHad.DeltaR(GenTauHad2), weight);
		combinedTaus = GenTauHad+GenTauHad2;
		findGenTauHad2 = true;

		for (unsigned int iDeepDiTaujet=0; iDeepDiTaujet<jet_pt->size(); iDeepDiTaujet++){
		  TLorentzVector DeepDiTaujet;
		  DeepDiTaujet.SetPtEtaPhiE(jet_pt->at(iDeepDiTaujet), jet_eta->at(iDeepDiTaujet), jet_phi->at(iDeepDiTaujet), jet_energy->at(iDeepDiTaujet));
		  if(DeepDiTaujet.DeltaR(combinedTaus) < smallestDRjettau){
		    // if(DeepDiTaujet.Pt() > 50{
			smallestDRjettau = DeepDiTaujet.DeltaR(combinedTaus);
			matchrecojetditau = true;
			deepvalue = DeepDiTauValue->at(iDeepDiTaujet);
			deepvalueMD = DeepDiTauValueMD->at(iDeepDiTaujet);
			matchedjet.SetPtEtaPhiE(jet_pt->at(iDeepDiTaujet), jet_eta->at(iDeepDiTaujet), jet_phi->at(iDeepDiTaujet), jet_energy->at(iDeepDiTaujet));	

			if(deepvalue > 0.4){
			  jetmass_fakelevela->Fill(matchedjet.M());
			}
			if(deepvalue > 0.85){
			  jetmass_fakelevelb->Fill(matchedjet.M());
			}
			if(deepvalue > 0.95){
			  jetmass_fakelevelc->Fill(matchedjet.M());
			}
			if (deepvalueMD > 0.75){
			  jetmass_md_fakelevela->Fill(matchedjet.M());
			}
			if(deepvalueMD > 0.85){
			  jetmass_md_fakelevelb->Fill(matchedjet.M());
			}
			if (deepvalueMD > 0.9){
			  jetmass_md_fakelevelc->Fill(matchedjet.M());
			}

			if (deepvalue > 0.1){
			  jetmass_deep01->Fill(matchedjet.M());
			}

			if (deepvalue > 0.2){
			  jetmass_deep02->Fill(matchedjet.M());
			}
			if (deepvalue > 0.3){
			  jetmass_deep03->Fill(matchedjet.M());
			}
			if (deepvalue > 0.4){
			  jetmass_deep04->Fill(matchedjet.M());
			}
			if (deepvalue > 0.5){
			  jetmass_deep05->Fill(matchedjet.M());
			}
			if (deepvalue > 0.6){
			  jetmass_deep06->Fill(matchedjet.M());
			}
			if (deepvalue > 0.7){
			  jetmass_deep07->Fill(matchedjet.M());
			  genmatched_mumujetMass_07->Fill((Mu1+Mu2+matchedjet).M());
			}
			if (deepvalue > 0.8){
			  jetmass_deep08->Fill(matchedjet.M());
			}
			if (deepvalue > 0.9){
			  jetmass_deep09->Fill(matchedjet.M());
			}

			if (deepvalueMD > 0.1){
			  jetmass_deepmd01->Fill(matchedjet.M());
			}
			if (deepvalueMD > 0.2){
			  jetmass_deepmd02->Fill(matchedjet.M());
			}
			if (deepvalueMD > 0.3){
			  jetmass_deepmd03->Fill(matchedjet.M());
			}
			if (deepvalueMD > 0.4){
			  jetmass_deepmd04->Fill(matchedjet.M());
			}
			if (deepvalueMD > 0.5){
			  jetmass_deepmd05->Fill(matchedjet.M());
			}
			if (deepvalueMD > 0.6){
			  jetmass_deepmd06->Fill(matchedjet.M());
			}
			if (deepvalueMD > 0.7){
			  jetmass_deepmd07->Fill(matchedjet.M());
			}
			if (deepvalueMD > 0.8){
			  jetmass_deepmd08->Fill(matchedjet.M());
			  std::cout << "jet index: " << jet_index->at(iDeepDiTaujet) << std::endl;
			}
			if (deepvalueMD > 0.9){
			  jetmass_deepmd09->Fill(matchedjet.M());
			}

			if (deepvalue > 0.8){
			  //			  std::cout << "deep value > 0.8" << std::endl;
			  highDDTVjetpt->Fill(matchedjet.Pt(),weight);
			  highDDTVtaudR->Fill(GenTauHad.DeltaR(GenTauHad2),weight);
			  highDDTV_genmatched_mumujetMass->Fill((Mu1+Mu2+matchedjet).M(), weight);
			}

			if (deepvalueMD > 0.7){
			  highDDTV_MDjetpt->Fill(matchedjet.Pt(),weight);
			  highDDTV_MDtaudR->Fill(GenTauHad.DeltaR(GenTauHad2),weight);
			  highDDTV_MD_genmatched_mumujetMass->Fill((Mu1+Mu2+matchedjet).M(), weight);
			}

			if (deepvalue < 0.5){
			  lowDDTVjetpt->Fill(matchedjet.Pt(), weight);
			  lowDDTVtaudR->Fill(GenTauHad.DeltaR(GenTauHad2),weight);
			}

			if (deepvalueMD < 0.4){
			  lowDDTV_MDjetpt->Fill(matchedjet.Pt(),weight);
			  lowDDTV_MDtaudR->Fill(GenTauHad.DeltaR(GenTauHad2),weight);
			}

			/*
			if (deepvalue > 0.8 && deepvalueMD < 0.4){
			  highDDTV_lowDDTV_MDjetpt->Fill(matchedjet.Pt(),weight);
			  highDDTV_lowDDTV_MDtaudR->Fill(GenTauHad.DeltaR(GenTauHad2),weight);
			}

			if (deepvalue < 0.5 && deepvalueMD > 0.7){
			  lowDDTV_highDDTV_MDjetpt->Fill(matchedjet.Pt(),weight);
			  lowDDTV_highDDTV_MDtaudR->Fill(GenTauHad.DeltaR(GenTauHad2),weight);
			}

			*/
			DDTV_matched_count++;
			DDTV_MD_matched_count++;

                        if (deepvalue > 0.1){
                          DDTV_01_count++;
                        }

                        if (deepvalue > 0.15){
                          DDTV_015_count++;
                        }

                        if (deepvalue > 0.2){
                          DDTV_02_count++;
                        }

                        if (deepvalue > 0.25){
                          DDTV_025_count++;
                        }

                        if(deepvalue > 0.3){
                          DDTV_03_count++;
                        }

                        if (deepvalue > 0.35){
                          DDTV_035_count++;
                        }

                        if (deepvalue > 0.4){
                          DDTV_04_count++;
                        }

                        if (deepvalue > 0.45){
                          DDTV_045_count++;
                        }


			if (deepvalue > 0.5){
                          DDTV_05_count++;
                        }
                        if (deepvalue > 0.55){
                          DDTV_055_count++;
                        }

                        if (deepvalue > 0.6){
                          DDTV_06_count++;
                        }
                        if (deepvalue > 0.65){
                          DDTV_065_count++;
                        }

                        if (deepvalue > 0.7){
                          DDTV_07_count++;
                        }

                        if (deepvalue > 0.75){
                          DDTV_075_count++;
                        }
                        if(deepvalue > 0.8){
                          DDTV_08_count++;
                        }
                        if (deepvalue > 0.85){
                          DDTV_085_count++;}

			if (deepvalue > 0.9){
                          DDTV_09_count++;
                        }

                        if (deepvalue > 0.95){
                          DDTV_095_count++;
                        }

                        if (deepvalueMD > 0.1){
                          DDTV_MD_01_count++;
                        }
                        if (deepvalueMD > 0.15){
                          DDTV_MD_015_count++;
                        }

                        if (deepvalueMD > 0.2){
                          DDTV_MD_02_count++;
                        }
                        if (deepvalueMD > 0.25){
                          DDTV_MD_025_count++;
                        }

			if(deepvalueMD > 0.3){
                          DDTV_MD_03_count++;
                        }

                        if (deepvalueMD > 0.35){
                          DDTV_MD_035_count++;
                        }

                        if (deepvalueMD > 0.4){
                          DDTV_MD_04_count++;
                        }
                        if (deepvalueMD > 0.45){
                          DDTV_MD_045_count++;
                        }

                        if (deepvalueMD > 0.5){
                          DDTV_MD_05_count++;
                        }
                        if (deepvalueMD > 0.55){
                          DDTV_MD_055_count++;
                        }
                        if (deepvalueMD > 0.6){
                          DDTV_MD_06_count++;
                        }

			if (deepvalueMD > 0.65){
                          DDTV_MD_065_count++;
                        }
                        if (deepvalueMD > 0.7){
                          DDTV_MD_07_count++;
                        }
                        if (deepvalueMD > 0.75){
                          DDTV_MD_075_count++;
                        }
                        if(deepvalueMD > 0.8){
                          DDTV_MD_08_count++;
                        }
                        if (deepvalueMD > 0.85){
                          DDTV_MD_085_count++;
                        }
                        if (deepvalueMD > 0.9){
                          DDTV_MD_09_count++;
                        }
                        if (deepvalueMD > 0.95){
                          DDTV_MD_095_count++;
                        }

		  } // smallest DRjet and two gen taus
		}// loop over DeepDiTau jets
		//	      } // smallest DR between taus
	    } // opp charge
	  } // if different tau
	} // second tau
	if(!findGenTauHad2) continue;
	else{
	  findgenTaugenTauPair = true;
	  break;
	} // end if find GenTauHad2
      
      } //first tau


      /*
      std::cout << "total matched: " << DDTV_matched_count << std::endl;
      std::cout << "ddtv 0.1: " << DDTV_01_count << std::endl;
      std::cout << "ddtv 0.2: " << DDTV_02_count << std::endl;
      std::cout << "ddtv 0.3: " << DDTV_03_count << std::endl;
      std::cout << "ddtv 0.4: " << DDTV_04_count << std::endl;
      std::cout << "ddtv 0.5: " << DDTV_05_count << std::endl;
      std::cout << "ddtv 0.6: " << DDTV_06_count << std::endl;
      std::cout << "ddtv 0.7: " << DDTV_07_count << std::endl;
      std::cout << "ddtv 0.8: " << DDTV_08_count << std::endl;
      std::cout << "ddtv 0.9: " << DDTV_09_count << std::endl;
      */

      /*
      // ---- prepare event weight info ----
      double weight = 1;
      if (isMC == true)
      {
          weight *= genEventWeight; 
      } // end if isMC == true
      */

      //      matchedDeepDiTauValue->Fill(deepvalue);
      //  if (findMu1 && findMu2 && findgenTaugenTauPair && matchrecojetditau){
      if (findgenTaugenTauPair && matchrecojetditau){
	gentauhadpairpt->Fill(combinedTaus.Pt(), weight);
	matchedrecojetpt->Fill(matchedjet.Pt(), weight);
	//dRgenTaugenTau->Fill(GenTauHad.DeltaR(GenTauHad2), weight);
	matchedDeepDiTauValue->Fill(deepvalue);
	matchedDeepDiTauValueMD->Fill(deepvalueMD);
	invMassgenMuMuTauHadTauHad->Fill((GenMu1+GenMu2+GenTauHad+GenTauHad2).M(), weight);
	genmatched_mumujetMass->Fill((Mu1+Mu2+matchedjet).M());

      }
      // ---- fill histograms ----
      if (findMu1 && findMu2 && findTauTauPair)
      {
          ptMu1Mu2->Fill((Mu1+Mu2).Pt(), weight);
          dRMu1Mu2->Fill(Mu1.DeltaR(Mu2), weight);
          invMassMu1Mu2->Fill((Mu1+Mu2).M(), weight);
          dRInvMassMu1Mu2->Fill(Mu1.DeltaR(Mu2), (Mu1+Mu2).M(), weight);

          mu1Iso->Fill(Mu1Iso, weight);
          mu2Iso->Fill(Mu2Iso, weight);

          mu1Pt->Fill(Mu1.Pt(), weight);
          mu1Eta->Fill(Mu1.Eta(), weight);
          mu1Phi->Fill(Mu1.Phi(), weight);

          mu2Pt->Fill(Mu2.Pt(), weight);
          mu2Eta->Fill(Mu2.Eta(), weight);
          mu2Phi->Fill(Mu2.Phi(), weight);

          ptTauTau->Fill((Tau1+Tau2).Pt(), weight);
          dRTauTau->Fill(Tau1.DeltaR(Tau2), weight);
          invMassTauTau->Fill((Tau1+Tau2).M(), weight);
          dRInvMassTauTau->Fill(Tau1.DeltaR(Tau2), (Tau1+Tau2).M(), weight);

          tauIsoMVA->Fill(Tau1Iso, weight);
          tauDecayMode->Fill(Tau1DM, weight);

          tau2IsoMVA->Fill(Tau2Iso, weight);
          tau2DecayMode->Fill(Tau2DM, weight);

          tauPt->Fill(Tau1.Pt(), weight);
          tauEta->Fill(Tau1.Eta(), weight);
          tauPhi->Fill(Tau1.Phi(), weight);
          tauMass->Fill(Tau1.M(), weight);

          tau2Pt->Fill(Tau2.Pt(), weight);
          tau2Eta->Fill(Tau2.Eta(), weight);
          tau2Phi->Fill(Tau2.Phi(), weight);
          tau2Mass->Fill(Tau2.M(), weight);

          dRMu1Tau->Fill(Mu1.DeltaR(Tau1), weight);
          dRMu1Tau2->Fill(Mu1.DeltaR(Tau2), weight);
          dRMu2Tau->Fill(Mu2.DeltaR(Tau1), weight);
          dRMu2Tau2->Fill(Mu2.DeltaR(Tau2), weight);

          ptMuMuTauHadTauHad->Fill((Mu1+Mu2+Tau1+Tau2).Pt(), weight);
          invMassMuMuTauHadTauHad->Fill((Mu1+Mu2+Tau1+Tau2).M(), weight);

	  //  matchedDeepDiTauValue->Fill(deepvalue, weight); 

          // ----- fill flat trees -----
          invMassMuMu = (Mu1+Mu2).M();
          visMassTauTau = (Tau1+Tau2).M();
          visMassMuMuTauTau = (Mu1+Mu2+Tau1+Tau2).M();

          deltaRMuMu = Mu1.DeltaR(Mu2);
          deltaRTauTau = Tau1.DeltaR(Tau2);

          Mu1Pt = Mu1.Pt();
          Mu1Eta = Mu1.Eta();

          Mu2Pt = Mu2.Pt();
          Mu2Eta = Mu2.Eta();

          Tau1Pt = Tau1.Pt();
          Tau1Eta = Tau1.Eta();
          Tau1DecayMode = Tau1DM;
          Tau1Isolation = Tau1Iso;

          Tau2Pt = Tau2.Pt();
          Tau2Eta = Tau2.Eta();
          Tau2DecayMode = Tau2DM;
          Tau2Isolation = Tau2Iso;

          eventWeight = weight/summedWeights;
          TreeMuMuTauTau->Fill();
      } // end if findMu1 && findMu2 && findTauTauPair
   }// end loop for events



   std::cout << "total matched: " << DDTV_matched_count << std::endl;
   std::cout << "ddtv 0.1 : " << DDTV_01_count << std::endl;
   std::cout << "ddtv 0.15: " << DDTV_015_count << std::endl;
   std::cout << "ddtv 0.2 : " << DDTV_02_count << std::endl;
   std::cout << "ddtv 0.25: " << DDTV_025_count << std::endl;
   std::cout << "ddtv 0.3 : " << DDTV_03_count << std::endl;
   std::cout << "ddtv 0.35: " << DDTV_035_count << std::endl;
   std::cout << "ddtv 0.4 : " << DDTV_04_count << std::endl;
   std::cout << "ddtv 0.45: " << DDTV_045_count << std::endl;
   std::cout << "ddtv 0.5 : " << DDTV_05_count << std::endl;
   std::cout << "ddtv 0.55: " << DDTV_055_count << std::endl;
   std::cout << "ddtv 0.6 : " << DDTV_06_count << std::endl;
   std::cout << "ddtv 0.65: " << DDTV_065_count << std::endl;
   std::cout << "ddtv 0.7 : " << DDTV_07_count << std::endl;
   std::cout << "ddtv 0.75: " << DDTV_075_count << std::endl;
   std::cout << "ddtv 0.8 : " << DDTV_08_count << std::endl;
   std::cout << "ddtv 0.85: " << DDTV_085_count <<std::endl;
   std::cout << "ddtv 0.9 : " << DDTV_09_count << std::endl;
   std::cout << "ddtv 0.95: " << DDTV_095_count << std::endl;

   DDTV_count_histo->SetBinContent(1,DDTV_matched_count);
   DDTV_count_histo->SetBinContent(2,DDTV_01_count);
   DDTV_count_histo->SetBinContent(3,DDTV_015_count);
   DDTV_count_histo->SetBinContent(4,DDTV_02_count);
   DDTV_count_histo->SetBinContent(5,DDTV_025_count);
   DDTV_count_histo->SetBinContent(6,DDTV_03_count);
   DDTV_count_histo->SetBinContent(7,DDTV_035_count);
   DDTV_count_histo->SetBinContent(8,DDTV_04_count);
   DDTV_count_histo->SetBinContent(9,DDTV_045_count);
   DDTV_count_histo->SetBinContent(10,DDTV_05_count);
   DDTV_count_histo->SetBinContent(11,DDTV_055_count);
   DDTV_count_histo->SetBinContent(12,DDTV_06_count);
   DDTV_count_histo->SetBinContent(13,DDTV_065_count);
   DDTV_count_histo->SetBinContent(14,DDTV_07_count);
   DDTV_count_histo->SetBinContent(15,DDTV_075_count);
   DDTV_count_histo->SetBinContent(16,DDTV_08_count);
   DDTV_count_histo->SetBinContent(17,DDTV_085_count);
   DDTV_count_histo->SetBinContent(18,DDTV_09_count);
   DDTV_count_histo->SetBinContent(19,DDTV_095_count);

   std::cout << "total matched: " << DDTV_MD_matched_count << std::endl;
   std::cout << "ddtv md 0.1 : " << DDTV_MD_01_count << std::endl;
   std::cout << "ddtv md 0.2 : " << DDTV_MD_02_count << std::endl;
   std::cout << "ddtv md 0.3 : " << DDTV_MD_03_count << std::endl;
   std::cout << "ddtv md 0.4 : " << DDTV_MD_04_count << std::endl;
   std::cout << "ddtv md 0.5 : " << DDTV_MD_05_count << std::endl;
   std::cout << "ddtv md 0.6 : " << DDTV_MD_06_count << std::endl;
   std::cout << "ddtv md 0.7 : " << DDTV_MD_07_count << std::endl;
   std::cout << "ddtv md 0.8 : " << DDTV_MD_08_count << std::endl;
   std::cout << "ddtv md 0.9 : " << DDTV_MD_09_count << std::endl;

   DDTV_MD_count_histo->SetBinContent(1,DDTV_MD_matched_count);
   DDTV_MD_count_histo->SetBinContent(2,DDTV_MD_01_count);
   DDTV_MD_count_histo->SetBinContent(3,DDTV_MD_015_count);
   DDTV_MD_count_histo->SetBinContent(4,DDTV_MD_02_count);
   DDTV_MD_count_histo->SetBinContent(5,DDTV_MD_025_count);
   DDTV_MD_count_histo->SetBinContent(6,DDTV_MD_03_count);
   DDTV_MD_count_histo->SetBinContent(7,DDTV_MD_035_count);
   DDTV_MD_count_histo->SetBinContent(8,DDTV_MD_04_count);
   DDTV_MD_count_histo->SetBinContent(9,DDTV_MD_045_count);
   DDTV_MD_count_histo->SetBinContent(10,DDTV_MD_05_count);
   DDTV_MD_count_histo->SetBinContent(11,DDTV_MD_055_count);
   DDTV_MD_count_histo->SetBinContent(12,DDTV_MD_06_count);
   DDTV_MD_count_histo->SetBinContent(13,DDTV_MD_065_count);
   DDTV_MD_count_histo->SetBinContent(14,DDTV_MD_07_count);
   DDTV_MD_count_histo->SetBinContent(15,DDTV_MD_075_count);
   DDTV_MD_count_histo->SetBinContent(16,DDTV_MD_08_count);
   DDTV_MD_count_histo->SetBinContent(17,DDTV_MD_085_count);
   DDTV_MD_count_histo->SetBinContent(18,DDTV_MD_09_count);
   DDTV_MD_count_histo->SetBinContent(19,DDTV_MD_095_count);



   std::cout << "vvloose Deep: " << DeepVsjetvvloose_count << std::endl;
   std::cout << "vloose Deep: " << DeepVsjetvloose_count << std::endl;
   std::cout << "loose Deep: " << DeepVsjetloose_count << std::endl;
   std::cout << "medium Deep: " << DeepVsjetmedium_count << std::endl;
   std::cout << "tight Deep: " << DeepVsjettight_count << std::endl;
   std::cout << "vtight Deep: " << DeepVsjetvtight_count << std::endl;


   outputFile->cd();

   int numberofhist = histColl.size();
   for(int i=0; i<numberofhist; i++){
     //       if (isMC) histColl[i]->Scale(lumiScale/summedWeights);
     histColl[i]->Write();
     //     outputFile->WriteObject(histColl[i],"OverWrite");
     //    histColl[i]->Write(),"OverWrite");
      } // end loop for writing all the histograms into the output file

   for(int j=0; j<numberofhist; j++){
       delete histColl[j];
   } // end loop for deleting all the histograms

   TreeMuMuTauTau->Write("TreeMuMuTauTau", TObject::kOverwrite);
   //TreeMuMuTauTau->WriteObject("TreeMuMuTauTau", "OverWrite");  
 outputFile->Close();
}
