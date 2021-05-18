#define MuMuTauHadAnalyzer_cxx
#include "MuMuTauHadAnalyzer.h"
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

void MuMuTauHadAnalyzer::Loop()
{
   TString outputfileName = createOutputFileName();
   TFile* outputFile = new TFile(outputfileName, "RECREATE");
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   if (nMaxEvents >= 0 && nMaxEvents  < nentries) nentries = nMaxEvents;
   cout << "We will run on " << nentries << " events" << endl;
   Long64_t nbytes = 0, nb = 0;

   bool matchRecGen  = true;

   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      if (jentry % 1000 == 0 && jentry > 0) cout << "*** Processing #Events: " << jentry << endl;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

      // ---- define varibles that will be used to be filled into histograms ---
      TLorentzVector Mu1;
      TLorentzVector Mu2;
      // TLorentzVector Mu3;
      TLorentzVector Tau;

      float Mu1Iso;
      float Mu2Iso;
      // float Mu3Iso;
      float TauIso;
      float TauDM;

      float TauMVAIso;
      float tauDeepVSjetrawscore;
      unsigned int indexMu1 = -1;
      unsigned int indexMu2 = -1;


      // ============================================================================

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

      if (!findMu1) continue;
      float smallestDR = 1.0; // dR cut between Mu1 and Mu2
      bool findMu2 = false;
      // ---- start loop on muon candidates for mu2 ----
      for (unsigned int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
	{
          if (iMuon == indexMu1) continue;
          if ((invertedMu2Iso == false && recoMuonIsolation->at(iMuon) > Mu2IsoThreshold) || (invertedMu2Iso == true && recoMuonIsolation->at(iMuon) < Mu2IsoThreshold)) continue;

          TLorentzVector Mu2Cand; // prepare this variable for dR(Mu1,Mu2) implementation
          Mu2Cand.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));

          if((Mu1.DeltaR(Mu2Cand) < smallestDR) && (recoMuonPDGId->at(indexMu1) == (-1) * recoMuonPDGId->at(iMuon)) && ((Mu1+Mu2Cand).M() > diMuonMassLowThreshold) && ((Mu1+Mu2Cand).M() < diMuonMassHighThreshold))
	    {
              Mu2.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
              Mu2Iso = recoMuonIsolation->at(iMuon);
              smallestDR = Mu1.DeltaR(Mu2);
              indexMu2 = iMuon;
              findMu2 = true;
	    } // end if pair candidates
	} // end loop for mu2
          
      if (!findMu2) continue;

      //  bool findMuTauPair = false;
      // ------- start loop on tau candidates -------
      bool findTauHad = false;   
      for (unsigned int iTau=0; iTau<recoTauPt->size(); iTau++)
	{
	  TLorentzVector TauCand;
	  TauCand.SetPtEtaPhiE(recoTauPt->at(iTau), recoTauEta->at(iTau), recoTauPhi->at(iTau), recoTauEnergy->at(iTau));

	  if (TauCand.DeltaR(Mu1) < 0.8 || TauCand.DeltaR(Mu2) < 0.8) continue;
	  if ((recoTauDecayMode->at(iTau) != tauDecayModeThreshold) && (tauDecayModeThreshold == 0 || tauDecayModeThreshold == 1 || tauDecayModeThreshold == 5 || tauDecayModeThreshold == 6 ||tauDecayModeThreshold == 10)) continue;
	   if(recoTauPt->at(iTau) < 30){
	  Tau.SetPtEtaPhiE(recoTauPt->at(iTau), recoTauEta->at(iTau), recoTauPhi->at(iTau), recoTauEnergy->at(iTau));
	  TauIso = deepTauID ? recoTauDeepVSjetraw->at(iTau) : recoTauIsoMVArawValue->at(iTau);
	  TauDM = recoTauDecayMode->at(iTau);

	  // float smallestDR = 0.8; // dR cut between Mu3 and tau
	  findTauHad = true;
	  } // pT under 30 requirement
	}// end loop for tau
	
      // ---- prepare event weight info ----
      double weight = 1;
      if (isMC == true)
      {
          weight *= genEventWeight; 
      } // end if isMC == true

      //      cout <<"find taud had: " << findTauHad << std::endl;
      // ---- fill histograms ----
      if (findMu1 && findMu2 && findTauHad)

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
	  	  
	  /*
	  tauIsoMVA->Fill(TauIso, weight);
	  recoTauDeepVSjetrawscore->Fill(TauIso, weight);
	  tauDecayMode->Fill(TauDM, weight);
	  */
          tauPt->Fill(Tau.Pt(), weight);
          tauEta->Fill(Tau.Eta(), weight);
          tauPhi->Fill(Tau.Phi(), weight);
          tauMass->Fill(Tau.M(), weight);


          dRMu1Tau->Fill(Mu1.DeltaR(Tau), weight);
          dRMu2Tau->Fill(Mu2.DeltaR(Tau), weight);


          // ----- fill flat trees -----
          invMassMuMu = (Mu1+Mu2).M();

	  // visMassMuMuTauTau = (Mu1+Mu2+Mu3+Tau).M();

          deltaRMuMu = Mu1.DeltaR(Mu2);


          Mu1Pt = Mu1.Pt();
          Mu1Eta = Mu1.Eta();

          Mu2Pt = Mu2.Pt();
          Mu2Eta = Mu2.Eta();

          Tau2Pt = Tau.Pt();
          Tau2Eta = Tau.Eta();
          Tau2DecayMode = TauDM;
          Tau2Isolation = TauIso;

          eventWeight = weight/summedWeights;
          TreeMuMuTauTau->Fill();


	  if (isMC && matchRecGen)
	    {
              TLorentzVector GenMu1;
              TLorentzVector GenMu2;
              //TLorentzVector GenMu3;
              //TLorentzVector GenTauMu;
              TLorentzVector GenTauHad;

              bool findMatchedRecGenMu1 = false;
              bool findMatchedRecGenMu2 = false;
	      //  bool findMatchedRecGenMu3 = false;
	      //              bool findMatchedRecGenTauMu = false;
              bool findMatchedRecGenTauHad = false;

	      //double GenTauHadVisiblePt = 0;
              unsigned int indexGenMu1 = -1;
              unsigned int indexGenMu2 = -1;


                  // --------- search for matched genMu1 for Mu1 --------------
	      if (genMuonPt->size()>0)
		{
                  double smallestDR = 0.15;
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
		} // gen Muonpt->Size>0


                  // --------- search for matched genTauHad for Tau --------------
	      if (genTauHadPt->size()>0)
		{
                  smallestDR = 0.15;
		  for (unsigned int iGenTauHad=0; iGenTauHad<genTauHadPt->size(); iGenTauHad++)
		    {
                      TLorentzVector GenTauHadCand;
                      GenTauHadCand.SetPtEtaPhiM(genTauHadPt->at(iGenTauHad), genTauHadEta->at(iGenTauHad), genTauHadPhi->at(iGenTauHad), genTauHadMass->at(iGenTauHad));
                      if (Tau.DeltaR(GenTauHadCand) <= smallestDR)
			{
                          smallestDR = Tau.DeltaR(GenTauHadCand);
                          findMatchedRecGenTauHad = true;
                          GenTauHad = GenTauHadCand;
			  //                          GenTauHadVisiblePt = genTauHadVisPt->at(iGenTauHad);
			} // end if Tau.DeltaR(GenTauHad) <= smallestDR
		    }// end for loop on GenTauHad
		} // end if size of genTauHad->Pt >0 
		 
	      if(findMatchedRecGenMu1 && findMatchedRecGenMu2 && findMatchedRecGenTauHad){
		genmuPt->Fill(GenMu1.Pt(), weight);
		genmuPt->Fill(GenMu2.Pt(),weight);
	
		genmu1Pt->Fill(GenMu1.Pt(), weight);
		genmu1Eta->Fill(GenMu1.Eta(), weight);
		genmu1Phi->Fill(GenMu1.Phi(), weight);
		genmu1Mass->Fill(GenMu1.M(), weight);
		mu1PtVSGenMu1Pt->Fill(Mu1.Pt(), GenMu1.Pt(), weight);
		mu1EtaVSGenMu1Eta->Fill(Mu1.Eta(), GenMu1.Eta(), weight);
		mu1PhiVSGenMu1Phi->Fill(Mu1.Phi(), GenMu1.Phi(), weight);

		genmu2Pt->Fill(GenMu2.Pt(), weight);
		genmu2Eta->Fill(GenMu2.Eta(), weight);
		genmu2Phi->Fill(GenMu2.Phi(), weight);
		genmu2Mass->Fill(GenMu2.M(), weight);
		mu2PtVSGenMu2Pt->Fill(Mu2.Pt(), GenMu2.Pt(), weight);
     		mu2EtaVSGenMu2Eta->Fill(Mu2.Eta(), GenMu2.Eta(), weight);
		mu2PhiVSGenMu2Phi->Fill(Mu2.Phi(), GenMu2.Phi(), weight);
	         
      	       	tauIsoMVA->Fill(TauIso, weight);
		recoTauDeepVSjetrawscore->Fill(TauIso, weight);
		tauDecayMode->Fill(TauDM, weight);
		

		gentauPt->Fill(GenTauHad.Pt(), weight);
		gentauEta->Fill(GenTauHad.Eta(), weight);
		gentauPhi->Fill(GenTauHad.Phi(), weight);
		gentauMass->Fill(GenTauHad.M(), weight);

		tauPtVSGenTauHadPt->Fill(Tau.Pt(), GenTauHad.Pt(), weight);
		tauEtaVSGenTauHadEta->Fill(Tau.Eta(), GenTauHad.Eta(), weight);
		tauPhiVSGenTauHadPhi->Fill(Tau.Phi(), GenTauHad.Phi(), weight);

		dRgenMu1genMu2->Fill(GenMu1.DeltaR(GenMu2), weight);

		dRgenMu1genTau->Fill(GenMu1.DeltaR(GenTauHad), weight);
	
		dRgenMu2genTau->Fill(GenMu2.DeltaR(GenTauHad), weight);
				
		ptgenMu1genMu2->Fill((GenMu1+GenMu2).Pt(), weight);

		dRInvMassgenMu1genMu2->Fill(GenMu1.DeltaR(GenMu2), (GenMu1+GenMu2).M(), weight);

	      } // if all gen reco matching satisfied
	    }//  end if MC & match genReco
	}// end if findMu1 && findMu2
   }// end loop for events
   
   outputFile->cd();

   int numberofhist = histColl.size();
   for(int i=0; i<numberofhist; i++){
       if (isMC) histColl[i]->Scale(lumiScale/summedWeights);
       histColl[i]->Write();
   } // end loop for writing all the histograms into the output file

   for(int j=0; j<numberofhist; j++){
       delete histColl[j];
   } // end loop for deleting all the histograms

   TreeMuMuTauTau->Write("TreeMuMuTauTau", TObject::kOverwrite);
   outputFile->Close();
}
