#define MuMuTauMuTauHadAnalyzer_cxx
#include "MuMuTauMuTauHadAnalyzer.h"
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

void MuMuTauMuTauHadAnalyzer::Loop()
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
      TLorentzVector Mu3;
      TLorentzVector Tau;

      float Mu1Iso;
      float Mu2Iso;
      float Mu3Iso;
      float TauIso;
      float TauDM;

      unsigned int indexMu1 = -1;
      unsigned int indexMu2 = -1;

      //      TLorentzVector genMu1;
      //TLorentzVector genMu2;
      //TLorentzVector genMu3;
      //TLorentzVector genTau;


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

      bool findMuTauPair = false;
      // ------- start loop on tau candidates -------
      for (unsigned int iTau=0; iTau<recoTauPt->size(); iTau++)
	{

	bool condTauMVARaw = tauMVAIsoRawORWP == true && recoTauIsoMVArawValue->at(iTau) > tauMVAIsoRawThreshold;
	bool condTauMVAWPVVLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VVLOOSE" && recoTauIsoMVAVVLoose->at(iTau)>0;
	bool condTauMVAWPVLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VLOOSE" && recoTauIsoMVAVLoose->at(iTau)>0;
	bool condTauMVAWPLoose = tauMVAIsoRawORWP == false && tauMVAIsoWP == "LOOSE" && recoTauIsoMVALoose->at(iTau)>0;
	bool condTauMVAWPMedium = tauMVAIsoRawORWP == false && tauMVAIsoWP == "MEDIUM" && recoTauIsoMVAMedium->at(iTau)>0;
	bool condTauMVAWPTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "TIGHT" && recoTauIsoMVATight->at(iTau)>0;
	bool condTauMVAWPVTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VTIGHT" && recoTauIsoMVAVTight->at(iTau)>0;
	bool condTauMVAWPVVTight = tauMVAIsoRawORWP == false && tauMVAIsoWP == "VVTIGHT" && recoTauIsoMVAVVTight->at(iTau)>0;

	bool passCondTauMVA = (condTauMVARaw || condTauMVAWPVVLoose || condTauMVAWPVLoose || condTauMVAWPLoose || condTauMVAWPMedium || condTauMVAWPTight || condTauMVAWPVTight || condTauMVAWPVVTight);

	if (!condTauMVAWPMedium) continue;
          TLorentzVector TauCand;
          TauCand.SetPtEtaPhiE(recoTauPt->at(iTau), recoTauEta->at(iTau), recoTauPhi->at(iTau), recoTauEnergy->at(iTau));
          
          if (TauCand.DeltaR(Mu1) < 0.8 || TauCand.DeltaR(Mu2) < 0.8) continue;
          Tau.SetPtEtaPhiE(recoTauPt->at(iTau), recoTauEta->at(iTau), recoTauPhi->at(iTau), recoTauEnergy->at(iTau));
          TauIso = recoTauIsoMVArawValue->at(iTau);
          TauDM = recoTauDecayMode->at(iTau);

          float smallestDR = 0.8; // dR cut between Mu3 and tau
          bool findMu3 = false;

          for (unsigned int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
	    {
              if (iMuon == indexMu1 || iMuon == indexMu2) continue;

              TLorentzVector Mu3Cand; // prepare this variable for dR(Mu3, tau) implementation
              Mu3Cand.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
              if ((Tau.DeltaR(Mu3Cand) < smallestDR) && (recoTauPDGId->at(iTau)/fabs(recoTauPDGId->at(iTau)) == (-1) * recoMuonPDGId->at(iMuon)/fabs(recoMuonPDGId->at(iMuon))) && ((Tau+Mu3Cand).M() < 60.0) && (Mu3Cand.DeltaR(Mu1) > 0.4) && (Mu3Cand.DeltaR(Mu2) > 0.4))
		{
		  Mu3.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
                  Mu3Iso = recoMuonIsolation->at(iMuon);
                  smallestDR = Tau.DeltaR(Mu3);
                  findMu3 = true;
		} // end if find mu3 with tau matched
	    } // end loop for mu3

          if (!findMu3) continue;
	  else{
	    findMuTauPair = true;
	    break;
          } // end if findMu3
	} // end loop for tau

      // ---- prepare event weight info ----
      double weight = 1;
      if (isMC == true)
      {
          weight *= genEventWeight; 
      } // end if isMC == true

      // ---- fill histograms ----
      if (findMu1 && findMu2 && findMuTauPair)
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

          ptMu3Tau->Fill((Mu3+Tau).Pt(), weight);
          dRMu3Tau->Fill(Mu3.DeltaR(Tau), weight);
          invMassMu3Tau->Fill((Mu3+Tau).M(), weight);
          dRInvMassMu3Tau->Fill(Mu3.DeltaR(Tau), (Mu3+Tau).M(), weight);

          mu3Iso->Fill(Mu3Iso, weight);
          tauIsoMVA->Fill(TauIso, weight);
          tauDecayMode->Fill(TauDM, weight);

          mu3Pt->Fill(Mu3.Pt(), weight);
          mu3Eta->Fill(Mu3.Eta(), weight);
          mu3Phi->Fill(Mu3.Phi(), weight);

          tauPt->Fill(Tau.Pt(), weight);
          tauEta->Fill(Tau.Eta(), weight);
          tauPhi->Fill(Tau.Phi(), weight);
          tauMass->Fill(Tau.M(), weight);

          dRMu1Mu3->Fill(Mu1.DeltaR(Mu3), weight);
          dRMu1Tau->Fill(Mu1.DeltaR(Tau), weight);
          dRMu2Mu3->Fill(Mu2.DeltaR(Mu3), weight);
          dRMu2Tau->Fill(Mu2.DeltaR(Tau), weight);

          ptMuMuTauMuTauHad->Fill((Mu1+Mu2+Mu3+Tau).Pt(), weight);
          invMassMuMuTauMuTauHad->Fill((Mu1+Mu2+Mu3+Tau).M(), weight);

	  if (isMC && matchRecGen)
	    {
              TLorentzVector GenMu1;
              TLorentzVector GenMu2;
              TLorentzVector GenMu3;
              TLorentzVector GenTauMu;
              TLorentzVector GenTauHad;

              bool findMatchedRecGenMu1 = false;
              bool findMatchedRecGenMu2 = false;
	      bool findMatchedRecGenMu3 = false;
              bool findMatchedRecGenTauMu = false;
              bool findMatchedRecGenTauHad = false;

	      double GenTauHadVisiblePt = 0;
              unsigned int indexGenMu1 = -1;
              unsigned int indexGenMu2 = -1;

              if (genMuonPt->size()>0)
		{
                  // --------- search for matched genMu1 for Mu1 --------------
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
		    } // end for loop on GenMu2

		 
                  // --------- search for matched genMu3 for Mu3 --------------
		  smallestDR = 0.15;
		  for (unsigned int iGenMu=0; iGenMu<genMuonPt->size(); iGenMu++)
		    {
		      TLorentzVector GenMuCand;
		      GenMuCand.SetPtEtaPhiM(genMuonPt->at(iGenMu), genMuonEta->at(iGenMu), genMuonPhi->at(iGenMu), genMuonMass->at(iGenMu));
		      if (Mu3.DeltaR(GenMuCand) <= smallestDR && iGenMu != indexGenMu1 && iGenMu != indexGenMu2)
			{
			  smallestDR = Mu3.DeltaR(GenMuCand);
			  findMatchedRecGenMu3 = true;
			  GenMu3 = GenMuCand;
			} // end if Mu3.DeltaR(GenMuCand) <= smallestDR && iGenMu != indexGenMu1 && iGenMu != indexGenMu2
		    } // end for loop on GenMu3
		} // end if genMuonPt->size()>0


	      if (genTauMuPt->size()>0)
		{
                  // --------- search for matched genTauMu for Mu3 --------------
                  double smallestDR = 0.15;
                  for (unsigned int iGenTauMu=0; iGenTauMu<genTauMuPt->size(); iGenTauMu++)
		    {
                      TLorentzVector GenTauMuCand;
                      GenTauMuCand.SetPtEtaPhiM(genTauMuPt->at(iGenTauMu), genTauMuEta->at(iGenTauMu), genTauMuPhi->at(iGenTauMu), genTauMuMass->at(iGenTauMu));
                      if (Mu3.DeltaR(GenTauMuCand) <= smallestDR)
			{
                          smallestDR = Mu3.DeltaR(GenTauMuCand);
                          findMatchedRecGenTauMu = true;
                          GenTauMu = GenTauMuCand;
			} // end if Mu3.DeltaR(GenTauMuCand) <= smallestDR
		    } // end for loop on GenTauMu
		} // end if genTauMuPt->size()>0

	      if (genTauHadPt->size()>0)
		{
                  // --------- search for matched genTauHad for Tau --------------
                  double smallestDR = 0.15;
		  for (unsigned int iGenTauHad=0; iGenTauHad<genTauHadPt->size(); iGenTauHad++)
		    {
                      TLorentzVector GenTauHadCand;
                      GenTauHadCand.SetPtEtaPhiM(genTauHadPt->at(iGenTauHad), genTauHadEta->at(iGenTauHad), genTauHadPhi->at(iGenTauHad), genTauHadMass->at(iGenTauHad));
                      if (Tau.DeltaR(GenTauHadCand) <= smallestDR)
			{
                          smallestDR = Tau.DeltaR(GenTauHadCand);
                          findMatchedRecGenTauHad = true;
                          GenTauHad = GenTauHadCand;
                          GenTauHadVisiblePt = genTauHadVisPt->at(iGenTauHad);
			} // end if Tau.DeltaR(GenTauHad) <= smallestDR
		    }// end for loop on GenTauHad
		} // end if genTauHadPt->size()>0

	      if(findMatchedRecGenMu1 && findMatchedRecGenMu2 && findMatchedRecGenTauMu && findMatchedRecGenTauHad){
		genmuPt->Fill(GenMu1.Pt(), weight);
		genmuPt->Fill(GenMu2.Pt(),weight);
		genmuPt->Fill(GenMu3.Pt(), weight);

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
	   
		genmu3Pt->Fill(GenMu3.Pt(), weight);
		genmu3Eta->Fill(GenMu3.Eta(), weight);
		genmu3Phi->Fill(GenMu3.Phi(), weight);
		genmu3Mass->Fill(GenMu3.M(), weight);
		mu3PtVSGenMu3Pt->Fill(Mu3.Pt(), GenMu3.Pt(), weight);                                                                   
		mu3EtaVSGenMu3Eta->Fill(Mu3.Eta(), GenMu3.Eta(), weight);                                                               
		mu3PhiVSGenMu3Phi->Fill(Mu3.Phi(), GenMu3.Phi(), weight); 

		gentauMuPt->Fill(GenTauMu.Pt(), weight);
		gentauMuEta->Fill(GenTauMu.Eta(), weight);
		gentauMuPhi->Fill(GenTauMu.Phi(), weight);
		gentauMuMass->Fill(GenTauMu.M(), weight);

		tauPtVSGenTauHadVisPt->Fill(Tau.Pt(), GenTauHadVisiblePt, weight);

		mu3PtVSGenMu3Pt->Fill(Mu3.Pt(), GenMu3.Pt(), weight);
		mu3EtaVSGenMu3Eta->Fill(Mu3.Eta(), GenMu3.Eta(), weight);
		mu3PhiVSGenMu3Phi->Fill(Mu3.Phi(), GenMu3.Phi(), weight);

		gentauPt->Fill(GenTauHad.Pt(), weight);
		gentauEta->Fill(GenTauHad.Eta(), weight);
		gentauPhi->Fill(GenTauHad.Phi(), weight);
		gentauMass->Fill(GenTauHad.M(), weight);
		recomatchedtauPt->Fill(Tau.Pt(), weight);
		tauPtVSGenTauHadPt->Fill(Tau.Pt(), GenTauHad.Pt(), weight);
		tauEtaVSGenTauHadEta->Fill(Tau.Eta(), GenTauHad.Eta(), weight);
		tauPhiVSGenTauHadPhi->Fill(Tau.Phi(), GenTauHad.Phi(), weight);

		dRgenMu1genMu2->Fill(GenMu1.DeltaR(GenMu2), weight);
		dRgenMu3genTau->Fill(GenMu3.DeltaR(GenTauHad), weight);

		dRgenMu1genMu3->Fill(GenMu1.DeltaR(GenMu3), weight);
		dRgenMu1genTau->Fill(GenMu1.DeltaR(GenTauHad), weight);
		dRgenMu2genMu3->Fill(GenMu2.DeltaR(GenMu3), weight);
		dRgenMu2genTau->Fill(GenMu2.DeltaR(GenTauHad), weight);
		dRgenTaugenTau->Fill(GenTauMu.DeltaR(GenTauHad), weight);

		invMassgenMu1genMu2->Fill((GenMu1+GenMu2).M(),weight);
		invMassgenMu3genTau->Fill((GenMu3+GenTauHad).M(), weight);
		invMassgenTaugenTau->Fill((GenTauMu+GenTauHad).M(), weight);

		invMassgenMuMuTauMuTauHad->Fill((GenMu1+GenMu2+GenMu3+GenTauHad).M(), weight);
				
		ptgenMu1genMu2->Fill((GenMu1+GenMu2).Pt(), weight);
		ptgenMu3genTau->Fill((GenMu3+GenTauHad).Pt(), weight);
		ptgenTaugenTau->Fill((GenTauMu+GenTauHad).Pt(), weight);

		ptgenMuMuTauMuTauHad->Fill((GenMu1+GenMu2+GenMu3+GenTauHad).Pt(), weight);
		
		dRInvMassgenMu1genMu2->Fill(GenMu1.DeltaR(GenMu2), (GenMu1+GenMu2).M(), weight);
		dRInvMassgenMu3genTau->Fill(GenMu3.DeltaR(GenTauHad), (GenMu3+GenTauHad).M(), weight);
		dRInvMassgenTaugenTau->Fill(GenTauMu.DeltaR(GenTauHad), (GenTauMu+GenTauHad).M(), weight);

	      } // if all gen reco matching satisfied
	    } //  end if MC & match genReco
	} // end if findMu1 && findMu2 && findMuTauPair
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

   outputFile->Close();
}
