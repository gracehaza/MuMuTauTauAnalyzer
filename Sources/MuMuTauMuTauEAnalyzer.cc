#define MuMuTauMuTauEAnalyzer_cxx
#include "MuMuTauMuTauEAnalyzer.h"
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

void MuMuTauMuTauEAnalyzer::Loop()
{
   TString outputfileName = createOutputFileName();
   TFile* outputFile = new TFile(outputfileName, "RECREATE");
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   if (nMaxEvents >= 0 && nMaxEvents  < nentries) nentries = nMaxEvents;
   cout << "We will run on " << nentries << " events" << endl;

   Long64_t nbytes = 0, nb = 0;

   bool matchRecGen = true;

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
      TLorentzVector Ele;

      float Mu1Iso;
      float Mu2Iso;
      float Mu3Iso;
      float EleIso;

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
      bool findMuElePair = false;

      // ------- start loop on electron candidates -------
      for (unsigned int iEle=0; iEle<recoElectronPt->size(); iEle++)
      {
          if ((invertedEle1Iso == false && recoElectronIsolation->at(iEle) > Ele1IsoThreshold) || (invertedEle1Iso == true && recoElectronIsolation->at(iEle) < Ele1IsoThreshold)) continue;
          TLorentzVector EleCand;
          EleCand.SetPtEtaPhiE(recoElectronPt->at(iEle), recoElectronEta->at(iEle), recoElectronPhi->at(iEle), recoElectronEnergy->at(iEle));

          if (EleCand.DeltaR(Mu1) < 0.4 || EleCand.DeltaR(Mu2) < 0.4) continue;
          Ele.SetPtEtaPhiE(recoElectronPt->at(iEle), recoElectronEta->at(iEle), recoElectronPhi->at(iEle), recoElectronEnergy->at(iEle));
          EleIso = recoElectronIsolation->at(iEle);

          float smallestDR = 1.0; // dR cut between Mu3 and electron
          bool findMu3 = false;

          for (unsigned int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
          {
              if (iMuon == indexMu1 || iMuon == indexMu2) continue;

              TLorentzVector Mu3Cand; // prepare this variable for dR(Mu3, electron) implementation
              Mu3Cand.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
              if ((Ele.DeltaR(Mu3Cand) < smallestDR) && (recoElectronPDGId->at(iEle)/fabs(recoElectronPDGId->at(iEle)) == (-1) * recoMuonPDGId->at(iMuon)/fabs(recoMuonPDGId->at(iMuon))) && ((Ele+Mu3Cand).M() < 60.0) && (Mu3Cand.DeltaR(Mu1) > 0.4) && (Mu3Cand.DeltaR(Mu2) > 0.4))
              {
                  Mu3.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
                  Mu3Iso = recoMuonIsolation->at(iMuon);
                  smallestDR = Ele.DeltaR(Mu3);
                  findMu3 = true;
              } // end if find mu3 with electron matched
          } // end loop for mu3

          if (!findMu3) continue;
          else{
              findMuElePair = true;
              break;
          } // end if findMu3
      } // end loop for electron

      // ---- prepare event weight info ----
      double weight = 1;
      if (isMC == true)
      {
          weight *= genEventWeight; 
      } // end if isMC == true

      // ---- fill histograms ----
      if (findMu1 && findMu2 && findMuElePair)
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

          ptMu3Ele->Fill((Mu3+Ele).Pt(), weight);
          dRMu3Ele->Fill(Mu3.DeltaR(Ele), weight);
          invMassMu3Ele->Fill((Mu3+Ele).M(), weight);
          dRInvMassMu3Ele->Fill(Mu3.DeltaR(Ele), (Mu3+Ele).M(), weight);

          mu3Iso->Fill(Mu3Iso, weight);
          ele1Iso->Fill(EleIso, weight);

          mu3Pt->Fill(Mu3.Pt(), weight);
          mu3Eta->Fill(Mu3.Eta(), weight);
          mu3Phi->Fill(Mu3.Phi(), weight);

          ele1Pt->Fill(Ele.Pt(), weight);
          ele1Eta->Fill(Ele.Eta(), weight);
          ele1Phi->Fill(Ele.Phi(), weight);

          dRMu1Mu3->Fill(Mu1.DeltaR(Mu3), weight);
          dRMu1Ele1->Fill(Mu1.DeltaR(Ele), weight);
          dRMu2Mu3->Fill(Mu2.DeltaR(Mu3), weight);
          dRMu2Ele1->Fill(Mu2.DeltaR(Ele), weight);

          ptMuMuTauMuTauEle->Fill((Mu1+Mu2+Mu3+Ele).Pt(), weight);
          invMassMuMuTauMuTauEle->Fill((Mu1+Mu2+Mu3+Ele).M(), weight);

          // ----- fill flat trees -----
          invMassMuMu = (Mu1+Mu2).M();
          visMassTauTau = (Mu3+Ele).M();
          visMassMuMuTauTau = (Mu1+Mu2+Mu3+Ele).M();

          deltaRMuMu = Mu1.DeltaR(Mu2);
          deltaRTauTau = Mu3.DeltaR(Ele);

          Mu1Pt = Mu1.Pt();
          Mu1Eta = Mu1.Eta();

          Mu2Pt = Mu2.Pt();
          Mu2Eta = Mu2.Eta();

          Tau1Pt = Mu3.Pt();
          Tau1Eta = Mu3.Eta();
          Tau1Isolation = Mu3Iso;

          Tau2Pt = Ele.Pt();
          Tau2Eta = Ele.Eta();
          Tau2Isolation = EleIso;

          eventWeight = weight/summedWeights;
          TreeMuMuTauTau->Fill();

          if (isMC && matchRecGen)
            {
              TLorentzVector GenMu1;
              TLorentzVector GenMu2;
              TLorentzVector GenMu3;
	      TLorentzVector GenEle;
              TLorentzVector GenTauEle;
              TLorentzVector GenTauMu;

              bool findMatchedRecGenMu1 = false;
              bool findMatchedRecGenMu2 = false;
              bool findMatchedRecGenMu3 = false;
	      bool findMatchedRecGenEle = false;
              bool findMatchedRecGenTauMu = false;
              bool findMatchedRecGenTauEle = false;

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
		} // if genMuonPt->size() > 0


	      // --------- search for matched genEle for Ele -------------- 
	      if(genElectronPt->size() > 0){
		double smallestDR = 0.15;
		for (unsigned int iGenEle=0; iGenEle<genElectronPt->size(); iGenEle++)
		  {
		    TLorentzVector GenEleCand;
		    GenEleCand.SetPtEtaPhiM(genElectronPt->at(iGenEle), genElectronEta->at(iGenEle), genElectronPhi->at(iGenEle), genElectronMass->at(iGenEle));
		    if (Ele.DeltaR(GenEleCand) <= smallestDR)
		      {
			smallestDR = Ele.DeltaR(GenEleCand);
			findMatchedRecGenEle = true;
			GenEle = GenEleCand;
		      } // end if Ele.DeltaR(GenMuCand) <= smallestDR && iGenMu != indexGenMu1 && iGenMu != indexGenMu2
		  } // end for loop on GenEle     
	      } // end  recoElectronPt->size() > 0

	      // --------- search for matched genTauMu for Mu3 --------------    
	      if (genTauMuPt->size()>0)
		{
		  smallestDR = 0.15;
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
		} // end if genTauMuPt->size() > 0

	      // --------- search for matched genTauEle for Ele --------------                                                                                         
	      if (genTauElePt->size()>0)
                {
		  smallestDR = 0.15;
		  for (unsigned int iGenTauEle=0; iGenTauEle<genTauElePt->size(); iGenTauEle++)
		    {
		      TLorentzVector GenTauEleCand;
		      GenTauEleCand.SetPtEtaPhiM(genTauElePt->at(iGenTauEle), genTauEleEta->at(iGenTauEle), genTauElePhi->at(iGenTauEle), genTauEleMass->at(iGenTauEle));
		      if (Ele.DeltaR(GenTauEleCand) <= smallestDR)
			{
			  smallestDR = Ele.DeltaR(GenTauEleCand);
			  findMatchedRecGenTauEle = true;
			  GenTauEle = GenTauEleCand;
			} // end if Ele.DeltaR(GenTauEleCand) <= smallestDR
		    } // end for loop on GenTauEle   
		} // end if genTauElePt->size() > 0


	      if(findMatchedRecGenMu1 && findMatchedRecGenMu2 && findMatchedRecGenTauEle && findMatchedRecGenTauMu){
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

		genmu3Pt->Fill(GenMu3.Pt(), weight);
                genmu3Eta->Fill(GenMu3.Eta(), weight);
                genmu3Phi->Fill(GenMu3.Phi(), weight);
                genmu3Mass->Fill(GenMu3.M(), weight);
                mu3PtVSGenMu3Pt->Fill(Mu3.Pt(), GenMu3.Pt(), weight);
                mu3EtaVSGenMu3Eta->Fill(Mu3.Eta(), GenMu3.Eta(), weight);
                mu3PhiVSGenMu3Phi->Fill(Mu3.Phi(), GenMu3.Phi(), weight);

                gentauMu1Pt->Fill(GenTauMu.Pt(), weight);
                gentauMu1Eta->Fill(GenTauMu.Eta(), weight);
                gentauMu1Phi->Fill(GenTauMu.Phi(), weight);
                gentauMu1Mass->Fill(GenTauMu.M(), weight);

		dRgenMu1genMu2->Fill(GenMu1.DeltaR(GenMu2), weight);
		dRgenMu1genMu3->Fill(GenMu1.DeltaR(GenMu3), weight);
		dRgenMu2genMu3->Fill(GenMu2.DeltaR(GenMu3), weight);
		dRgenTaugenTau->Fill(GenTauMu.DeltaR(GenTauEle), weight);

		invMassgenMu1genMu2->Fill((GenMu1+GenMu2).M(),weight);
		invMassgenTaugenTau->Fill((GenTauMu+GenTauEle).M(), weight);
                invMassgenMuMuTauMuTauEle->Fill((GenMu1+GenMu2+GenTauMu+GenTauEle).M(), weight);

                ptgenMu1genMu2->Fill((GenMu1+GenMu2).Pt(), weight);
		ptgenTaugenTau->Fill((GenTauMu+GenTauEle).Pt(), weight);
                ptgenMuMuTauMuTauHad->Fill((GenMu1+GenMu2+GenTauMu+GenTauEle).Pt(), weight);

		dRInvMassgenMu1genMu2->Fill(GenMu1.DeltaR(GenMu2), (GenMu1+GenMu2).M(), weight);
		dRInvMassgenTaugenTau->Fill(GenTauMu.DeltaR(GenTauEle), (GenTauMu+GenTauEle).M(), weight);

	      } // end if findMatchedRecGenMu1 && findMatchedRecGenMu2 && findMatchedRecGenTauEle && findMatchedRecGenTauMu 
	    } //end ifMC and match RecGen
      } // end if findMu1 && findMu2 && findMuElePair
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
