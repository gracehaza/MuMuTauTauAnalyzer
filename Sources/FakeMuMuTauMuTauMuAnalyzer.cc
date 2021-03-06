#define FakeMuMuTauMuTauMuAnalyzer_cxx
#include "FakeMuMuTauMuTauMuAnalyzer.h"
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

void FakeMuMuTauMuTauMuAnalyzer::Loop()
{
   TString outputfileName = createOutputFileName();
   TFile* outputFile = new TFile(outputfileName, "RECREATE");
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   if (nMaxEvents >= 0 && nMaxEvents  < nentries) nentries = nMaxEvents;
   cout << "We will run on " << nentries << " events" << endl;

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
      TLorentzVector Mu3;
      TLorentzVector Mu4;

      float Mu1Iso;
      float Mu2Iso;
      float Mu3Iso;
      float Mu4Iso;

      unsigned int indexMu1 = -1;
      unsigned int indexMu2 = -1;
      unsigned int indexMu4 = -1;
      // =============================================================================

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
      float dRCut = 0.3; // dR cut between Mu1 and Mu2
      float highestPt = 0;
      float invMassLowThre = 60.0;
      float invMassHighThre = 120.0;
      bool findMu2 = false;

      // ---- start loop on muon candidates for mu2 ----
      for (unsigned int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
      {
          if (iMuon == indexMu1) continue;
          if (recoMuonIsolation->at(iMuon) > 0.25) continue;

          TLorentzVector Mu2Cand; // prepare this variable for dR(Mu1,Mu2) implementation
          Mu2Cand.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
          if((Mu1.DeltaR(Mu2Cand) > dRCut) && (Mu2Cand.Pt() > highestPt) 
                  && ((Mu1+Mu2Cand).M() > invMassLowThre) && ((Mu1+Mu2Cand).M() < invMassHighThre)
                  && (recoMuonPDGId->at(indexMu1) == (-1) * recoMuonPDGId->at(iMuon)))
          {
              Mu2.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
              Mu2Iso = recoMuonIsolation->at(iMuon);
              highestPt = Mu2Cand.Pt();
              findMu2 = true;
          } // end if pair candidates
      } // end loop for mu2

      if (!findMu2) continue;
      bool findMuMuPair = false;

      // ---- search for an additional muon pair for fake rate study ----
      for (unsigned int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
      {
          if (iMuon == indexMu1 || iMuon == indexMu2) continue;
          if ((invertedMu4Iso == false && recoMuonIsolation->at(iMuon) > Mu4IsoThreshold) || (invertedMu4Iso == true && recoMuonIsolation->at(iMuon) < Mu4IsoThreshold)) continue;
          
          TLorentzVector Mu4Cand;
          Mu4Cand.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));

          if (Mu4Cand.DeltaR(Mu1) < 0.4 || Mu4Cand.DeltaR(Mu2) < 0.4) continue;
          Mu4.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
          Mu4Iso = recoMuonIsolation->at(iMuon);
          indexMu4 = iMuon;

          float smallestDR = 4.0; // dR cut between Mu3 and Mu4
          bool findMu3 = false;

          for (unsigned int iMuon3=0; iMuon3<recoMuonPt->size(); iMuon3++)
          {
              if (iMuon3 == indexMu1 || iMuon3 == indexMu2 || iMuon3 == indexMu4) continue;

              TLorentzVector Mu3Cand; // prepare this variable for dR(Mu3, Mu4) implementation
              Mu3Cand.SetPtEtaPhiE(recoMuonPt->at(iMuon3), recoMuonEta->at(iMuon3), recoMuonPhi->at(iMuon3), recoMuonEnergy->at(iMuon3));
              if ((Mu4.DeltaR(Mu3Cand) < smallestDR) && (recoMuonPDGId->at(iMuon3) == (-1) * recoMuonPDGId->at(iMuon)) && ((Mu4+Mu3Cand).M() < 60.0) && (Mu3Cand.DeltaR(Mu1) > 0.4) && (Mu3Cand.DeltaR(Mu2) > 0.4))
              {
                  Mu3.SetPtEtaPhiE(recoMuonPt->at(iMuon3), recoMuonEta->at(iMuon3), recoMuonPhi->at(iMuon3), recoMuonEnergy->at(iMuon3));
                  Mu3Iso = recoMuonIsolation->at(iMuon3);
                  smallestDR = Mu4.DeltaR(Mu3);
                  findMu3 = true;
              } // end if find mu3 with mu4 matched
          } // end loop for mu3

          if (!findMu3) continue;
          else{
              findMuMuPair = true;
              break;
          } // end if findMu3
      } // end loop for Mu4

      // ---- prepare event weight info ----
      double weight = 1;
      if (isMC == true)
      {
          weight *= genEventWeight; 
      } // end if isMC == true

      // ---- fill histograms ----
      if (findMu1 && findMu2 && findMuMuPair)
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

          ptMu3Mu4->Fill((Mu3+Mu4).Pt(), weight);
          dRMu3Mu4->Fill(Mu3.DeltaR(Mu4), weight);
          invMassMu3Mu4->Fill((Mu3+Mu4).M(), weight);
          dRInvMassMu3Mu4->Fill(Mu3.DeltaR(Mu4), (Mu3+Mu4).M(), weight);

          mu3Iso->Fill(Mu3Iso, weight);
          mu4Iso->Fill(Mu4Iso, weight);

          mu3Pt->Fill(Mu3.Pt(), weight);
          mu3Eta->Fill(Mu3.Eta(), weight);
          mu3Phi->Fill(Mu3.Phi(), weight);

          mu4Pt->Fill(Mu4.Pt(), weight);
          mu4Eta->Fill(Mu4.Eta(), weight);
          mu4Phi->Fill(Mu4.Phi(), weight);

          dRMu1Mu3->Fill(Mu1.DeltaR(Mu3), weight);
          dRMu1Mu4->Fill(Mu1.DeltaR(Mu4), weight);
          dRMu2Mu3->Fill(Mu2.DeltaR(Mu3), weight);
          dRMu2Mu4->Fill(Mu2.DeltaR(Mu4), weight);

          ptMuMuTauMuTauMu->Fill((Mu1+Mu2+Mu3+Mu4).Pt(), weight);
          invMassMuMuTauMuTauMu->Fill((Mu1+Mu2+Mu3+Mu4).M(), weight);
      } // end if findMu1 && findMu2 && findMuMuPair
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
