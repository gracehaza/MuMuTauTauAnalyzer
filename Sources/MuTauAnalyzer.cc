#define MuTauAnalyzer_cxx
#include "MuTauAnalyzer.h"
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

void MuTauAnalyzer::Loop()
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

      // ---- prepare for the vector of matched muon pairs and muon-tau pairs ---
      vector<TLorentzVector> Mu1s;
      vector<TLorentzVector> Mu2s;
      vector<TLorentzVector> Mu3s;
      vector<TLorentzVector> Taus;

      vector<float> Mu1Iso;
      vector<float> Mu2Iso;
      vector<float> Mu3Iso;
      vector<float> TauIso;
      vector<float> TauDM;
      vector<float> TauMediumIsoDisc;

      Mu1s.clear();
      Mu2s.clear();
      Mu3s.clear();
      Taus.clear();

      Mu1Iso.clear();
      Mu2Iso.clear();
      Mu3Iso.clear();
      TauIso.clear();
      TauDM.clear();
      TauMediumIsoDisc.clear();
      // ========================================================================

      // ---- these vectors containing the rank of each matched muon to avoid double counting ---
      vector<int> indexMu1s;
      vector<int> indexMu2s;
      vector<int> indexMu3s;

      indexMu1s.clear();
      indexMu2s.clear();
      indexMu3s.clear();
      // =============================================================================

      // ---- these vectors containing the muons and taus that are not matched into pairs --- 
      vector<TLorentzVector> unMatchedMus;
      vector<TLorentzVector> unMatchedTaus;

      vector<float> unMatchedMuonIso;
      vector<float> unMatchedTauIso;
      vector<float> unMatchedTauDM;

      unMatchedMus.clear();
      unMatchedTaus.clear();

      unMatchedMuonIso.clear();
      unMatchedTauIso.clear();
      unMatchedTauDM.clear();
      // ============================================================================

      // ---- define varibles that will be used to be pushed into the above vectors ---
      TLorentzVector Mu1;
      TLorentzVector Mu2;
      TLorentzVector Mu3;
      TLorentzVector Tau;
      TLorentzVector unMatchedMu;
      // ============================================================================

      // ---- start loop on muon candidates ----
      for (unsigned int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
      {
          if (indexMu2s.size() > 0) 
          {
              std::vector<int>::iterator iter = std::find(indexMu2s.begin(), indexMu2s.end(), iMuon);
              if (iter != indexMu2s.end()) continue;
          } // end if there is any matched Mu2 candidiate

          //cout << "******** Mu1 index: " << iMuon << endl;
          if (recoMuonIsolation->at(iMuon) > 0.25) continue;
          Mu1.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
          float smallestDR = 1.0; // dR cut between Mu1 and Mu2
          bool findMu2 = false;
          int indexMu2 = 0;

          for (unsigned int iMuon2=iMuon+1; iMuon2<recoMuonPt->size(); iMuon2++)
          {
              std::vector<int>::iterator iter2 = std::find(indexMu2s.begin(), indexMu2s.end(), iMuon2);
              if (iter2 != indexMu2s.end()) continue;

              //cout << "******** Mu2 index: " << iMuon2 << endl;
              if (recoMuonIsolation->at(iMuon2) > 0.25) continue;
              TLorentzVector Mu2Cand; // prepare this variable for dR(Mu1,Mu2) implementation
              Mu2Cand.SetPtEtaPhiE(recoMuonPt->at(iMuon2), recoMuonEta->at(iMuon2), recoMuonPhi->at(iMuon2), recoMuonEnergy->at(iMuon2));
              if((Mu1.DeltaR(Mu2Cand) < smallestDR) && (recoMuonPDGId->at(iMuon) == (-1) * recoMuonPDGId->at(iMuon2)))
              {
                  Mu2.SetPtEtaPhiE(recoMuonPt->at(iMuon2), recoMuonEta->at(iMuon2), recoMuonPhi->at(iMuon2), recoMuonEnergy->at(iMuon2));
                  smallestDR = Mu1.DeltaR(Mu2);
                  findMu2 = true;
                  indexMu2 = iMuon2;
              } // end if pair candidates
          } // end loop for mu2
          
          if (findMu2 == true)
          {
              Mu1s.push_back(Mu1);
              Mu2s.push_back(Mu2);

              indexMu1s.push_back(iMuon);
              indexMu2s.push_back(indexMu2);

              Mu1Iso.push_back(recoMuonIsolation->at(iMuon));
              Mu2Iso.push_back(recoMuonIsolation->at(indexMu2));
          } // end if findMu2 

      } // end loop for mu1

      // ------- start loop on tau candidates -------
      for (unsigned int iTau=0; iTau<recoTauPt->size(); iTau++)
      {
          if (recoTauIsoMVAMedium->at(iTau) <= 0) continue;
          Tau.SetPtEtaPhiE(recoTauPt->at(iTau), recoTauEta->at(iTau), recoTauPhi->at(iTau), recoTauEnergy->at(iTau));
          float smallestDR = 0.8; // dR cut between Mu3 and tau
          bool findMu3 = false;
          int indexMu3 = 0;

          for (unsigned int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
          {
              std::vector<int>::iterator iter1 = std::find(indexMu1s.begin(), indexMu1s.end(), iMuon);
              std::vector<int>::iterator iter2 = std::find(indexMu2s.begin(), indexMu2s.end(), iMuon);
              if (iter1 != indexMu1s.end() || iter2 != indexMu2s.end()) continue;

              //cout << "******** Mu3 index: " << iMuon << endl;
              TLorentzVector Mu3Cand; // prepare this variable for dR(Mu3, tau) implementation
              Mu3Cand.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
              if ((Tau.DeltaR(Mu3Cand) < smallestDR) && (recoTauPDGId->at(iTau)/fabs(recoTauPDGId->at(iTau)) == (-1) * recoMuonPDGId->at(iMuon)/fabs(recoMuonPDGId->at(iMuon))))
              {
                  Mu3.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
                  smallestDR = Tau.DeltaR(Mu3);
                  findMu3 = true;
                  indexMu3 = iMuon;
              } // end if find mu3 with tau matched
          } // end loop for mu3

          if (findMu3 == true)
          {
              Mu3s.push_back(Mu3);
              Taus.push_back(Tau);

              indexMu3s.push_back(indexMu3);

              Mu3Iso.push_back(recoMuonIsolation->at(indexMu3));
              TauIso.push_back(recoTauIsoMVArawValue->at(iTau));
              TauMediumIsoDisc.push_back(recoTauIsoMVAMedium->at(iTau));
              TauDM.push_back(recoTauDecayMode->at(iTau));
          } // end if findMu3

          else{
              unMatchedTaus.push_back(Tau);
              unMatchedTauIso.push_back(recoTauIsoMVArawValue->at(iTau));
              unMatchedTauDM.push_back(recoTauDecayMode->at(iTau));
          } // end else findMu3
      } // end loop for tau

      // ---- search for unMatched muon candidates ----
      for (unsigned int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
      {
          std::vector<int>::iterator iter1 = std::find(indexMu1s.begin(), indexMu1s.end(), iMuon);
          std::vector<int>::iterator iter2 = std::find(indexMu2s.begin(), indexMu2s.end(), iMuon);
          std::vector<int>::iterator iter3 = std::find(indexMu3s.begin(), indexMu3s.end(), iMuon);

          if (iter1 == indexMu1s.end() && iter2 == indexMu2s.end() && iter3 == indexMu3s.end())
          {
              unMatchedMu.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
              unMatchedMus.push_back(unMatchedMu);
              unMatchedMuonIso.push_back(recoMuonIsolation->at(iMuon));
          } // end if find unMatched Mu
      } // end loop for unMatched muon candidates

      // ---- prepare event weight info ----
      double weight = 1;
      if (isMC == true)
      {
          weight *= genEventWeight; 
      } // end if isMC == true

      // ---- fill histograms ----
      nMatchedMuPair->Fill(Mu1s.size(), weight);
      nMatchedMuTauPair->Fill(Taus.size(), weight);
      nUnMatchedMu->Fill(unMatchedMus.size(), weight);
      nUnMatchedTau->Fill(unMatchedTaus.size(), weight);
      nMatchedMuPairNMatchedMuTauPair->Fill(Mu1s.size(), Taus.size(), weight);

      if (Mu1s.size() >0 && Taus.size() >0)
      {
          // --- filling histograms of four-body of mu-mu-mu-tau ---
          for (unsigned int iMuon=0; iMuon<Mu1s.size(); iMuon++)
          {
              Mu1 = Mu1s.at(iMuon);
              Mu2 = Mu2s.at(iMuon);
              TLorentzVector Mu1Mu2 = Mu1 + Mu2;
              bool passDR = false; // dR between mu-mu pair and mu-tau pair

              for (unsigned int iTau=0; iTau<Taus.size(); iTau++)
              {
                  Mu3 = Mu3s.at(iTau);
                  Tau = Taus.at(iTau);
                  TLorentzVector MuTau = Mu3 + Tau;

                  if (Mu1.DeltaR(Mu3) > 0.4 && Mu2.DeltaR(Mu3) > 0.4 && Mu1.DeltaR(Tau) > 0.8 && Mu2.DeltaR(Tau) > 0.8)
                  {
                      passDR = true;

                      dRMuTau->Fill(Mu3.DeltaR(Tau), weight);
                      invMassMuTau->Fill(MuTau.M(), weight);
                      dRInvMassMuTau->Fill(Mu3.DeltaR(Tau), MuTau.M(), weight);

                      Mu3IsoMuTauPair->Fill(Mu3Iso.at(iTau), weight);
                      TauIsoMVAMuTauPair->Fill(TauIso.at(iTau), weight);
                      tauDecayMode->Fill(TauDM.at(iTau), weight);

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
                      break;
                  } // end if dR between mu-mu pair and mu-tau pair
              } // end loop for mu-tau pairs
              
              if (passDR == true)
              {
                  dRMuMu->Fill(Mu1.DeltaR(Mu2), weight);
                  invMassMuMu->Fill(Mu1Mu2.M(), weight);
                  dRInvMassMuMu->Fill(Mu1.DeltaR(Mu2), Mu1Mu2.M(), weight);

                  Mu1IsoMuMuPair->Fill(Mu1Iso.at(iMuon), weight);
                  Mu2IsoMuMuPair->Fill(Mu2Iso.at(iMuon), weight);

                  mu1Pt->Fill(Mu1.Pt(), weight);
                  mu1Eta->Fill(Mu1.Eta(), weight);
                  mu1Phi->Fill(Mu1.Phi(), weight);

                  mu2Pt->Fill(Mu2.Pt(), weight);
                  mu2Eta->Fill(Mu2.Eta(), weight);
                  mu2Phi->Fill(Mu2.Phi(), weight);
                  break;
              } // end if passDR between mu-mu pair and mu-tau pair
          } // end loop for mu-mu pairs
      } // end if mu-mu & mu-tau pairs

      for (unsigned int iMuon=0; iMuon<unMatchedMus.size(); iMuon++)
      {
          unMatchedMuIso->Fill(unMatchedMuonIso.at(iMuon), weight);
      } // end loop for unMatched muons

      for (unsigned int iTau=0; iTau<unMatchedTaus.size(); iTau++)
      {
          unMatchedTauIsoMVA->Fill(unMatchedTauIso.at(iTau), weight);
          unMatchedTauDecayMode->Fill(unMatchedTauDM.at(iTau), weight);
      } // end loop for unMatched taus

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

   outputFile->Write();
   outputFile->Close();
}
