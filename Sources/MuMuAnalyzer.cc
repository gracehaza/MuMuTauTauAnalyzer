#define MuMuAnalyzer_cxx
#include "MuMuAnalyzer.h"
#include "RoccoR.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TRandom.h>
#include <iomanip>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>
using namespace std;

void MuMuAnalyzer::Loop()
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

      // ---- prepare for the vector of matched muon pairs ---
      vector<TLorentzVector> Mu1s;
      vector<TLorentzVector> Mu2s;

      vector<int> Mu1Charge;
      vector<int> Mu2Charge;

      vector<int> Mu1NTrackerLayers;
      vector<int> Mu2NTrackerLayers;

      vector<float> Mu1Iso;
      vector<float> Mu2Iso;

      Mu1s.clear();
      Mu2s.clear();
      Mu1Iso.clear();
      Mu2Iso.clear();
      Mu1Charge.clear();
      Mu2Charge.clear();
      Mu1NTrackerLayers.clear();
      Mu2NTrackerLayers.clear();
      // ========================================================================

      // ---- these vectors containing the rank of each matched muon to avoid double counting ---
      vector<int> indexMu1s;
      vector<int> indexMu2s;

      indexMu1s.clear();
      indexMu2s.clear();
      // =============================================================================

      // ---- these vectors containing the muons that are not matched into pairs --- 
      vector<TLorentzVector> unMatchedMus;
      vector<float> unMatchedMuonIso;

      unMatchedMus.clear();
      unMatchedMuonIso.clear();
      // ============================================================================

      // ---- define varibles that will be used to be pushed into the above vectors ---
      TLorentzVector Mu1;
      TLorentzVector Mu2;
      TLorentzVector unMatchedMu;
      // ============================================================================

      // ---- define the variable for rochester correction ----
      RoccoR rc("text/RoccoR2017.txt");
      double dataRochesterSF = 1;
      double mcRochesterSF = 1;

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

          float dRCut = 0.3; // dR cut between Mu1 and Mu2
          float highestPt = 0;
          float invMassLowThre = 60.0;
          float invMassHighThre = 120.0;
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
              if((Mu1.DeltaR(Mu2Cand) > dRCut) && (Mu2Cand.Pt() > highestPt) 
                      && ((Mu1+Mu2Cand).M() > invMassLowThre) && ((Mu1+Mu2Cand).M() < invMassHighThre)
                      && (recoMuonPDGId->at(iMuon) == (-1) * recoMuonPDGId->at(iMuon2)))
              {
                  Mu2.SetPtEtaPhiE(recoMuonPt->at(iMuon2), recoMuonEta->at(iMuon2), recoMuonPhi->at(iMuon2), recoMuonEnergy->at(iMuon2));
                  highestPt = Mu2Cand.Pt();
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

              Mu1Charge.push_back(recoMuonPDGId->at(iMuon)/fabs(recoMuonPDGId->at(iMuon)));
              Mu2Charge.push_back(recoMuonPDGId->at(indexMu2)/fabs(recoMuonPDGId->at(indexMu2)));

              Mu1NTrackerLayers.push_back(recoMuonNTrackerLayers->at(iMuon));
              Mu2NTrackerLayers.push_back(recoMuonNTrackerLayers->at(indexMu2));

              Mu1Iso.push_back(recoMuonIsolation->at(iMuon));
              Mu2Iso.push_back(recoMuonIsolation->at(indexMu2));
              break;
          } // end if findMu2 
      } // end loop for mu1

      // ---- search for unMatched muon candidates ----
      for (unsigned int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
      {
          std::vector<int>::iterator iter1 = std::find(indexMu1s.begin(), indexMu1s.end(), iMuon);
          std::vector<int>::iterator iter2 = std::find(indexMu2s.begin(), indexMu2s.end(), iMuon);

          if (iter1 == indexMu1s.end() && iter2 == indexMu2s.end())
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
      nUnMatchedMu->Fill(unMatchedMus.size(), weight);
      
      if (Mu1s.size() >0)
      {
          // --- filling histograms of mu-mu ---
          for (unsigned int iMuon=0; iMuon<Mu1s.size(); iMuon++)
          {
              // ---- apply rochester correction ----
              if (isMC == false)
              {
                  dataRochesterSF = rc.kScaleDT(Mu1Charge.at(iMuon), Mu1s.at(iMuon).Pt(), Mu1s.at(iMuon).Eta(), Mu1s.at(iMuon).Phi(), 0, 0);
                  Mu1.SetPtEtaPhiE(Mu1s.at(iMuon).Pt()*dataRochesterSF, Mu1s.at(iMuon).Eta(), Mu1s.at(iMuon).Phi(), Mu1s.at(iMuon).E());

                  dataRochesterSF = rc.kScaleDT(Mu2Charge.at(iMuon), Mu2s.at(iMuon).Pt(), Mu2s.at(iMuon).Eta(), Mu2s.at(iMuon).Phi(), 0, 0);
                  Mu2.SetPtEtaPhiE(Mu2s.at(iMuon).Pt()*dataRochesterSF, Mu2s.at(iMuon).Eta(), Mu2s.at(iMuon).Phi(), Mu2s.at(iMuon).E());
              } // end if isMC == false

              else{
                  double rng = gRandom->Rndm();
                  mcRochesterSF = rc.kSmearMC(Mu1Charge.at(iMuon), Mu1s.at(iMuon).Pt(), Mu1s.at(iMuon).Eta(), Mu1s.at(iMuon).Phi(), Mu1NTrackerLayers.at(iMuon), rng, 0, 0);
                  Mu1.SetPtEtaPhiE(Mu1s.at(iMuon).Pt()*mcRochesterSF, Mu1s.at(iMuon).Eta(), Mu1s.at(iMuon).Phi(), Mu1s.at(iMuon).E());

                  rng = gRandom->Rndm();
                  mcRochesterSF = rc.kSmearMC(Mu2Charge.at(iMuon), Mu2s.at(iMuon).Pt(), Mu2s.at(iMuon).Eta(), Mu2s.at(iMuon).Phi(), Mu2NTrackerLayers.at(iMuon), rng, 0, 0);
                  Mu2.SetPtEtaPhiE(Mu2s.at(iMuon).Pt()*mcRochesterSF, Mu2s.at(iMuon).Eta(), Mu2s.at(iMuon).Phi(), Mu2s.at(iMuon).E());
              } // end isMC == true

              TLorentzVector Mu1Mu2 = Mu1 + Mu2;

              dRMu1Mu2->Fill(Mu1.DeltaR(Mu2), weight);
              invMassMu1Mu2->Fill(Mu1Mu2.M(), weight);
              ptMu1Mu2->Fill(Mu1Mu2.Pt(), weight);
              dRInvMassMu1Mu2->Fill(Mu1.DeltaR(Mu2), Mu1Mu2.M(), weight);

              mu1Iso->Fill(Mu1Iso.at(iMuon), weight);
              mu2Iso->Fill(Mu2Iso.at(iMuon), weight);

              mu1Pt->Fill(Mu1.Pt(), weight);
              mu1Eta->Fill(Mu1.Eta(), weight);
              mu1Phi->Fill(Mu1.Phi(), weight);

              mu2Pt->Fill(Mu2.Pt(), weight);
              mu2Eta->Fill(Mu2.Eta(), weight);
              mu2Phi->Fill(Mu2.Phi(), weight);
          } // end loop for mu-mu pairs
      } // end if mu-mu pairs

      for (unsigned int iMuon=0; iMuon<unMatchedMus.size(); iMuon++)
      {
          unMatchedMuPt->Fill(unMatchedMus.at(iMuon).Pt(), weight);
          unMatchedMuEta->Fill(unMatchedMus.at(iMuon).Eta(), weight);
          unMatchedMuPhi->Fill(unMatchedMus.at(iMuon).Phi(), weight);
          unMatchedMuIso->Fill(unMatchedMuonIso.at(iMuon), weight);
      } // end loop for unMatched muons

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
