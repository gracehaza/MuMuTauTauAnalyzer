#define MuMuAnalyzer_cxx
#include "MuMuAnalyzer.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <sys/time.h>
#include <iomanip>
#include <TLorentzVector.h>
#include <math.h>
using namespace std;

void MuMuAnalyzer::Loop()
{
   TString outputfileName = createOutputFileName();
   TFile* outputFile = new TFile(outputfileName, "RECREATE");
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   struct timeval t0;
   gettimeofday(&t0, 0);
   int mess_every_n =  std::min(1000000LL, nentries/10);

   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      if (jentry % mess_every_n == 0 && jentry > 0){
          timeval t1;
          gettimeofday(&t1, 0);
          double dt = (t1.tv_sec - t0.tv_sec) + 1.e-6*(t1.tv_usec - t0.tv_usec);
          dt *= double(nentries  - jentry) / mess_every_n;
          int dt_s = int(dt + 0.5);
          int dt_h = int(dt_s) / 3600; 
          dt_s -= dt_h * 3600;
          int dt_m = dt_s / 60;
          dt_s -= dt_m *60;
          cout << TString::Format("%4.1f%%", (100. * jentry) / nentries)
              << "\t" << std::setw(11) << jentry << " / " << nentries
              << "\t " << std::setw(4) << int(dt / mess_every_n * 1.e6 + 0.5) << " us/event"
              << "\t Remaining time for this dataset loop: " 
              << std::setw(2) << dt_h << " h "
              << std::setw(2) << dt_m << " min "
              << std::setw(2) << dt_s << " s"
              << "\r" << std::flush;
          t0 = t1;
      } // end if for timer settings

      //cout << "*** iEntry: " << jentry << endl;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

      // ---- prepare for the vector of matched muon pairs ---
      vector<TLorentzVector> Mu1s;
      vector<TLorentzVector> Mu2s;

      vector<float> Mu1Iso;
      vector<float> Mu2Iso;

      Mu1s.clear();
      Mu2s.clear();
      Mu1Iso.clear();
      Mu2Iso.clear();

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

      // ---- start loop on muon candidates ----
      for (int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
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

          for (int iMuon2=iMuon+1; iMuon2<recoMuonPt->size(); iMuon2++)
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

              Mu1Iso.push_back(recoMuonIsolation->at(iMuon));
              Mu2Iso.push_back(recoMuonIsolation->at(indexMu2));
          } // end if findMu2 

      } // end loop for mu1

      // ---- search for unMatched muon candidates ----
      for (int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
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

      // ---- fill histograms ----
      nMatchedMuPair->Fill(Mu1s.size());
      nUnMatchedMu->Fill(unMatchedMus.size());
      
      if (Mu1s.size() >0)
      {
          // --- filling histograms of mu-mu ---
          for (int iMuon=0; iMuon<Mu1s.size(); iMuon++)
          {
              Mu1 = Mu1s.at(iMuon);
              Mu2 = Mu2s.at(iMuon);
              TLorentzVector Mu1Mu2 = Mu1 + Mu2;

              dRMuMu->Fill(Mu1.DeltaR(Mu2));
              invMassMuMu->Fill(Mu1Mu2.M());
              dRInvMassMuMu->Fill(Mu1.DeltaR(Mu2), Mu1Mu2.M());

              Mu1IsoMuMuPair->Fill(Mu1Iso.at(iMuon));
              Mu2IsoMuMuPair->Fill(Mu2Iso.at(iMuon));

              mu1Pt->Fill(Mu1.Pt());
              mu1Eta->Fill(Mu1.Eta());
              mu1Phi->Fill(Mu1.Phi());

              mu2Pt->Fill(Mu2.Pt());
              mu2Eta->Fill(Mu2.Eta());
              mu2Phi->Fill(Mu2.Phi());
          } // end loop for mu-mu pairs
      } // end if mu-mu pairs

      for (int iMuon=0; iMuon<unMatchedMus.size(); iMuon++)
      {
          unMatchedMuPt->Fill(unMatchedMus.at(iMuon).Pt());
          unMatchedMuEta->Fill(unMatchedMus.at(iMuon).Eta());
          unMatchedMuPhi->Fill(unMatchedMus.at(iMuon).Phi());
          unMatchedMuIso->Fill(unMatchedMuonIso.at(iMuon));
      } // end loop for unMatched muons

   }// end loop for events

   outputFile->cd();

   int numberofhist = histColl.size();
   for(int i=0; i<numberofhist; i++){
       if (!isData) histColl[i]->Scale(lumiScale/summedWeights);
       histColl[i]->Write();
   } // end loop for writing all the histograms into the output file

   for(int j=0; j<numberofhist; j++){
       delete histColl[j];
   } // end loop for deleting all the histograms

   outputFile->Write();
   outputFile->Close();
}