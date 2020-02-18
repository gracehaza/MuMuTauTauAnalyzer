#define MuMuTauMuTauMuAnalyzer_cxx
#include "MuMuTauMuTauMuAnalyzer.h"
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

void MuMuTauMuTauMuAnalyzer::Loop()
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

      // ---- prepare for the vector of matched muon pairs and muon-muon pairs ---
      vector<TLorentzVector> Mu1s;
      vector<TLorentzVector> Mu2s;
      vector<TLorentzVector> Mu3s;
      vector<TLorentzVector> Mu4s;
      vector<TLorentzVector> genMu1s;
      vector<TLorentzVector> genMu2s;
      vector<TLorentzVector> genMu3s;
      vector<TLorentzVector> genMu4s;

      vector<float> Mu1Iso;
      vector<float> Mu2Iso;
      vector<float> Mu3Iso;
      vector<float> Mu4Iso;

      Mu1s.clear();
      Mu2s.clear();
      Mu3s.clear();
      Mu4s.clear();

      genMu1s.clear();
      genMu2s.clear();
      genMu3s.clear();
      genMu4s.clear();

      Mu1Iso.clear();
      Mu2Iso.clear();
      Mu3Iso.clear();
      Mu4Iso.clear();
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
      TLorentzVector Mu3;
      TLorentzVector Mu4;
      TLorentzVector unMatchedMu;

      TLorentzVector genMu1;
      TLorentzVector genMu2;
      TLorentzVector genMu3;
      TLorentzVector genMu4;

      // ============================================================================

      for (unsigned int igenMuon=0; igenMuon<genMuonPt->size(); igenMuon++){
        genMu1.SetPtEtaPhiE(genMuonPt->at(igenMuon), genMuonEta->at(igenMuon), genMuonPhi->at(igenMuon), genMuonMass->at(igenMuon));
        genMu1s.push_back(genMu1);
      }

      /*      for (unsigned int igenMuon2=0; igenMuon2<genMuon2Pt->size(); igenMuon2++){
        genMu2.SetPtEtaPhiE(genMuon2Pt->at(igenMuon2), genMuon2Eta->at(igenMuon2), genMuon2Phi->at(igenMuon2), genMuon2Energy->at(igenMuon2));
        genMu2s.push_back(genMu2);
      }
      for (unsigned int igenMuon3=0; igenMuon3<genMuon3Pt->size(); igenMuon3++){
        genMu3.SetPtEtaPhiE(genMuon3Pt->at(igenMuon3), genMuon3Eta->at(igenMuon3), genMuon3Phi->at(igenMuon3), genMuon3Energy->at(igenMuon3));
        genMu3s.push_back(genMu3);
      }

      for (unsigned int igenMuon4=0; igenMuon4<genMuon4Pt->size(); igenMuon4++){
        genMu4.SetPtEtaPhiE(genMuon4Pt->at(igenMuon4), genMuon4Eta->at(igenMuon4), genMuon4Phi->at(igenMuon4), genMuon4Energy->at(igenMuon4));
        genMu4s.push_back(genMu4);
      }
      */
      // ---- start loop on muon candidates ----
      for (unsigned int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
      {
          if (indexMu2s.size() > 0) 
          {
              std::vector<int>::iterator iter = std::find(indexMu2s.begin(), indexMu2s.end(), iMuon);
              if (iter != indexMu2s.end()) continue;
          } // end if there is any matched Mu2 candidiate

          if (recoMuonIsolation->at(iMuon) > 0.25) continue;
          Mu1.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
          float smallestDR = 1.0; // dR cut between Mu1 and Mu2
          bool findMu2 = false;
          int indexMu2 = 0;

          for (unsigned int iMuon2=iMuon+1; iMuon2<recoMuonPt->size(); iMuon2++)
          {
              std::vector<int>::iterator iter2 = std::find(indexMu2s.begin(), indexMu2s.end(), iMuon2);
              if (iter2 != indexMu2s.end()) continue;

              if ((invertedMu2Iso == false && recoMuonIsolation->at(iMuon2) > Mu2IsoThreshold) || (invertedMu2Iso == true && recoMuonIsolation->at(iMuon2) < Mu2IsoThreshold)) continue;
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
              break; // only select one pair of Mu1Mu2 and give the others to Mu3Mu4 candidates, otherwise some events may be lost (multiple Mu1Mu2 pair but no Mu3Mu4 pair due to Mu1Mu2 vector takes potential Mu3Mu4 pair!)
          } // end if findMu2 
      } // end loop for mu1

      // ------- start loop on the second muon pair candidates -------
      for (unsigned int iMuon4=0; iMuon4<recoMuonPt->size(); iMuon4++)
      {
          std::vector<int>::iterator iter1 = std::find(indexMu1s.begin(), indexMu1s.end(), iMuon4);
          std::vector<int>::iterator iter2 = std::find(indexMu2s.begin(), indexMu2s.end(), iMuon4);
          if (iter1 != indexMu1s.end() || iter2 != indexMu2s.end()) continue;

          if ((invertedMu4Iso == false && recoMuonIsolation->at(iMuon4) > Mu4IsoThreshold) || (invertedMu4Iso == true && recoMuonIsolation->at(iMuon4) < Mu4IsoThreshold)) continue;
          Mu4.SetPtEtaPhiE(recoMuonPt->at(iMuon4), recoMuonEta->at(iMuon4), recoMuonPhi->at(iMuon4), recoMuonEnergy->at(iMuon4));
          float smallestDR = 0.8; // dR cut between Mu3 and Mu4
          bool findMu3 = false;
          int indexMu3 = 0;

          for (unsigned int iMuon3=iMuon4+1; iMuon3<recoMuonPt->size(); iMuon3++)
          {
              std::vector<int>::iterator iter3 = std::find(indexMu1s.begin(), indexMu1s.end(), iMuon3);
              std::vector<int>::iterator iter4 = std::find(indexMu2s.begin(), indexMu2s.end(), iMuon3);
              if (iter3 != indexMu1s.end() || iter4 != indexMu2s.end()) continue;

              TLorentzVector Mu3Cand; // prepare this variable for dR(Mu3, Mu4) implementation
              Mu3Cand.SetPtEtaPhiE(recoMuonPt->at(iMuon3), recoMuonEta->at(iMuon3), recoMuonPhi->at(iMuon3), recoMuonEnergy->at(iMuon3));
              if ((Mu4.DeltaR(Mu3Cand) < smallestDR) && (recoMuonPDGId->at(iMuon3) == (-1) * recoMuonPDGId->at(iMuon4)))
              {
                  Mu3.SetPtEtaPhiE(recoMuonPt->at(iMuon3), recoMuonEta->at(iMuon3), recoMuonPhi->at(iMuon3), recoMuonEnergy->at(iMuon3));
                  smallestDR = Mu4.DeltaR(Mu3);
                  findMu3 = true;
                  indexMu3 = iMuon3;
              } // end if find mu3 with mu4 matched
          } // end loop for mu3

          if (findMu3 == true)
          {
              Mu3s.push_back(Mu3);
              Mu4s.push_back(Mu4);

              Mu3Iso.push_back(recoMuonIsolation->at(indexMu3));
              Mu4Iso.push_back(recoMuonIsolation->at(iMuon4));
          } // end if findMu3

          else{
              unMatchedMus.push_back(Mu4);
              unMatchedMuonIso.push_back(recoMuonIsolation->at(iMuon4));
          } // end else findMu3
      } // end loop for Mu4

      // ---- prepare event weight info ----
      double weight = 1;
      if (isMC == true)
      {
          weight *= genEventWeight; 
      } // end if isMC == true

      if(genMu1s.size() >0){
        for (unsigned int igenMuon=0; igenMuon<genMu1s.size(); igenMuon++)
	   {
            genMu1 = genMu1s.at(igenMuon);
            genmu1Pt->Fill(genMu1.Pt(), weight);
	    /* if(genMu2s.size() >0){
              for (unsigned int igenMuon2=0; igenMuon2<genMu2s.size(); igenMuon2++)
                {
                  TLorentzVector genMu1genMu2 = genMu1 + genMu2;
                  genMu2 = genMu2s.at(igenMuon2);
                  genmu2Pt->Fill(genMu2.Pt(), weight);
                  dRgenMu1genMu2->Fill(genMu1.DeltaR(genMu2), weight);
                  dRInvMassgenMu1genMu2->Fill(genMu1.DeltaR(genMu2), genMu1genMu2.M(), weight);
                  if(genMu3s.size() >0){
                    for (unsigned int igenMuon3=0; igenMuon3<genMu3s.size(); igenMuon3++)
                      {
                        genMu3 = genMu3s.at(igenMuon3);
                        genele1Pt->Fill(genMu3.Pt(), weight);
                        if(genMu4s.size() >0){
                          for (unsigned int igenMuon4=0; igenMuon4<genMu4s.size(); igenMuon4++)
                            {
                              TLorentzVector genMuon3genMuon4 = genMu3 + genMu4;
                              genMu4 = genMu4s.at(igenMuon4);
                              genele2Pt->Fill(genMu4.Pt(), weight);
                              dRgenMu3genMu4->Fill(genMu3.DeltaR(genMu4), weight);
                              dRInvMassgenMu3genMu4->Fill(genMu3.DeltaR(genMu4), genMuon3genMuon4.M(), weight);
                            }
                        }
                      }
                  }

                }
            } */
	   }
      }


      // ---- fill histograms ----
      nMatchedMuPair->Fill(Mu1s.size(), weight);
      nMatchedMuMuPair->Fill(Mu4s.size(), weight);
      nUnMatchedMu->Fill(unMatchedMus.size(), weight);
      nMatchedMuPairNMatchedMuMuPair->Fill(Mu1s.size(), Mu4s.size(), weight);

      if (Mu1s.size() >0 && Mu4s.size() >0)
      {
          // --- filling histograms of four-body of mu-mu-mu-mu ---
          for (unsigned int iMuon=0; iMuon<Mu1s.size(); iMuon++)
          {
              Mu1 = Mu1s.at(iMuon);
              Mu2 = Mu2s.at(iMuon);
              TLorentzVector Mu1Mu2 = Mu1 + Mu2;
              bool passDR = false; // dR between mu-mu pair and mu-mu pair

              for (unsigned int iMuon4=0; iMuon4<Mu4s.size(); iMuon4++)
              {
                  Mu3 = Mu3s.at(iMuon4);
                  Mu4 = Mu4s.at(iMuon4);
                  TLorentzVector Mu3Mu4 = Mu3 + Mu4;
                  TLorentzVector MuMuMuMu = Mu1Mu2 + Mu3Mu4;

                  if (Mu1.DeltaR(Mu3) > 0.4 && Mu2.DeltaR(Mu3) > 0.4 && Mu1.DeltaR(Mu4) > 0.8 && Mu2.DeltaR(Mu4) > 0.8)
                  {
                      passDR = true;

                      ptMu3Mu4->Fill(Mu3Mu4.Pt(), weight);
                      dRMu3Mu4->Fill(Mu3.DeltaR(Mu4), weight);
                      invMassMu3Mu4->Fill(Mu3Mu4.M(), weight);
                      dRInvMassMu3Mu4->Fill(Mu3.DeltaR(Mu4), Mu3Mu4.M(), weight);

                      mu3Iso->Fill(Mu3Iso.at(iMuon4), weight);
                      mu4Iso->Fill(Mu4Iso.at(iMuon4), weight);

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

                      ptMuMuTauMuTauMu->Fill(MuMuMuMu.Pt(), weight);
                      invMassMuMuTauMuTauMu->Fill(MuMuMuMu.M(), weight);
                      break;
                  } // end if dR between mu-mu pair and mu3-mu4 pair
              } // end loop for mu3-mu4 pairs
              
              if (passDR == true)
              {
                  ptMu1Mu2->Fill(Mu1Mu2.Pt(), weight);
                  dRMu1Mu2->Fill(Mu1.DeltaR(Mu2), weight);
                  invMassMu1Mu2->Fill(Mu1Mu2.M(), weight);
                  dRInvMassMu1Mu2->Fill(Mu1.DeltaR(Mu2), Mu1Mu2.M(), weight);

                  mu1Iso->Fill(Mu1Iso.at(iMuon), weight);
                  mu2Iso->Fill(Mu2Iso.at(iMuon), weight);

                  mu1Pt->Fill(Mu1.Pt(), weight);
                  mu1Eta->Fill(Mu1.Eta(), weight);
                  mu1Phi->Fill(Mu1.Phi(), weight);

                  mu2Pt->Fill(Mu2.Pt(), weight);
                  mu2Eta->Fill(Mu2.Eta(), weight);
                  mu2Phi->Fill(Mu2.Phi(), weight);
                  break;
              } // end if passDR between mu-mu pair and mu3-mu4 pair
          } // end loop for mu-mu pairs
      } // end if mu-mu & mu3-mu4 pairs

      for (unsigned int iMuon=0; iMuon<unMatchedMus.size(); iMuon++)
      {
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
