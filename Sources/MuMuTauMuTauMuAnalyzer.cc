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

   bool matchRecGen = true;

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

      bool findMuMuPair = false;
      // ------- start loop on the second muon pair candidates -------
      for (unsigned int iMuon=0; iMuon<recoMuonPt->size(); iMuon++)
      {
          if (iMuon == indexMu1 || iMuon == indexMu2) continue;
          if ((invertedMu4Iso == false && recoMuonIsolation->at(iMuon) > Mu4IsoThreshold) || (invertedMu4Iso == true && recoMuonIsolation->at(iMuon) < Mu4IsoThreshold)) continue;
          
          TLorentzVector Mu4Cand;
          Mu4Cand.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));

          if (Mu4Cand.DeltaR(Mu1) < 0.4 || Mu4Cand.DeltaR(Mu2) < 0.4) continue;
          Mu4.SetPtEtaPhiE(recoMuonPt->at(iMuon), recoMuonEta->at(iMuon), recoMuonPhi->at(iMuon), recoMuonEnergy->at(iMuon));
          Mu4Iso = recoMuonIsolation->at(iMuon);

          float smallestDR = 1.0; // dR cut between Mu3 and Mu4
          bool findMu3 = false;

          for (unsigned int iMuon3=0; iMuon3<recoMuonPt->size(); iMuon3++)
          {
              if (iMuon3 == indexMu1 || iMuon3 == indexMu2 || iMuon3 == iMuon) continue;

              TLorentzVector Mu3Cand; // prepare this variable for dR(Mu3, Mu4) implementation
              Mu3Cand.SetPtEtaPhiE(recoMuonPt->at(iMuon3), recoMuonEta->at(iMuon3), recoMuonPhi->at(iMuon3), recoMuonEnergy->at(iMuon3));

              if ((Mu4.DeltaR(Mu3Cand) < smallestDR) && (recoMuonPDGId->at(iMuon3) == (-1) * recoMuonPDGId->at(iMuon)) && ((Mu4+Mu3Cand).M() < 60.0) && (Mu3Cand.DeltaR(Mu1) > 0.4) && (Mu3Cand.DeltaR(Mu2) > 0.4) && (recoMuonIsolation->at(iMuon3) < Mu3IsoThreshold))

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


          // ----- fill flat trees -----                                                                                                                                
          invMassMuMu = (Mu1+Mu2).M();
          visMassTauTau = (Mu3+Mu4).M();
          visMassMuMuTauTau = (Mu1+Mu2+Mu3+Mu4).M();

          deltaRMuMu = Mu1.DeltaR(Mu2);
          deltaRTauTau = Mu3.DeltaR(Mu4);

          Mu1Pt = Mu1.Pt();
          Mu1Eta = Mu1.Eta();

          Mu2Pt = Mu2.Pt();
          Mu2Eta = Mu2.Eta();

          Tau1Pt = Mu3.Pt();
          Tau1Eta = Mu3.Eta();
          Tau1Isolation = Mu3Iso;

          Tau2Pt = Mu4.Pt();
          Tau2Eta = Mu4.Eta();
          Tau2Isolation = Mu4Iso;

          eventWeight = weight/summedWeights;
          TreeMuMuTauTau->Fill();


	  if (isMC && matchRecGen)
	    {
	      TLorentzVector GenMu1;
              TLorentzVector GenMu2;
              TLorentzVector GenMu3;
              TLorentzVector GenMu4;
              TLorentzVector GenTauMu1;
              TLorentzVector GenTauMu2;

              bool findMatchedRecGenMu1 = false;
              bool findMatchedRecGenMu2 = false;
              bool findMatchedRecGenMu3 = false;
              bool findMatchedRecGenMu4 = false;
              bool findMatchedRecGenTauMu = false;
              bool findMatchedRecGenTauMu2 = false;

              unsigned int indexGenMu1 = -1;
              unsigned int indexGenMu2 = -1;
              unsigned int indexGenMu3 = -1;
	      //              unsigned int indexGenMu4 = -1;
              unsigned int indexGenTauMu1 = -1;

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
                      if (iGenMu == indexGenMu1) continue;
                      TLorentzVector GenMu2Cand;
                      GenMu2Cand.SetPtEtaPhiM(genMuonPt->at(iGenMu), genMuonEta->at(iGenMu), genMuonPhi->at(iGenMu), genMuonMass->at(iGenMu));
                      if (Mu2.DeltaR(GenMu2Cand) <= smallestDR && iGenMu != indexGenMu1)
                        {
                          smallestDR = Mu2.DeltaR(GenMu2Cand);
                          findMatchedRecGenMu2 = true;
                          GenMu2 = GenMu2Cand;
                          indexGenMu2 = iGenMu;
                        } // end if Mu2.DeltaR(GenMuCand) <= smallestDR && iGenMu != indexGenMu1
		    } // end for loop on GenMu2

	      // --------- search for matched genMu3 for Mu3 --------------      
		  smallestDR = 0.15;
		  for (unsigned int iGenMu=0; iGenMu<genMuonPt->size(); iGenMu++)
		    {
		      TLorentzVector GenMu3Cand;
		      GenMu3Cand.SetPtEtaPhiM(genMuonPt->at(iGenMu), genMuonEta->at(iGenMu), genMuonPhi->at(iGenMu), genMuonMass->at(iGenMu));
		      if (Mu3.DeltaR(GenMu3Cand) <= smallestDR && iGenMu != indexGenMu1 && iGenMu != indexGenMu2)
			{
			  smallestDR = Mu3.DeltaR(GenMu3Cand);
			  findMatchedRecGenMu3 = true;
			  GenMu3 = GenMu3Cand;
			  indexGenMu3 = iGenMu;
			} // end if Mu3.DeltaR(GenMuCand) <= smallestDR && iGenMu != indexGenMu1 && iGenMu != indexGenMu2
		    } // looping over GenMu for mu 3

              // --------- search for matched genMu4 for Mu4 --------------  
		  smallestDR = 0.15;
		  for (unsigned int iGenMu=0; iGenMu<genMuonPt->size(); iGenMu++)
		    {
		      if (iGenMu == indexGenMu1) continue;
		      if (iGenMu == indexGenMu2) continue;
		      if (iGenMu == indexGenMu3) continue;
		      TLorentzVector GenMu4Cand;
		      GenMu4Cand.SetPtEtaPhiM(genMuonPt->at(iGenMu), genMuonEta->at(iGenMu), genMuonPhi->at(iGenMu), genMuonMass->at(iGenMu));
		      if (Mu4.DeltaR(GenMu4Cand) <= smallestDR)
			{
			  smallestDR = Mu4.DeltaR(GenMu4Cand);
			  findMatchedRecGenMu4 = true;
			  GenMu4 = GenMu4Cand;
			} // if mu4.deltaR(genMuCand) <= smallestDR)
		    }// end for loop on GenMu                                   
		} // gen muon pt-> size() > 0


	      // --------- search for matched genTauMu for Mu3 --------------  
	      if (genTauMuPt->size()>0)
                {
		  double smallestDR = 0.15;
                  for (unsigned int iGenTauMu=0; iGenTauMu<genTauMuPt->size(); iGenTauMu++)
                    {
                      TLorentzVector GenTauMuCand;
                      GenTauMuCand.SetPtEtaPhiM(genTauMuPt->at(iGenTauMu), genTauMuEta->at(iGenTauMu), genTauMuPhi->at(iGenTauMu), genTauMuMass->at(iGenTauMu));
                      if (Mu3.DeltaR(GenTauMuCand) <= smallestDR)
                        {
                          smallestDR = Mu3.DeltaR(GenTauMuCand);
                          findMatchedRecGenTauMu = true;
                          GenTauMu1 = GenTauMuCand;
                          indexGenTauMu1 = iGenTauMu;
                        } // end if Mu3.DeltaR(GenTauMuCand) <= smallestDR
		    } // end for loop on GenTauMu
		} // end if genTauMuPt->size()>0    

	      // --------- search for matched genTauMu2 for Mu4 --------------  
	      if (genTauMuPt->size()>1)
                {
		  double smallestDR = 0.15;
                  for (unsigned int iGenTauMu=0; iGenTauMu<genTauMuPt->size(); iGenTauMu++)
                    {
                      if (iGenTauMu == indexGenTauMu1) continue;
                      TLorentzVector GenTauMu2Cand;
                      GenTauMu2Cand.SetPtEtaPhiM(genTauMuPt->at(iGenTauMu), genTauMuEta->at(iGenTauMu), genTauMuPhi->at(iGenTauMu), genTauMuMass->at(iGenTauMu));
                      if (Mu4.DeltaR(GenTauMu2Cand) <= smallestDR && iGenTauMu != indexGenTauMu1)
                        {
                          smallestDR = Mu4.DeltaR(GenTauMu2Cand);
                          findMatchedRecGenTauMu2 = true;
                          GenTauMu2 = GenTauMu2Cand;
                        } // end if Mu4.DeltaR(GenTauMuCand) <= smallestDR
		    } // end for loop on GenTauMu
		} // end if genTauMuPt->size()>0   

	      if(findMatchedRecGenMu1 && findMatchedRecGenMu2 && findMatchedRecGenMu3 && findMatchedRecGenMu4 && findMatchedRecGenTauMu && findMatchedRecGenTauMu2){
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

		genmu4Pt->Fill(GenMu4.Pt(), weight);
		genmu4Eta->Fill(GenMu4.Eta(), weight);
		genmu4Phi->Fill(GenMu4.Phi(), weight);
		genmu4Mass->Fill(GenMu4.M(), weight);
		mu4PtVSGenMu4Pt->Fill(Mu4.Pt(), GenMu4.Pt(), weight);
		mu4EtaVSGenMu4Eta->Fill(Mu4.Eta(), GenMu4.Eta(), weight);
		mu4PhiVSGenMu4Phi->Fill(Mu4.Phi(), GenMu4.Phi(), weight);

		gentauMu1Pt->Fill(GenTauMu1.Pt(), weight);
		gentauMu1Eta->Fill(GenTauMu1.Eta(), weight);
		gentauMu1Phi->Fill(GenTauMu1.Phi(), weight);
		gentauMu1Mass->Fill(GenTauMu1.M(), weight);

		gentauMu2Pt->Fill(GenTauMu2.Pt(), weight);
		gentauMu2Eta->Fill(GenTauMu2.Eta(), weight);
		gentauMu2Phi->Fill(GenTauMu2.Phi(), weight);
		gentauMu2Mass->Fill(GenTauMu2.M(), weight);

		dRgenMu1genMu2->Fill(GenMu1.DeltaR(GenMu2), weight);
		dRgenMu3genTauMu1->Fill(GenMu3.DeltaR(GenTauMu1), weight);
		dRgenMu4genTauMu2->Fill(GenMu4.DeltaR(GenTauMu2), weight);
		dRgenTaugenTau->Fill(GenTauMu1.DeltaR(GenTauMu2), weight);
		dRgenMu3genMu4->Fill(GenMu3.DeltaR(GenMu4), weight);

		invMassgenMu1genMu2->Fill((GenMu1+GenMu2).M(),weight);
		invMassgenTauMugenTauMu->Fill((GenTauMu1+GenTauMu2).M(), weight);
		ptgenMu1genMu2->Fill((GenMu1+GenMu2).Pt(), weight);
		ptgenTaugenTau->Fill((GenTauMu1+GenTauMu2).Pt(), weight);
		ptgenMuMuTauMuTauMu->Fill((GenMu1+GenMu2+GenTauMu1+GenTauMu2).Pt(), weight);
		dRInvMassgenMu1genMu2->Fill(GenMu1.DeltaR(GenMu2), (GenMu1+GenMu2).M(), weight);
		dRInvMassgenTaugenTau->Fill(GenTauMu1.DeltaR(GenTauMu2), (GenTauMu1+GenTauMu2).M(), weight);

	      } // if find rec gen for 4 muons
	    } // isMC && matchRecGen
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

   TreeMuMuTauTau->Write("TreeMuMuTauTau", TObject::kOverwrite);
   outputFile->Close();
}
