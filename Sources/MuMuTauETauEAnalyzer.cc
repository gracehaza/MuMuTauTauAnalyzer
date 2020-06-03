#define MuMuTauETauEAnalyzer_cxx
#include "MuMuTauETauEAnalyzer.h"
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

void MuMuTauETauEAnalyzer::Loop()
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
      TLorentzVector Ele1;
      TLorentzVector Ele2;

      float Mu1Iso;
      float Mu2Iso;
      float Ele1Iso;
      float Ele2Iso;

      unsigned int indexMu1 = -1;

      TLorentzVector genMu;

      bool matchRecGen = true;
      // ============================================================================

      // ---- start loop on muon candidates for mu1 ----
      // ---- prepare event weight info ----                                                                                                                                                                           
      //      double weight = 1;

      //      bool findgenMu = false;                                                                                                                                                                                  

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
              findMu2 = true;
          } // end if pair candidates
      } // end loop for mu2
          
      if (!findMu2) continue;

      bool findMuElePair = false;
      // ------- start loop on electron candidates -------
      for (unsigned int iEle=0; iEle<recoElectronPt->size(); iEle++)
      {
	double Ele1IsoThresholdpT;
	if(abs(recoElectronEta->at(iEle)) > 1.479){
	  Ele1IsoThresholdpT = Ele1IsoThresholdEndcap + 0.963/(recoElectronPt->at(iEle));}
	else{
	 Ele1IsoThresholdpT = Ele1IsoThresholdBarrel + 0.506/(recoElectronPt->at(iEle));} 
	//	double Ele1IsoThresholdpT = 0.3;

          if ((invertedEle1Iso == false && recoElectronIsolation->at(iEle) > Ele1IsoThresholdpT) || (invertedEle1Iso == true && recoElectronIsolation->at(iEle) < Ele1IsoThresholdpT)) continue;
          TLorentzVector Ele1Cand;
          Ele1Cand.SetPtEtaPhiE(recoElectronPt->at(iEle), recoElectronEta->at(iEle), recoElectronPhi->at(iEle), recoElectronEnergy->at(iEle));

          if (Ele1Cand.DeltaR(Mu1) < 0.4 || Ele1Cand.DeltaR(Mu2) < 0.4) continue;
          Ele1.SetPtEtaPhiE(recoElectronPt->at(iEle), recoElectronEta->at(iEle), recoElectronPhi->at(iEle), recoElectronEnergy->at(iEle));
          Ele1Iso = recoElectronIsolation->at(iEle);
          
          float smallestDR = 1.0; // dR cut between Ele1 and Ele2
          bool findEle2 = false;

          for (unsigned int iEle2=iEle+1; iEle2<recoElectronPt->size(); iEle2++)
          {
	    double Ele2IsoThresholdpT;// = 0.3;	   
	    if(abs(recoElectronEta->at(iEle)) > 1.479){
	      Ele2IsoThresholdpT = Ele1IsoThresholdEndcap + 0.963/(recoElectronPt->at(iEle));}
	    else{
	      Ele2IsoThresholdpT = Ele1IsoThresholdBarrel + 0.506/(recoElectronPt->at(iEle));}
              TLorentzVector Ele2Cand; // prepare this variable for dR(Ele1, Ele2) implementation
              Ele2Cand.SetPtEtaPhiE(recoElectronPt->at(iEle2), recoElectronEta->at(iEle2), recoElectronPhi->at(iEle2), recoElectronEnergy->at(iEle2));
              if ((Ele1.DeltaR(Ele2Cand) < smallestDR) && (recoElectronPDGId->at(iEle) == (-1) * recoElectronPDGId->at(iEle2)) && ((Ele1+Ele2Cand).M() < 60.0) && (Ele2Cand.DeltaR(Mu1) > 0.4) && (Ele2Cand.DeltaR(Mu2) > 0.4) && (recoElectronIsolation->at(iEle2) < Ele2IsoThresholdpT))
              {
                  Ele2.SetPtEtaPhiE(recoElectronPt->at(iEle2), recoElectronEta->at(iEle2), recoElectronPhi->at(iEle2), recoElectronEnergy->at(iEle2));
                  Ele2Iso = recoElectronIsolation->at(iEle2);
                  smallestDR = Ele1.DeltaR(Ele2);
                  findEle2 = true;
              } // end if find ele2 with electron matched
          } // end loop for ele2

          if (!findEle2) continue;
          else{
              findMuElePair = true;
              break;
          } // end if findEle2
      } // end loop for electron

      
      // ---- prepare event weight info ----
      double weight = 1;

      //      bool findgenMu = false;
      if (isMC == true)
      {
          weight *= genEventWeight; 
      }// end if isMC == true


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

          ptEleEle->Fill((Ele1+Ele2).Pt(), weight);
          dREleEle->Fill(Ele1.DeltaR(Ele2), weight);
          invMassEleEle->Fill((Ele1+Ele2).M(), weight);
          dRInvMassEleEle->Fill(Ele1.DeltaR(Ele2), (Ele1+Ele2).M(), weight);

          ele1Iso->Fill(Ele1Iso, weight);
          ele2Iso->Fill(Ele2Iso, weight);

          ele1Pt->Fill(Ele1.Pt(), weight);
          ele1Eta->Fill(Ele1.Eta(), weight);
          ele1Phi->Fill(Ele1.Phi(), weight);

          ele2Pt->Fill(Ele2.Pt(), weight);
          ele2Eta->Fill(Ele2.Eta(), weight);
          ele2Phi->Fill(Ele2.Phi(), weight);

          dRMu1Ele1->Fill(Mu1.DeltaR(Ele1), weight);
          dRMu1Ele2->Fill(Mu1.DeltaR(Ele2), weight);
          dRMu2Ele1->Fill(Mu2.DeltaR(Ele1), weight);
          dRMu2Ele2->Fill(Mu2.DeltaR(Ele2), weight);

          ptMuMuTauEleTauEle->Fill((Mu1+Mu2+Ele1+Ele2).Pt(), weight);
          invMassMuMuTauEleTauEle->Fill((Mu1+Mu2+Ele1+Ele2).M(), weight);
	  if (isMC && matchRecGen)
            {
              TLorentzVector GenMu1;
              TLorentzVector GenMu2;
              TLorentzVector GenEle1;
	      TLorentzVector GenEle2;
              TLorentzVector GenTauEle1;
              TLorentzVector GenTauEle2;

              bool findMatchedRecGenMu1 = false;
              bool findMatchedRecGenMu2 = false;
              bool findMatchedRecGenEle = false;
	      bool findMatchedRecGenEle2 = false;    
	      bool findMatchedRecGenTauEle = false;
              bool findMatchedRecGenTauEle2 = false;

              //double GenTauHadVisiblePt = 0;
              unsigned int indexGenMu1 = -1;
              unsigned int indexGenMu2 = -1;
	      unsigned int indexGenEle = -1;
	      unsigned int indexGenEle2 = -1;
	      unsigned int indexGenTauEle = -1;
	      unsigned int indexGenTauEle2 = -1;

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
		}
		  // --------- search for matched genEle for Ele --------------                                                                                      
                  smallestDR = 0.15;
                  for (unsigned int iGenEle=0; iGenEle<genElectronPt->size(); iGenEle++)
                    {
                      TLorentzVector GenEleCand;
                      GenEleCand.SetPtEtaPhiM(genElectronPt->at(iGenEle), genElectronEta->at(iGenEle), genElectronPhi->at(iGenEle), genElectronMass->at(iGenEle));
                      if (Ele1.DeltaR(GenEleCand) <= smallestDR)// && iGenEle != indexGenEle1 && iGenEle != indexGenEle2)                                             
                        {
                          smallestDR = Ele1.DeltaR(GenEleCand);
                          findMatchedRecGenEle = true;
                          GenEle1 = GenEleCand;
			  indexGenEle = iGenEle;
                        } // end if Ele.DeltaR(GenMuCand) <= smallestDR && iGenMu != indexGenMu1 && iGenMu != indexGenMu2                                            
                    } // end for loop on GenEle                                                                                                                      
		  // } // end if genMuonPt->size()>0

	      // --------- search for matched genEle2 for Ele2 --------------                                                                                                       
	      smallestDR = 0.15;
	      for (unsigned int iGenEle=0; iGenEle<genElectronPt->size(); iGenEle++)
		{
		  if (iGenEle == indexGenEle) continue;
		  TLorentzVector GenEle2Cand;
		  GenEle2Cand.SetPtEtaPhiM(genElectronPt->at(iGenEle), genElectronEta->at(iGenEle), genElectronPhi->at(iGenEle), genElectronMass->at(iGenEle));
		  if (Ele2.DeltaR(GenEle2Cand) <= smallestDR && iGenEle != indexGenEle)                                                              
		    {
		      smallestDR = Ele2.DeltaR(GenEle2Cand);
		      findMatchedRecGenEle2 = true;
		      GenEle2 = GenEle2Cand;
		    } // end if Ele.DeltaR(GenMuCand) <= smallestDR && iGenMu != indexGenMu1 && iGenMu != indexGenMu2                                                             
		}// end for loop on GenEle                                                                                                                                       
	      // end if genMuonPt->size()>0  


	      if (genTauElePt->size()>0)
                {
                  // --------- search for matched genTauEle for Ele --------------                                                                                                 
                  double smallestDR = 0.15;
                  for (unsigned int iGenTauEle=0; iGenTauEle<genTauElePt->size(); iGenTauEle++)
                    {
                      TLorentzVector GenTauEleCand;
                      GenTauEleCand.SetPtEtaPhiM(genTauElePt->at(iGenTauEle), genTauEleEta->at(iGenTauEle), genTauElePhi->at(iGenTauEle), genTauEleMass->at(iGenTauEle));
                      if (Ele1.DeltaR(GenTauEleCand) <= smallestDR)
                        {
                          smallestDR = Ele1.DeltaR(GenTauEleCand);
                          findMatchedRecGenTauEle = true;
                          GenTauEle1 = GenTauEleCand;
			  indexGenTauEle = iGenTauEle;
                        } // end if Ele.DeltaR(GenTauEleCand) <= smallestDR         
                    } // end for loop on GenTauEle
		} // end if genTauElePt->size()>0  

	      if (genTauElePt->size()>1)
                {
                  // --------- search for matched genTauEle for Ele --------------                                                                        
                  double smallestDR = 0.15;
                  for (unsigned int iGenTauEle=0; iGenTauEle<genTauElePt->size(); iGenTauEle++)
                    {
		      if (iGenTauEle == indexGenTauEle) continue;
                      TLorentzVector GenTauEle2Cand;
                      GenTauEle2Cand.SetPtEtaPhiM(genTauElePt->at(iGenTauEle), genTauEleEta->at(iGenTauEle), genTauElePhi->at(iGenTauEle), genTauEleMass->at(iGenTauEle));
                      if (Ele2.DeltaR(GenTauEle2Cand) <= smallestDR && iGenTauEle != indexGenTauEle)
                        {
                          smallestDR = Ele2.DeltaR(GenTauEle2Cand);
                          findMatchedRecGenTauEle2 = true;
                          GenTauEle2 = GenTauEle2Cand;
			} // end if Ele.DeltaR(GenTauEleCand) <= smallestDR                                                                               
                    } // end for loop on GenTauEle                                                                                                        
                } // end if genTauElePt->size()>0   

              if(findMatchedRecGenMu1 && findMatchedRecGenMu2 && findMatchedRecGenEle && findMatchedRecGenEle2 && findMatchedRecGenTauEle && findMatchedRecGenTauEle2){

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

		genele1Pt->Fill(GenEle1.Pt(), weight);
                genele1Eta->Fill(GenEle1.Eta(), weight);
                genele1Phi->Fill(GenEle1.Phi(), weight);
                genele1Mass->Fill(GenEle1.M(), weight);
                ele1PtVSGenEle1Pt->Fill(Ele1.Pt(), GenEle1.Pt(), weight);
                ele1EtaVSGenEle1Eta->Fill(Ele1.Eta(), GenEle1.Eta(), weight);
                ele1PhiVSGenEle1Phi->Fill(Ele1.Phi(), GenEle1.Phi(), weight);

                genele2Pt->Fill(GenEle2.Pt(), weight);
                genele2Eta->Fill(GenEle2.Eta(), weight);
                genele2Phi->Fill(GenEle2.Phi(), weight);
                genele2Mass->Fill(GenEle2.M(), weight);
                ele2PtVSGenEle2Pt->Fill(Ele2.Pt(), GenEle2.Pt(), weight);
                ele2EtaVSGenEle2Eta->Fill(Ele2.Eta(), GenEle2.Eta(), weight);
                ele2PhiVSGenEle2Phi->Fill(Ele2.Phi(), GenEle2.Phi(), weight);

		gentauEle1Pt->Fill(GenTauEle1.Pt(), weight);
                gentauEle1Eta->Fill(GenTauEle1.Eta(), weight);
                gentauEle1Phi->Fill(GenTauEle1.Phi(), weight);
                gentauEle1Mass->Fill(GenTauEle1.M(), weight);

		gentauEle2Pt->Fill(GenTauEle2.Pt(), weight);
                gentauEle2Eta->Fill(GenTauEle2.Eta(), weight);
                gentauEle2Phi->Fill(GenTauEle2.Phi(), weight);
                gentauEle2Mass->Fill(GenTauEle2.M(), weight);

		dRgenMu1genMu2->Fill(GenMu1.DeltaR(GenMu2), weight);
		dRgenMu1genTauEle2->Fill(GenMu1.DeltaR(GenTauEle2), weight);
		dRgenMu2genTauEle1->Fill(GenMu2.DeltaR(GenTauEle1), weight);
                dRgenMu2genTauEle2->Fill(GenMu2.DeltaR(GenTauEle2), weight);
                dRgenTaugenTau->Fill(GenTauEle1.DeltaR(GenTauEle2), weight);
		dRgenEle1genEle2->Fill(GenEle1.DeltaR(GenEle2), weight);

		invMassgenMu1genMu2->Fill((GenMu1+GenMu2).M(),weight);
		invMassgenTauElegenTauEle->Fill((GenTauEle1+GenTauEle2).M(), weight);

		ptgenMu1genMu2->Fill((GenMu1+GenMu2).Pt(), weight);
		ptgenTaugenTau->Fill((GenTauEle1+GenTauEle2).Pt(), weight);

		ptgenMuMuTauEleTauEle->Fill((GenMu1+GenMu2+GenTauEle1+GenTauEle2).Pt(), weight);

		dRInvMassgenMu1genMu2->Fill(GenMu1.DeltaR(GenMu2), (GenMu1+GenMu2).M(), weight);
        
                dRInvMassgenTaugenTau->Fill(GenTauEle1.DeltaR(GenTauEle2), (GenTauEle1+GenTauEle2).M(), weight);
	


	      } // if all gen reco matching satisfied
	    } //  end if MC and matchGenReco
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

   outputFile->Close();
}
