//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 18 14:19:21 2019 by ROOT version 6.10/08
// from TTree objectTree/objectTree
// found on file: MuMuTauTauTreelization.root
//////////////////////////////////////////////////////////

#ifndef MuTauAnalyzer_h
#define MuTauAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "Histomutau.h"
#include <TString.h>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

class MuTauAnalyzer : public Histomutau {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *recoMuonPt;
   vector<float>   *recoMuonEta;
   vector<float>   *recoMuonPhi;
   vector<float>   *recoMuonEnergy;
   vector<int>     *recoMuonPDGId;
   vector<float>   *recoMuonIsolation;
   vector<float>   *recoTauPt;
   vector<float>   *recoTauEta;
   vector<float>   *recoTauPhi;
   vector<float>   *recoTauEnergy;
   vector<int>     *recoTauPDGId;
   vector<float>   *recoTauDecayModeFinding;
   vector<float>   *recoTauIsoMVArawValue;
   vector<float>   *recoTauIsoMVAVVLoose;
   vector<float>   *recoTauIsoMVAVLoose;
   vector<float>   *recoTauIsoMVALoose;
   vector<float>   *recoTauIsoMVAMedium;
   vector<float>   *recoTauIsoMVATight;
   vector<float>   *recoTauIsoMVAVTight;
   vector<float>   *recoTauIsoMVAVVTight;
   vector<float>   *recoJetPt;
   vector<float>   *recoJetEta;
   vector<float>   *recoJetPhi;
   vector<float>   *recoJetEnergy;
   vector<float>   *recoJetCSV;
   vector<float>   *recoMET;
   vector<float>   *recoMETPhi;
   Int_t           recoNPrimaryVertex;
   Int_t           recoNPU;
   Int_t           trueNInteraction;
   Float_t         genEventWeight;

   // List of branches
   TBranch        *b_recoMuonPt;   //!
   TBranch        *b_recoMuonEta;   //!
   TBranch        *b_recoMuonPhi;   //!
   TBranch        *b_recoMuonEnergy;   //!
   TBranch        *b_recoMuonPDGId;   //!
   TBranch        *b_recoMuonIsolation;   //!
   TBranch        *b_recoTauPt;   //!
   TBranch        *b_recoTauEta;   //!
   TBranch        *b_recoTauPhi;   //!
   TBranch        *b_recoTauEnergy;   //!
   TBranch        *b_recoTauPDGId;   //!
   TBranch        *b_recoTauDecayModeFinding;   //!
   TBranch        *b_recoTauIsoMVArawValue;   //!
   TBranch        *b_recoTauIsoMVAVVLoose;   //!
   TBranch        *b_recoTauIsoMVAVLoose;   //!
   TBranch        *b_recoTauIsoMVALoose;   //!
   TBranch        *b_recoTauIsoMVAMedium;   //!
   TBranch        *b_recoTauIsoMVATight;   //!
   TBranch        *b_recoTauIsoMVAVTight;   //!
   TBranch        *b_recoTauIsoMVAVVTight;   //!
   TBranch        *b_recoJetPt;   //!
   TBranch        *b_recoJetEta;   //!
   TBranch        *b_recoJetPhi;   //!
   TBranch        *b_recoJetEnergy;   //!
   TBranch        *b_recoJetCSV;   //!
   TBranch        *b_recoMET;   //!
   TBranch        *b_recoMETPhi;   //!
   TBranch        *b_recoNPrimaryVertex;   //!
   TBranch        *b_recoNPU;   //!
   TBranch        *b_trueNInteraction;   //!
   TBranch        *b_genEventWeight;   //!

   TString fileName;
   float lumiScale;
   float summedWeights; // these two factors contribute to the MC normalization

   MuTauAnalyzer(TString fileName_, float lumiScale_, float summedWeights_ = 1.0);
   string createOutputFileName();
   virtual ~MuTauAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init();
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   bool isData;
};

#endif

#ifdef MuTauAnalyzer_cxx
MuTauAnalyzer::MuTauAnalyzer(TString fileName_, float lumiScale_, float summedWeights_) : Histomutau() 
{
    fileName = fileName_;
    lumiScale = lumiScale_;
    summedWeights = summedWeights_;

    TChain *chain = new TChain("", "");
    TString treePath = fileName + ".root/MuMuTauMuTauHadAnalyzer/objectTree";
    chain->Add(treePath);
    fChain = chain;
    Init();
    isData = ((fileName.Index("Data") >= 0) || (fileName.Index("data") >= 0));
}

string MuTauAnalyzer::createOutputFileName()
{
    ostringstream outputName;
    fileName.ReplaceAll("Tree/","");
    outputName << "Histogram/";
    outputName << fileName;
    outputName << "_histogram";
    outputName << ".root";
    return outputName.str();
}

MuTauAnalyzer::~MuTauAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MuTauAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MuTauAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MuTauAnalyzer::Init()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   recoMuonPt = 0;
   recoMuonEta = 0;
   recoMuonPhi = 0;
   recoMuonEnergy = 0;
   recoMuonPDGId = 0;
   recoMuonIsolation = 0;
   recoTauPt = 0;
   recoTauEta = 0;
   recoTauPhi = 0;
   recoTauEnergy = 0;
   recoTauPDGId = 0;
   recoTauDecayModeFinding = 0;
   recoTauIsoMVArawValue = 0;
   recoTauIsoMVAVVLoose = 0;
   recoTauIsoMVAVLoose = 0;
   recoTauIsoMVALoose = 0;
   recoTauIsoMVAMedium = 0;
   recoTauIsoMVATight = 0;
   recoTauIsoMVAVTight = 0;
   recoTauIsoMVAVVTight = 0;
   recoJetPt = 0;
   recoJetEta = 0;
   recoJetPhi = 0;
   recoJetEnergy = 0;
   recoJetCSV = 0;
   recoMET = 0;
   recoMETPhi = 0;
   // Set branch addresses and branch pointers
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("recoMuonPt", &recoMuonPt, &b_recoMuonPt);
   fChain->SetBranchAddress("recoMuonEta", &recoMuonEta, &b_recoMuonEta);
   fChain->SetBranchAddress("recoMuonPhi", &recoMuonPhi, &b_recoMuonPhi);
   fChain->SetBranchAddress("recoMuonEnergy", &recoMuonEnergy, &b_recoMuonEnergy);
   fChain->SetBranchAddress("recoMuonPDGId", &recoMuonPDGId, &b_recoMuonPDGId);
   fChain->SetBranchAddress("recoMuonIsolation", &recoMuonIsolation, &b_recoMuonIsolation);
   fChain->SetBranchAddress("recoTauPt", &recoTauPt, &b_recoTauPt);
   fChain->SetBranchAddress("recoTauEta", &recoTauEta, &b_recoTauEta);
   fChain->SetBranchAddress("recoTauPhi", &recoTauPhi, &b_recoTauPhi);
   fChain->SetBranchAddress("recoTauEnergy", &recoTauEnergy, &b_recoTauEnergy);
   fChain->SetBranchAddress("recoTauPDGId", &recoTauPDGId, &b_recoTauPDGId);
   fChain->SetBranchAddress("recoTauDecayModeFinding", &recoTauDecayModeFinding, &b_recoTauDecayModeFinding);
   fChain->SetBranchAddress("recoTauIsoMVArawValue", &recoTauIsoMVArawValue, &b_recoTauIsoMVArawValue);
   fChain->SetBranchAddress("recoTauIsoMVAVVLoose", &recoTauIsoMVAVVLoose, &b_recoTauIsoMVAVVLoose);
   fChain->SetBranchAddress("recoTauIsoMVAVLoose", &recoTauIsoMVAVLoose, &b_recoTauIsoMVAVLoose);
   fChain->SetBranchAddress("recoTauIsoMVALoose", &recoTauIsoMVALoose, &b_recoTauIsoMVALoose);
   fChain->SetBranchAddress("recoTauIsoMVAMedium", &recoTauIsoMVAMedium, &b_recoTauIsoMVAMedium);
   fChain->SetBranchAddress("recoTauIsoMVATight", &recoTauIsoMVATight, &b_recoTauIsoMVATight);
   fChain->SetBranchAddress("recoTauIsoMVAVTight", &recoTauIsoMVAVTight, &b_recoTauIsoMVAVTight);
   fChain->SetBranchAddress("recoTauIsoMVAVVTight", &recoTauIsoMVAVVTight, &b_recoTauIsoMVAVVTight);
   fChain->SetBranchAddress("recoJetPt", &recoJetPt, &b_recoJetPt);
   fChain->SetBranchAddress("recoJetEta", &recoJetEta, &b_recoJetEta);
   fChain->SetBranchAddress("recoJetPhi", &recoJetPhi, &b_recoJetPhi);
   fChain->SetBranchAddress("recoJetEnergy", &recoJetEnergy, &b_recoJetEnergy);
   fChain->SetBranchAddress("recoJetCSV", &recoJetCSV, &b_recoJetCSV);
   fChain->SetBranchAddress("recoMET", &recoMET, &b_recoMET);
   fChain->SetBranchAddress("recoMETPhi", &recoMETPhi, &b_recoMETPhi);
   fChain->SetBranchAddress("recoNPrimaryVertex", &recoNPrimaryVertex, &b_recoNPrimaryVertex);
   fChain->SetBranchAddress("recoNPU", &recoNPU, &b_recoNPU);
   fChain->SetBranchAddress("trueNInteraction", &trueNInteraction, &b_trueNInteraction);
   fChain->SetBranchAddress("genEventWeight", &genEventWeight, &b_genEventWeight);
   Notify();
}

Bool_t MuTauAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MuTauAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MuTauAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MuTauAnalyzer_cxx