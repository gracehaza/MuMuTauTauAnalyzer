//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov  8 17:35:56 2019 by ROOT version 6.10/08
// from TTree objectTree/objectTree
// found on file: ZTauMuTauHadTreelization_100.root
//////////////////////////////////////////////////////////

#ifndef ZTauMuTauHadAnalyzer_h
#define ZTauMuTauHadAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "HistoZmutau.h"

class ZTauMuTauHadAnalyzer : public HistoZmutau {
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
   vector<float>   *recoMuonDXY;
   vector<float>   *recoMuonDZ;
   vector<float>   *recoTauPt;
   vector<float>   *recoTauEta;
   vector<float>   *recoTauPhi;
   vector<float>   *recoTauEnergy;
   vector<int>     *recoTauPDGId;
   vector<float>   *recoTauDecayMode;
   vector<float>   *recoTauIsoMVArawValue;
   vector<float>   *recoTauIsoMVAVVLoose;
   vector<float>   *recoTauIsoMVAVLoose;
   vector<float>   *recoTauIsoMVALoose;
   vector<float>   *recoTauIsoMVAMedium;
   vector<float>   *recoTauIsoMVATight;
   vector<float>   *recoTauIsoMVAVTight;
   vector<float>   *recoTauIsoMVAVVTight;
   vector<float>   *recoTauAntiMuMVALoose;
   vector<float>   *recoTauAntiMuMVATight;
   vector<float>   *recoTauAntiEleMVArawValue;
   vector<float>   *recoTauAntiEleMVAVLoose;
   vector<float>   *recoTauAntiEleMVALoose;
   vector<float>   *recoTauAntiEleMVAMedium;
   vector<float>   *recoTauAntiEleMVATight;
   vector<float>   *recoTauAntiEleMVAVTight;
   vector<float>   *recoJetPt;
   vector<float>   *recoJetEta;
   vector<float>   *recoJetPhi;
   vector<float>   *recoJetEnergy;
   vector<float>   *recoJetCSV;
   vector<float>   *recoMET;
   vector<float>   *recoMETPhi;
   vector<float>   *recoMETPx;
   vector<float>   *recoMETPy;
   Int_t           recoNPrimaryVertex;
   vector<float>   *genMuonPt;
   vector<float>   *genMuonEta;
   vector<float>   *genMuonPhi;
   vector<float>   *genMuonMass;
   vector<int>     *genMuonPDGId;
   vector<int>     *genMuonMotherPDGId;
   vector<float>   *genTauMuPt;
   vector<float>   *genTauMuEta;
   vector<float>   *genTauMuPhi;
   vector<float>   *genTauMuMass;
   vector<int>     *genTauMuPDGId;
   vector<int>     *genTauMuMotherPDGId;
   vector<float>   *genTauMuVisPt;
   vector<float>   *genTauMuVisMass;
   vector<float>   *genTauHadPt;
   vector<float>   *genTauHadEta;
   vector<float>   *genTauHadPhi;
   vector<float>   *genTauHadMass;
   vector<int>     *genTauHadPDGId;
   vector<int>     *genTauHadMotherPDGId;
   vector<float>   *genTauHadVisPt;
   vector<float>   *genTauHadVisMass;
   vector<int>     *genTauHadNPionZero;
   vector<int>     *genTauHadNChargedHadrons;
   Int_t           recoNPU;
   Int_t           trueNInteraction;
   Float_t         genEventWeight;
   Int_t           eventID;

   // List of branches
   TBranch        *b_recoMuonPt;   //!
   TBranch        *b_recoMuonEta;   //!
   TBranch        *b_recoMuonPhi;   //!
   TBranch        *b_recoMuonEnergy;   //!
   TBranch        *b_recoMuonPDGId;   //!
   TBranch        *b_recoMuonIsolation;   //!
   TBranch        *b_recoMuonDXY;   //!
   TBranch        *b_recoMuonDZ;   //!
   TBranch        *b_recoTauPt;   //!
   TBranch        *b_recoTauEta;   //!
   TBranch        *b_recoTauPhi;   //!
   TBranch        *b_recoTauEnergy;   //!
   TBranch        *b_recoTauPDGId;   //!
   TBranch        *b_recoTauDecayMode;   //!
   TBranch        *b_recoTauIsoMVArawValue;   //!
   TBranch        *b_recoTauIsoMVAVVLoose;   //!
   TBranch        *b_recoTauIsoMVAVLoose;   //!
   TBranch        *b_recoTauIsoMVALoose;   //!
   TBranch        *b_recoTauIsoMVAMedium;   //!
   TBranch        *b_recoTauIsoMVATight;   //!
   TBranch        *b_recoTauIsoMVAVTight;   //!
   TBranch        *b_recoTauIsoMVAVVTight;   //!
   TBranch        *b_recoTauAntiMuMVALoose;   //!
   TBranch        *b_recoTauAntiMuMVATight;   //!
   TBranch        *b_recoTauAntiEleMVArawValue;   //!
   TBranch        *b_recoTauAntiEleMVAVLoose;   //!
   TBranch        *b_recoTauAntiEleMVALoose;   //!
   TBranch        *b_recoTauAntiEleMVAMedium;   //!
   TBranch        *b_recoTauAntiEleMVATight;   //!
   TBranch        *b_recoTauAntiEleMVAVTight;   //!
   TBranch        *b_recoJetPt;   //!
   TBranch        *b_recoJetEta;   //!
   TBranch        *b_recoJetPhi;   //!
   TBranch        *b_recoJetEnergy;   //!
   TBranch        *b_recoJetCSV;   //!
   TBranch        *b_recoMET;   //!
   TBranch        *b_recoMETPhi;   //!
   TBranch        *b_recoMETPx;   //!
   TBranch        *b_recoMETPy;   //!
   TBranch        *b_recoNPrimaryVertex;   //!
   TBranch        *b_genMuonPt;   //!
   TBranch        *b_genMuonEta;   //!
   TBranch        *b_genMuonPhi;   //!
   TBranch        *b_genMuonMass;   //!
   TBranch        *b_genMuonPDGId;   //!
   TBranch        *b_genMuonMotherPDGId;   //!
   TBranch        *b_genTauMuPt;   //!
   TBranch        *b_genTauMuEta;   //!
   TBranch        *b_genTauMuPhi;   //!
   TBranch        *b_genTauMuMass;   //!
   TBranch        *b_genTauMuPDGId;   //!
   TBranch        *b_genTauMuMotherPDGId;   //!
   TBranch        *b_genTauMuVisPt;   //!
   TBranch        *b_genTauMuVisMass;   //!
   TBranch        *b_genTauHadPt;   //!
   TBranch        *b_genTauHadEta;   //!
   TBranch        *b_genTauHadPhi;   //!
   TBranch        *b_genTauHadMass;   //!
   TBranch        *b_genTauHadPDGId;   //!
   TBranch        *b_genTauHadMotherPDGId;   //!
   TBranch        *b_genTauHadVisPt;   //!
   TBranch        *b_genTauHadVisMass;   //!
   TBranch        *b_genTauHadNPionZero;   //!
   TBranch        *b_genTauHadNChargedHadrons;   //!
   TBranch        *b_recoNPU;   //!
   TBranch        *b_trueNInteraction;   //!
   TBranch        *b_genEventWeight;   //!
   TBranch        *b_eventID;   //!

   TString fileName;
   TString outputDir;
   Long_t  nMaxEvents;
   float lumiScale;
   float summedWeights; // these two factors contribute to the MC normalization
   bool isMC;
   bool invertedMu1Iso;
   float Mu1IsoThreshold;
   bool tauMVAIsoRawORWP;
   double tauMVAIsoRawThreshold;
   TString tauMVAIsoWP;
   bool signSameOROpposite;
   float mTMuMetLowThreshold;
   float mTMuMetHighThreshold;
   bool invertedPzetaCut;
   float pzetaThreshold;
   float tauPtLowThreshold;
   float tauPtHighThreshold;
   TString tauAntiMuDisc;
   TString tauAntiEleDisc;
   TString doWhatSample;

   ZTauMuTauHadAnalyzer(TString fileName_, TString outputDir_, float lumiScale_, float summedWeights_ = 1.0, Long_t nMaxEvents_ = 0, bool isMC_ = false, bool invertedMu1Iso_ = false, float Mu1IsoThreshold_ = 0.25, bool tauMVAIsoRawORWP_ = false, double tauMVAIsoRawThreshold_ = -0.5, TString tauMVAIsoWP_ = "MEDIUM", bool signSameOROpposite_ = false, float mTMuMetLowThreshold_ = 0, float mTMuMetHighThreshold_ = 160.0, bool invertedPzetaCut_ = false, float pzetaThreshold_ = -125.0, float tauPtLowThreshold_ = 10.0, float tauPtHighThreshold_ = 10000.0, TString tauAntiMuDisc_ = "TIGHT", TString tauAntiEleDisc_ = "LOOSE", TString doWhatSample_ = "ZTT");
   string createOutputFileName();
   virtual ~ZTauMuTauHadAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init();
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ZTauMuTauHadAnalyzer_cxx
ZTauMuTauHadAnalyzer::ZTauMuTauHadAnalyzer(TString fileName_, TString outputDir_, float lumiScale_, float summedWeights_, Long_t nMaxEvents_, bool isMC_, bool invertedMu1Iso_, float Mu1IsoThreshold_, bool tauMVAIsoRawORWP_, double tauMVAIsoRawThreshold_, TString tauMVAIsoWP_, bool signSameOROpposite_, float mTMuMetLowThreshold_, float mTMuMetHighThreshold_, bool invertedPzetaCut_, float pzetaThreshold_, float tauPtLowThreshold_, float tauPtHighThreshold_, TString tauAntiMuDisc_, TString tauAntiEleDisc_, TString doWhatSample_) : HistoZmutau() 
{
    fileName = fileName_;
    outputDir = outputDir_;
    lumiScale = lumiScale_;
    summedWeights = summedWeights_;
    nMaxEvents = nMaxEvents_;
    isMC = isMC_;
    invertedMu1Iso = invertedMu1Iso_;
    Mu1IsoThreshold = Mu1IsoThreshold_;
    tauMVAIsoRawORWP = tauMVAIsoRawORWP_;
    tauMVAIsoRawThreshold = tauMVAIsoRawThreshold_;
    tauMVAIsoWP = tauMVAIsoWP_;
    signSameOROpposite = signSameOROpposite_;
    mTMuMetLowThreshold = mTMuMetLowThreshold_;
    mTMuMetHighThreshold = mTMuMetHighThreshold_;
    invertedPzetaCut = invertedPzetaCut_;
    pzetaThreshold = pzetaThreshold_;
    tauPtLowThreshold = tauPtLowThreshold_;
    tauPtHighThreshold = tauPtHighThreshold_;
    tauAntiMuDisc = tauAntiMuDisc_;
    tauAntiEleDisc = tauAntiEleDisc_;
    doWhatSample = doWhatSample_;

    //--- Create output directory if necessary ---
    if (nMaxEvents > 0) {
        outputDir.Remove(TString::kTrailing, '/');
        outputDir += TString::Format("_%ldevts/", nMaxEvents);
        cout << "Output directory has been changed to " << outputDir << endl;
    }

    TString command = "mkdir -p " + outputDir;
    system(command);

    TChain *chain = new TChain("", "");
    TString treePath = fileName + "/ZTauMuTauHadAnalyzer/objectTree";
    chain->Add(treePath);
    fChain = chain;
    Init();
}

string ZTauMuTauHadAnalyzer::createOutputFileName()
{
    ostringstream outputName;
    fileName.Replace(0, fileName.Last('/'), "");
    fileName.ReplaceAll(".root","");
    outputName << outputDir;
    outputName << "/";
    outputName << fileName;
    outputName << "_histogram";
    outputName << ".root";
    return outputName.str();
}

ZTauMuTauHadAnalyzer::~ZTauMuTauHadAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ZTauMuTauHadAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ZTauMuTauHadAnalyzer::LoadTree(Long64_t entry)
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

void ZTauMuTauHadAnalyzer::Init()
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
   recoMuonDXY = 0;
   recoMuonDZ = 0;
   recoTauPt = 0;
   recoTauEta = 0;
   recoTauPhi = 0;
   recoTauEnergy = 0;
   recoTauPDGId = 0;
   recoTauDecayMode = 0;
   recoTauIsoMVArawValue = 0;
   recoTauIsoMVAVVLoose = 0;
   recoTauIsoMVAVLoose = 0;
   recoTauIsoMVALoose = 0;
   recoTauIsoMVAMedium = 0;
   recoTauIsoMVATight = 0;
   recoTauIsoMVAVTight = 0;
   recoTauIsoMVAVVTight = 0;
   recoTauAntiMuMVALoose = 0;
   recoTauAntiMuMVATight = 0;
   recoTauAntiEleMVArawValue = 0;
   recoTauAntiEleMVAVLoose = 0;
   recoTauAntiEleMVALoose = 0;
   recoTauAntiEleMVAMedium = 0;
   recoTauAntiEleMVATight = 0;
   recoTauAntiEleMVAVTight = 0;
   recoJetPt = 0;
   recoJetEta = 0;
   recoJetPhi = 0;
   recoJetEnergy = 0;
   recoJetCSV = 0;
   recoMET = 0;
   recoMETPhi = 0;
   recoMETPx = 0;
   recoMETPy = 0;
   genMuonPt = 0;
   genMuonEta = 0;
   genMuonPhi = 0;
   genMuonMass = 0;
   genMuonPDGId = 0;
   genMuonMotherPDGId = 0;
   genTauMuPt = 0;
   genTauMuEta = 0;
   genTauMuPhi = 0;
   genTauMuMass = 0;
   genTauMuPDGId = 0;
   genTauMuMotherPDGId = 0;
   genTauMuVisPt = 0;
   genTauMuVisMass = 0;
   genTauHadPt = 0;
   genTauHadEta = 0;
   genTauHadPhi = 0;
   genTauHadMass = 0;
   genTauHadPDGId = 0;
   genTauHadMotherPDGId = 0;
   genTauHadVisPt = 0;
   genTauHadVisMass = 0;
   genTauHadNPionZero = 0;
   genTauHadNChargedHadrons = 0;
   // Set branch addresses and branch pointers
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("recoMuonPt", &recoMuonPt, &b_recoMuonPt);
   fChain->SetBranchAddress("recoMuonEta", &recoMuonEta, &b_recoMuonEta);
   fChain->SetBranchAddress("recoMuonPhi", &recoMuonPhi, &b_recoMuonPhi);
   fChain->SetBranchAddress("recoMuonEnergy", &recoMuonEnergy, &b_recoMuonEnergy);
   fChain->SetBranchAddress("recoMuonPDGId", &recoMuonPDGId, &b_recoMuonPDGId);
   fChain->SetBranchAddress("recoMuonIsolation", &recoMuonIsolation, &b_recoMuonIsolation);
   fChain->SetBranchAddress("recoMuonDXY", &recoMuonDXY, &b_recoMuonDXY);
   fChain->SetBranchAddress("recoMuonDZ", &recoMuonDZ, &b_recoMuonDZ);
   fChain->SetBranchAddress("recoTauPt", &recoTauPt, &b_recoTauPt);
   fChain->SetBranchAddress("recoTauEta", &recoTauEta, &b_recoTauEta);
   fChain->SetBranchAddress("recoTauPhi", &recoTauPhi, &b_recoTauPhi);
   fChain->SetBranchAddress("recoTauEnergy", &recoTauEnergy, &b_recoTauEnergy);
   fChain->SetBranchAddress("recoTauPDGId", &recoTauPDGId, &b_recoTauPDGId);
   fChain->SetBranchAddress("recoTauDecayMode", &recoTauDecayMode, &b_recoTauDecayMode);
   fChain->SetBranchAddress("recoTauIsoMVArawValue", &recoTauIsoMVArawValue, &b_recoTauIsoMVArawValue);
   fChain->SetBranchAddress("recoTauIsoMVAVVLoose", &recoTauIsoMVAVVLoose, &b_recoTauIsoMVAVVLoose);
   fChain->SetBranchAddress("recoTauIsoMVAVLoose", &recoTauIsoMVAVLoose, &b_recoTauIsoMVAVLoose);
   fChain->SetBranchAddress("recoTauIsoMVALoose", &recoTauIsoMVALoose, &b_recoTauIsoMVALoose);
   fChain->SetBranchAddress("recoTauIsoMVAMedium", &recoTauIsoMVAMedium, &b_recoTauIsoMVAMedium);
   fChain->SetBranchAddress("recoTauIsoMVATight", &recoTauIsoMVATight, &b_recoTauIsoMVATight);
   fChain->SetBranchAddress("recoTauIsoMVAVTight", &recoTauIsoMVAVTight, &b_recoTauIsoMVAVTight);
   fChain->SetBranchAddress("recoTauIsoMVAVVTight", &recoTauIsoMVAVVTight, &b_recoTauIsoMVAVVTight);
   fChain->SetBranchAddress("recoTauAntiMuMVALoose", &recoTauAntiMuMVALoose, &b_recoTauAntiMuMVALoose);
   fChain->SetBranchAddress("recoTauAntiMuMVATight", &recoTauAntiMuMVATight, &b_recoTauAntiMuMVATight);
   fChain->SetBranchAddress("recoTauAntiEleMVArawValue", &recoTauAntiEleMVArawValue, &b_recoTauAntiEleMVArawValue);
   fChain->SetBranchAddress("recoTauAntiEleMVAVLoose", &recoTauAntiEleMVAVLoose, &b_recoTauAntiEleMVAVLoose);
   fChain->SetBranchAddress("recoTauAntiEleMVALoose", &recoTauAntiEleMVALoose, &b_recoTauAntiEleMVALoose);
   fChain->SetBranchAddress("recoTauAntiEleMVAMedium", &recoTauAntiEleMVAMedium, &b_recoTauAntiEleMVAMedium);
   fChain->SetBranchAddress("recoTauAntiEleMVATight", &recoTauAntiEleMVATight, &b_recoTauAntiEleMVATight);
   fChain->SetBranchAddress("recoTauAntiEleMVAVTight", &recoTauAntiEleMVAVTight, &b_recoTauAntiEleMVAVTight);
   fChain->SetBranchAddress("recoJetPt", &recoJetPt, &b_recoJetPt);
   fChain->SetBranchAddress("recoJetEta", &recoJetEta, &b_recoJetEta);
   fChain->SetBranchAddress("recoJetPhi", &recoJetPhi, &b_recoJetPhi);
   fChain->SetBranchAddress("recoJetEnergy", &recoJetEnergy, &b_recoJetEnergy);
   fChain->SetBranchAddress("recoJetCSV", &recoJetCSV, &b_recoJetCSV);
   fChain->SetBranchAddress("recoMET", &recoMET, &b_recoMET);
   fChain->SetBranchAddress("recoMETPhi", &recoMETPhi, &b_recoMETPhi);
   fChain->SetBranchAddress("recoMETPx", &recoMETPx, &b_recoMETPx);
   fChain->SetBranchAddress("recoMETPy", &recoMETPy, &b_recoMETPy);
   fChain->SetBranchAddress("recoNPrimaryVertex", &recoNPrimaryVertex, &b_recoNPrimaryVertex);
   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   if (isMC) 
   {
       fChain->SetBranchAddress("genMuonPt", &genMuonPt, &b_genMuonPt);
       fChain->SetBranchAddress("genMuonEta", &genMuonEta, &b_genMuonEta);
       fChain->SetBranchAddress("genMuonPhi", &genMuonPhi, &b_genMuonPhi);
       fChain->SetBranchAddress("genMuonMass", &genMuonMass, &b_genMuonMass);
       fChain->SetBranchAddress("genMuonPDGId", &genMuonPDGId, &b_genMuonPDGId);
       fChain->SetBranchAddress("genMuonMotherPDGId", &genMuonMotherPDGId, &b_genMuonMotherPDGId);
       fChain->SetBranchAddress("genTauMuPt", &genTauMuPt, &b_genTauMuPt);
       fChain->SetBranchAddress("genTauMuEta", &genTauMuEta, &b_genTauMuEta);
       fChain->SetBranchAddress("genTauMuPhi", &genTauMuPhi, &b_genTauMuPhi);
       fChain->SetBranchAddress("genTauMuMass", &genTauMuMass, &b_genTauMuMass);
       fChain->SetBranchAddress("genTauMuPDGId", &genTauMuPDGId, &b_genTauMuPDGId);
       fChain->SetBranchAddress("genTauMuMotherPDGId", &genTauMuMotherPDGId, &b_genTauMuMotherPDGId);
       fChain->SetBranchAddress("genTauMuVisPt", &genTauMuVisPt, &b_genTauMuVisPt);
       fChain->SetBranchAddress("genTauMuVisMass", &genTauMuVisMass, &b_genTauMuVisMass);
       fChain->SetBranchAddress("genTauHadPt", &genTauHadPt, &b_genTauHadPt);
       fChain->SetBranchAddress("genTauHadEta", &genTauHadEta, &b_genTauHadEta);
       fChain->SetBranchAddress("genTauHadPhi", &genTauHadPhi, &b_genTauHadPhi);
       fChain->SetBranchAddress("genTauHadMass", &genTauHadMass, &b_genTauHadMass);
       fChain->SetBranchAddress("genTauHadPDGId", &genTauHadPDGId, &b_genTauHadPDGId);
       fChain->SetBranchAddress("genTauHadMotherPDGId", &genTauHadMotherPDGId, &b_genTauHadMotherPDGId);
       fChain->SetBranchAddress("genTauHadVisPt", &genTauHadVisPt, &b_genTauHadVisPt);
       fChain->SetBranchAddress("genTauHadVisMass", &genTauHadVisMass, &b_genTauHadVisMass);
       fChain->SetBranchAddress("genTauHadNPionZero", &genTauHadNPionZero, &b_genTauHadNPionZero);
       fChain->SetBranchAddress("genTauHadNChargedHadrons", &genTauHadNChargedHadrons, &b_genTauHadNChargedHadrons);
       fChain->SetBranchAddress("recoNPU", &recoNPU, &b_recoNPU);
       fChain->SetBranchAddress("trueNInteraction", &trueNInteraction, &b_trueNInteraction);
       fChain->SetBranchAddress("genEventWeight", &genEventWeight, &b_genEventWeight);
   } // end if isMC
   Notify();
}

Bool_t ZTauMuTauHadAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ZTauMuTauHadAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ZTauMuTauHadAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ZTauMuTauHadAnalyzer_cxx
