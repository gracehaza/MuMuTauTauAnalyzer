#ifndef _Histomutau_h_
#define _Histomutau_h_

#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include <vector>
#include <string>
#include "TArray.h"
using namespace std;

class Histomutau{

    public:

        vector<TH1*> histColl;

        TH1D* newTH1D(string, string, int, double*);
        TH1D* newTH1D(string, string, int, double, double);
        TH1D* newTH1D(string, string, vector<double>&);

        TH2D* newTH2D(string, string, string, int, double*, int, double*);
        TH2D* newTH2D(string, string, string, int, double*, int, double, double);
        TH2D* newTH2D(string, string, string, int, double, double, int, double*);
        TH2D* newTH2D(string, string, string, int, double, double, int, double, double);

        Histomutau();
        ~Histomutau();

        TH1D* dRMu1Mu2;
        TH1D* dRMu3Mu4;
        TH1D* dRMu3Ele;
        TH1D* dREleEle;
        TH1D* dRMu3Tau;
        TH1D* dREleTau;
        TH1D* dRTauTau;

        TH1D* invMassMu1Mu2;
        TH1D* invMassMu3Mu4;
        TH1D* invMassMu3Ele;
        TH1D* invMassEleEle;
        TH1D* invMassMu3Tau;
        TH1D* invMassEleTau;
        TH1D* invMassTauTau;

        TH1D* invMassMuMuTauMuTauMu;
        TH1D* invMassMuMuTauMuTauEle;
        TH1D* invMassMuMuTauEleTauEle;
        TH1D* invMassMuMuTauMuTauHad;
        TH1D* invMassMuMuTauEleTauHad;
        TH1D* invMassMuMuTauHadTauHad;

        TH1D* ptMu1Mu2;
        TH1D* ptMu3Mu4;
        TH1D* ptMu3Ele;
        TH1D* ptEleEle;
        TH1D* ptMu3Tau;
        TH1D* ptEleTau;
        TH1D* ptTauTau;

        TH1D* ptMuMuTauMuTauMu;
        TH1D* ptMuMuTauMuTauEle;
        TH1D* ptMuMuTauEleTauEle;
        TH1D* ptMuMuTauMuTauHad;
        TH1D* ptMuMuTauEleTauHad;
        TH1D* ptMuMuTauHadTauHad;

        TH1D* mu1Iso;
        TH1D* mu2Iso;
        TH1D* mu3Iso;
        TH1D* mu4Iso;
        TH1D* ele1Iso;
        TH1D* ele2Iso;
        TH1D* tauIsoMVA;
        TH1D* tau2IsoMVA;

        TH1D* mu1Pt;
        TH1D* mu1Eta;
        TH1D* mu1Phi;
        
        TH1D* mu2Pt;
        TH1D* mu2Eta;
        TH1D* mu2Phi;

        TH1D* mu3Pt;
        TH1D* mu3Eta;
        TH1D* mu3Phi;

        TH1D* mu4Pt;
        TH1D* mu4Eta;
        TH1D* mu4Phi;

        TH1D* ele1Pt;
        TH1D* ele1Eta;
        TH1D* ele1Phi;

        TH1D* ele2Pt;
        TH1D* ele2Eta;
        TH1D* ele2Phi;

        TH1D* tauPt;
        TH1D* tauEta;
        TH1D* tauPhi;
        TH1D* tauMass;
        TH1D* tauDecayMode;

        TH1D* tau2Pt;
        TH1D* tau2Eta;
        TH1D* tau2Phi;
        TH1D* tau2Mass;
        TH1D* tau2DecayMode;

        TH1D* dRMu1Mu3;
        TH1D* dRMu1Mu4;
        TH1D* dRMu1Ele1;
        TH1D* dRMu1Ele2;
        TH1D* dRMu1Tau;
        TH1D* dRMu1Tau2;

        TH1D* dRMu2Mu3;
        TH1D* dRMu2Mu4;
        TH1D* dRMu2Ele1;
        TH1D* dRMu2Ele2;
        TH1D* dRMu2Tau;
        TH1D* dRMu2Tau2;

        TH2D* dRInvMassMu1Mu2;
        TH2D* dRInvMassMu3Mu4;
        TH2D* dRInvMassMu3Ele;
        TH2D* dRInvMassEleEle;
        TH2D* dRInvMassMu3Tau;
        TH2D* dRInvMassEleTau;
        TH2D* dRInvMassTauTau;

	TH1D* dRgenMu1genMu2;
        TH1D* dRgenMu3genMu4;
        TH1D* dRgenMu3genEle;
        TH1D* dRgenElegenEle;
        TH1D* dRgenMu3genTau;
        TH1D* dRgenElegenTau;
        TH1D* dRgenTaugenTau;

        TH1D* invMassgenMu1genMu2;
        TH1D* invMassgenMu3genMu4;
        TH1D* invMassgenMu3genEle;
        TH1D* invMassgenElegenEle;
        TH1D* invMassgenMu3genTau;
        TH1D* invMassgenElegenTau;
        TH1D* invMassgenTaugenTau;

	TH1D* invMassgenMuMuTauMuTauMu;
        TH1D* invMassgenMuMuTauMuTauEle;
        TH1D* invMassgenMuMuTauEleTauEle;
        TH1D* invMassgenMuMuTauMuTauHad;
        TH1D* invMassgenMuMuTauEleTauHad;
        TH1D* invMassgenMuMuTauHadTauHad;

        TH1D* ptgenMu1genMu2;
        TH1D* ptgenMu3genMu4;
        TH1D* ptgenMu3genEle;
        TH1D* ptgenElegenEle;
        TH1D* ptgenMu3genTau;
        TH1D* ptgenElegenTau;
        TH1D* ptgenTaugenTau;

        TH1D* ptgenMuMuTauMuTauMu;
        TH1D* ptgenMuMuTauMuTauEle;
        TH1D* ptgenMuMuTauEleTauEle;
        TH1D* ptgenMuMuTauMuTauHad;
        TH1D* ptgenMuMuTauEleTauHad;
        TH1D* ptgenMuMuTauHadTauHad;

	TH1D* genmuPt;

        TH1D* genmu1Pt;
        TH1D* genmu1Eta;
        TH1D* genmu1Phi;
	TH1D* genmu1Mass;

	TH1D* genmu2Pt;
        TH1D* genmu2Eta;
        TH1D* genmu2Phi;
	TH1D* genmu2Mass;

        TH1D* genmu3Pt;
        TH1D* genmu3Eta;
        TH1D* genmu3Phi;
	TH1D* genmu3Mass;

        TH1D* genmu4Pt;
        TH1D* genmu4Eta;
        TH1D* genmu4Phi;
	TH1D* genmu4Mass;

        TH1D* genele1Pt;
        TH1D* genele1Eta;
        TH1D* genele1Phi;
	TH1D* genele1Mass;

        TH1D* genele2Pt;
        TH1D* genele2Eta;
        TH1D* genele2Phi;
	TH1D* genele2Mass;

	TH1D* gentauPt;
        TH1D* gentauEta;
        TH1D* gentauPhi;
        TH1D* gentauMass;

	TH1D* gentauElePt;
        TH1D* gentauEleEta;
        TH1D* gentauElePhi;
        TH1D* gentauEleMass;

	TH1D* gentauMuPt;
        TH1D* gentauMuEta;
        TH1D* gentauMuPhi;
        TH1D* gentauMuMass;
	
	TH1D* recomatchedtauPt;

        TH1D* gentau2Pt;
        TH1D* gentau2Eta;
        TH1D* gentau2Phi;
        TH1D* gentau2Mass;

        TH1D* dRgenMu1genMu3;
        TH1D* dRgenMu1genMu4;
        TH1D* dRgenMu1genEle1;
        TH1D* dRgenMu1genEle2;
        TH1D* dRgenMu1genTau;
        TH1D* dRgenMu1genTau2;

	TH1D* dRgenMu2genMu3;
        TH1D* dRgenMu2genMu4;
        TH1D* dRgenMu2genEle1;
        TH1D* dRgenMu2genEle2;
        TH1D* dRgenMu2genTau;
        TH1D* dRgenMu2genTau2;

        TH2D* dRInvMassgenMu1genMu2;
        TH2D* dRInvMassgenMu3genMu4;
        TH2D* dRInvMassgenMu3genEle;
        TH2D* dRInvMassgenElegenEle;
        TH2D* dRInvMassgenMu3genTau;
        TH2D* dRInvMassgenElegenTau;
        TH2D* dRInvMassgenTaugenTau;

	TH2D* mu1PtVSGenMu1Pt;
        TH2D* mu1EtaVSGenMu1Eta;
        TH2D* mu1PhiVSGenMu1Phi;

	TH2D* mu2PtVSGenMu2Pt;
        TH2D* mu2EtaVSGenMu2Eta;
        TH2D* mu2PhiVSGenMu2Phi;

	TH2D* mu3PtVSGenMu3Pt;
	TH2D* mu3EtaVSGenMu3Eta;
	TH2D* mu3PhiVSGenMu3Phi;

	TH2D* mu3PtVSGenTauMuPt;

	TH2D* tauPtVSGenTauHadPt;
        TH2D* tauEtaVSGenTauHadEta;
        TH2D* tauPhiVSGenTauHadPhi;
        TH2D* tauPtVSGenTauHadVisPt;

	TH2D* elePtVSGenElePt;
	TH2D* elePtVSGenTauElePt;
	TH2D* eleEtaVSGenEleEta;
	TH2D* elePhiVSGenElePhi;
	TH2D* elePtVSGenTauHadVisPt;


};

#endif
