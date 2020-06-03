#include <TH1D.h>
#include "Histomutau.h"
#include <TString.h>
#include <TAxis.h>
#define PI 3.14159265359
using namespace std;

// ---------- customize the axis settings of TH1D and TH2D ----------------
TH1D* Histomutau::newTH1D(string name, string xTitle, int nBins, double *xBins){
    TH1D* hist = new TH1D(name.c_str(), name.c_str(), nBins, xBins);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle("# Events");
    histColl.push_back(hist);
    return hist;
}

TH1D* Histomutau::newTH1D(string name, string xTitle, vector<double>& xBinsVect)
{
    int nBins = xBinsVect.size()-1;
    double *xBins = new double[xBinsVect.size()];
    std::copy(xBinsVect.begin(), xBinsVect.end(), xBins);
    TH1D* hist = new TH1D(name.c_str(), name.c_str(), nBins, xBins);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle("# Events");
    delete [] xBins;
    histColl.push_back(hist);
    return hist;
}

TH1D* Histomutau::newTH1D(string name, string xTitle, int nBins, double xLow, double xUp){
    TH1D* hist = new TH1D(name.c_str(), name.c_str(), nBins, xLow, xUp);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle("# Events");
    hist->SetOption("HIST");
    histColl.push_back(hist);
    return hist;
}

TH2D* Histomutau::newTH2D(string name, string xTitle, string yTitle, int nBinsX, double *xBins, int nBinsY, double *yBinsY){
    TH2D* hist = new TH2D(name.c_str(), name.c_str(), nBinsX, xBins, nBinsY, yBinsY);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());
    hist->GetZaxis()->SetTitle("# Events");
    histColl.push_back(hist);
    return hist;
}

TH2D* Histomutau::newTH2D(string name, string xTitle, string yTitle, int nBinsX, double *xBins, int nBinsY, double yLow, double yUp){
    TH2D* hist = new TH2D(name.c_str(), name.c_str(), nBinsX, xBins, nBinsY, yLow, yUp);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());
    hist->GetZaxis()->SetTitle("# Events");
    histColl.push_back(hist);
    return hist;
}

TH2D* Histomutau::newTH2D(string name, string xTitle, string yTitle, int nBinsX, double xLow, double xUp, int nBinsY, double *yBins){
    TH2D* hist = new TH2D(name.c_str(), name.c_str(), nBinsX, xLow, xUp, nBinsY, yBins);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());
    hist->GetZaxis()->SetTitle("# Events");
    histColl.push_back(hist);
    return hist;
}

TH2D* Histomutau::newTH2D(string name, string xTitle, string yTitle, int nBinsX, double xLow, double xUp, int nBinsY, double yLow, double yUp){
    TH2D* hist = new TH2D(name.c_str(), name.c_str(), nBinsX, xLow, xUp, nBinsY, yLow, yUp);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());
    hist->GetZaxis()->SetTitle("# Events");
    hist->SetOption("HIST");
    histColl.push_back(hist);
    return hist;
}

// --------------- customize the binning of output histograms -----------------
Histomutau::Histomutau(){

  double Mu1PtBin [] = {3, 5,10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150};
  double Mu2PtBin [] = {3, 5,10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150};
    double Mu3PtBin [] = {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56, 66, 78, 92, 108, 200}; 
    double Mu4PtBin [] = {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56, 66, 78, 92, 108, 200}; 
    double Ele1PtBin [] = {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56, 66, 78, 92, 108, 200}; 
    double Ele2PtBin [] = {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56, 66, 78, 92, 108, 200}; 
    double Mu1Mu2PtBin [] = {0, 5, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 75, 80, 86, 92, 99, 107, 116, 126, 137, 149, 162, 176, 200};
    //double tauPtBin [] = {8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 66, 74, 84, 96, 110, 126, 144, 164, 186, 210};
    double tauPtBin [] = {8, 10, 15, 20, 25, 30, 35, 40, 45, 50}; 
   double Mu3TauPtBin [] = {0, 5, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 75, 80, 86, 92, 99, 107, 116, 126, 137, 149, 162, 176, 200};
   double Mu1Mu2Mu3TauPtBin [] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250}; 
   //{0, 5, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 75, 80, 86, 92, 99, 107, 116, 126, 137, 149, 162, 176, 200};



    double genMuPtBin [] = {3, 18, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93};
    //, 97, 101, 105, 110, 115, 120, 126, 133, 141, 150, 160, 171, 183, 196, 210, 225, 241, 258, 276, 300};
    double genMu1PtBin [] = {3, 5,10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150};
    //, 97, 101, 105, 110, 115, 120, 126, 133, 141, 150, 160, 171, 183, 196, 210, 225, 241, 258, 276, 300};
    double genMu2PtBin [] = {3, 5,10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150};
    //, 97, 101, 105, 110, 115, 120, 126, 133, 141, 150, 160, 171, 183, 196, 210, 225, 241, 258, 276, 300};
    double genMu3PtBin []= {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56};
    //, 66, 78, 92, 108, 200};
    double genMu4PtBin []= {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56, 66, 78, 92, 108, 200};
    double genEle1PtBin [] = {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56, 66, 78, 92, 108, 600};
    double genEle2PtBin [] = {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56, 66, 78, 92, 108, 600};
    double genTauEle1PtBin [] = {8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 66, 74, 84, 96, 110, 126, 144, 164, 186, 210};
    double genTauMu1PtBin [] = {8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 66, 74, 84, 96, 110, 126, 144, 164, 186, 210};
    double genTauMu2PtBin [] = {8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 66, 74, 84, 96, 110, 126, 144, 164, 186, 210};



    int NBinsMu1Pt = sizeof(Mu1PtBin)/sizeof(Mu1PtBin[0])-1;
    int NBinsMu2Pt = sizeof(Mu2PtBin)/sizeof(Mu2PtBin[0])-1;
    int NBinsMu3Pt = sizeof(Mu3PtBin)/sizeof(Mu3PtBin[0])-1;
    int NBinsMu4Pt = sizeof(Mu4PtBin)/sizeof(Mu4PtBin[0])-1;
    int NBinsEle1Pt = sizeof(Ele1PtBin)/sizeof(Ele1PtBin[0])-1;
    int NBinsEle2Pt = sizeof(Ele2PtBin)/sizeof(Ele2PtBin[0])-1;
    int NBinsMu1Mu2Pt = sizeof(Mu1Mu2PtBin)/sizeof(Mu1Mu2PtBin[0])-1;
    int NBinsTauPt = sizeof(tauPtBin)/sizeof(tauPtBin[0])-1;
    int NBinsMu3TauPt = sizeof(Mu3TauPtBin)/sizeof(Mu3TauPtBin[0])-1;
    int NBinsMu1Mu2Mu3TauPt = sizeof(Mu1Mu2Mu3TauPtBin)/sizeof(Mu1Mu2Mu3TauPtBin[0])-1;

    int NBinsgenMuPt = sizeof(genMuPtBin)/sizeof(genMuPtBin[0])-1;
    int NBinsgenMu1Pt = sizeof(genMu1PtBin)/sizeof(genMu1PtBin[0])-1;
    int NBinsgenMu2Pt = sizeof(genMu2PtBin)/sizeof(genMu2PtBin[0])-1;
    int NBinsgenMu3Pt = sizeof(genMu3PtBin)/sizeof(genMu3PtBin[0])-1;
    int NBinsgenMu4Pt = sizeof(genMu4PtBin)/sizeof(genMu4PtBin[0])-1;
    int NBinsgenEle1Pt = sizeof(genEle1PtBin)/sizeof(genEle1PtBin[0])-1;
    int NBinsgenEle2Pt= sizeof(genEle2PtBin)/sizeof(genEle2PtBin[0])-1;
    int NBinsgenTauEle1Pt = sizeof(genTauEle1PtBin)/sizeof(genTauEle1PtBin[0])-1;
    int NBinsgenTauMu1Pt = sizeof(genTauMu1PtBin)/sizeof(genTauMu1PtBin[0])-1;
    int NBinsgenTauMu2Pt = sizeof(genTauMu2PtBin)/sizeof(genTauMu2PtBin[0])-1;

    dRMu1Mu2 = newTH1D("dRMu1Mu2", "#Delta R(#mu_{1}#mu_{2})", 40, 0, 1.5);
    dRMu3Mu4 = newTH1D("dRMu3Mu4", "#Delta R(#mu_{3}#mu_{4})", 25, 0, 1.0);
    dRMu3Ele = newTH1D("dRMu3Ele", "#Delta R(#mu_{3}e)", 25, 0, 1.0);
    dREleEle = newTH1D("dREleEle", "#Delta R(ee)", 25, 0, 1.0);
    dRMu3Tau = newTH1D("dRMu3Tau", "#Delta R(#mu_{3}#tau)", 25, 0, 1.0);
    dREleTau = newTH1D("dREleTau", "#Delta R(e#tau)", 25, 0, 1.0);
    dRTauTau = newTH1D("dRTauTau", "#Delta R(#tau#tau)", 25, 0, 1.0);

    invMassMu1Mu2 = newTH1D("invMassMu1Mu2", "M(#mu_{1}#mu_{2})[GeV]", 100, 0, 100);
    invMassMu3Mu4 = newTH1D("invMassMu3Mu4", "M(#mu_{3}#mu_{4})[GeV]", 100, 0, 100);
    invMassMu3Ele = newTH1D("invMassMu3Ele", "M(#mu_{3}e)[GeV]", 100, 0, 100);
    invMassEleEle = newTH1D("invMassEleEle", "M(ee)[GeV]", 100, 0, 100);
    invMassMu3Tau = newTH1D("invMassMu3Tau", "M(#mu_{3}#tau)[GeV]", 100, 0, 100);
    invMassEleTau = newTH1D("invMassEleTau", "M(e#tau)[GeV]", 100, 0, 100);
    invMassTauTau = newTH1D("invMassTauTau", "M(#tau#tau)[GeV]", 100, 0, 100);

    invMassMuMuTauMuTauMu = newTH1D("invMassMuMuTauMuTauMu", "M(#mu#mu#mu#mu)[GeV]", 100, 20, 300);
    invMassMuMuTauMuTauEle = newTH1D("invMassMuMuTauMuTauEle", "M(#mu#mu#mue)[GeV]", 100, 20, 300);
    invMassMuMuTauEleTauEle = newTH1D("invMassMuMuTauEleTauEle", "M(#mu#muee)[GeV]", 100, 20, 300);
    invMassMuMuTauMuTauHad = newTH1D("invMassMuMuTauMuTauHad", "M(#mu#mu#mu#tau)[GeV]", 100, 20, 300);
    invMassMuMuTauEleTauHad = newTH1D("invMassMuMuTauEleTauHad", "M(#mu#mue#tau)[GeV]", 100, 20, 300);
    invMassMuMuTauHadTauHad = newTH1D("invMassMuMuTauHadTauHad", "M(#mu#mu#tau#tau)[GeV]", 100, 20, 300);

    ptMu1Mu2 = newTH1D("ptMu1Mu2", "p_{T}(#mu_{1}#mu_{2})[GeV]", NBinsMu1Mu2Pt, Mu1Mu2PtBin);
    ptMu3Mu4 = newTH1D("ptMu3Mu4", "p_{T}(#mu_{3}#mu_{4})[GeV]", NBinsMu3TauPt, Mu3TauPtBin);
    ptMu3Ele = newTH1D("ptMu3Ele", "p_{T}(#mu_{3}e)[GeV]", NBinsMu3TauPt, Mu3TauPtBin);
    ptEleEle = newTH1D("ptEleEle", "p_{T}(ee)[GeV]", NBinsMu3TauPt, Mu3TauPtBin);
    ptMu3Tau = newTH1D("ptMu3Tau", "p_{T}(#mu_{3}#tau)[GeV]", NBinsMu3TauPt, Mu3TauPtBin);
    ptEleTau = newTH1D("ptEleTau", "p_{T}(e#tau)[GeV]", NBinsMu3TauPt, Mu3TauPtBin);
    ptTauTau = newTH1D("ptTauTau", "p_{T}(#tau#tau)[GeV]", NBinsMu3TauPt, Mu3TauPtBin);

    ptMuMuTauMuTauMu = newTH1D("ptMuMuTauMuTauMu", "p_{T}(#mu#mu#mu#mu)[GeV]", NBinsMu1Mu2Mu3TauPt, Mu1Mu2Mu3TauPtBin);
    ptMuMuTauMuTauEle = newTH1D("ptMuMuTauMuTauEle", "p_{T}(#mu#mu#mue)[GeV]", NBinsMu1Mu2Mu3TauPt, Mu1Mu2Mu3TauPtBin);
    ptMuMuTauEleTauEle = newTH1D("ptMuMuTauEleTauEle", "p_{T}(#mu#muee)[GeV]", NBinsMu1Mu2Mu3TauPt, Mu1Mu2Mu3TauPtBin);
    ptMuMuTauMuTauHad = newTH1D("ptMuMuTauMuTauHad", "p_{T}(#mu#mu#mu#tau)[GeV]", NBinsMu1Mu2Mu3TauPt, Mu1Mu2Mu3TauPtBin);
    ptMuMuTauEleTauHad = newTH1D("ptMuMuTauEleTauHad", "p_{T}(#mu#mue#tau)[GeV]", NBinsMu1Mu2Mu3TauPt, Mu1Mu2Mu3TauPtBin);
    ptMuMuTauHadTauHad = newTH1D("ptMuMuTauHadTauHad", "p_{T}(#mu#mu#tau#tau)[GeV]", NBinsMu1Mu2Mu3TauPt, Mu1Mu2Mu3TauPtBin);

    mu1Iso = newTH1D("mu1Iso", "#mu_{1}^{iso}", 20, 0, 0.25);
    mu2Iso = newTH1D("mu2Iso", "#mu_{2}^{iso}", 20, 0, 0.25);
    mu3Iso = newTH1D("mu3Iso", "#mu_{3}^{iso}", 20, 0, 20);
    mu4Iso = newTH1D("mu4Iso", "#mu_{4}^{iso}", 20, 0, 20);
    ele1Iso = newTH1D("ele1Iso", "e_{1}^{iso}", 20, 0, 20);
    ele2Iso = newTH1D("ele2Iso", "e_{2}^{iso}", 20, 0, 20);
    tauIsoMVA = newTH1D("tauIsoMVA", "#tau^{iso}", 20, -1, 1);
    tau2IsoMVA = newTH1D("tau2IsoMVA", "#tau_{2}^{iso}", 20, -1, 1);

    /*
    unMatchedEleIso = newTH1D("unMatchedEleIso", "e^{iso}(unMatch)", 200, 0, 200);
    unMatchedTauIsoMVA = newTH1D("unMatchedTauIsoMVA", "#tau^{iso}(unMatch)", 20, -1, 1);
    unMatchedTauDecayMode = newTH1D("unMatchedTauDecayMode", "DM(unMatch #tau)", 11, 0, 11);

    unMatchedMuPt = newTH1D("unMatchedMuPt", "p_{T}(#mu_{unMatched}) [GeV]", 50, 3, 203);
    unMatchedMuEta = newTH1D("unMatchedMuEta", "#eta(#mu_{unMatched})", 20, -2.1, 2.1);
    unMatchedMuPhi = newTH1D("unMatchedMuPhi", "#phi(#mu_{unMatched})", 20, -2.5, 2.5);

    unMatchedElePt = newTH1D("unMatchedElePt", "p_{T}(e_{unMatched}) [GeV]", 50, 3, 203);
    unMatchedEleEta = newTH1D("unMatchedEleEta", "#eta(e_{unMatched})", 20, -2.1, 2.1);
    unMatchedElePhi = newTH1D("unMatchedElePhi", "#phi(e_{unMatched})", 20, -2.5, 2.5);

    unMatchedTauPt = newTH1D("unMatchedTauPt", "p_{T}(#tau_{unMatched}) [GeV]", 50, 3, 203);
    unMatchedTauEta = newTH1D("unMatchedTauEta", "#eta(#tau_{unMatched})", 20, -2.1, 2.1);
    unMatchedTauPhi = newTH1D("unMatchedTauPhi", "#phi(#tau_{unMatched})", 20, -2.5, 2.5);
    */
    mu1Pt = newTH1D("mu1Pt", "p_{T}(#mu_{1}) [GeV]", NBinsMu1Pt, Mu1PtBin);
    mu1Eta = newTH1D("mu1Eta", "#eta(#mu_{1})", 20, -2.1, 2.1);
    mu1Phi = newTH1D("mu1Phi", "#phi(#mu_{1})", 20, -3.1, 3.1);

    mu2Pt = newTH1D("mu2Pt", "p_{T}(#mu_{2}) [GeV]", NBinsMu2Pt, Mu2PtBin);
    mu2Eta = newTH1D("mu2Eta", "#eta(#mu_{2})", 20, -2.1, 2.1);
    mu2Phi = newTH1D("mu2Phi", "#phi(#mu_{2})", 20, -3.1, 3.1);

    mu3Pt = newTH1D("mu3Pt", "p_{T}(#mu_{3}) [GeV]", NBinsMu3Pt, Mu3PtBin);
    mu3Eta = newTH1D("mu3Eta", "#eta(#mu_{3})", 20, -2.1, 2.1);
    mu3Phi = newTH1D("mu3Phi", "#phi(#mu_{3})", 20, -3.1, 3.1);

    mu4Pt = newTH1D("mu4Pt", "p_{T}(#mu_{4}) [GeV]", NBinsMu4Pt, Mu4PtBin);
    mu4Eta = newTH1D("mu4Eta", "#eta(#mu_{4})", 20, -2.1, 2.1);
    mu4Phi = newTH1D("mu4Phi", "#phi(#mu_{4})", 20, -3.1, 3.1);

    ele1Pt = newTH1D("ele1Pt", "p_{T}(e_{1}) [GeV]", NBinsEle1Pt, Ele1PtBin);
    ele1Eta = newTH1D("ele1Eta", "#eta(e_{1})", 20, -2.1, 2.1);
    ele1Phi = newTH1D("ele1Phi", "#phi(e_{1})", 20, -3.1, 3.1);

    ele2Pt = newTH1D("ele2Pt", "p_{T}(e_{2}) [GeV]", NBinsEle2Pt, Ele2PtBin);
    ele2Eta = newTH1D("ele2Eta", "#eta(e_{2})", 20, -2.1, 2.1);
    ele2Phi = newTH1D("ele2Phi", "#phi(e_{2})", 20, -3.1, 3.1);

    tauPt = newTH1D("tauPt", "p_{T}(#tau) [GeV]", NBinsTauPt, tauPtBin);
    tauEta = newTH1D("tauEta", "#eta(#tau)", 20, -2.1, 2.1);
    tauPhi = newTH1D("tauPhi", "#phi(#tau)", 20, -3.1, 3.1);
    tauMass = newTH1D("tauMass", "M(#tau) [GeV]", 10, 0, 5);
    tauDecayMode = newTH1D("tauDecayMode", "DecayMode(#tau)", 11, 0, 11);

    tau2Pt = newTH1D("tau2Pt", "p_{T}(#tau_{2}) [GeV]", NBinsTauPt, tauPtBin);
    tau2Eta = newTH1D("tau2Eta", "#eta(#tau_{2})", 20, -2.1, 2.1);
    tau2Phi = newTH1D("tau2Phi", "#phi(#tau_{2})", 20, -3.1, 3.1);
    tau2Mass = newTH1D("tau2Mass", "M(#tau_{2}) [GeV]", 10, 0, 5);
    tau2DecayMode = newTH1D("tau2DecayMode", "DecayMode(#tau_{2})", 11, 0, 11);

    dRMu1Mu3 = newTH1D("dRMu1Mu3", "#Delta R(#mu_{1}#mu_{3})", 25, 0, 5);
    dRMu1Mu4 = newTH1D("dRMu1Mu4", "#Delta R(#mu_{1}#mu_{4})", 25, 0, 5);
    dRMu1Ele1 = newTH1D("dRMu1Ele1", "#Delta R(#mu_{1}e_{1})", 25, 0, 5);
    dRMu1Ele2 = newTH1D("dRMu1Ele2", "#Delta R(#mu_{1}e_{2})", 25, 0, 5);
    dRMu1Tau = newTH1D("dRMu1Tau", "#Delta R(#mu_{1}#tau)", 25, 0, 5);
    dRMu1Tau2 = newTH1D("dRMu1Tau2", "#Delta R(#mu_{1}#tau_{2})", 25, 0, 5);

    dRMu2Mu3 = newTH1D("dRMu2Mu3", "#Delta R(#mu_{2}#mu_{3})", 25, 0, 5);
    dRMu2Mu4 = newTH1D("dRMu2Mu4", "#Delta R(#mu_{2}#mu_{4})", 25, 0, 5);
    dRMu2Ele1 = newTH1D("dRMu2Ele1", "#Delta R(#mu_{2}e_{1})", 25, 0, 5);
    dRMu2Ele2 = newTH1D("dRMu2Ele2", "#Delta R(#mu_{2}e_{2})", 25, 0, 5);
    dRMu2Tau = newTH1D("dRMu2Tau", "#Delta R(#mu_{2}#tau)", 25, 0, 5);
    dRMu2Tau2 = newTH1D("dRMu2Tau2", "#Delta R(#mu_{2}#tau_{2})", 25, 0, 5);

    dRInvMassMu1Mu2 = newTH2D("dRInvMassMu1Mu2", "#Delta R(#mu_{1}#mu_{2})", "M(#mu_{1}#mu_{2})[GeV]", 25, 0, 1, 100, 0, 100);
    dRInvMassMu3Mu4 = newTH2D("dRInvMassMu3Mu4", "#Delta R(#mu_{3}#mu_{4})", "M(#mu_{3}#mu_{4})[GeV]", 25, 0, 1, 100, 0, 100);
    dRInvMassMu3Ele = newTH2D("dRInvMassMu3Ele", "#Delta R(#mu_{3}e)", "M(#mu_{3}e)[GeV]", 25, 0, 1, 100, 0, 100);
    dRInvMassEleEle = newTH2D("dRInvMassEleEle", "#Delta R(ee)", "M(ee)[GeV]", 25, 0, 1, 100, 0, 6);
    dRInvMassMu3Tau = newTH2D("dRInvMassMu3Tau", "#Delta R(#mu_{3}#tau)", "M(#mu_{3}#tau)[GeV]", 25, 0, 1, 100, 0, 100);
    dRInvMassEleTau = newTH2D("dRInvMassEleTau", "#Delta R(e#tau)", "M(e#tau)[GeV]", 25, 0, 1, 100, 0, 100);
    dRInvMassTauTau = newTH2D("dRInvMassTauTau", "#Delta R(#tau#tau)", "M(#tau#tau)[GeV]", 25, 0, 1, 100, 0, 100);
    
    /*    nMatchedMuPairNMatchedMuMuPair = newTH2D("nMatchedMuPairNMatchedMuMuPair", "N(#mu_{1}#mu_{2})", "N(#mu_{3}#mu_{4})", 4, 0, 4, 4, 0, 4);
    nMatchedMuPairNMatchedMuElePair = newTH2D("nMatchedMuPairNMatchedMuElePair", "N(#mu_{1}#mu_{2})", "N(#mu_{3}e)", 4, 0, 4, 4, 0, 4);
    nMatchedMuPairNMatchedEleElePair = newTH2D("nMatchedMuPairNMatchedEleElePair", "N(#mu_{1}#mu_{2})", "N(ee)", 4, 0, 4, 4, 0, 4);
    nMatchedMuPairNMatchedMuTauPair = newTH2D("nMatchedMuPairNMatchedMuTauPair", "N(#mu_{1}#mu_{2})", "N(#mu_{3}#tau)", 4, 0, 4, 4, 0, 4);
    nMatchedMuPairNMatchedEleTauPair = newTH2D("nMatchedMuPairNMatchedEleTauPair", "N(#mu_{1}#mu_{2})", "N(e#tau)", 4, 0, 4, 4, 0, 4);
    nMatchedMuPairNMatchedTauTauPair = newTH2D("nMatchedMuPairNMatchedTauTauPair", "N(#mu_{1}#mu_{2})", "N(#tau#tau)", 4, 0, 4, 4, 0, 4);
    */
   
    invMassgenMu1genMu2 = newTH1D("invMassgenMu1genMu2", "M(#mu_{1}#mu_{2})[GeV]", 100, 9.5, 10.5);
    invMassgenMu3genMu4 = newTH1D("invMassgenMu3genMu4", "M(#mu_{3}#mu_{4})[GeV]", 100, 9.5, 10.5);
    invMassgenMu3genEle = newTH1D("invMassgenMu3genEle", "M(#mu_{3}e)[GeV]", 100, 4.5, 5.5);
    invMassgenElegenEle = newTH1D("invMassgenElegenEle", "M(ee)[GeV]", 100, 4.5, 5.5);
    invMassgenMu3genTau = newTH1D("invMassgenMu3genTau", "M(#mu_{3}#tau)[GeV]", 100, 0, 11);
    invMassgenElegenTau = newTH1D("invMassgenElegenTau", "M(e#tau)[GeV]", 100, 4.5, 5.5);
    invMassgenTaugenTau = newTH1D("invMassgenTaugenTau", "M(#tau#tau)[GeV]", 100, 9.5, 10.5);
    invMassgenTauElegenTauEle = newTH1D("invMassgenTauElegenTauEle", "M(#tau_{e}#tau_{e}",100,0,10);
    invMassgenTauMugenTauMu = newTH1D("invMassgenTauMugenTauMu", "M(#tau#tau)[GeV]", 100, 0, 100);

    invMassgenMuMuTauMuTauMu = newTH1D("invMassgenMuMuTauMuTauMu", "M(#mu#mu#mu#mu)[GeV]", 100, 20, 300);
    invMassgenMuMuTauMuTauEle = newTH1D("invMassgenMuMuTauMuTauEle", "M(#mu#mu#mue)[GeV]", 100, 20, 300);
    invMassgenMuMuTauEleTauEle = newTH1D("invMassgenMuMuTauEleTauEle", "M(#mu#muee)[GeV]", 100, 20, 300);
    invMassgenMuMuTauMuTauHad = newTH1D("invMassgenMuMuTauMuTauHad", "M(#mu#mu#mu#tau)[GeV]", 100, 20, 300);
    invMassgenMuMuTauEleTauHad = newTH1D("invMassgenMuMuTauEleTauHad", "M(#mu#mue#tau)[GeV]", 100, 20, 300);
    invMassgenMuMuTauHadTauHad = newTH1D("invMassgenMuMuTauHadTauHad", "M(#mu#mu#tau#tau)[GeV]", 100, 20, 300);

    ptgenMu1genMu2 = newTH1D("ptgenMu1genMu2", "p_{T}(#mu_{1}#mu_{2})[GeV]", NBinsMu1Mu2Pt, Mu1Mu2PtBin);
    ptgenMu3genMu4 = newTH1D("ptgenMu3genMu4", "p_{T}(#mu_{3}#mu_{4})[GeV]", NBinsMu3TauPt, Mu3TauPtBin);
    ptgenMu3genEle = newTH1D("ptgenMu3genEle", "p_{T}(#mu_{3}e)[GeV]", NBinsMu3TauPt, Mu3TauPtBin);
    ptgenElegenEle = newTH1D("ptgenElegenEle", "p_{T}(ee)[GeV]", NBinsMu3TauPt, Mu3TauPtBin);
    ptgenMu3genTau = newTH1D("ptgenMu3genTau", "p_{T}(#mu_{3}#tau)[GeV]", NBinsMu3TauPt, Mu3TauPtBin);
    ptgenElegenTau = newTH1D("ptgenElegenTau", "p_{T}(e#tau)[GeV]", NBinsMu3TauPt, Mu3TauPtBin);
    ptgenTaugenTau = newTH1D("ptgenTaugenTau", "p_{T}(#tau#tau)[GeV]", NBinsMu3TauPt, Mu3TauPtBin);

    ptgenMuMuTauMuTauMu = newTH1D("ptgenMuMuTauMuTauMu", "p_{T}(#mu#mu#mu#mu)[GeV]", NBinsMu1Mu2Mu3TauPt, Mu1Mu2Mu3TauPtBin);
    ptgenMuMuTauMuTauEle = newTH1D("ptgenMuMuTauMuTauEle", "p_{T}(#mu#mu#mue)[GeV]", NBinsMu1Mu2Mu3TauPt, Mu1Mu2Mu3TauPtBin);
    ptgenMuMuTauEleTauEle = newTH1D("ptgenMuMuTauEleTauEle", "p_{T}(#mu#muee)[GeV]", NBinsMu1Mu2Mu3TauPt, Mu1Mu2Mu3TauPtBin);
    ptgenMuMuTauMuTauHad = newTH1D("ptgenMuMuTauMuTauHad", "p_{T}(#mu#mu#tau#tau)[GeV]", NBinsMu1Mu2Mu3TauPt, Mu1Mu2Mu3TauPtBin);
    ptgenMuMuTauEleTauHad = newTH1D("ptgenMuMuTauEleTauHad", "p_{T}(#mu#mue#tau)[GeV]", NBinsMu1Mu2Mu3TauPt, Mu1Mu2Mu3TauPtBin);
    ptgenMuMuTauHadTauHad = newTH1D("ptgenMuMuTauHadTauHad", "p_{T}(#mu#mu#tau#tau)[GeV]", NBinsMu1Mu2Mu3TauPt, Mu1Mu2Mu3TauPtBin);


 genmuPt = newTH1D("genmuPt", "p_{T}(#mu) [GeV]", NBinsgenMuPt, genMuPtBin);

    genmu1Eta = newTH1D("genmu1Eta", "#eta(#mu_{1})", 20, -2.5, 2.5);
    genmu1Pt = newTH1D("genmu1Pt", "p_{T}(#mu_{1}) [GeV]", NBinsMu1Pt, Mu1PtBin);
    genmu1Phi = newTH1D("genmu1Phi", "#phi(#mu_{1})", 20, -3.1, 3.1);
    genmu1Mass = newTH1D("genmu1Mass", "mass(#mu_{1})",20, 0.08, 0.15);

    genmu2Pt = newTH1D("genmu2Pt", "p_{T}(#mu_{2}) [GeV]", NBinsMu2Pt, Mu2PtBin);
    genmu2Eta = newTH1D("genmu2Eta", "#eta(#mu_{2})", 20, -2.5, 2.5);
    genmu2Phi = newTH1D("genmu2Phi", "#phi(#mu_{2})", 20, -3.1, 3.1);
    genmu2Mass = newTH1D("genmu2Mass", "mass(#mu_{2})",20, 0.08, 0.15);

    genmu3Pt = newTH1D("genmu3Pt", "p_{T}(#mu_{3}) [GeV]", NBinsMu3Pt, Mu3PtBin);
    genmu3Eta = newTH1D("genmu3Eta", "#eta(#mu_{3})", 20, -2.5, 2.5);
    genmu3Phi = newTH1D("genmu3Phi", "#phi(#mu_{3})", 20, -3.1, 3.1);
    genmu3Mass = newTH1D("genmu3Mass", "mass(#mu_{3})", 20, 0.08, 0.15);

    genmu4Pt = newTH1D("genmu4Pt", "p_{T}(#mu_{4}) [GeV]", NBinsMu4Pt, Mu4PtBin);
    genmu4Eta = newTH1D("genmu4Eta", "#eta(#mu_{4})", 20, -2.5, 2.5);
    genmu4Phi = newTH1D("genmu4Phi", "#phi(#mu_{4})", 20, -3.1, 3.1);
    genmu4Mass  = newTH1D("genmu4Mass", "mass(#mu_{4})", 20, 0.08, 0.15);

    genele1Eta = newTH1D("genele1Eta", "#eta(e_{1})", 20, -2.5, 2.5);
    genele1Pt = newTH1D("genele1Pt", "p_{T}(e_{1}) [GeV]", NBinsgenEle1Pt,genEle1PtBin);
    genele1Phi = newTH1D("genele1Phi", "#phi(e_{1})", 20, -3.1, 3.1);
    genele1Mass = newTH1D("genele1Mass", "mass(e_{1})",20, 0.00, 0.7);

    genele2Eta = newTH1D("genele2Eta", "#eta(e_{2})", 20, -2.5, 2.5);
    genele2Pt = newTH1D("genele2Pt", "p_{T}(e_{2}) [GeV]", NBinsgenEle2Pt,genEle2PtBin);
    genele2Phi = newTH1D("genele2Phi", "#phi(e_{2})", 20, -3.1, 3.1);
    genele2Mass = newTH1D("genele2Mass", "mass(e_{2})",20, 0.00, 0.7);


    recomatchedtauPt = newTH1D("recomatchedtauPt","p_{T} [GeV]", NBinsTauPt,tauPtBin);
    gentauPt = newTH1D("gentauPt", "p_{T}(#tau) [GeV]", NBinsTauPt, tauPtBin);
    gentauEta = newTH1D("gentauEta", "#eta(#tau)", 20, -2.5, 2.5);
    gentauPhi = newTH1D("gentauPhi", "#phi(#tau)", 20, -3.1, 3.1);
    gentauMass = newTH1D("gentauMass", "M(#tau) [GeV]", 10, 0, 5);

    gentauMu1Pt = newTH1D("gentauMu1Pt", "p_{T}(#tau) [GeV]", NBinsTauPt, tauPtBin);
    gentauMu1Eta = newTH1D("gentauMu1Eta", "#eta(#tau)", 20, -2.5, 2.5);
    gentauMu1Phi = newTH1D("gentauMu1Phi", "#phi(#tau)", 20, -3.1, 3.1);
    gentauMu1Mass = newTH1D("gentauMu1Mass", "M(tau) [GeV]", 10, 0, 5);

    gentauMu2Pt = newTH1D("gentauMu2Pt", "p_{T}(#tau) [GeV]", NBinsTauPt, tauPtBin);
    gentauMu2Eta = newTH1D("gentauMu2Eta", "#eta(#tau)", 20, -2.5, 2.5);
    gentauMu2Phi = newTH1D("gentauMu2Phi", "#phi(#tau)", 20, -3.1, 3.1);
    gentauMu2Mass = newTH1D("gentauMu2Mass", "M(tau) [GeV]", 10, 0, 5);

    gentau2Pt = newTH1D("gentau2Pt", "p_{T}(gem#tau_{2}) [GeV]", NBinsTauPt, tauPtBin);
    gentau2Eta = newTH1D("gentau2Eta", "#eta(#tau_{2})", 20, -2.5, 2.5);
    gentau2Phi = newTH1D("gentau2Phi", "#phi(#tau_{2})", 20, -3.1, 3.1);
    gentau2Mass = newTH1D("gentau2Mass", "M(#tau_{2}) [GeV]", 10, 0, 5);

    gentauEle1Pt = newTH1D("gentauEle1Pt", "p_{T}(#tau_{e}) [GeV]", NBinsgenTauEle1Pt, genTauEle1PtBin);
    gentauEle1Eta = newTH1D("gentauEle1Eta", "#eta(#tau_{e})", 20, -2.5, 2.5);
    gentauEle1Phi = newTH1D("gentauEle1Phi", "#phi(#tau_{e})", 20, -3.1, 3.1);
    gentauEle1Mass = newTH1D("gentauEle1Mass", "M(tau_{e}) [GeV]", 10, 0, 5);

    gentauEle2Pt = newTH1D("gentauEle2Pt", "p_{T}(#tau_{e2}) [GeV]", NBinsgenTauEle1Pt, genTauEle1PtBin);
    gentauEle2Eta = newTH1D("gentauEle2Eta", "#eta(#tau_{e2})", 20, -2.5, 2.5);
    gentauEle2Phi = newTH1D("gentauEle2Phi", "#phi(#tau_{e2})", 20, -3.1, 3.1);
    gentauEle2Mass = newTH1D("gentauEle2Mass", "M(tau_{e2}) [GeV]", 10, 0, 5);

    //gentauEle2Pt = newTH1D("gentauEle2Pt", "p_{T}(#mu_{1}) [GeV]", NBinsgenTauEle1Pt, genTauEle1PtBin);
    //gentauMu2Pt = newTH1D("gentauMu2Pt", "p_{T}(#mu_{1}) [GeV]", NBinsgenTauMu2Pt, genTauMu2PtBin);

    dRgenMu1genMu2 = newTH1D("dRgenMu1genMu2", "#Delta R(#mu_{1}#mu_{2})", 40, 0, 1.5);
    dRgenMu3genMu4 = newTH1D("dRgenMu3genMu4", "#Delta R(#mu_{3}#mu_{4})", 25, 0, 1.0);
    dRgenMu3genEle = newTH1D("dRgenMu3genEle", "#Delta R(gen#mu_{3}gen e)", 25, 0, 1.0);
    //    dRgenElegenEle = newTH1D("dRgenEle1genEle2", "#Delta_R(e_{1}e_{2}", 25, 0, 1.0);
    dRgenMu3genTau = newTH1D("dRgenMu3genTau", "#Delta R(gen#mu_{3}gen#tau)", 25, 0, 1.0);
    dRgenElegenTau = newTH1D("dRgenElegenTau", "#Delta R(genegen#tau)", 25, 0, 1.0);
    dRgenTaugenTau = newTH1D("dRgenTaugenTau", "#Delta R(gen#taugen#tau)", 25, 0, 1.0);

    dRgenMu1genMu3 = newTH1D("dRgenMu1genMu3", "#Delta R(#mu_{1}#mu_{3})", 25, 0, 5);
    dRgenMu1genMu4 = newTH1D("dRgenMu1genMu4", "#Delta R(#mu_{1}#mu_{4})", 25, 0, 5);
    dRgenMu1genEle1 = newTH1D("dRgenMu1genEle1", "#Delta R(#mu_{1}e_{1})", 25, 0, 5);
    dRgenMu1genEle2 = newTH1D("dRgenMu1genEle2", "#Delta R(#mu_{1}e_{2})", 25, 0, 5);
    dRgenMu1genTau = newTH1D("dRgenMu1genTau", "#Delta R(#mu_{1}#tau)", 25, 0, 5);
    dRgenMu3genTauMu1 = newTH1D("dRgenMu3genTauMu1", "#Delta R(#mu_{3}#tau_{1})", 25, 0, 5);
    dRgenMu4genTauMu2 = newTH1D("dRgenMu4genTauMu2", "#Delta R(#mu_{4}#tau_{2})", 25, 0, 5);

    dRgenMu2genMu3 = newTH1D("dRgenMu2genMu3", "#Delta R(#mu_{2}#mu_{3})", 25, 0, 5);
    dRgenMu2genMu4 = newTH1D("dRgenMu2genMu4", "#Delta R(#mu_{2}#mu_{4})", 25, 0, 5);
    dRgenMu2genEle1 = newTH1D("dRgenMu2genEle1", "#Delta R(#mu_{2}e_{1})", 25, 0, 5);
    dRgenMu2genEle2 = newTH1D("dRgenMu2genEle2", "#Delta R(#mu_{2}e_{2})", 25, 0, 5);
    dRgenMu2genTau = newTH1D("dRgenMu2genTau", "#Delta R(#mu_{2}#tau)", 25, 0, 5);
    dRgenMu2genTau2 = newTH1D("dRgenMu2genTau2", "#Delta R(#mu_{2}#tau_{2})", 25, 0, 5);

    dRgenMu1genTauEle1= newTH1D("dRgenMu1genTauEle1","Delta R(#mu_{1}#tau_{e1}",25, 0, 5);
    dRgenMu1genTauEle2 = newTH1D("dRgenMu1genTauEle2","Delta R(#mu_{1}#tau_{e2}",25, 0, 5);
    dRgenMu2genTauEle1 = newTH1D("dRgenMu2genTauEle1","Delta R(#mu_{2}#tau_{e1}",25, 0, 5);
    dRgenMu2genTauEle2 = newTH1D("dRgenMu2genTauEle2","Delta R(#mu_{2}#tau_{e2}",25, 0, 5);
    dRgenEle1genEle2 = newTH1D("dRgenEle1genEle2", "Delta R(e_{1}e_{2})",25, 0, 1);



       //  dRgenTauEle1genTauEle2 = newTH1D("dRgenTauEle1genTauEle2", "#Delta_R(#tau#tau)", 25, 0, 1.0);
    //dRgenTauMu1genTauMu2 =  newTH1D("dRgenTauMu1genTauMu2", "#Delta_R(#tau#tau)", 25, 0, 1.0);
  
    //    dRInvMassgenEle1genEle2 = newTH2D("dRInvMassgenEle1genEle2", "#Delta R(ee)", "M(ee)[GeV]", 25, 0, 0.5, 100, 0, 6);
    dRInvMassgenMu1genMu2 = newTH2D("dRInvMassgenMu1genMu2", "#Delta R(#mu#mu)", "M(#mu_{1}#mu_{2})[GeV]", 25, 0, 0.5, 100, 0, 6);
    dRInvMassgenMu3genMu4 = newTH2D("dRInvMassgenMu3genMu4", "#Delta R(#mu#mu)", "M(#mu_{3}#mu_{4})[GeV]", 25, 0, 0.5, 100, 0, 6); 
    dRInvMassgenMu3genEle = newTH2D("dRInvMassgenMu3genEle", "#Delta R(#mu_{3}e)", "M(#mu_{3}e)[GeV]", 25, 0, 1, 100, 0, 100);
    dRInvMassgenElegenEle = newTH2D("dRInvMassgenElegenEle", "#Delta R(ee)", "M(ee)[GeV]", 25, 0, 1, 100, 0, 6);
    dRInvMassgenMu3genTau = newTH2D("dRInvMassgenMu3genTau", "#Delta R(#mu_{3}#tau)", "M(#mu_{3}#tau)[GeV]", 25, 0, 1, 100, 0, 100);
    dRInvMassgenElegenTau = newTH2D("dRInvMassgenElegenTau", "#Delta R(e#tau)", "M(e#tau)[GeV]", 25, 0, 1, 100, 0, 100);
    dRInvMassgenTaugenTau = newTH2D("dRInvMassgenTaugenTau", "#Delta R(#tau#tau)", "M(#tau#tau)[GeV]", 25, 0, 1, 100, 0, 100);

   
      //  dRInvMassgenTauEle1genTauEle2 = newTH2D("dRInvMassgenTauEle1genTauEle2", "#Delta R(#tau#tau)", "M(#tau#tau)[GeV]", 25, 0, 0.5, 100, 0, 6);
    //dRInvMassgenTauMu1genTauMu2 = newTH2D("dRInvMassgenTauMu1genTauMu2", "#Delta R(#tau#tau)", "M(#tau#tau)[GeV]", 25, 0, 0.5, 100, 0, 6);

    mu1PtVSGenMu1Pt = newTH2D("mu1PtVSGenMu1Pt", "p_{T}^{rec}(#mu_{1})[GeV]", "p_{T}^{gen}(#mu_{1})[GeV]", NBinsMu1Pt, Mu1PtBin, NBinsMu1Pt, Mu1PtBin);
    mu1EtaVSGenMu1Eta = newTH2D("mu1EtaVSGenMu1Eta", "#eta^{rec}(#mu_{1})", "#eta^{gen}(#mu_{1})", 20, -2.1, 2.1, 20, -2.1, 2.1);
    mu1PhiVSGenMu1Phi = newTH2D("mu1PhiVSGenMu1Phi", "#phi^{rec}(#mu_{1})", "#phi^{gen}(#mu_{1})", 20, -3.1, 3.1, 20, -3.1, 3.1);

    mu2PtVSGenMu2Pt = newTH2D("mu2PtVSGenMu2Pt", "p_{T}^{rec}(#mu_{2})[GeV]", "p_{T}^{gen}(#mu_{2})[GeV]", NBinsMu2Pt, Mu2PtBin, NBinsMu2Pt, Mu2PtBin);
    mu2EtaVSGenMu2Eta = newTH2D("mu2EtaVSGenMu2Eta", "#eta^{rec}(#mu_{2})", "#eta^{gen}(#mu_{2})", 20, -2.1, 2.1, 20, -2.1, 2.1);
    mu2PhiVSGenMu2Phi = newTH2D("mu2PhiVSGenMu2Phi", "#phi^{rec}(#mu_{2})", "#phi^{gen}(#mu_{2})", 20, -3.1, 3.1, 20, -3.1, 3.1);

    mu3PtVSGenMu3Pt = newTH2D("mu3PtVSGenMu3Pt", "p_{T}^{rec}(#mu_{3})[GeV]", "p_{T}^{gen}(#mu_{3})[GeV]", NBinsMu3Pt, Mu3PtBin, NBinsMu3Pt, Mu3PtBin);
    mu3PtVSGenTauMu1Pt = newTH2D("mu3PtVSGenTauMu1Pt", "p_{T}^{rec}(#mu_{3})[GeV]", "p_{T}^{gen}(#tau_{#mu})[GeV]", NBinsMu3Pt, Mu3PtBin, NBinsMu3Pt, Mu3PtBin);
    mu3EtaVSGenMu3Eta = newTH2D("mu3EtaVSGenMu3Eta", "#eta^{rec}(#mu_{3})", "#eta^{gen}(#mu_{3})", 20, -2.1, 2.1, 20, -2.1, 2.1);
    mu3PhiVSGenMu3Phi = newTH2D("mu3PhiVSGenMu3Phi", "#phi^{rec}(#mu_{3})", "#phi^{gen}(#mu_{3})", 20, -3.1, 3.1, 20, -3.1, 3.1);

    mu4PtVSGenMu4Pt = newTH2D("mu4PtVSGenMu4Pt", "p_{T}^{rec}(#mu_{4})[GeV]", "p_{T}^{gen}(#mu_{4})[GeV]", NBinsMu3Pt, Mu3PtBin, NBinsMu3Pt, Mu3PtBin);
    mu4PtVSGenTauMu2Pt = newTH2D("mu4PtVSGenTauMu2Pt", "p_{T}^{rec}(#mu_{4})[GeV]", "p_{T}^{gen(#tau_{#mu})[GeV]", NBinsMu3Pt, Mu3PtBin, NBinsMu3Pt, Mu3PtBin);
    mu4EtaVSGenMu4Eta = newTH2D("mu4EtaVSGenMu4Eta", "#eta^{rec}(#mu_{4})", "#eta^{gen}(#mu_{4})", 20, -2.1, 2.1, 20, -2.1, 2.1);
    mu4PhiVSGenMu4Phi = newTH2D("mu4PhiVSGenMu4Phi", "#phi^{rec}(#mu_4})", "#phi^{gen}(#mu_{4})", 20, -3.1, 3.1, 20, -3.1, 3.1);

    ele1PtVSGenEle1Pt = newTH2D("ele1PtVSGenEle1Pt", "p_{T}^{rec}(e)[GeV]", "p_{T}^{gen}(e)[GeV]", NBinsEle1Pt, Ele1PtBin, NBinsEle1Pt, Ele1PtBin);
    ele1PtVSGenTauEle1Pt = newTH2D("ele1PtVSGenTauEle1Pt", "p_{T}^{rec}(e_)[GeV]", "p_{T}^{gen}(#tau_{e})[GeV]", NBinsEle1Pt, Ele1PtBin, NBinsEle1Pt, Ele1PtBin);
    ele1EtaVSGenEle1Eta = newTH2D("ele1EtaVSGenEle1Eta", "#eta^{rec}(e)", "#eta^{gen}(e)", 20, -2.1, 2.1, 20, -2.1, 2.1);
    ele1PhiVSGenEle1Phi = newTH2D("ele1PhiVSGenEle1Phi", "#phi^{rec}(e)", "#phi^{gen}(e)", 20, -3.1, 3.1, 20, -3.1, 3.1);

    ele2PtVSGenEle2Pt = newTH2D("ele2PtVSGenEle2Pt", "p_{T}^{rec}(e)[GeV]", "p_{T}^{gen}(e)[GeV]", NBinsEle1Pt, Ele1PtBin, NBinsEle1Pt, Ele1PtBin);
    ele2PtVSGenTauEle2Pt = newTH2D("ele2PtVSGenTauEle2Pt", "p_{T}^{rec}(e_)[GeV]", "p_{T}^{gen}(#tau_{e})[GeV]", NBinsEle1Pt, Ele1PtBin, NBinsEle1Pt, Ele1PtBin);
    ele2EtaVSGenEle2Eta = newTH2D("ele2EtaVSGenEle2Eta", "#eta^{rec}(e)", "#eta^{gen}(e)", 20, -2.1, 2.1, 20, -2.1, 2.1);
    ele2PhiVSGenEle2Phi = newTH2D("ele2PhiVSGenEle2Phi", "#phi^{rec}(e)", "#phi^{gen}(e)", 20, -3.1, 3.1, 20, -3.1, 3.1);



    tauPtVSGenTauHadPt = newTH2D("tauPtVSGenTauHadPt", "p_{T}^{rec}(#tau)[GeV]", "p_{T}^{gen}(#tau_{h})[GeV]", NBinsTauPt, tauPtBin, NBinsTauPt, tauPtBin);
    tauEtaVSGenTauHadEta = newTH2D("tauEtaVSGenTauHadEta", "#eta^{rec}(#tau)", "#eta^{gen}(#tau_{h})", 20, -2.1, 2.1, 20, -2.1, 2.1);
    tauPhiVSGenTauHadPhi = newTH2D("tauPhiVSGenTauHadPhi", "#phi^{rec}(#tau)", "#phi^{gen}(#tau_{h})", 20, -3.1, 3.1, 20, -3.1, 3.1);
    tauPtVSGenTauHadVisPt = newTH2D("tauPtVSGenTauHadVisPt", "p_{T}^{rec}(#tau)[GeV]", "p_{T}^{gen}(#tau_{h}^{vis})[GeV]", NBinsTauPt, tauPtBin, NBinsTauPt, tauPtBin);

   
}

Histomutau::~Histomutau()
{}
