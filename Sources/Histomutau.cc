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

    double Mu1PtBin [] = {3, 18, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93, 97, 101, 105, 110, 115, 120, 126, 133, 141, 150, 160, 171, 183, 196, 210, 225, 241, 258, 276, 300};
    double Mu2PtBin [] = {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 75, 80, 86, 92, 99, 107, 116, 126, 137, 149, 162, 176, 200};
    double Mu3PtBin [] = {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56, 66, 78, 92, 108, 200}; 
    double Mu4PtBin [] = {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56, 66, 78, 92, 108, 200}; 
    double Ele1PtBin [] = {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56, 66, 78, 92, 108, 200}; 
    double Ele2PtBin [] = {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56, 66, 78, 92, 108, 200}; 
    double Mu1Mu2PtBin [] = {0, 5, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 75, 80, 86, 92, 99, 107, 116, 126, 137, 149, 162, 176, 200};
    double tauPtBin [] = {8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 66, 74, 84, 96, 110, 126, 144, 164, 186, 210};
    double Mu3TauPtBin [] = {0, 5, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 75, 80, 86, 92, 99, 107, 116, 126, 137, 149, 162, 176, 200};
    double Mu1Mu2Mu3TauPtBin [] = {0, 5, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 75, 80, 86, 92, 99, 107, 116, 126, 137, 149, 162, 176, 200};



    double genMu1PtBin [] = {3, 18, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93, 97, 101, 105, 110, 115, 120, 126, 133, 141, 150, 160, 171, 183, 196, 210, 225, 241, 258, 276, 300};
    double genMu2PtBin [] = {3, 18, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93, 97, 101, 105, 110, 115, 120, 126, 133, 141, 150, 160, 171, 183, 196, 210, 225, 241, 258, 276, 300};
    //    double genEle1PtBin [] = {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56, 66, 78, 92, 108, 200};
    //double genEle2PtBin [] = {3, 10, 14, 18, 22, 26, 30, 34, 38, 42, 48, 56, 66, 78, 92, 108, 200};


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

    int NBinsgenMu1Pt = sizeof(genMu1PtBin)/sizeof(genMu1PtBin[0])-1;
    int NBinsgenMu2Pt = sizeof(genMu2PtBin)/sizeof(genMu2PtBin[0])-1;
    // int NBinsgenEle1Pt = sizeof(genEle1PtBin)/sizeof(genEle1PtBin[0])-1;
    //int NBinsgenEle2Pt= sizeof(genEle2PtBin)/sizeof(genEle2PtBin[0])-1;

    nMatchedMuPair = newTH1D("nMatchedMuPair", "N(#mu_{1}#mu_{2})", 5, 0, 5);
    nMatchedMuMuPair = newTH1D("nMatchedMuMuPair", "N(#mu_{3}#mu_{4})", 5, 0, 5);
    nMatchedMuElePair = newTH1D("nMatchedMuElePair", "N(#mu_{3}e)", 5, 0, 5);
    nMatchedEleElePair = newTH1D("nMatchedEleElePair", "N(ee)", 5, 0, 5);
    nMatchedMuTauPair = newTH1D("nMatchedMuTauPair", "N(#mu_{3}#tau)", 5, 0, 5);
    nMatchedEleTauPair = newTH1D("nMatchedEleTauPair", "N(e#tau)", 5, 0, 5);
    nMatchedTauTauPair = newTH1D("nMatchedTauTauPair", "N(#tau#tau)", 5, 0, 5);
    nUnMatchedMu = newTH1D("nUnMatchedMu", "N_{#mu}(unMatch)", 5, 0, 5);
    nUnMatchedEle = newTH1D("nUnMatchedEle", "N_{e}(unMatch)", 5, 0, 5);
    nUnMatchedTau = newTH1D("nUnMatchedTau", "N_{#tau}(unMatch)", 5, 0, 5);

    dRMu1Mu2 = newTH1D("dRMu1Mu2", "#Delta R(#mu_{1}#mu_{2})", 25, 0, 1.0);
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

    unMatchedMuIso = newTH1D("unMatchedMuIso", "#mu^{iso}(unMatch)", 200, 0, 200);
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

    mu1Pt = newTH1D("mu1Pt", "p_{T}(#mu_{1}) [GeV]", NBinsMu1Pt, Mu1PtBin);
    mu1Eta = newTH1D("mu1Eta", "#eta(#mu_{1})", 20, -2.1, 2.1);
    mu1Phi = newTH1D("mu1Phi", "#phi(#mu_{1})", 20, -2.5, 2.5);

    mu2Pt = newTH1D("mu2Pt", "p_{T}(#mu_{2}) [GeV]", NBinsMu2Pt, Mu2PtBin);
    mu2Eta = newTH1D("mu2Eta", "#eta(#mu_{2})", 20, -2.1, 2.1);
    mu2Phi = newTH1D("mu2Phi", "#phi(#mu_{2})", 20, -2.5, 2.5);

    mu3Pt = newTH1D("mu3Pt", "p_{T}(#mu_{3}) [GeV]", NBinsMu3Pt, Mu3PtBin);
    mu3Eta = newTH1D("mu3Eta", "#eta(#mu_{3})", 20, -2.1, 2.1);
    mu3Phi = newTH1D("mu3Phi", "#phi(#mu_{3})", 20, -2.5, 2.5);

    mu4Pt = newTH1D("mu4Pt", "p_{T}(#mu_{4}) [GeV]", NBinsMu4Pt, Mu4PtBin);
    mu4Eta = newTH1D("mu4Eta", "#eta(#mu_{4})", 20, -2.1, 2.1);
    mu4Phi = newTH1D("mu4Phi", "#phi(#mu_{4})", 20, -2.5, 2.5);

    ele1Pt = newTH1D("ele1Pt", "p_{T}(e_{1}) [GeV]", NBinsEle1Pt, Ele1PtBin);
    ele1Eta = newTH1D("ele1Eta", "#eta(e_{1})", 20, -2.1, 2.1);
    ele1Phi = newTH1D("ele1Phi", "#phi(e_{1})", 20, -2.5, 2.5);

    ele2Pt = newTH1D("ele2Pt", "p_{T}(e_{2}) [GeV]", NBinsEle2Pt, Ele2PtBin);
    ele2Eta = newTH1D("ele2Eta", "#eta(e_{2})", 20, -2.1, 2.1);
    ele2Phi = newTH1D("ele2Phi", "#phi(e_{2})", 20, -2.5, 2.5);

    tauPt = newTH1D("tauPt", "p_{T}(#tau) [GeV]", NBinsTauPt, tauPtBin);
    tauEta = newTH1D("tauEta", "#eta(#tau)", 20, -2.1, 2.1);
    tauPhi = newTH1D("tauPhi", "#phi(#tau)", 20, -2.5, 2.5);
    tauMass = newTH1D("tauMass", "M(#tau) [GeV]", 10, 0, 5);
    tauDecayMode = newTH1D("tauDecayMode", "DecayMode(#tau)", 11, 0, 11);

    tau2Pt = newTH1D("tau2Pt", "p_{T}(#tau_{2}) [GeV]", NBinsTauPt, tauPtBin);
    tau2Eta = newTH1D("tau2Eta", "#eta(#tau_{2})", 20, -2.1, 2.1);
    tau2Phi = newTH1D("tau2Phi", "#phi(#tau_{2})", 20, -2.5, 2.5);
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
    dRInvMassEleEle = newTH2D("dRInvMassEleEle", "#Delta R(ee)", "M(ee)[GeV]", 25, 0, 1, 100, 0, 100);
    dRInvMassMu3Tau = newTH2D("dRInvMassMu3Tau", "#Delta R(#mu_{3}#tau)", "M(#mu_{3}#tau)[GeV]", 25, 0, 1, 100, 0, 100);
    dRInvMassEleTau = newTH2D("dRInvMassEleTau", "#Delta R(e#tau)", "M(e#tau)[GeV]", 25, 0, 1, 100, 0, 100);
    dRInvMassTauTau = newTH2D("dRInvMassTauTau", "#Delta R(#tau#tau)", "M(#tau#tau)[GeV]", 25, 0, 1, 100, 0, 100);
    
    nMatchedMuPairNMatchedMuMuPair = newTH2D("nMatchedMuPairNMatchedMuMuPair", "N(#mu_{1}#mu_{2})", "N(#mu_{3}#mu_{4})", 4, 0, 4, 4, 0, 4);
    nMatchedMuPairNMatchedMuElePair = newTH2D("nMatchedMuPairNMatchedMuElePair", "N(#mu_{1}#mu_{2})", "N(#mu_{3}e)", 4, 0, 4, 4, 0, 4);
    nMatchedMuPairNMatchedEleElePair = newTH2D("nMatchedMuPairNMatchedEleElePair", "N(#mu_{1}#mu_{2})", "N(ee)", 4, 0, 4, 4, 0, 4);
    nMatchedMuPairNMatchedMuTauPair = newTH2D("nMatchedMuPairNMatchedMuTauPair", "N(#mu_{1}#mu_{2})", "N(#mu_{3}#tau)", 4, 0, 4, 4, 0, 4);
    nMatchedMuPairNMatchedEleTauPair = newTH2D("nMatchedMuPairNMatchedEleTauPair", "N(#mu_{1}#mu_{2})", "N(e#tau)", 4, 0, 4, 4, 0, 4);
    nMatchedMuPairNMatchedTauTauPair = newTH2D("nMatchedMuPairNMatchedTauTauPair", "N(#mu_{1}#mu_{2})", "N(#tau#tau)", 4, 0, 4, 4, 0, 4);

    //    if (isMC){
    genmu1Pt = newTH1D("genmu1Pt", "p_{T}(#mu_{1}) [GeV]", NBinsgenMu1Pt, genMu1PtBin);
    genmu2Pt = newTH1D("genmu2Pt", "p_{T}(#mu_{1}) [GeV]", NBinsgenMu2Pt, genMu2PtBin);

    //genele1Pt = newTH1D("genele1Pt", "p_{T}(e_{1}) [GeV]", NBinsgenEle1Pt, genEle1PtBin);
    //genele2Pt = newTH1D("genele2Pt", "p_{T}(e_{2}) [GeV]", NBinsgenEle2Pt, genEle2PtBin);

    dRgenMu1genMu2 = newTH1D("dRgenMu1genMu2", "#Delta R(#mu_{1}#mu_{2})", 25, 0, 1.0);
    //    dRgenEle1genEle2 = newTH1D("dRgenEle1genEle2", "#Delta_R(e_{1}e_{2}", 25, 0, 1.0);
     
 //mu1Eta = newTH1D("mu1Eta", "#eta(#mu_{1})", 20, -2.1, 2.1);
      //mu1Phi = newTH1D("mu1Phi", "#phi(#mu_{1})", 20, -2.5, 2.5);

    // }


}

Histomutau::~Histomutau()
{}
