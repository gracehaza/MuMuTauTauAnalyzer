#include <fstream>
#include <TString.h>
#include "ArgParser.h"
#include "ConfigArg.h"
#include "MuMuTauMuTauMuAnalyzer.h"
#include "lumiana.h"

using namespace std;

int main(int argc, char **argv)
{
    //--- Load configuration ---
    ConfigArg cfg;

    TString inputFile      = cfg.getS("inputFile");
    TString outputDir      = cfg.getS("outputDir");
    TString doWhat         = cfg.getS("doWhat", "DYJETS");
    Long_t maxEvents       = cfg.getL("maxEvents", -1);
    double lumi            = cfg.getD("lumi", 1);
    bool invertedMu2Iso    = cfg.getB("invertedMu2Iso", false);
    bool invertedMu4Iso    = cfg.getB("invertedMu4Iso", false);
    double Mu2IsoThreshold = cfg.getD("Mu2IsoThreshold", 0.25);
    double Mu3IsoThreshold = cfg.getD("Mu3IsoThreshold", 0.25);   
    double Mu4IsoThreshold = cfg.getD("Mu4IsoThreshold", 0.25);
    double diMuonMassLowThreshold  = cfg.getD("diMuonMassLowThreshold", 0);
    double diMuonMassHighThreshold = cfg.getD("diMuonMassHighThreshold", 25.0);

    //--- Parse the arguments -----------------------------------------------------
    if (argc > 1)
    {
        for (int i = 1; i < argc; ++i)
        {
            TString currentArg = argv[i];
            
            //--- possible options ---
            if (currentArg.BeginsWith("inputFile="))
            {
                getArg(currentArg, inputFile);
            }

            else if (currentArg.BeginsWith("outputDir="))
            {
                getArg(currentArg, outputDir);
            }

            else if (currentArg.BeginsWith("doWhat="))
            {
                getArg(currentArg, doWhat);
            }

            else if (currentArg.BeginsWith("maxEvents="))
            {
                getArg(currentArg, maxEvents);
            }

            else if (currentArg.BeginsWith("lumi="))
            {
                getArg(currentArg, lumi);
            }

            else if (currentArg.BeginsWith("invertedMu2Iso="))
            {
                getArg(currentArg, invertedMu2Iso);
            }

            else if (currentArg.BeginsWith("invertedMu4Iso="))
            {
                getArg(currentArg, invertedMu4Iso);
            }

            else if (currentArg.BeginsWith("Mu2IsoThreshold="))
            {
                getArg(currentArg, Mu2IsoThreshold);
            }

	    else if (currentArg.BeginsWith("Mu3IsoThreshold="))
	      {
                getArg(currentArg, Mu3IsoThreshold);
	      }

            else if (currentArg.BeginsWith("Mu4IsoThreshold="))
            {
                getArg(currentArg, Mu4IsoThreshold);
            }

	    else if (currentArg.BeginsWith("diMuonMassLowThreshold="))
            {
                getArg(currentArg, diMuonMassLowThreshold);
            }

            else if (currentArg.BeginsWith("diMuonMassHighThreshold="))
            {
                getArg(currentArg, diMuonMassHighThreshold);
            }
        } // end for loop in argc
    } // end if argc > 1
    
    doWhat.ToUpper();

    //------------------------- data histograms production -------------------
    if (doWhat == "DATA" || doWhat == "ALL")
    {
        if (inputFile.EndsWith(".root"))
        {
            MuMuTauMuTauMuAnalyzer DataHist(inputFile, outputDir, 1, 1, maxEvents, false, invertedMu2Iso, invertedMu4Iso, Mu2IsoThreshold, Mu4IsoThreshold, diMuonMassLowThreshold, diMuonMassHighThreshold);
            DataHist.Loop();
        } // end if inputFile.EndsWith(".root")
        
        else{
            ifstream finTree;
            finTree.open(inputFile);
            string fileName;
            while (getline(finTree, fileName))
            {
                MuMuTauMuTauMuAnalyzer DataHist(fileName, outputDir, 1, 1, maxEvents, false, invertedMu2Iso, invertedMu4Iso, Mu2IsoThreshold, Mu4IsoThreshold, diMuonMassLowThreshold, diMuonMassHighThreshold);
                DataHist.Loop();
            } // end while loop on file list 
        } // end else inputFile.EndsWith(".root")
    } // end if data

    // ----- This variable is to be used for MC reweighting -----
    float summedWeights = 0;

    // ------------------------ MC histogram production ------------------
    if (doWhat == "DYJETS" || doWhat == "ALL")
    {
        if (inputFile.EndsWith(".root"))
        {
            lumiana DYJetsLumi(inputFile);
            summedWeights = DYJetsLumi.Loop();
            MuMuTauMuTauMuAnalyzer DYJetsHist(inputFile, outputDir, lumi*2075.14*3*1000, summedWeights, maxEvents, true, invertedMu2Iso, invertedMu4Iso, Mu2IsoThreshold, Mu4IsoThreshold, diMuonMassLowThreshold, diMuonMassHighThreshold);
            DYJetsHist.Loop();
        } // end if inputFile.EndsWith(".root")

        else{
            ifstream finLumi;
            finLumi.open(inputFile);
            string fileName;
            while (getline(finLumi, fileName))
            {
                lumiana DYJetsLumi(fileName);
                summedWeights += DYJetsLumi.Loop();
            }// end while loop for weight sum

            ifstream finTree;
            finTree.open(inputFile);
            while (getline(finTree, fileName))
            {
                MuMuTauMuTauMuAnalyzer DYJetsHist(fileName, outputDir, lumi*2075.14*3*1000, summedWeights, maxEvents, true, invertedMu2Iso, invertedMu4Iso, Mu2IsoThreshold, Mu4IsoThreshold, diMuonMassLowThreshold, diMuonMassHighThreshold);
                DYJetsHist.Loop();
            } // end while loop on input file list
        } // end else
    } // end if DYJets

    // --- always need to reinitialize the weight parameter for new sample -----
    summedWeights = 0;


    if (doWhat == "WJETS" || doWhat == "ALL")
    {
        if (inputFile.EndsWith(".root"))
        {
            lumiana WJetsLumi(inputFile);
            summedWeights = WJetsLumi.Loop();
            MuMuTauMuTauMuAnalyzer WJetsHist(inputFile, outputDir, lumi*61526.7*1000, summedWeights, maxEvents, true, invertedMu2Iso, invertedMu4Iso, Mu2IsoThreshold, Mu4IsoThreshold, diMuonMassLowThreshold, diMuonMassHighThreshold);
            WJetsHist.Loop();
        } // end if inputFile.EndsWith(".root")

        else{
            ifstream finLumi;
            finLumi.open(inputFile);
            string fileName;
            while (getline(finLumi, fileName))
            {
                lumiana WJetsLumi(fileName);
                summedWeights += WJetsLumi.Loop();
            }// end while loop for weight sum

            ifstream finTree;
            finTree.open(inputFile);
            while (getline(finTree, fileName))
            {
                MuMuTauMuTauMuAnalyzer WJetsHist(fileName, outputDir, lumi*61526.7*1000, summedWeights, maxEvents, true, invertedMu2Iso, invertedMu4Iso, Mu2IsoThreshold, Mu4IsoThreshold, diMuonMassLowThreshold, diMuonMassHighThreshold);
                WJetsHist.Loop();
            } // end while loop on input file list
        } // end else
    } // end if WJets

    // --- always need to reinitialize the weight parameter for new sample -----
    summedWeights = 0;


    if (doWhat == "H125AA5" || doWhat == "ALL")
    {
        if (inputFile.EndsWith(".root"))
        {
            lumiana H125AA5Lumi(inputFile);
            summedWeights = H125AA5Lumi.Loop();
            MuMuTauMuTauMuAnalyzer H125AA5Hist(inputFile, outputDir, 1, 1, maxEvents, true, invertedMu2Iso, invertedMu4Iso, Mu2IsoThreshold, Mu3IsoThreshold, Mu4IsoThreshold, diMuonMassLowThreshold, diMuonMassHighThreshold);
            H125AA5Hist.Loop();
        } // end if inputFile.EndsWith(".root")
        
        else{
            ifstream finLumi;
            finLumi.open(inputFile);
            string fileName;
            while (getline(finLumi, fileName))
            {
                lumiana H125AA5Lumi(fileName);
                summedWeights += H125AA5Lumi.Loop();
            } // end while loop for weight sum

            ifstream finTree;
            finTree.open(inputFile);
            while (getline(finTree, fileName))
            {
	      MuMuTauMuTauMuAnalyzer H125AA5Hist(fileName, outputDir, 1, 1, maxEvents, true, invertedMu2Iso, invertedMu4Iso, Mu2IsoThreshold, Mu3IsoThreshold, Mu4IsoThreshold, diMuonMassLowThreshold, diMuonMassHighThreshold);
                H125AA5Hist.Loop();
            } // end while loop on input file list
        } // end else 
    } // end if H125AA5 signal

    return 0;
}
