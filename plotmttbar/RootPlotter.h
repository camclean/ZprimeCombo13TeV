#ifndef ROOTPLOTTER_H
#define ROOTPLOTTER_H

#include <TH1.h>
#include <THStack.h>
#include <TH2.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TPostScript.h>
#include <TLatex.h>

class RootPlotter
{

public:

    RootPlotter();
    ~RootPlotter();

    // Open the Root files
    void OpenRootFiles(TObjArray* filenames, const char* cyclename = 0);

    // Plotters
    void PlotHistos(const char* psfilename);
    void PlotYields(TH1* hist);

    // Utilities
    void ShapeNormalize(TH1* hist);
    void ShapeNormalize(THStack* stack);
    TH1* MakePtFractionHisto(TH2D* in2dhist);
    TH1* MakeJetShapeMeanHist(TH2D* in2dhist);
    TH1* MakeMatchingEfficiencyHisto(TH2D* In2DHist);
    TH1* GetSliceHisto(TH2D* In2DHist, Int_t bin);
    TObjArray* FindSubdirs(TDirectory* dir);
    void DrawPageNum(TCanvas* can, Int_t num);
    TObjArray* GetPtBalanceHists();
    TObjArray* GetResolutionHists();
    TH1* StackHistograms(TString subdirname, TString histname);
    TH1* BuildBackground(TString subdirname, TString histname);
    TPostScript* MakeNewPsFile(const char* psfilename);
    void MakeSubTitle(TH1* hist);
    Bool_t ShouldBeStacked(const char* name);
    void Cosmetics(TH1* hist, Int_t isample);
    void CopyStyle(TH1& h1, TH1* h2);
    void ApplyWeight(TH1* hist, Int_t isample);

    // Setters  
    void SetShapeNorm(Bool_t flag = true){bShapeNorm = flag;}
    void SetPortraitMode(Bool_t flag = true){bPortrait = flag;}
    void SetDrawEntries(Bool_t flag = true){bDrawEntries = flag;}
    void PerformFit(Bool_t flag = true){bFitPtBalanceHists = flag;}
    void SetJetShapesPerSlice(Bool_t flag = true){bJetShapesPerSlice = flag;}
    void SetSubstractBkgd(Bool_t flag = true){bSubstractBkgd = flag;}
    void SetPlotRatio(Bool_t flag=true){bPlotRatio = flag;}
    void SetDrawLumi(Bool_t flag=true){bDrawLumi = flag;}
    
    void SetSampleNames(TObjArray* SampleNames);
    void SetHistColors(TArrayI colors);
    void SetHistMarkers(TArrayI markers);
    void SetSamplesToStack(TObjArray* names);
    void SetSamplesWeight(TArrayF weights);
    
private:

    TFile** fRootfiles;         // the Rootfiles 
    Int_t fNumOfSamples;        // how many chains should be plotted?
    
    TObjArray* fSampleNames;    // the sample names

    TObjArray* fSamplesToStack; // the samples that should be stacked
    Int_t fNumOfSamplesToStack; // number of samples to stack

    TArrayF fSamplesWeight;     // weights of the chains
    TArrayI fSampleColors;      // colors of the chains
    TArrayI fSampleMarkers;     // markers of the chains

    Bool_t  bShapeNorm;         // use shape normalization
    Bool_t  bPortrait;          // portrait or landscape mode
    Bool_t  bDrawEntries;       // display the number of entries 
    Bool_t  bDrawLumi;          // display the lumi information 
    Bool_t  bFitPtBalanceHists; // perform fit of Pt-balance histograms?
    Bool_t  bJetShapesPerSlice; // perform fit of Pt-balance histograms?
    Bool_t  bSubstractBkgd;     // substract background from the first chain?
    Bool_t  bPlotRatio;         // should a ratio be plotted?
  
};

#endif //  __JETANAPLOTTER_H
