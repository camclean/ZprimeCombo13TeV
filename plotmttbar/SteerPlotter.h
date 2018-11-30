#ifndef STEERPLOTTER__
#define STEERPLOTTER__
#include <cstdlib>
#include <TObjArray.h>
#include <TArrayI.h>
#include <TArrayF.h>
#include <TString.h>
#include "BaseSteer.h"


class SteerPlotter: public BaseSteer {

public:
    SteerPlotter();
    virtual ~SteerPlotter();
  
    virtual void Print(Option_t* option="") const;
  
    void SetShapeNorm(Bool_t flag);
    Bool_t GetShapeNorm();

    void SetLumiNorm(Bool_t flag);
    Bool_t GetLumiNorm();
   
    void SetRatioPlot(Bool_t flag);
    Bool_t GetRatioPlot();

    void SetZScoreInRatio(Bool_t flag);
    Bool_t GetZScoreInRatio();

    void SetPortrait(Bool_t flag);
    Bool_t GetPortrait();

    void SetFitPtBalanceHists(Bool_t flag);
    Bool_t GetFitPtBalanceHists();

    void SetJetShapesPerSlice(Bool_t flag);
    Bool_t GetJetShapesPerSlice();

    void SetDrawEntries(Bool_t flag);
    Bool_t GetDrawEntries();

    void SetSampleNames(const char* in);
    TObjArray* GetSampleNames();

    void SetInputFiles(const char* in);
    TObjArray* GetInputFiles();

    void SetOutputPsFile(const char* in);
    const char* GetOutputPsFile();

    void SetCycleName(const char* in);
    const char* GetCycleName();

    void SetLegStrings(const char* in);
    TObjArray* GetLegStrings();

    void SetScaleSysUnc(const char* in);
    TObjArray* GetScaleSysUnc();

    void SetSysUncWeight(const char* in);
    TArrayF GetSysUncWeight();

    void SetHistColors(const char* in);
    TArrayI GetHistColors();
   
    void SetHistMarkers(const char* in);
    TArrayI GetHistMarkers();

    void SetSamplesToStack(const char* in);
    TObjArray* GetSamplesToStack();

    void SetSamplesWeight(const char* in);
    TArrayF GetSamplesWeight();

    void SetSamplesUnc(const char* in);
    TArrayF GetSamplesUnc();
    
    void SetSubstractBkgd(Bool_t flag);
    Bool_t GetSubstractBkgd();

    void SetDrawLumi(Bool_t flag);
    Bool_t GetDrawLumi();

    void SetForPrelim(Bool_t flag);
    Bool_t GetForPrelim();

    void SetForPublication(Bool_t flag);
    Bool_t GetForPublication();

    void SetDrawLegend(Bool_t flag);
    Bool_t GetDrawLegend();

    void SetLumi(Float_t lumi);
    Float_t GetLumi();

    void SetSysError(Float_t err);
    Float_t GetSysError();

    void SetDoCumulative(Bool_t flag);
    Bool_t GetDoCumulative();

    void SetSingleEPS(Bool_t flag);
    Bool_t GetSingleEPS();

    void SetIgnoreEmptyBins(Bool_t flag);
    Bool_t GetIgnoreEmptyBins();

    void SetPlotThetaFile(Bool_t flag);
    Bool_t GetPlotThetaFile();

    void SetLogy(Bool_t flag);
    Bool_t GetLogy();

private:

    Bool_t    bShapeNorm;         // Shape normalization?
    Bool_t    bLumiNorm;          // Lumi normalization?
    Bool_t    bRatioPlot;         // plot ratios
    Bool_t    bZScoreInRatio;     // plot z-score instead of usual ratio
    Bool_t    bLogy;              // plot y-axis on log scale?
    Bool_t    bPortrait;          // portrait or landscape 
    Bool_t    bDrawEntries;       // draw the number of entries?
    Bool_t    bDrawLumi;          // draw the lumi information?
    Bool_t    bForPrelim;         // write "CMS Preliminary"
    Bool_t    bForPublication;    // write "CMS"
    Bool_t    bDrawLegend;        // draw the legend everywhere?
    Bool_t    bFitPtBalanceHists; // fit Pt-balance histograms?
    Bool_t    bJetShapesPerSlice; // plot each slide of the jet shape histograms?
    Bool_t    bDoCumulative;      // do cumulative distributions instead of normal plots
    Bool_t    bSingleEPS;         // make one EPS file for each histogram
    Bool_t    bPlotThetaFile;     // take input histograms from theta file
    Bool_t    bIgnoreEmptyBins;   // ignore empty bins in the ratio
    Int_t     fNumOfSamples;      // how many analysis samples should be plotted
    Float_t   fLumi;              // integrated luminosity of sample
    Float_t   fSysError;          // systematic error on normalisation
    TObjArray  fSampleNames;      // all sample name

    TObjArray fSamplesToStack;    // name of samples that should be stacked on top of each other
    Int_t     fNumOfSamplesToStack; // how many samples should be stacked
    Bool_t    bSubstractBkgd;     // substract all background samples?

    TArrayF   fSamplesWeight;     // weights for the different samples
    TArrayF   fSamplesUnc;        // uncertainty on the normalisation of various samples
    
    TObjArray fInputFiles;        // input filenames
    TString   fCycleName;         // name of the cycle, can be used as prefix for the filenames
    TString   fOutputPsFile;      // name of the resulting ps file

    TArrayI   fHistColors;        // the histogram colors
    TArrayI   fHistMarkers;       // the histogram markers

    TObjArray fLegStrings;        // legend entries

    TObjArray fScaleSysUnc;       // names of systematic uncertainties to scale 
    TArrayF fSysUncWeight;        // scaling factor for systematic uncertainties
 
    ClassDef(SteerPlotter,0)      // steering class for the SFrame Plotter

};

#endif
