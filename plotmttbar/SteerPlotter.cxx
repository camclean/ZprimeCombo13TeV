#include "SteerPlotter.h"

#include <iostream>
#include <iomanip>
#include <TObjString.h>

using namespace std;

ClassImp(SteerPlotter)

SteerPlotter::SteerPlotter()
{
  // Here you have to set the default values which will be used,
  // if this steering class is required but not given in the steering file
   
   fNumOfSamples = 0;  
   fOutputPsFile = "Plots/test.ps";
   bSubstractBkgd = false;
   bDrawEntries = false;
   bDrawLumi = false;
   bDrawLegend = true;
   bSingleEPS = false;
   bIgnoreEmptyBins = false;
   bDoCumulative = false;
   fNumOfSamplesToStack = 0;
   fLumi = 0;
   fSysError = -1.;
   bPlotThetaFile = false;
   bLogy = false;
   bZScoreInRatio = false;
}

SteerPlotter::~SteerPlotter()
{
   fLegStrings.Delete();
}

void SteerPlotter::Print(Option_t* opt) const
{
  // Prints all settings of the steering  

  // First: perform some sanity checks
  if (fNumOfSamples != fInputFiles.GetEntries()){
    cout << "Error: Number of input files is not the same as the number of samples to be plotted." << endl;
    exit(3);
  }
  if (fNumOfSamples != fSamplesWeight.GetSize()){
    cout << "Error: Number of given weights is not the same as number of samples." << endl;
    exit(3);
  }
  if (fSamplesUnc.GetSize()!=0){
    if (fNumOfSamples != fSamplesUnc.GetSize()){
      cout << "Error: Number of given uncertainties for normalisation of samples is not the same as number of samples." << endl;
      exit(3);
    }
  }
  if (fNumOfSamples != fHistColors.GetSize()){
    cout << "Error: Number of colours given is not the same as number of samples." << endl;
    exit(3);
  }
  if (fNumOfSamples != fHistMarkers.GetSize()){
    cout << "Error: Number of markers given is not the same as number of samples." << endl;
    exit(3);
  }
  if (fNumOfSamples < fSamplesToStack.GetEntries() ){
    cout << "Error: Number of samples is smaller than the number of samples to stack." << endl;
    exit(3);
  }


  cout << endl;
  cout << "-------------------------------------------------------- SteerPlotter " << opt << "-----------------------------------------------" << endl;

  cout << "Number of analysis samples to be plotted: " << fNumOfSamples << endl;
  if (fCycleName.Length()>0){
    cout << "Cylcle Name (used as prefix for root files): " << fCycleName << endl;
  }
  for (Int_t i=0; i<fNumOfSamples; ++i){
    TString name(((TObjString*) fSampleNames.At(i))->GetName() );
    cout << "File of sample " << i << ":  " << setw(25) << ((TObjString*)fInputFiles.At(i))->GetName()
	 << "   name =  " << setw(15) << name
    	 << "   color = " << setw(4) << fHistColors.At(i) 
	 << "   marker = " << setw(4) << fHistMarkers.At(i) 
	 << "   with weight " << fSamplesWeight.At(i);
    if (fSamplesUnc.GetSize()!=0){
      cout << "  and uncertainty of " << fSamplesUnc.At(i)*100 << "%" << endl;
    } else {
      cout << endl;
    }
  }
  cout << "Output Ps File:                " << fOutputPsFile << endl;
  cout << endl;
  if (fNumOfSamplesToStack>0){
      cout << "These samples will be stacked:" << endl;
      for (Int_t i=0; i<fNumOfSamplesToStack; ++i){
          TString name(((TObjString*) fSamplesToStack.At(i))->GetName() );
          cout << "    Name of sample " << i << " in stack :      " << setw(15) << name << endl;
      }
      
  } else {
      cout << "No stacking will be plotted." << endl;
  }
  if (bPlotThetaFile){
    if (fScaleSysUnc.GetEntries()>0){
      cout << endl;
      if (fScaleSysUnc.GetEntries() != fSysUncWeight.GetSize()){
	cout << "Error: inconsistent number of entries in fScaleSysUnc and fSysUncWeight! Arrays must have same length, please correct steering." << endl;
	exit(3);
      }
      cout << "Systematic errors will be scaled with these factors:" << endl;
      for (Int_t i=0; i<fScaleSysUnc.GetEntries(); ++i){
          TString name(((TObjString*) fScaleSysUnc.At(i))->GetName() );
          cout << "    Name of systematic unc " << i << " = " << setw(15) << name << " scale with factor : " << fSysUncWeight.At(i) << endl;
      }
      
    }
  }

  if (bSubstractBkgd){
    cout << "Background will be substracted from sample 0: " << (((TObjString*) fSampleNames.At(0))->GetName() ) << endl;
  } else {
    cout << "No background substraction" << endl;
  }
  cout << (bRatioPlot? "Ratios will be plotted." : "No ratios will be plotted") << endl;
  if (bRatioPlot){
    cout << (bZScoreInRatio? "Z-Score will be plotted in bottom pad." : "Usual data/MC ratios will be plotted in bottom pad.") << endl;
  }
  cout << (bDrawEntries? "Number of histogram entries will be plotted." : "Number of histogram entries will not be plotted") << endl;
  cout << (bDrawLumi? "Lumi inforamtion will be plotted." : "Lumi inforamtion will not be plotted") << endl;
  cout << "Integrated luminosity = " << fLumi << " fb-1" << endl;
  if (fSysError>0){
    cout << "Normalisation error of " << fSysError*100 << "% will be drawn." << endl;
  } else {
    cout << "No normalisation error will be drawn." << endl;
  }

  if (bDrawLegend) cout << "Legend will be plotted everywhere." << endl;
  else cout << "Legend will be plotted on first plot only" << endl;

  if (bShapeNorm) cout << "Shape normalization" << endl;
  else cout << "No shape normalization" << endl;

  if (bDoCumulative)  cout << "Cumulative distributions will be plotted." << endl;
  else cout << "Normal distributions will be plotted" << endl;

  if (bSingleEPS) cout << "Creating one eps file per histogram." << endl;
  else cout << "Creating one ps file with all histograms for each histogram collection." << endl;

  if (bIgnoreEmptyBins) cout << "Empty bins will not be plotted in the ratio." << endl;
  else cout << "Empty bins will have infinite error in the ratio." << endl;
  
  if (bLumiNorm) cout << "Luminosity normalization" << endl; 
  else cout << "No lumi normalization" << endl;
  
  if (bLumiNorm) cout << "Luminosity normalization" << endl;
  else cout << "No lumi normalization" << endl;

  if (bPortrait) cout << "Setting the page to portrait mode" << endl;
  else cout << "Setting the page to landscape mode" << endl;

  if (bPlotThetaFile) cout << "Creating plots from one input theta file." << endl;
  else cout << "Using standard SFrame output for plots." << endl;
  cout << "--------------------------------------------------------------------------------------------------------------------" << endl;

}

void SteerPlotter::SetShapeNorm(Bool_t flag){bShapeNorm = flag;}
Bool_t SteerPlotter::GetShapeNorm(){return bShapeNorm;}

void SteerPlotter::SetLumiNorm(Bool_t flag){bLumiNorm = flag;}
Bool_t SteerPlotter::GetLumiNorm(){return bLumiNorm;}

void SteerPlotter::SetRatioPlot(Bool_t flag){bRatioPlot = flag;}
Bool_t SteerPlotter::GetRatioPlot(){return bRatioPlot;}

void SteerPlotter::SetZScoreInRatio(Bool_t flag){bZScoreInRatio = flag;}
Bool_t SteerPlotter::GetZScoreInRatio(){return bZScoreInRatio;}

void SteerPlotter::SetPortrait(Bool_t flag){bPortrait = flag;}
Bool_t SteerPlotter::GetPortrait(){return bPortrait;}

void SteerPlotter::SetFitPtBalanceHists(Bool_t flag){bFitPtBalanceHists = flag;}
Bool_t SteerPlotter::GetFitPtBalanceHists(){return bFitPtBalanceHists;}

void SteerPlotter::SetDrawEntries(Bool_t flag){bDrawEntries = flag;}
Bool_t SteerPlotter::GetDrawEntries(){return bDrawEntries;}

void SteerPlotter::SetDrawLumi(Bool_t flag){bDrawLumi = flag;}
Bool_t SteerPlotter::GetDrawLumi(){return bDrawLumi;}

void SteerPlotter::SetForPrelim(Bool_t flag){bForPrelim = flag;}
Bool_t SteerPlotter::GetForPrelim(){return bForPrelim;}

void SteerPlotter::SetForPublication(Bool_t flag){bForPublication = flag;}
Bool_t SteerPlotter::GetForPublication(){return bForPublication;}

void SteerPlotter::SetDrawLegend(Bool_t flag){bDrawLegend = flag;}
Bool_t SteerPlotter::GetDrawLegend(){return bDrawLegend;}

void SteerPlotter::SetDoCumulative(Bool_t flag){bDoCumulative = flag;}
Bool_t SteerPlotter::GetDoCumulative(){return bDoCumulative;}

void SteerPlotter::SetSingleEPS(Bool_t flag){bSingleEPS = flag;}
Bool_t SteerPlotter::GetSingleEPS(){return bSingleEPS;}

void SteerPlotter::SetIgnoreEmptyBins(Bool_t flag){bIgnoreEmptyBins = flag;}
Bool_t SteerPlotter::GetIgnoreEmptyBins(){return bIgnoreEmptyBins;}

void SteerPlotter::SetPlotThetaFile(Bool_t flag){bPlotThetaFile = flag;}
Bool_t SteerPlotter::GetPlotThetaFile(){return bPlotThetaFile;}

void SteerPlotter::SetJetShapesPerSlice(Bool_t flag){bJetShapesPerSlice = flag;}
Bool_t SteerPlotter::GetJetShapesPerSlice(){return bJetShapesPerSlice;}

void SteerPlotter::SetLogy(Bool_t flag){bLogy = flag;}
Bool_t SteerPlotter::GetLogy(){return bLogy;}

void SteerPlotter::SetLumi(Float_t lumi){fLumi = lumi;}
Float_t SteerPlotter::GetLumi(){return fLumi;}

void SteerPlotter::SetSysError(Float_t err){fSysError = err;}
Float_t SteerPlotter::GetSysError(){return fSysError;}

void SteerPlotter::SetSampleNames(const char* in) {
    this->SplitString(in,",",&fSampleNames);
    fNumOfSamples = fSampleNames.GetEntries();
}
TObjArray* SteerPlotter::GetSampleNames(){return &fSampleNames;}

void SteerPlotter::SetInputFiles(const char* in){ this->SplitString(in,",",&fInputFiles);}
TObjArray* SteerPlotter::GetInputFiles() {return &fInputFiles;}

void SteerPlotter::SetOutputPsFile(const char* in) {fOutputPsFile = in;}
const char* SteerPlotter::GetOutputPsFile() {return fOutputPsFile.Data();}

void SteerPlotter::SetCycleName(const char* in) {fCycleName = in;}
const char* SteerPlotter::GetCycleName() {return fCycleName.Data();}

void SteerPlotter::SetLegStrings(const char* in){this->SplitString(in,",",&fLegStrings);}
TObjArray* SteerPlotter::GetLegStrings() {return &fLegStrings;}

void SteerPlotter::SetHistColors(const char* in){this->StringToArray(in, fHistColors);}
TArrayI SteerPlotter::GetHistColors(){return fHistColors;}

void SteerPlotter::SetHistMarkers(const char* in){this->StringToArray(in, fHistMarkers);}
TArrayI SteerPlotter::GetHistMarkers(){return fHistMarkers;}

void SteerPlotter::SetSamplesToStack(const char* in){ this->SplitString(in,",",&fSamplesToStack); fNumOfSamplesToStack = fSamplesToStack.GetEntries();}
TObjArray* SteerPlotter::GetSamplesToStack(){return &fSamplesToStack;}

void SteerPlotter::SetScaleSysUnc(const char* in){ this->SplitString(in,",",&fScaleSysUnc);}
TObjArray* SteerPlotter::GetScaleSysUnc(){return &fScaleSysUnc;}

void SteerPlotter::SetSysUncWeight(const char* in){ this->StringToArray(in, fSysUncWeight);}
TArrayF SteerPlotter::GetSysUncWeight(){return fSysUncWeight;}

void SteerPlotter::SetSamplesWeight(const char* in){ this->StringToArray(in, fSamplesWeight);}
TArrayF SteerPlotter::GetSamplesWeight(){return fSamplesWeight;}

void SteerPlotter::SetSamplesUnc(const char* in){ this->StringToArray(in, fSamplesUnc);}
TArrayF SteerPlotter::GetSamplesUnc(){return fSamplesUnc;}

void SteerPlotter::SetSubstractBkgd(Bool_t flag){ bSubstractBkgd = flag; }
Bool_t SteerPlotter::GetSubstractBkgd(){ return bSubstractBkgd; }
