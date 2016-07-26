// Dear emacs, this is -*- c++ -*-
//______________________________________________________
//
// Plots.cxx
//
// All the output histograms of the SFrame analysis (also other interesting 
// information - can be included at a later time) are being saved into a 
// root file as genuine root objects.
//
// This program can be used to combine, sum and plot all histograms (data, 
// background, signal...) and save them to a postscript file. The program 
// can be steered with a steer file, usually called Plots.steer. 
//
// Authors    : Roman Kogler
// Created    : 2012
// Last update: 
//          by: 

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include <TROOT.h>
#include <TH1.h>
#include <TRint.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TObjString.h>

#include "SteerParser.h"
#include "SteerPlotter.h"
#include "FileParser.h"
#include "SHist.h"
#include "SPlotter.h"


using namespace std;

// Initialize the root framework
TROOT Plots("Plots","SFrame Plots");

int main(int argc, char** argv)
{

  // ________ parse the command line ___________
  TString filename;
  if (argc < 3) { // Check the value of argc. If not enough parameters have been passed, inform user and exit.
    cout << "Not enough arguments provided. Usage:\n"; // Inform the user of how to use the program
    cout << argv[0] << " -f <steerfile>\n";
    exit(0);
  } else { // we got enough parameters
    cout << "Executing: " << argv[0] << " ";
    for (int i = 1; i < argc; i++) { // iterate over argv[] to get the parameters
      std::cout << argv[i] << " ";
      if (i+1 != argc){ // Check that we haven't finished parsing already
	TString arg(argv[i]);
	if (arg == "-f") {
	  filename = argv[i+1];
	} else {
	  std::cout << "Not enough or invalid arguments, please try again.\n";
	  exit(0);
	}
      }
    }
    cout << endl;
  }

  TRint theApp("App",&argc,argv); 
  gROOT->Reset();

  // Stops ROOT keeping references to all histograms
  TH1::AddDirectory(kFALSE);
  // Run in batch mode by default - this supresses the canvas -> speed increase
  gROOT->SetBatch(kTRUE);   


  // ____________ process the steering ______________
  
  SteerParser parser;
  parser.ParseFile(filename);

  SteerPlotter* steerfile = (SteerPlotter*) parser.GetSteer(SteerPlotter::Class());
  //steerfile->Print();

  TString CycleName          = steerfile->GetCycleName();
  TObjArray* InputFilenames  = steerfile->GetInputFiles();
  TString PsFilename         = steerfile->GetOutputPsFile();
  Bool_t ShapeNorm           = steerfile->GetShapeNorm();
  Bool_t RatioPlot           = steerfile->GetRatioPlot();
  Bool_t zscore              = steerfile->GetZScoreInRatio();
  Bool_t PortraitMode        = steerfile->GetPortrait();
  Bool_t DrawEntries         = steerfile->GetDrawEntries();
  //Bool_t PtBalanceFitOpt     = steerfile->GetFitPtBalanceHists();
  //Bool_t SubstractBkgd       = steerfile->GetSubstractBkgd();
  //Bool_t JetShapesPerSlice   = steerfile->GetJetShapesPerSlice();
  TObjArray* SampleNames     = steerfile->GetSampleNames();
  TObjArray* SamplesToStack  = steerfile->GetSamplesToStack(); 
  TArrayI HistColors         = steerfile->GetHistColors();
  TArrayI HistMarkers        = steerfile->GetHistMarkers();
  TArrayF SamplesWeight      = steerfile->GetSamplesWeight();
  TArrayF SamplesUnc         = steerfile->GetSamplesUnc();
  Bool_t DrawLumi            = steerfile->GetDrawLumi();
  Bool_t ForPrelim           = steerfile->GetForPrelim();
  Bool_t ForPublication      = steerfile->GetForPublication();
  Float_t Lumi               = steerfile->GetLumi();
  Float_t SysErr             = steerfile->GetSysError();
  Bool_t DrawLegend          = steerfile->GetDrawLegend();
  Bool_t DoCumulative        = steerfile->GetDoCumulative();
  Bool_t SingleEPS           = steerfile->GetSingleEPS();
  Bool_t ThetaFile           = steerfile->GetPlotThetaFile();
  Bool_t Logy                = steerfile->GetLogy();
  Bool_t IgnoreEmptyBins     = steerfile->GetIgnoreEmptyBins();

  TObjArray* ScaleSysUnc     = steerfile->GetScaleSysUnc();
  TArrayF  SysUncWeight      = steerfile->GetSysUncWeight();
  

  // _______________ loop over files and get all histograms ______________

  FileParser fp;
  //fp.SetDebug();
  fp.SetDoCumulative(DoCumulative);  
  vector<TObjArray*> harr;
  vector<TObjArray*> harr_sys;

  if (!ThetaFile){ // standard case
    for (int i=0; i<InputFilenames->GetEntries(); ++i){
      TString file = ((TObjString*)InputFilenames->At(i))->GetString();
      TString legname = ((TObjString*)SampleNames->At(i))->GetString();
      fp.OpenFile(file, CycleName);
      fp.BrowseFile();
      float unc = 0.;
      if (SamplesUnc.GetSize()!=0) unc = SamplesUnc[i];
      fp.SetInfo(legname, SamplesWeight[i], HistColors[i], HistMarkers[i], unc);
      fp.CloseFile();
      harr.push_back( fp.GetHists() );
    }

  } else { // use theta input for plots
    fp.OpenThetaFile(CycleName);
    for (int i=0; i<InputFilenames->GetEntries(); ++i){
      TString sample = ((TObjString*)InputFilenames->At(i))->GetString();
      TString legname = ((TObjString*)SampleNames->At(i))->GetString();
      fp.BrowseThetaFile(sample);
      float unc = 0.;
      if (SamplesUnc.GetSize()!=0) unc = SamplesUnc[i];
      fp.SetInfo(legname, SamplesWeight[i], HistColors[i], HistMarkers[i], unc);
      harr.push_back( new TObjArray(*fp.GetHists()) );
      if ((fp.GetShapeSys())->GetEntries()>0){
	harr_sys.push_back( new TObjArray(*fp.GetShapeSys()) );
      }
      fp.Clear();
    }
    fp.CloseFile();
    //exit(4);
  }

  // _______________ stack histograms ______________
  
  SPlotter pl;
  //pl.SetDebug();
  if (!ThetaFile){
    pl.DoStacking(harr, SamplesToStack);
  } else {
    pl.DoStacking(harr, SamplesToStack, true);
  }

  // ____________ set up the plotter ______________

  pl.SetShapeNorm(ShapeNorm);
  pl.SetPortraitMode(PortraitMode);
  pl.SetDrawEntries(DrawEntries);
  pl.SetPlotRatio(RatioPlot);
  pl.SetZScoreInRatio(zscore);
  pl.SetDrawLumi(DrawLumi); 
  pl.SetDrawLegend(DrawLegend);
  pl.SetPsFilename(PsFilename);
  pl.SetLumi(Lumi);
  pl.SetNormError(SysErr);
  pl.SetSingleEPSMode(SingleEPS);
  pl.SetForPublication(ForPublication);
  pl.SetForPrelim(ForPrelim);
  pl.SetLogy(Logy);
  pl.SetIgnoreEmptyBins(IgnoreEmptyBins);

  if (ThetaFile){
    pl.SetShapeSysHists(harr_sys);

    pl.SetScaleSysUnc(ScaleSysUnc);
    pl.SetSysUncWeight(SysUncWeight);
  }
  
  // _______________ do the plotting ______________
  
  pl.ProcessAndPlot(harr);

  cout << "\nDone processing all plots. Wrote files to: " << endl;
  cout << PsFilename << endl << endl;

  // Done! Exit Root
  gSystem->Exit(0);

  return 0;

}
