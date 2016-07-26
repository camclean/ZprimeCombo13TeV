#include "RootPlotter.h"

#include <iostream>
#include <iomanip>

#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TStyle.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TClass.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TText.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TProfile.h>
#include <TF1.h>
#include <THStack.h>
#include <TFile.h>

#include "tdrstyle.C"


using namespace std;

/********************************************************
 *
 * RootPlotter
 *
 * A plotting tool using ROOT's functionality. 
 * 
 *******************************************************/


RootPlotter::RootPlotter()
{
  // some default settings
  fRootfiles = NULL;
  fNumOfSamples = 0;
  bShapeNorm = true;
  fSampleNames = NULL;
  bPortrait = true;
  bDrawEntries = false;
  bFitPtBalanceHists = false;
  bJetShapesPerSlice = false;
  bSubstractBkgd = false;
  bPlotRatio = false;
  bDrawLumi = false;

}

RootPlotter::~RootPlotter()
{
  // destructor: clean up
  if (fRootfiles){
    delete[] fRootfiles;
  }
  
}

void RootPlotter::OpenRootFiles(TObjArray* filenames, const char* cyclename)
{
  // open the root files with names given in the TObjArray

  cout << endl << "----------------------------------- RootPlotter -----------------------------------" << endl;
  fRootfiles = new TFile*[fNumOfSamples];
  fNumOfSamples = filenames->GetEntries();
  for (Int_t i=0; i<fNumOfSamples; ++i){
    TString name(((TObjString*)filenames->At(i))->GetName());

    if (cyclename){
      TString Prefix(cyclename);
      Prefix.Append(".");
      name.Prepend(Prefix);
    }

    // check if name consists of a wildcard, if so use hadd to combine histograms
    if (name.Contains("*")){
      TString target(name);
      target.ReplaceAll("*","");
      TString command = "hadd -f " + target + " " + name;
      gSystem->Exec(command);
      name = target;
    }

    //cout << "Opening file with name " << name << "..." << endl;
    fRootfiles[i] = new TFile(name, "READ");
    //cout << "... success! pointer = " << fRootfiles[i] << endl;
    //cout << "name = " << fRootfiles[i]->GetName() << endl;
    //fRootfiles[i]->ls();
    
    //cout << " is open? " << fRootfiles[i]->IsOpen() << endl;
    
    if (!fRootfiles[i]->IsOpen()) {
      cout << endl << "RootPlotter: File " << name << " does not exist!!!" << endl;
      exit(EXIT_FAILURE);
    } else { // success!
      cout << "RootPlotter: Successfully opened file " << name << endl;
    }

  }
  cout << "-----------------------------------------------------------------------------------" << endl << endl;
  return;

}

void RootPlotter::SetSampleNames(TObjArray* SampleNames)
{
  // set the names of the chains which should be plotted
  
  fSampleNames = SampleNames;
  fNumOfSamples = fSampleNames->GetEntries();
}

void RootPlotter::SetSamplesToStack(TObjArray* names)
{
  // set the names of the chains which should be stacked
  // and plotted in one histogram
  
  fSamplesToStack = names;
  fNumOfSamplesToStack = fSamplesToStack->GetEntries();
}

void RootPlotter::SetSamplesWeight(TArrayF weights)
{
  // set the names of the chains which should be summed
  // and plotted in one histogram
  
  fSamplesWeight = weights;
}

void RootPlotter::SetHistColors(TArrayI colors)
{
  // set the histogram colors
  
  fSampleColors = colors;
}

void RootPlotter::SetHistMarkers(TArrayI markers)
{
  // set the histogram markers
  
  fSampleMarkers = markers;
}

TPostScript* RootPlotter::MakeNewPsFile(const char* psfilename)
{
  // create a new ps file with name psfilename
    TPostScript *PSFile = NULL;
  if (bPortrait){
    PSFile = new TPostScript(psfilename, 111);  // ps output
    PSFile->Range(20.0, 30.0);
  } else {
    PSFile = new TPostScript(psfilename, 112);  // ps output
    PSFile->Range(27.0, 18.0);
  }

  return PSFile;

}

Bool_t RootPlotter::ShouldBeStacked(const char* name)
{
  // check if the sample with name 'name' is in the list of 
  // samples to stack

  TString inname(name);
  for (Int_t i=0; i<fNumOfSamplesToStack; ++i){
    TString sname(fSamplesToStack->At(i)->GetName());
    // easy case:
    if (inname == sname) return true;

    // difficult case: name contains directory structure
    TObjArray* arr = inname.Tokenize("/");
    TString last;
    for (Int_t it=0; it<arr->GetEntries(); ++it){
      TString tok(arr->At(it)->GetName());
      if (tok.Contains(".root")) last = tok;
    }
    if (last.Contains(sname)) return true;
  }

  return false;

}


void RootPlotter::PlotHistos(const char* psfilename)
{

  // loops through all directories found in the 'data' directory
  // and writes the contents (histograms, hopefully) to a ps file
  // looks for every histogram for the corresponding MC reconstructed
  // histogram and plots it onto the same canvas

  static Int_t iPageNum = 0;
  Int_t CanWidth;
  Int_t CanHeight;
  if (bPortrait){
    CanWidth = 600;
    CanHeight = 830;
  } else {
    CanWidth =  800;
    CanHeight = 600;
  }
 
 
  // find all subdirectories (former histogram collections) in the first file
  TDirectory* firstdir = (TDirectory*) fRootfiles[0];
  TObjArray* dirs = FindSubdirs(firstdir);

  // general appearance and style
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);

  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);

  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");

  gStyle->UseCurrentStyle();


    
  // set up the canvas
  TCanvas *can = new TCanvas("canvas","Control Plots", CanWidth, CanHeight);

  Float_t yplot = 0.3;
  Float_t yratio = 0.15;


                                                //  coordinates:
                                                //  
  // set up the coordinates of the two pads:    //  y6 +-------------+
  Float_t y1, y2, y3, y4, y5, y6;               //     |             |
  y6 = 0.97;                                    //     |     pad1    |
  y5 = y6-yplot;                                //  y5 |-------------|
  y4 = y5-yratio;                               //     |     rp1     |
  y3 = 0.49;                                    //  y4 +-------------+
  y2 = y3-yplot;                                //  
  y1 = y2-yratio;                               //  y3 +-------------+
  Float_t x1, x2;                               //     |             |
  x1 = 0.01;                                    //     |     pad2    |
  x2 = 0.99;                                    //  y2 |-------------|
                                                //     |     rp2     |
                                                //  y1 +-------------+
                                                //     x1            x2



  TPad *pad1=0; TPad *rp1=0; TPad *pad2=0; TPad *rp2=0;
  if (bPlotRatio){
      
    pad1 = new TPad("pad1", "Control Plots 1", x1, y5, x2, y6);
    rp1 = new TPad("rp1", "Ratio1", x1, y4, x2, y5);
    
    pad2 = new TPad("pad2", "Control Plots 2", x1, y2, x2, y3);
    rp2 = new TPad("rp2", "Ratio2", x1, y1, x2, y2);
  
  } else {

    pad1 = new TPad("pad1", "Control Plots 1", x1, y4, x2, y6);
    pad2 = new TPad("pad2", "Control Plots 2", x1, y1, x2, y3);

  }
 
  // set margins for portrait mode
  if (bPortrait){

    pad1->SetTopMargin(0.07); pad1->SetBottomMargin(0.13);  pad1->SetLeftMargin(0.19); pad1->SetRightMargin(0.05);
    pad2->SetTopMargin(0.07); pad2->SetBottomMargin(0.13);  pad2->SetLeftMargin(0.19); pad2->SetRightMargin(0.05);
    
    if (bPlotRatio){
      pad1->SetTopMargin(0.07); pad1->SetBottomMargin(0.0);  pad1->SetLeftMargin(0.19); pad1->SetRightMargin(0.05);
      pad2->SetTopMargin(0.07); pad2->SetBottomMargin(0.0);  pad2->SetLeftMargin(0.19); pad2->SetRightMargin(0.05);
      rp1->SetTopMargin(0.0);    rp1->SetBottomMargin(0.33);  rp1->SetLeftMargin(0.19);  rp1->SetRightMargin(0.05);
      rp2->SetTopMargin(0.0);    rp2->SetBottomMargin(0.33);  rp2->SetLeftMargin(0.19);  rp2->SetRightMargin(0.05);
    }
	    
  // margins for landscape
  } else {

    pad1->SetTopMargin(0.02); pad1->SetBottomMargin(0.0);  pad1->SetLeftMargin(0.13); pad1->SetRightMargin(0.05);        
    pad2->SetTopMargin(0.02); pad2->SetBottomMargin(0.0);  pad2->SetLeftMargin(0.13); pad2->SetRightMargin(0.05);
    
    if (bPlotRatio){
      rp1->SetTopMargin(0.0);    rp1->SetBottomMargin(0.35);  rp1->SetLeftMargin(0.13);  rp1->SetRightMargin(0.05);
      rp2->SetTopMargin(0.0);    rp2->SetBottomMargin(0.35);  rp2->SetLeftMargin(0.13);  rp2->SetRightMargin(0.05);
    }
  }
  
  pad1->Draw(); 
  pad2->Draw();
  
  if (bPlotRatio){
    rp1->Draw();
    rp2->Draw();
  }
 
  Int_t ihist = 0; // histogram index
  TObjArray* OneDHistArray = new TObjArray();
  TH1* OneDHist = NULL;
  THStack* StackHist = NULL;

  TObjArray* TwoDHistArray = new TObjArray();
  TH2D* TwoDHist = NULL;

  TH1*  BkgdHist = NULL; // Background Histogram

  // array which histograms are already in the stacked one, so that they are not plotted twice
  Bool_t IsInStack[fNumOfSamples];
  for (Int_t i=0; i<fNumOfSamples; ++i) IsInStack[i] = false;
  
  // collector for cleaning up 
  TObjArray* Collector = new TObjArray();
  
  // the postscript file
  TPostScript *PSFile = NULL;

  // loop over all directories and plot the histograms
  for (Int_t i=0; i<dirs->GetEntries(); ++i){

    Int_t NumHistosPerObject = 0; // how many histograms per object should be plotted?
                                  // e.g. for a 2d-histo plot slices of it

    TString dirname(((TObjString*) dirs->At(i))->GetString());
    TString text(dirname.Data());

    text.Prepend("Plotting all histograms in directory ");

    cout << "\n+---------------------- histograms with ratios -----------------------+" << endl;
    cout <<   "| " << setw(60)<< text                                    << "        |" << endl;
    cout <<   "+---------------------------------------------------------------------+" << endl;
    ihist = 0;

    // create a new ps file
    if (PSFile){
      PSFile->Close();
      delete PSFile;
      PSFile = NULL;
    }
    iPageNum = 0;
    TString filename(psfilename);
    filename.ReplaceAll(".ps","");
    filename.Append("_");
    filename.Append(dirname);
    filename.Append(".ps");
    PSFile = MakeNewPsFile(filename);
    
    // make a page title
    can->cd();
    //TPaveText* pagetitle = new TPaveText(0.2,0.97,0.8,1.0, "NDC");
    //pagetitle->SetBorderSize(0);
    //pagetitle->SetFillColor(0);
    //pagetitle->AddText(((TObjString*) dirs->At(i))->GetString());
    //pagetitle->SetTextFont(42);
    //pagetitle->SetTextColor(kBlack);
    //pagetitle->SetTextSize(0.05);
    //pagetitle->Draw();

    firstdir->cd();
    gDirectory->Cd(((TObjString*) dirs->At(i))->GetString());

    Float_t TotMax = 1;
    
    // loop over all histograms in the directory and plot them
    TKey *key;
    TIter nextkey( gDirectory->GetListOfKeys() );
    while ( (key = (TKey*)nextkey())) {

      TotMax = 1;

      TObject *obj = key->ReadObj();

      Int_t Nhists = fNumOfSamples; 
      Double_t TotMax = 0.;

      if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {

	// first histogram found
	OneDHist = (TH1*) obj;
	TString histname(OneDHist->GetName());
	TString oldchainname(((TObjString*) fSampleNames->At(0))->GetName());

	if (OneDHist->GetMaximum()>TotMax){
	  TotMax = OneDHist->GetMaximum();
	}

	// cosmetics
	Cosmetics(OneDHist, 0);
	ApplyWeight(OneDHist, 0);

	if (OneDHist->GetMaximum()>TotMax) TotMax = OneDHist->GetMaximum();

	// if there is a lumi yield plot, plot only the first chain
	if (dirname.Contains("Yield")){
	  TString NewName = TString::Format("%s;1", OneDHist->GetName());
	  OneDHist = (TH1*) gDirectory->Get(NewName);
	  PSFile->NewPage();
	  PlotYields(OneDHist);
	  OneDHistArray->Clear();
	  break;
	}

	// build and substract background if requested
	if (bSubstractBkgd){
	  BkgdHist = BuildBackground(dirname, histname);
	  Collector->Add(BkgdHist);
	  OneDHist->Add(BkgdHist, -1);
	}	       

	if (ShouldBeStacked(firstdir->GetName())){
	  if (StackHist==NULL){
	    StackHist = new THStack(histname, "stack");
	  } 
	  Nhists--;
	  StackHist->Add(OneDHist);
	} else { 
	  OneDHistArray->Add(OneDHist);
	}
	
	if (obj->IsA()->InheritsFrom( TH2D::Class() ) ) {
	  TwoDHist = (TH2D*) OneDHist;
	  TwoDHist->SetTitle(firstdir->GetName());
	  TwoDHistArray->Add(TwoDHist);
	}


	// find the same histogram for the other samples
	for (Int_t ichain=1; ichain<fNumOfSamples; ++ichain){

	  TFile* thisdir = fRootfiles[ichain];	  
	  thisdir->cd();
	  gDirectory->Cd(((TObjString*) dirs->At(i))->GetString());

 
	  // get the right histogram
	  OneDHist = NULL;
	  TObject* obj = gDirectory->Get(histname);
	  if (!obj){
	    cerr << "Could not find histogram with name " << histname << " in file " << thisdir->GetName() 
		 << " and directory " << ((TObjString*) dirs->At(i))->GetString() << " - please check." << endl;
	    exit(EXIT_FAILURE);
	  }

	  if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
	    OneDHist = (TH1*) obj;
	    Cosmetics(OneDHist, ichain);
	    ApplyWeight(OneDHist, ichain);
	    if (OneDHist->GetMaximum()>TotMax) TotMax = OneDHist->GetMaximum();

	  } else {
	    cerr << "Found an object with name " << histname << " in file " << thisdir->GetName() 
		 << " and directory " << ((TObjString*) dirs->At(i))->GetString() << "." << endl;
	    cerr << "But it's not a TH1! Please correct the error." << endl;
	    exit(EXIT_FAILURE);
	  }


	  // if the requested sample name is in the list of samples to stack, add it to the stack
	  // chains to create on histogram
	  TString chainname( thisdir->GetName() );
	  if (ShouldBeStacked(chainname)){
	    if (StackHist==NULL){
	      StackHist = new THStack(histname, "stack");
	    }
	    Nhists--;
	    StackHist->Add(OneDHist);
	  } else {	    
	    OneDHistArray->Add(OneDHist);
	  }

	  if (obj->IsA()->InheritsFrom( "TH2D" ) ) {
	    TwoDHist = (TH2D*) OneDHist;
	    TwoDHist->SetTitle(chainname);
	    TwoDHistArray->Add(TwoDHist);
	  }
	}

	NumHistosPerObject = 1;

	
	// check if the 2d-histograms need 'special care'	
	if (obj->IsA()->InheritsFrom( "TH2D" ) ) {
	  if (histname.Contains("PtFraction")) {
	    OneDHistArray->Clear();
	    for (Int_t j=0; j<TwoDHistArray->GetEntries(); ++j){
	      TH2D* hist = (TH2D*) TwoDHistArray->At(j);
	      TH1* fractionhist = MakePtFractionHisto(hist);
	      OneDHistArray->Add(fractionhist);
	    }
	    NumHistosPerObject = 1;
	    //	  } else if ( histname.Contains("JetShape") ){
	    //	    NumHistosPerObject = OneDHist->GetNbinsX()+1;
	    //	    IsJetShape = true;    
	  } else {
	    NumHistosPerObject = OneDHist->GetNbinsX();
	  }	  
	}
      }  else { // do nothing for an unknown object
	continue;
      }
      

      // do nothing for empty histograms
      /*
      cout << "checkin entries..." << endl;
      if (OneDHist->GetEntries() == 0) {
	cout << "Entries are 0! Histname = " << OneDHist->GetName() << endl;
	TwoDHistArray->Clear();
	OneDHistArray->Clear();
	continue; 
      }
      */

      // skip the ptbalance calibration histograms -> plot them separately
      TString histname(OneDHist->GetName());
      if ( histname.Contains("PtBalanceFit")  ){
	TwoDHistArray->Clear();
	OneDHistArray->Clear();
	continue;
      }

      // skip denominators of efficiency histograms
      if ( histname.Contains("EffBot")  ){
	TwoDHistArray->Clear();
	OneDHistArray->Clear();
	continue;
      }


      // loop over the histograms that should be plotted per canvas
      for (Int_t jobjhist=0; jobjhist<NumHistosPerObject; ++jobjhist){
	
	//	if( IsJetShape && !bJetShapesPerSlice) jobjhist = NumHistosPerObject-1; // exception for the jetshapes
	
	// get 1-dimensional histograms out of the 2d ones
	if (NumHistosPerObject>1){
	    
	  // reset the 1-d hists, to be filled with slices
	  OneDHistArray->Clear();
	  if (StackHist){
	    delete StackHist;
	    StackHist=NULL;
	  }
	  Nhists=fNumOfSamples;

	  // create the slices
	  for (Int_t j=0; j<TwoDHistArray->GetEntries(); ++j){
	    TH2D* hist = (TH2D*) TwoDHistArray->At(j);
	    TString hname = hist->GetName();
	      
	    OneDHist = GetSliceHisto(hist, jobjhist+1);
	    if (OneDHist->GetMaximum()>TotMax) TotMax = OneDHist->GetMaximum();

	    TString htitle = hist->GetTitle();
	    if (ShouldBeStacked(htitle)){
	      if (StackHist==NULL){
		StackHist = new THStack(histname, "stack");
	      }
	      Nhists--;
	      StackHist->Add(OneDHist);
	    } else {
	      //Collector->Add(OneDHist);
	      OneDHistArray->Add(OneDHist);
	    }	  
	    	    
	    if (Nhists != TwoDHistArray->GetEntries()){
	      cout << "Could not convert 2d histogram into 1d one. HistName = " << ((TH2D*) TwoDHistArray->At(0))->GetName() << "  Aborting." << endl;
	      exit(EXIT_FAILURE);
	    }
	  }
	}

	if ((Nhists+fNumOfSamplesToStack) != fNumOfSamples){
	  cout << "For histogram " << ((TH2D*) TwoDHistArray->At(0))->GetName() 
	       << ", found " << OneDHistArray->GetEntries() << " different samples, requested were " << fNumOfSamples << endl;
	  cout << "Skipping histogram." << endl;
	  continue;
	}
	
	if (StackHist){

	  //invert the ordering of all histograms in the stack
	  THStack *StackHist2 = new THStack("stack2","stack2");
	  TList* histlist = StackHist->GetHists();
	  for(int i=fNumOfSamplesToStack-1; i>=0; --i){
	    StackHist2->Add( (TH1*)histlist->At(i));
	  }
	  
	  StackHist = NULL;
	  StackHist = StackHist2;

	  if (StackHist->GetMaximum()>TotMax) TotMax = StackHist->GetMaximum();
	}

	TH1* FirstHist = (TH1*) OneDHistArray->At(0);
	TString histname(FirstHist->GetName());
	
	cout << "ihist = " << ihist << " histname = " << histname << endl;

	// set up the canvas
	if (ihist%2==0){

	  if (ihist != 0){
	    ++iPageNum;
	    DrawPageNum(can, iPageNum);	
	  }	

	  can->Update();
	  PSFile->NewPage();

	  // clean up	 
	  if (ihist!=0){
	    can->Clear("D");
	    for (Int_t i=0;i<Collector->GetEntries();++i){
	      TObject* obj = Collector->At(i);
	      delete obj;
	    }
	    Collector->Clear();
	  }	  

	}

	switch (ihist%2){
	case 0:  pad1->cd(); break;
	case 1:  pad2->cd(); break;
	default: break;	
	}
     
	gPad->SetLogx(0);
	gPad->SetLogy(0);

	/*------------------------------------------------------------------\
	| DRAW THE HISTOGRAMS AS WELL AS THE RATIOS                         |
	\------------------------------------------------------------------*/

	if (StackHist){
	  TList* list = StackHist->GetHists();
	  TH1* hist = (TH1*) list->At(list->GetEntries()-1);
	  if (hist->GetMaximum()>TotMax){
	    TotMax = hist->GetMaximum();	
	  }
	}


	//////////////////////////////////////
	// draw the histograms              //
	//////////////////////////////////////
	      	
	Int_t NOneDHists = OneDHistArray->GetEntries();

	if (bShapeNorm){
	  ShapeNormalize(FirstHist);
	  for (Int_t jhist=1; jhist<NOneDHists; ++jhist){
	    ShapeNormalize((TH1*) OneDHistArray->At(jhist));	    
	  }
	  // shape normalisation of the stack (overall factor for the sum)
	  if (StackHist){
	    ShapeNormalize(StackHist);	    
	  }
	}	  

	// if the first hist has no entries, then take the second one
	if (FirstHist->GetEntries()==0){
	  FirstHist = (TH1*) OneDHistArray->At(1);
	}
	
	// store the subtitle if it was set in GetSliceHistos
	TString oldtitle = FirstHist->GetTitle();
	TString subtitle;
	Bool_t PlotSubtitle = false;
	if (oldtitle.Contains("SUBTITLE") ){
	  subtitle = oldtitle;
	  subtitle.ReplaceAll("SUBTITLE:","");
	  PlotSubtitle = true;
	}

	// set the title, logscales and maxima
	Ssiz_t pos2 = oldtitle.First("(");
	TString title(oldtitle(0,pos2));
	FirstHist->SetTitle(title);
	histname.ReplaceAll("(",", "); histname.ReplaceAll("_Data",""); histname.ReplaceAll(")","");
	FirstHist->SetTitle("");
	//FirstHist->SetTitle(histname);

	Float_t MaxScale = 1.2;
	Float_t LogScale = 20.0;
	Bool_t  IsLogPlot = false;
	Bool_t IsEffHist = false;
	TString CheckEff(FirstHist->GetYaxis()->GetTitle());
	if (CheckEff.Contains("#epsilon")) IsEffHist=true;

	if (bDrawEntries){
	  MaxScale += fNumOfSamples / 10.;
	  LogScale *= fNumOfSamples;
	}

	if (histname.EndsWith("_lxy")){
	  IsLogPlot = true;
	  gPad->SetLogx(1); 
	  gPad->SetLogy(1);
	  if (TotMax>0){
	    FirstHist->SetMaximum(LogScale*TotMax);
	  }
	  if (!bShapeNorm){
	    //FirstHist->SetMinimum( 0.3 );
	  }
	} else if (histname.EndsWith("_ly")) {
	  IsLogPlot = true;
	  gPad->SetLogy(1);
	  if (TotMax>0){
	    FirstHist->SetMaximum(LogScale*TotMax);
	  }
 	  if (!bShapeNorm){
	    //FirstHist->SetMinimum( 0.3 );
	  }

	// no log-y
	} else { 
	  
	  if (histname.EndsWith("_lx")){
	    gPad->SetLogx(1);
	  }
	  
	  FirstHist->SetMinimum(0.001);
	  if (TotMax>0){
	    if (FirstHist->GetMaximum() != 0){
	      FirstHist->SetMaximum(MaxScale*TotMax);
	    } else {
	      FirstHist->SetMaximum(TotMax);
	    }
	  }
	}

	// set range for resolution plots
	if (FirstHist->InheritsFrom(TProfile::Class())){
	  FirstHist->SetMinimum(-0.2); 
	  FirstHist->SetMaximum(0.2);
	}
	
	// efficiency with efficiency close to 1: redefine minimum
	if (IsEffHist){
	  if (FirstHist->GetMaximum() > 1.1){
	    FirstHist->SetMinimum(0.76);
	    FirstHist->SetMaximum(1.23);
	  }
	}

	// draw the histograms
	if (FirstHist->GetMarkerStyle() < 2){
	  FirstHist->Draw("HIST");
	} else {
	  FirstHist->Draw("P");
	}
	
	// draw background if requested
	/*
	if (bSubstractBkgd && BkgdHist){
	  BkgdHist->SetLineColor(kGray);
	  BkgdHist->SetFillColor(kGray);
	  BkgdHist->SetFillStyle(3002);
	  BkgdHist->Draw("Histsame");
	}
	*/	

	// now draw the stack if it exists
	if (StackHist){
	  StackHist->Draw("hist,same");
	}

	for (Int_t mhist=1; mhist<OneDHistArray->GetEntries(); ++mhist){
	  TH1* hist = (TH1*) OneDHistArray->At(mhist);

	  // check entries in case of log plot
	  if (IsLogPlot){
	    if (hist->GetEntries()<1) continue;
	  }
	  
	  if (hist->GetMarkerStyle() == 0){
	    hist->Draw("HISTsame");
	  } else {
	    hist->Draw("P same");
	  }
	}

	if (FirstHist->GetMarkerStyle() < 2){
	  FirstHist->Draw("HISTsame");
	} else {
	  FirstHist->Draw("P same");
	}


	/*        	
	TPaveText *p = (TPaveText*) (gPad->GetListOfPrimitives()->FindObject("title"));
	if (p){
	  p->SetLineColor(kWhite);
	  p->Draw("same");
	}
	*/

	gPad->RedrawAxis();


	////////////////////// 
	// draw a legend    //
	//////////////////////
	Float_t yfrac = 0.06;
	if (!bPlotRatio) yfrac = 0.05;
	Float_t top = 0.92;
	Float_t ysize = yfrac*fNumOfSamples;
	Float_t xleft = 0.65;
	Float_t xright = 0.92;
	if (!bPortrait){
	  top = 0.99;
	  ysize = 0.07*fNumOfSamples;
	  xleft = 0.72;
	  xright = 0.96;
	}
	
	TLegend *leg = new TLegend(xleft,top-ysize,xright,top, NULL,"brNDC");
	Collector->Add(leg);
	leg->SetFillColor(0);
	leg->SetLineColor(1);
	leg->SetBorderSize(0);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	for (Int_t i=0; i<fNumOfSamples; ++i){
	  TString legname = TString::Format("leg_entry_%i",i);
	  TString legtitle = ((TObjString*) fSampleNames->At(i))->GetName(); 
	  legtitle.ReplaceAll("SPACE", " ");
	  legtitle.ReplaceAll("[", "{");
	  legtitle.ReplaceAll("]", "}");
	  TLegendEntry* entry = NULL;
	  if(fSampleMarkers.At(i)!=0) {
	    entry = leg->AddEntry(legname, legtitle, "lpf");
	    entry->SetLineColor(fSampleColors.At(i));
	  } else if (fSampleMarkers.At(i)==0) {
	    entry = leg->AddEntry(legname, legtitle, "f");
	  } else if (fSampleMarkers.At(i)<0) {
	    entry = leg->AddEntry(legname, legtitle, "l");
	  }
	  entry->SetLineWidth(1);
	  if (fSampleMarkers.At(i)>0){
	    entry->SetMarkerColor(fSampleColors.At(i));
	    entry->SetMarkerStyle(fSampleMarkers.At(i));
	    entry->SetMarkerSize(0.8);
	    //entry->SetLineWidth(2);
	  } else {
	    entry->SetMarkerStyle(0);
	    entry->SetMarkerSize(0);
	    entry->SetMarkerColor(fSampleColors.At(i));
	    //entry->SetLineWidth(2);

	    // only line
	    if (fSampleMarkers.At(i)<0){
	      entry->SetLineWidth(2);
	      if (fSampleMarkers.At(i)==-1) entry->SetLineStyle(kSolid);
	      if (fSampleMarkers.At(i)==-2) entry->SetLineStyle(kDotted);
	      if (fSampleMarkers.At(i)==-3) entry->SetLineStyle(kDashDotted);
	      if (fSampleMarkers.At(i)==-4) entry->SetLineStyle(kDashDotted);
	    // fill
	    } else {
	      entry->SetFillColor(fSampleColors.At(i));
	      entry->SetFillStyle(1001);
	    }
	  }
	  entry->SetTextAlign(12);
	  //entry->SetTextColor(fSampleColors.At(i));
	}
	leg->Draw();
        
	// subtitle, if given
	if ( PlotSubtitle ){
	  TPaveText* text;
	  if (bDrawEntries){
	    text = new TPaveText(xleft,top-ysize-0.1, xright-0.03, top-ysize-0.03, "NDC");
	  } else {
	    text = new TPaveText(0.23, 0.88, xleft-0.03, 0.94, "NDC");
	  }
	  Collector->Add(text);
	  text->SetBorderSize(0);
	  text->SetFillColor(0);
	  text->AddText(subtitle);
	  text->Draw();
	  text->Paint();
	}


	///////////////////
	// draw entries
	//////////////////
	if (bDrawEntries){
	  
	  TPaveText* NumEntries = new TPaveText(xleft-0.4, top-ysize, xleft-0.02, top, "brNDC");
	  Collector->Add(NumEntries);
	  NumEntries->SetFillColor(kWhite);
	  NumEntries->SetBorderSize(2);
	  NumEntries->SetShadowColor(kWhite);
	  // bug: (TODO!)
	  // use OneDHistArray->GetEntries() instead of fNumOfSamples
	  // question: how to get the numbers from the stack?
	  for (Int_t jhist=0; jhist<fNumOfSamples; ++jhist){
	    TH1* hist = (TH1*) OneDHistArray->At(jhist);
	    // number of entries
	    Double_t HistEntries;
	    if (!bShapeNorm){
	      HistEntries = hist->GetSumOfWeights();
	    } else {
	      HistEntries = hist->GetEntries();
	    }
	    Int_t Num = static_cast<Int_t>(HistEntries);
	    Float_t HistMean = static_cast<Float_t>(hist->GetMean());
	    TString entlegtitle;
	    entlegtitle.Form("entries: %i  mean: %5.3f", Num, HistMean);
	    TText* text = NumEntries->AddText(entlegtitle);
	    text->SetTextSize(0.03);
	    text->SetTextAlign(12);
	    text->SetTextColor(fSampleColors.At(jhist));
	  }
	  NumEntries->Draw();
	}

	if(bDrawLumi){
	  TLatex *text1 = new TLatex(3.570061,23.08044,"CMS Preliminary");
	  text1->SetNDC();
	  text1->SetTextAlign(13);
	  text1->SetX(0.22);
	  text1->SetY(0.918);
	  text1->SetTextFont(42);
	  text1->SetTextSizePixels(24);
	  text1->Draw();
	  
	  TLatex *text2 = new TLatex(3.570061,23.08044,"5.2 fb^{-1} at #sqrt{s} = 8 TeV");
	  text2->SetNDC();
	  text2->SetTextAlign(13);
	  text2->SetX(0.22);
	  text2->SetY(0.87);
	  text2->SetTextFont(42);
	  text2->SetTextSizePixels(24);
	  text2->Draw();
	}

	////////////////////////
	// draw the ratio     //
	////////////////////////
	if (bPlotRatio){

	  switch (ihist%2){
	  case 0:  rp1->cd(); break;
	  case 1:  rp2->cd(); break;
	  default: break;
	  }
	  
	  gPad->SetLogx(0);
	  gPad->SetLogy(0);
	
	  // set up ratio 
	  TObjArray* RatioHistArray = new TObjArray();

	  // special treatment for TProfiles
	  TH1D* tempRatioHist = NULL;

	  for ( Int_t iChain=1; iChain<NOneDHists; iChain++ ){

	    if (FirstHist->InheritsFrom(TProfile::Class())){
	      tempRatioHist = ((TProfile*)FirstHist)->ProjectionX();
	      Collector->Add(tempRatioHist);
	      TH1D* denom = (TProfile*) OneDHistArray->At(iChain);
	    
	      tempRatioHist->Divide( denom );	    
	      tempRatioHist->GetYaxis()->SetRangeUser(0.3, 1.7);
	      CopyStyle((*tempRatioHist), denom);
	      RatioHistArray->Add(tempRatioHist);

	    } else {
	      tempRatioHist = (TH1D*) FirstHist->Clone();
	      Collector->Add(tempRatioHist);
	      TH1D* denom = (TH1D*) OneDHistArray->At(iChain);

	      // divide
	      tempRatioHist->Divide(denom);
	      tempRatioHist->GetYaxis()->SetRangeUser(0.3, 1.7);
	      CopyStyle((*tempRatioHist), denom);
	      RatioHistArray->Add(tempRatioHist);
	    
	      if (IsEffHist){
		tempRatioHist->GetYaxis()->SetRangeUser(0.955, 1.045);
		tempRatioHist->GetYaxis()->SetNdivisions(505);
	      }

	    }
	  }

	  Int_t NRatioPlots = NOneDHists;

	  // divide by the stack
	  if (StackHist){

	    tempRatioHist = (TH1D*) FirstHist->Clone();
	    Collector->Add(tempRatioHist);
	  
	    TObjArray* arr = StackHist->GetStack();
	    // the last element is the sum
	    TH1D* denom = (TH1D*) arr->At(arr->GetEntries()-1);
	  
	    // first: one histogram for the MC statistical error
	    TH1D* MCstat = (TH1D*) FirstHist->Clone();
	    Collector->Add(MCstat);

	    for (Int_t ibin=1;ibin<denom->GetNbinsX()+1; ++ibin){
	      Double_t val = denom->GetBinContent(ibin);
	      Double_t err = denom->GetBinError(ibin);
	      MCstat->SetBinContent(ibin,  1.0);
	      MCstat->SetBinError(ibin,  err/val);
	    }
	    MCstat->SetMarkerStyle(0);
	    MCstat->SetMarkerSize(0);
	    MCstat->SetLineColor(kGray+1);
	    MCstat->SetFillColor(kGray+1);
	    RatioHistArray->Add(MCstat);
	    ++NRatioPlots;

	    // now divide and add the result
	    tempRatioHist->Divide(denom);
	    tempRatioHist->GetYaxis()->SetRangeUser(0.3, 1.7);
	    Cosmetics(tempRatioHist, 0);
	    tempRatioHist->SetMarkerSize(0.7);
	    RatioHistArray->Add(tempRatioHist);
	    ++NRatioPlots;

	  }
	  //exit(3);


	  // position of the line
	  Float_t LinePos = 1.0;

	  tempRatioHist = (TH1D*)RatioHistArray->At(0);
	
	  TString ratiotitle = ((TObjString*) fSampleNames->At(0))->GetName();
	  ratiotitle.Append(" / MC");	  
	  ((TH1D*)RatioHistArray->At(0))->GetYaxis()->SetTitle(ratiotitle);
 	
	  if (histname.EndsWith("_lxy") || histname.EndsWith("_lx")){
	    gPad->SetLogx(1);
	  }


	  // cosmetics for portrait mode 
	  if (bPortrait){
	    tempRatioHist->SetTitle("");
	    tempRatioHist->SetTitleOffset(1.1, "X");
	    tempRatioHist->SetTitleOffset(0.43, "Y");
	    tempRatioHist->SetLabelOffset(0.02, "X");
	    tempRatioHist->SetLabelOffset(0.02, "Y");
	  
	    tempRatioHist->GetXaxis()->SetLabelSize(0.13);
	    tempRatioHist->GetXaxis()->SetTickLength(0.07);
	    tempRatioHist->GetXaxis()->SetTitleSize(0.13);
	    tempRatioHist->GetXaxis()->SetTitleOffset(1.3);
	    tempRatioHist->GetXaxis()->SetLabelFont(42);
	    tempRatioHist->GetXaxis()->SetTitleFont(42);
	  
	    tempRatioHist->GetYaxis()->CenterTitle();
	    tempRatioHist->GetYaxis()->SetTitleSize(0.13);
	    tempRatioHist->GetYaxis()->SetLabelSize(0.12);
	    tempRatioHist->GetYaxis()->SetNdivisions(210);
	    tempRatioHist->GetYaxis()->SetTickLength(0.02);
	    tempRatioHist->GetYaxis()->SetLabelFont(42);
	    tempRatioHist->GetYaxis()->SetLabelOffset(0.011);

	    if (histname.Contains("ElecCalib")){
	      tempRatioHist->GetYaxis()->SetNdivisions(205);
	    }	  

	    // cosmetics for landscape mode 
	  } else {

	    tempRatioHist->SetTitle("");
	    tempRatioHist->SetTitleOffset(1.1, "X");
	    tempRatioHist->SetTitleOffset(0.5, "Y");
	    tempRatioHist->SetLabelOffset(0.02, "X");
	    tempRatioHist->SetLabelOffset(0.01, "Y");
	  
	    tempRatioHist->GetXaxis()->SetLabelSize(0.14);
	    tempRatioHist->GetXaxis()->SetTickLength(0.07);
	    tempRatioHist->GetXaxis()->SetTitleSize(0.15);
	    tempRatioHist->GetXaxis()->SetLabelFont(42);
	  
	    tempRatioHist->GetYaxis()->CenterTitle();
	    tempRatioHist->GetYaxis()->SetTitleSize(0.11);
	    tempRatioHist->GetYaxis()->SetLabelSize(0.12);
	    tempRatioHist->GetYaxis()->SetNdivisions(505);
	    tempRatioHist->GetYaxis()->SetTickLength(0.03);
	    tempRatioHist->GetYaxis()->SetLabelFont(42);
	  
	  }


	  Double_t FirstBinEdge = tempRatioHist->GetXaxis()->GetXmin();
	  Double_t LastBinEdge = tempRatioHist->GetXaxis()->GetXmax();
	  TLine* line = new TLine(FirstBinEdge,LinePos,LastBinEdge,LinePos);
	  line->SetLineColor(kGray+2);
	  Collector->Add(line);

	  // draw them
	  if (FirstHist->GetEntries() != 0){
	    tempRatioHist->Draw("P");
	    line->Draw("same");

	    // more lines for the electron calibration plot
	    if (histname.Contains("ElecCalib")){
	      TLine* line1 = new TLine(FirstBinEdge,1.005,LastBinEdge,1.005);
	      line1->SetLineStyle(kDotted);
	      Collector->Add(line1);
	      TLine* line2 = new TLine(FirstBinEdge,0.995,LastBinEdge,0.995);
	      line2->SetLineStyle(kDotted);
	      Collector->Add(line2);
	      TLine* line3 = new TLine(FirstBinEdge,1.01,LastBinEdge,1.01);
	      line3->SetLineStyle(kDashed);
	      Collector->Add(line3);
	      TLine* line4 = new TLine(FirstBinEdge,0.99,LastBinEdge,0.99);
	      line4->SetLineStyle(kDashed);
	      Collector->Add(line4);
	      line1->Draw("same");
	      line2->Draw("same");
	      line3->Draw("same");
	      line4->Draw("same");
	    }

	    if (IsEffHist){
	      TLine* line3 = new TLine(FirstBinEdge,1.01,LastBinEdge,1.01);
	      line3->SetLineColor(kGreen);
	      line3->SetLineStyle(kDashed);
	      Collector->Add(line3);
	      TLine* line4 = new TLine(FirstBinEdge,0.99,LastBinEdge,0.99);
	      line4->SetLineColor(kGreen);
	      line4->SetLineStyle(kDashed);
	      Collector->Add(line4);
	      line3->Draw("same");
	      line4->Draw("same");
	    }


	    for ( Int_t iChain = 1; iChain<NRatioPlots ; iChain++ ){
	      tempRatioHist = (TH1D*)RatioHistArray->At(iChain-1);
	      if (tempRatioHist->GetMarkerStyle()>1){
		tempRatioHist->Draw("Psame");
	      } else {
		if (tempRatioHist->GetFillColor()<2){
		  tempRatioHist->Draw("HISTsame");
		} else {
		  tempRatioHist->DrawCopy("E2same");
		  tempRatioHist->SetFillColor(0);
		  tempRatioHist->DrawCopy("HISTsame");
		}
	      }
	    }
  
	  }
	  gPad->RedrawAxis();

	  // clean up
	  RatioHistArray->Clear();
	  delete RatioHistArray;

	} // end "if (bRatioPlot)
	
	++ihist;
	
	// clean up
	OneDHistArray->Clear();
	StackHist=NULL;
	//delete BkgdHist;

      } // loop over histograms that should be plotted per object

      TwoDHistArray->Clear();

    } // loop over histograms


  } // loop over directories


  ///////////////////////////////////////
  // PLOT THE PT-BALANCE HISTOGRAMS    //
  ///////////////////////////////////////

  /*
  if (bFitPtBalanceHists){

    cout << "\n+---------------------- calibration histograms -----------------------+" << endl;
    cout <<   "|              Plotting Pt-balance histograms with fits               |" << endl;
    cout <<   "+---------------------------------------------------------------------+" << endl;
     
    // create a new ps file
    if (PSFile){
      PSFile->Close();
      delete PSFile;
      PSFile = NULL;
    }
    iPageNum = 0;
    TString filename(psfilename);
    filename.ReplaceAll(".ps","");
    filename.Append("_PtBalance");
    filename.Append(".ps");
    PSFile = MakeNewPsFile(filename);

    // get all ptbalance plots and hand them over to the fitting/plotting program
    TObjArray* PtBalanceHists = GetPtBalanceHists();
    
    FitCalibHistos fitter;
    fitter.SetLegEntries(fLegEntries);
    fitter.SetHistColors(fChainColors);
    fitter.SetHistMarkers(fChainMarkers);
    fitter.SetPortraitMode(bPortrait);
    //fitter.SetFitFunction(FitCalibHistos::GAUSS);
    //fitter.SetFitFunction(FitCalibHistos::ASYMMGAUSS);
    fitter.SetFitFunction(FitCalibHistos::STUDENT);
    fitter.FitAndPlot(PSFile, PtBalanceHists);
  }
  */

  // clean up and finish
  PSFile->Close();
  delete PSFile;
  delete OneDHistArray;
  delete TwoDHistArray;

  TString outtext(psfilename);
  outtext.ReplaceAll(".ps","_");
  outtext.Append("<HistManager>.ps");
  
  cout << endl << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl;
  cout << "<> Wrote all histograms to ps-file " << outtext << endl;
  cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endl << endl;

}

void RootPlotter::ApplyWeight(TH1* hist, Int_t isample)
{
  // apply weight as given in the array fSamplesWeight

  if (isample >= fNumOfSamples){
    cerr << "\n\n Cannot apply weight for histogram index " << isample << ", because only " << fNumOfSamples << " samples are given." << endl;
    return;
  }

  hist->Scale(fSamplesWeight.At(isample));

  return;
}

void RootPlotter::CopyStyle(TH1& h1, TH1* h2)
{
  // copy the style from hist2 to hist1
  h1.SetMarkerStyle(h2->GetMarkerStyle());
  h1.SetMarkerSize(h2->GetMarkerSize());
  h1.SetMarkerColor(h2->GetMarkerColor());
  h1.SetLineWidth(h2->GetLineWidth());   
  h1.SetLineStyle(h2->GetLineStyle());
  h1.SetLineColor(h2->GetLineColor());
  h1.SetFillColor(h2->GetFillColor());
  return;
}

void RootPlotter::Cosmetics(TH1* hist, Int_t isample)
{
  // apply cosmetics

  if (isample >= fNumOfSamples){
    cerr << "\n\n Cannot set cosmetics for histogram index " << isample << ", because only " << fNumOfSamples << " samples are given." << endl;
    return;
  }
  hist->SetLineColor(fSampleColors.At(isample));

  
  if (fSampleMarkers.At(isample) > 1 ){
    hist->SetMarkerStyle(fSampleMarkers.At(isample));
    hist->SetMarkerColor(fSampleColors.At(isample));
    hist->SetMarkerSize(0.8);
  } else {
    hist->SetMarkerStyle(0);
    hist->SetMarkerSize(0);
    hist->SetMarkerColor(fSampleColors.At(isample));
    hist->SetLineWidth(2);   
    //    if(fSampleMarkers.At(isample) >=0 ){
    //      hist->SetLineColor(fSampleColors.At(isample));
    //      hist->SetLineWidth(1); 
    //    }
  }

  // histogram is transparent if marker < 0  
  if (fSampleMarkers.At(isample) < 0 ){
    // change line style
    if (fSampleMarkers.At(isample)==-1) hist->SetLineStyle(kSolid);
    if (fSampleMarkers.At(isample)==-2) hist->SetLineStyle(kDotted);
    if (fSampleMarkers.At(isample)==-3) hist->SetLineStyle(kDashDotted);
    if (fSampleMarkers.At(isample)==-4) hist->SetLineStyle(kDashDotted);

  } else {
    hist->SetFillColor(fSampleColors.At(isample));
  }

  // set Y-axis title
  hist->GetYaxis()->SetTitle("Entries");  
  
  // set X-axis title
  hist->GetXaxis()->SetTitle(hist->GetTitle()); 

  if (bShapeNorm) {
    hist->GetYaxis()->SetTitle("#DeltaN/N");
  }

  
  
  // portrait mode
  if (bPortrait){
    
    hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetLabelSize(0.06);
    //hist->GetYaxis()->SetNdivisions(1005);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelOffset(0.01);
    
    hist->GetYaxis()->SetTitleOffset(0.8);
    hist->GetYaxis()->SetTickLength(0.02);
    
    // move the y-axis title depending on the number of entries
    //if (hist->GetMaximum() < 0.05) hist->GetYaxis()->SetTitleOffset(0.9);	
    //if (hist->GetMaximum() > 10000 && hist->GetMaximum() < 100000) hist->GetYaxis()->SetTitleOffset(1.1);	
    
    
    // landscape mode
  } 
  else {
    
    hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetLabelSize(0.06);
    //hist->GetYaxis()->SetNdivisions(1005);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelOffset(0.01);
    
    hist->GetYaxis()->SetTitleOffset(0.95);
  }
  
  if (!bPlotRatio){
    hist->SetTitle("");
    
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetLabelOffset(0.008);
    hist->GetXaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleFont(42);
    
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.045);
    //hist->GetYaxis()->SetNdivisions(210);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelOffset(0.011);
  }
  

  return;

}

void RootPlotter::ShapeNormalize(TH1* hist)
{
  // normalize the histogram according to its area:
  // Integral hist dx = 1
  
  Int_t beg, end;
  Float_t area;

  beg=1; // integrate from bin 1 
  end=hist->GetNbinsX();
  area=hist->Integral(beg,end);
  if (area>0){
    hist->Scale(1./area);
  }
  
}

void RootPlotter::ShapeNormalize(THStack* stack)
{
  // normalize the stack according to its area:
  // Integral stack dx = 1, meaning that all 
  // histograms in the stack get an overall normalisation 
  // factor
  
  Int_t beg, end;
  Float_t area = 0.;

  TList* list = stack->GetHists();
  for (Int_t i=0; i<list->GetEntries();++i){
    TH1* hist = (TH1*) list->At(i);
    beg=1; // integrate from bin 1 
    end=hist->GetNbinsX();
    area += hist->Integral(beg,end);
  }

  // now normalise
  for (Int_t i=0; i<list->GetEntries();++i){
    TH1* hist = (TH1*) list->At(i);
    if (area>0){
      hist->Scale(1./area);
    }
  }
  return;
}


TH1* RootPlotter::MakePtFractionHisto(TH2D* In2DHist)
{
  // make a 'Pt fraction' histogram out of the 2d-histogram
  // this is done via averaging over the y-axis

  if (In2DHist==NULL)
    return NULL;

  // for every x-bin calculate the average pt-fraction
  TString histname = In2DHist->GetName();
  histname.Append("_1d");
  TH1D* finalhist = In2DHist->ProjectionX(histname, 0, -1, "e");

  for (Int_t ibin=1; ibin<In2DHist->GetNbinsX()+1; ++ibin){
    TH1D* slice = In2DHist->ProjectionY("temp", ibin, ibin, "e");

    Double_t meanfraction = slice->GetMean();
    Double_t meanerror = slice->GetMeanError();

    finalhist->SetBinContent(ibin, meanfraction);
    finalhist->SetBinError(ibin, meanerror);
    delete slice;
  }

  TString ytitle( In2DHist->GetYaxis()->GetTitle() );
  finalhist->GetYaxis()->SetTitle( ytitle );

  return finalhist;

}

TH1* RootPlotter::MakeJetShapeMeanHist(TH2D* In2DHist){

  // exception for empty histogram
  if (In2DHist==NULL)  return NULL;
  
  // calculate the mean jet shape as a function of the radius (bin)
  TString histname = In2DHist->GetName();
  histname.Append("_1d");
  TH1D* finalhist = In2DHist->ProjectionX(histname, 0, -1, "e");

  for (Int_t ibin=1; ibin<In2DHist->GetNbinsX()+1; ++ibin){
    TH1D* slice = In2DHist->ProjectionY("temp", ibin, ibin, "e");

    Double_t meanfraction = slice->GetMean();
    Double_t meanerror = slice->GetMeanError();

    finalhist->SetBinContent(ibin, meanfraction);
    finalhist->SetBinError(ibin, meanerror);
    delete slice;
  }

  TString ytitle( In2DHist->GetYaxis()->GetTitle() );
  finalhist->GetYaxis()->SetTitle( ytitle );

  return finalhist;

}

TH1* RootPlotter::MakeMatchingEfficiencyHisto(TH2D* In2DHist)
{
  // make a 'jet matching efficiency' histogram
  
  if (In2DHist == NULL)
    return NULL;

  TString histname = In2DHist->GetName();
  histname.Append("_1d");
  cout << endl;
  TH1D* finalhist = In2DHist->ProjectionX(histname, 0, -1, "e");
  
  for (int ibin=1; ibin<In2DHist->GetNbinsX()+1; ++ibin){
    TH1D* slice = In2DHist->ProjectionY("temp", ibin, ibin, "e");
    Double_t allentries = slice->Integral();
    Double_t matched = allentries - slice->GetBinContent(1);
    
    if (allentries){
      finalhist->SetBinContent(ibin, matched/allentries);
      finalhist->SetBinError(ibin, TMath::Sqrt(matched*(1.0-matched/allentries))/allentries);
    }
    delete slice;
  }

  finalhist->SetMinimum(0.4);
  finalhist->GetYaxis()->SetTitle("Jet Matching Efficiency");
  return finalhist;

}


TH1* RootPlotter::GetSliceHisto(TH2D* In2DHist, Int_t bin)
{
  // get a slice of the 2d histogram

  if (In2DHist==NULL)
    return NULL;

  TString histname = In2DHist->GetTitle();

  TString postfix;
  postfix.Form("_Bin_%i",bin);
  TObjArray* toktitle = histname.Tokenize("_");
  TString BinVariable;
  Bool_t found = false;
  for (Int_t i=0; i<toktitle->GetEntries(); ++i){
    TString str( ((TObjString*) toktitle->At(i))->GetString() );
    if (str.Contains("Bin") ){
      BinVariable = str(0, str.First("Bin"));
      found = true;
    } else if ( str.Contains("bin") ) {
      BinVariable = str(0, str.First("bin"));
      found = true;
    }
    
    // append the lxy, lx, ly again
    if (str.Contains("lx")) postfix.Append("_lx");
    if (str.Contains("ly")) postfix.Append("_ly");
    if (str.Contains("lxy")) postfix.Append("_lxy");	

  }
  TAxis* axis = In2DHist->GetXaxis();
  TString low = TString::Format("%f", axis->GetBinLowEdge(bin));
  if (low(low.First(".")+1, 1)=="0"){
    low = low(0, low.First(".")); // cut out everything behind the comma
  } else {
    low = low(0, low.First(".")+2); // keep first digit after the comma
  }
  low.ReplaceAll(" ", "");
  TString up = TString::Format("%f", axis->GetBinUpEdge(bin));
  if (up(up.First(".")+1, 1)=="0"){
    up = up(0, up.First("."));  // cut out everything behind the comma
  } else {
    up = up(0, up.First(".")+2);  // keep first digit after the comma
  }
  up.ReplaceAll(" ", "");
  TString subtitle = TString::Format("SUBTITLE:%s < %s < %s", low.Data(), BinVariable.Data(), up.Data());
  histname.Append(postfix);
  TH1D* slice = In2DHist->ProjectionY(histname, bin, bin, "e");
  if (found){
    slice->SetTitle(subtitle);
  } else {
    slice->SetTitle(histname);
  }
  return slice;
 
}


TObjArray* RootPlotter::FindSubdirs(TDirectory* dir)
{
  // find all subdirectories (former histogram collections) in the input directory
  // return a TObjArray with the names of the subdirectories 
  dir->cd();

  TObjArray* subdirnames = new TObjArray();
  TKey *key;
  TIter nextkey( gDirectory->GetListOfKeys() );
  while ( (key = (TKey*)nextkey())) {
    TObject *obj = key->ReadObj();
    if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {    // found a subdirectory! 
      TString dirname(((TDirectory*) obj)->GetName());
      subdirnames->Add(new TObjString(dirname));
    }
  }
  return subdirnames;
}


TObjArray* RootPlotter::GetPtBalanceHists()
{
  // get the pt-balance histograms that should be fitted and plotted accordingly
  // they are stored in a TObjArray, whereas every entry is an array of 
  // histograms that belong together (e.g. data and different mcs)

  /*
  TObjArray* PtBalanceHists = new TObjArray();
  TObjArray* TwoDHistArray = NULL;
    
  if (!fChainDirs)
    FindDirectories();

  if (fChainDirs->GetEntries() == 0){
    cout << "No directories available. Exiting." << endl;
    exit(EXIT_FAILURE);
  }

  // find all subdirectories (former histmanagers) in the first directory
  TDirectory* firstdir = (TDirectory*) fChainDirs->At(0);
  TObjArray* dirs = FindSubdirs(firstdir);

  // loop over all directories and find the right histograms
  for (Int_t i=0; i<dirs->GetEntries(); ++i){

    firstdir->cd();
    TString dirname(((TObjString*) dirs->At(i))->GetString());
    gDirectory->Cd(dirname.Data());

  
    // loop over all histograms in the directory
    TKey *key;
    TIter nextkey( gDirectory->GetListOfKeys() );
    while ( (key = (TKey*)nextkey())) {
      TObject *obj = key->ReadObj();

      // set the array to 0 = no histograms found yet
      TwoDHistArray = NULL;  
      
      // get the histograms for all chains
      if ( obj->IsA()->InheritsFrom( "TH2" ) ) {

	// first histogram found
	TH2* TwoDHist = (TH2*) obj;
	TString histname(TwoDHist->GetName());

	// see if it's one of the needed histograms
	if ( ! histname.Contains("PtBalanceFit") ){
	  delete TwoDHist;
	  continue;
	}	
	

	TString oldchainname(((TObjString*) fChainNames->At(0))->GetName());

	// plot the sum if it's at the first position
	if (bPlotSum && oldchainname.Contains("SUM")){
	  TString histtitle = TwoDHist->GetTitle();
	  Ssiz_t pos = histtitle.First("(");
	  TString histending(histtitle(0,pos));
	  TwoDHist = (TH2*) SumHistograms(dirname, histending);
	}

	// histograms found! 	   
	TString histtitle = TwoDHist->GetTitle();
	TwoDHist->SetTitle(dirname);
	TwoDHistArray = new TObjArray();
	TwoDHistArray->Add(TwoDHist);

	// find the same histogram for the other chains
	Int_t ichaindir=0;
	for (Int_t ichain=1; ichain<fNumOfSamples; ++ichain){
	  ++ichaindir;
	  Ssiz_t pos = histtitle.First("(");
	  TString histending(histtitle(0,pos));


	  // if the requested chain name contains SUM, sum the corresponding
	  // chains to create on histogram
	  TString chainname( ((TObjString*) fChainNames->At(ichain))->GetName() );
	  if (bPlotSum && chainname.Contains("SUM")){

	    TwoDHist = (TH2*) SumHistograms(dirname, histending);
	    --ichaindir;
	    
	    if (!TwoDHist){
	      cout << "Could not sum histograms " << histending << " for directory name = " << dirname << ". Aborting." << endl;
	      exit(EXIT_FAILURE);
	    }
	    
	  
	  // don't build a sum, get the right histogram from the other chains
	  } else {
	    
	    TDirectory* thisdir = (TDirectory*) fChainDirs->At(ichaindir);
	    thisdir->cd();
	    gDirectory->Cd(((TObjString*) dirs->At(i))->GetString());

	    // loop over all histograms and get the right one 
	    TwoDHist = NULL;
	    TKey *dirkey;
	    TIter nextdirkey( gDirectory->GetListOfKeys() );
	    while ( (dirkey = (TKey*) nextdirkey())) {
	      TObject* dirobj = dirkey->ReadObj();
	      if ( dirobj->IsA()->InheritsFrom( "TH2" ) ) {
		TH2* dirhist = (TH2*) dirobj;
		TString dirhisttitle = dirhist->GetTitle();
		Ssiz_t pos2 = dirhisttitle.First("(");
		TString dirhistending(dirhisttitle(0,pos2));				    
		if (dirhistending == histending){
		  TwoDHist = (TH2*) dirobj;
		} else {
		  delete dirobj;
		}
	      }
	    }

	    if (!TwoDHist){
	      cout << "Could not find histogram " << histname << " in chain " << thisdir->GetName() << " and subdirectory " << ((TObjString*) dirs->At(i))->GetString() << endl;
	      cout << "Aborting" << endl;
	      exit(EXIT_FAILURE);
	    }
	    
	  }
	  
	  TwoDHist->SetTitle(dirname);
	  TwoDHistArray->Add(TwoDHist);
	  	  
	}

      }  else { // do nothing for an unknown object
	delete obj;
	continue;
      }
      
      if (TwoDHistArray != 0){
	PtBalanceHists->Add(TwoDHistArray);
      }


    } // end loop over all histograms in the directory
      
  } // end loop over all directories
  

  return PtBalanceHists;
  */
  return NULL;

}


TObjArray* RootPlotter::GetResolutionHists()
{
  // get the resolution histograms that should be fitted and plotted accordingly
  // they are stored in a TObjArray, whereas every entry is an array of 
  // histograms that belong together (e.g. data and different mcs)

  /*
  TObjArray* ResolutionHists = new TObjArray();
  TObjArray* TwoDHistArray = NULL;
    
  if (!fChainDirs)
    FindDirectories();

  if (fChainDirs->GetEntries() == 0){
    cout << "No directories available. Exiting." << endl;
    exit(EXIT_FAILURE);
  }

  // find all subdirectories (former histmanagers) in the first directory
  TDirectory* firstdir = (TDirectory*) fChainDirs->At(0);
  TObjArray* dirs = FindSubdirs(firstdir);

  // loop over all directories and find the right histograms
  for (Int_t i=0; i<dirs->GetEntries(); ++i){

    firstdir->cd();
    TString dirname(((TObjString*) dirs->At(i))->GetString());
    gDirectory->Cd(dirname);
    
    // loop over all histograms in the directory
    TKey *key;
    TIter nextkey( gDirectory->GetListOfKeys() );
    while ( (key = (TKey*)nextkey())) {
      TObject *obj = key->ReadObj();

      // set the array to 0 = no histograms found yet
      TwoDHistArray = NULL;  
      
      // get the histograms for all chains
      if ( obj->IsA()->InheritsFrom( "TH2" ) ) {

	// first histogram found
	TH2* TwoDHist = (TH2*) obj;
	TString histname(TwoDHist->GetName());

	// see if it's one of the needed histograms
	if ( ! (histname.Contains("Resolution") && histname.Contains("Bins")) ){
	  continue;
	}
	if (histname.Contains("Resolution") && histname.Contains("JetMatchingEfficiency"))
	  continue;

	TString oldchainname(((TObjString*) fChainNames->At(0))->GetName());
	
	// plot the sum if it's at the first position
	if (bPlotSum && oldchainname.Contains("SUM")){
	  TString histtitle = TwoDHist->GetTitle();
	  Ssiz_t pos = histtitle.First("(");
	  TString histending(histtitle(0,pos));
	  TwoDHist = (TH2*) SumHistograms(dirname, histending);
	}

	// histograms found! 
	TwoDHistArray = new TObjArray();
	TwoDHistArray->Add(TwoDHist);


	// find the same histogram for the other chains
	Int_t ichaindir=0;
	for (Int_t ichain=1; ichain<fNumOfSamples; ++ichain){
	  ++ichaindir;
	  TString histtitle = TwoDHist->GetTitle();
	  Ssiz_t pos = histtitle.First("(");
	  TString histending(histtitle(0,pos));

	  // if the requested chain name contains SUM, sum the corresponding
	  // chains to create on histogram
	  TString chainname( ((TObjString*) fChainNames->At(ichain))->GetName() );
	  if (bPlotSum && chainname.Contains("SUM")){

	    TwoDHist = (TH2*) SumHistograms(dirname, histending);
	    --ichaindir;
	    
	    if (!TwoDHist){
	      cout << "Could not sum histograms " << histending << " for directory name = " << dirname << ". Aborting." << endl;
	      exit(EXIT_FAILURE);
	    }
	   
	  
	  // don't build a sum, get the right histogram from the other chains
	  } else {

	    
	    TDirectory* thisdir = (TDirectory*) fChainDirs->At(ichain);
	    thisdir->cd();
	    gDirectory->Cd(((TObjString*) dirs->At(i))->GetString());
	    
	    // loop over all histograms and get the right one 
	    TwoDHist = NULL;
	    TKey *dirkey;
	    TIter nextdirkey( gDirectory->GetListOfKeys() );
	    while ( (dirkey = (TKey*) nextdirkey())) {
	      TObject* dirobj = dirkey->ReadObj();
	      if ( dirobj->IsA()->InheritsFrom( "TH2" ) ) {
		TH2* dirhist = (TH2*) dirobj;
		TString dirhisttitle = dirhist->GetTitle();
		Ssiz_t pos2 = dirhisttitle.First("(");
		TString dirhistending(dirhisttitle(0,pos2));				    
		if (dirhistending == histending){
		  TwoDHist = (TH2*) dirobj;
		} else {
		  delete dirhist;
		}
	      } else {
		delete dirobj;
	      }
	    }

	    if (!TwoDHist){
	      cout << "Could not find histogram " << histname << " in chain " << thisdir->GetName() << " and subdirectory " << ((TObjString*) dirs->At(i))->GetString() << endl;
	      cout << "Aborting" << endl;
	      exit(EXIT_FAILURE);
	    }
	    
	  }
	 
	  TwoDHistArray->Add(TwoDHist);
	  	  
	}

      }  else { // do nothing for an unknown object
	continue;
      }
      
      if (TwoDHistArray != 0){
	ResolutionHists->Add(TwoDHistArray);
      }


    } // end loop over all histograms in the directory
      
  } // end loop over all directories

  return ResolutionHists;
  */
  return NULL;
}

void RootPlotter::MakeSubTitle(TH1* hist)
{
  // see if a calo wheel is in the hist's title, 
  // make a subtitle out of it

  Bool_t found = false;
  
  TString name = hist->GetName();
  TObjArray* parts = name.Tokenize("_");
  TString wheel;
  for (Int_t iword=0; iword<parts->GetEntries();++iword){
    TString word = ((TObjString*) parts->At(iword))->GetString();
    if (word.Contains("CB") || 
	word.Contains("IF") || 
	word.Contains("FB") || 
	word.Contains("OF") || 
	word.Contains("BBE") ){
      wheel.Append(word);
      found = true;
    }
  }
  if (found){
    TString subtitle = TString::Format("SUBTITLE:%s", wheel.Data());
    hist->SetTitle(subtitle);
  }
}


void RootPlotter::PlotYields(TH1* hist)
{
  // plot the lumi yield

  static TCanvas* can = new TCanvas("LumiYieldCanvas","Lumi Yield Plot",600,800);
  can->Divide(1,2);
  can->cd(1);

  hist->SetTitleOffset(1.2,"Y");
  hist->SetTitleOffset(0.9,"X"); 
  hist->SetTitle("");
  hist->GetXaxis()->SetTitle("Run Number");
  hist->SetLineWidth(2);
  hist->SetMinimum(0.0);
  hist->SetMaximum(1.2*hist->GetMaximum());
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.8);

  // fit it 
  Float_t min = hist->GetXaxis()->GetXmin();
  Float_t max = hist->GetXaxis()->GetXmax();
  TF1* func = new TF1("lumi_yield", "[0]", min, max);
  func->SetLineWidth(2);
  func->SetLineColor(kAzure-2);
  hist->Fit(func, "N");
  hist->Draw("P ");
  
  Float_t mean = func->GetParameter(0);
  Float_t err = mean*0.025;

  TPave* errbar = new TPave(min, mean-err, max, mean+err);
  errbar->SetBorderSize(0);
  errbar->SetFillColor(kAzure-9);
  errbar->SetFillStyle(3013);
  errbar->Draw();
  func->Draw("same");

  TPaveText* pave = new TPaveText(0.45, 0.15, 0.85, 0.23, "brNDC");
  pave->SetFillColor(kWhite);
  pave->SetBorderSize(0);
  TString Result = TString::Format("< %s > = %4.0f #pm %2.0f", hist->GetYaxis()->GetTitle(), mean, err);
  pave->AddText(Result);
  pave->SetTextSize(0.05);
  pave->SetTextFont(22);
  pave->SetTextAlign(12);
  pave->SetTextColor(kAzure-2);
  pave->Draw();

  hist->Draw("P same");

  return;

}

void RootPlotter::DrawPageNum(TCanvas* can, Int_t num)
{
  // draw a pagenumber in the bottom right corner

  can->cd();
  TPaveText* text; 
  TString s;
  s.Form("%i",num);
  if (bPortrait){
    text = new TPaveText(0.93,0.00, 0.97, 0.03, "NDC");
  } else {
    text = new TPaveText(0.03,0.00, 0.06, 0.03, "NDC");
  }
  text->SetBorderSize(0);
  text->SetFillColor(0);
  text->AddText(s.Data());
  text->Draw("same");
  
}



TH1* RootPlotter::BuildBackground(TString subdirname, TString histname)
{
  // sum histograms from all background samples, i.e. all samples
  // that have the string Bkgd in their name

  TH1* Sum = NULL;

  /*
  TDirectory* origdir = gDirectory;
  fRootfile.cd();
  if (!fRootfile.Cd("JetsAtHighQ2")){
    cout << "Could not find analysis 'JetsAtHighQ2' in file " << fRootfile.GetName() << ". Abort!" << endl;
    exit(EXIT_FAILURE);
  }

  TKey *chain;
  TIter nextchain( gDirectory->GetListOfKeys() );
  while ( (chain = (TKey*) nextchain())) {

    // move into the next directory
    TObject* chaindir = chain->ReadObj();
    if (chaindir->IsA()->InheritsFrom( TDirectory::Class() )){
      ((TDirectory*) chaindir)->cd();
    } else {
      continue;
    }
    
    // check if it's a background chain (by its name)
    TString dirname(((TDirectory*) chaindir)->GetName());
    if (!dirname.Contains("Bkgd")) continue;

    // go into the histogram manager's directory
    gDirectory->Cd( subdirname.Data() );

    // loop over all histograms and get the right one 
    TKey *dirkey;
    TIter nextdirkey( gDirectory->GetListOfKeys() );
    while ( (dirkey = (TKey*) nextdirkey())) {
      TObject* dirobj = dirkey->ReadObj();
      if ( dirobj->IsA()->InheritsFrom( "TH1" ) ) {
	TH1* dirhist = (TH1*) dirobj;
	TString dirhisttitle = dirhist->GetTitle();
	Ssiz_t pos2 = dirhisttitle.First("(");
	TString dirhistending(dirhisttitle(0,pos2));				    
	if (dirhistending == histname){

	  // found the right object, build sum
	  if (!Sum){
	    Sum = (TH1*) dirhist->Clone();
	  } else {
	    Sum->Add(dirhist);
	  }

	} else {
	  delete dirhist;
	}

      } else {
	delete dirobj;
      }

    } // end loop over keys
    
  } // end loop over all chains

  // let's return to the original directory
  origdir->cd();
  */

  return Sum;

}
