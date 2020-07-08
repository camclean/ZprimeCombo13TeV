#include <iostream>
#include <iomanip>

#include <TObjArray.h>
#include <TObjString.h>
#include <TStyle.h>
#include <TFile.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include "Math/QuantFuncMathCore.h"
#include <TColor.h>
#include "SPlotter.h"


using namespace std;

SPlotter::SPlotter()
{
  m_can = NULL;
  m_ps  = NULL;
  m_ps_name = "default.ps";

  m_pad1 = NULL;
  m_pad2 = NULL;

  m_rp1_top = NULL;
  m_rp1  = NULL;
  m_rp2_top = NULL;
  m_rp2  = NULL;

  m_page       = 0;
  m_lumi       = 0;
  m_syserr     = -1;
  debug        = false;
  bShapeNorm   = false;
  bPortrait    = true;
  bDrawEntries = false;
  bDrawLumi    = true;
  bForPrelim   = false;
  bForPublication = false;
  bDrawLegend  = true;
  bPlotRatio   = false;
  bZScoreInRatio = false;
  bSingleEPS   = false;
  need_update  = true;
  bPlotLogy    = false;
  bIgnoreEmptyBins = false;
  bPubStyleErrors = true;

  m_printout = false;

}

SPlotter::~SPlotter()
{
  Cleanup();
}

void SPlotter::SetPsFilename(TString name)
{
  if (!name.EndsWith(".ps")){
    cerr << "SPlotter::SetPsFilename, given filename: " << name
	 << " does not end with .ps, intention? Please correct steering." << endl;
    exit(EXIT_FAILURE);
  }
  m_ps_name = name;
  
}

void SPlotter::DoStacking(vector<TObjArray*>& hists, TObjArray* StackNames, bool rename)
{
  if (hists.size()==0){
    cerr << "SPlotter::DoStacking: Empty array of histograms. Aborting." << endl;
    exit(EXIT_FAILURE);
  }

  if (!StackNames){ // trivial case: do nothing
    return;
  }

  // loop over all histogram arrays
  int narr = hists.size();
  for (int i=narr-1; i>=0; --i){
    TObjArray* ha = hists[i];
    if (ha->GetEntries()<1) continue;    
    TString proc = ((SHist*)ha->At(0))->GetProcessName();
    if (debug) cout << "SPlotter::DoStacking, hist array = " << i 
		    << " process name = " << proc << endl;

    // loop over all stack-names
    for (int j=0; j<StackNames->GetEntries(); ++j){
      TString sname = ((TObjString*)StackNames->At(j))->GetString();
      if (debug) cout << " stack name = " << sname << endl;
      if (proc.Contains(sname)){
	if (debug) cout << " -> found match, stacking this array." << endl;
	StackHists(hists, i, rename);
	break;
      }
    }
  }

  if (debug) cout << "SPlotter::DoStacking: Done." << endl;

  return;

}

void SPlotter::StackHists(std::vector<TObjArray*>& hists, int index, bool rename)
{
  // stack histograms at position 'index' with an existing array of stacks
  // in hists
  // if the stacks don't exist, they are created and added to the array
  
  // get the stack (create a new one if it doesn't exist yet)
  TObjArray* stacks = GetStacks(hists, index);

  // add the histograms at 'index' to the stack
  for (int i=0; i<stacks->GetEntries(); ++i){
    SHist* stack = (SHist*)stacks->At(i);
    SHist* hist = (SHist*)hists[index]->At(i);
    if (!stack || !hist){
      cerr << "SPlotter::StackHists: stack or hist at position " << i 
	   << " does not exist! Abort." << endl;
      cerr << "index of hists = " << hists.size() << " histograms = " << hists[index]->GetEntries() << endl;
      exit(EXIT_FAILURE);
    }
    // sanity check: compare names
    TString stackname = stack->GetStack()->GetName();
    TString histname = hist->GetHist()->GetName();
    if (!stackname.Contains(histname)){
      cerr << "SPlotter::StackHists: incompatible histograms at position " << i 
	   << ", stackname = " << stackname << " histname = " << histname 
	   << ". Prefer to exit because of consistency." << endl;
      exit(EXIT_FAILURE);
    }
    // still here? do the stackin'!
    hist->GetHist()->SetFillColor(hist->GetHist()->GetLineColor());
    hist->GetHist()->SetFillStyle(1001);
    if (rename){
      TString pname = hist->GetProcessName();
      hist->SetName(histname + "__" + pname);
    }
    stack->GetStack()->Add(hist->GetHist());
    hist->SetIsUsedInStack(true);
    hist->SetDoDraw(false);
    stack->SetUnc(hist->GetUnc(), stack->GetStack()->GetHists()->GetSize()-1);
    if (debug) cout << "stacking hist " << histname << " on " << stackname << " for process " << hist->GetProcessName()
		    << " (dir = " << stack->GetDir() << ")" << endl;
    if (i==0){
      cout << "stacking process " << hist->GetProcessName() << " with weight " << hist->GetWeight() << " and uncertainty " << hist->GetUnc() << endl;
    }
  }  
  
  return;

}

TObjArray* SPlotter::GetStacks(std::vector<TObjArray*>& hists, int index)
{
  // get the array of stacks from the input hists
  // if it doesn't exist, a new array will be created if index>0
  // and the hists at position 'index' will be used as blue-print

  // try to find a stack in the array
  TObjArray* arr = NULL;
  int narr = hists.size();
  for (int i=0; i<narr; ++i){
    if (hists[i]->GetEntries()==0){
      cerr << "SPlotter::GetStacks: Got no histograms in array " << i 
	   << " unexpected behaviour - abort." << endl;
      exit(EXIT_FAILURE);
    }

    arr = hists[i];
    SHist* sh = (SHist*) arr->At(0);
    if (sh->IsStack()){
      if (debug) cout << "SPlotter::GetStacks: Found stack at position " << i << endl;
      return arr;
    }
  }

  // no stack found, create a new array with THStacks -> use position 'index'
  // in the array as a blue-print
  if (index>-1){
    if (debug) cout << "SPlotter::GetStacks: Creating new array of THStacks " 
		    << "using position " << index << " as blueprint." << endl;
    if (index>narr){
      cerr << "SPlotter::GetStacks: Can not create an array of stacks from array"
	   << " index " << index << ", since size is only " << hists.size() 
	   << ". Unexpected behaviour - abort." << endl;
      exit(EXIT_FAILURE);
    }

    arr = new TObjArray();
    for (int i=0; i<hists[index]->GetEntries();++i){      
      TString hname = ((SHist*)hists[index]->At(i))->GetHist()->GetName();
      TString name = hname + "_stack";
      THStack* st = new THStack(name, "");
      SHist* sh = new SHist(st);
      sh->SetDir(((SHist*)hists[index]->At(i))->GetDir());
      sh->SetProcessName("SM");
      sh->SetDoDraw(true);
      arr->Add(sh);
      if (debug) cout << "SPlotter::GetStacks: Adding stack with name " << name
		      << " in directory " << sh->GetDir() << " at position " << i <<  endl;
    }

    hists.push_back(arr);
    if (debug) cout << "SPlotter::GetStacks: Added stack array to collection. " 
		    << "New size = " << hists.size() << ", old = " << narr << endl;
  }

  return arr;

}

void SPlotter::CopyStyle(TH1& h1, TH1* h2)
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

void SPlotter::SetupGlobalStyle()
{
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

}

void SPlotter::Cleanup()
{
  // do what the name suggests

  ClosePostscript();
  if (m_can){
    delete m_can;
    m_can = NULL;
  }
}

void SPlotter::SetupCanvas()
{
  // set up a canvas, different possibilities 
  // to take into account portrait/landscape 
  // and ratio/no ratio plots

  Int_t CanWidth;
  Int_t CanHeight;
  if (bPortrait){
    CanWidth = 600;
    CanHeight = 830;
  } else {
    CanWidth =  800;
    CanHeight = 600;
  }

  // set up the canvas
  m_can = new TCanvas("canvas","Control Plots", CanWidth, CanHeight);

  Float_t yplot = 0.3;
  Float_t yratio = 0.17;

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


      
  m_rp1_top = new TPad("pad1", "Control Plots 1", x1, y5, x2, y6);
  m_rp1 = new TPad("rp1", "Ratio1", x1, y4, x2, y5);
  
  m_rp2_top = new TPad("pad2", "Control Plots 2", x1, y2, x2, y3);
  m_rp2 = new TPad("rp2", "Ratio2", x1, y1, x2, y2);
  

  m_pad1 = new TPad("pad1", "Control Plots 1", x1, y4, x2, y6);
  m_pad2 = new TPad("pad2", "Control Plots 2", x1, y1, x2, y3);
  
  // set margins for portrait mode
  if (bPortrait){
    
    m_pad1->SetTopMargin(0.05); m_pad1->SetBottomMargin(0.16);  m_pad1->SetLeftMargin(0.19); m_pad1->SetRightMargin(0.05);
    m_pad2->SetTopMargin(0.05); m_pad2->SetBottomMargin(0.16);  m_pad2->SetLeftMargin(0.19); m_pad2->SetRightMargin(0.05);
    
    m_rp1_top->SetTopMargin(0.065); m_rp1_top->SetBottomMargin(0.0);  m_rp1_top->SetLeftMargin(0.19); m_rp1_top->SetRightMargin(0.05);
    m_rp2_top->SetTopMargin(0.065); m_rp2_top->SetBottomMargin(0.0);  m_rp2_top->SetLeftMargin(0.19); m_rp2_top->SetRightMargin(0.05);
    m_rp1->SetTopMargin(0.0);    m_rp1->SetBottomMargin(0.35);  m_rp1->SetLeftMargin(0.19);  m_rp1->SetRightMargin(0.05);
    m_rp2->SetTopMargin(0.0);    m_rp2->SetBottomMargin(0.35);  m_rp2->SetLeftMargin(0.19);  m_rp2->SetRightMargin(0.05);    
    
      // margins for landscape
  } else {
    
    m_rp1_top->SetTopMargin(0.065); m_rp1_top->SetBottomMargin(0.0);  m_rp1_top->SetLeftMargin(0.13); m_rp1_top->SetRightMargin(0.05);        
    m_rp2_top->SetTopMargin(0.065); m_rp2_top->SetBottomMargin(0.0);  m_rp2_top->SetLeftMargin(0.13); m_rp2_top->SetRightMargin(0.05);
    
    if (bPlotRatio){
	m_rp1->SetTopMargin(0.0);    m_rp1->SetBottomMargin(0.35);  m_rp1->SetLeftMargin(0.13);  m_rp1->SetRightMargin(0.05);
	m_rp2->SetTopMargin(0.0);    m_rp2->SetBottomMargin(0.35);  m_rp2->SetLeftMargin(0.13);  m_rp2->SetRightMargin(0.05);
    }
  }

  
  if (debug){
    m_rp1_top->SetFillColor(kYellow);
    m_rp2_top->SetFillColor(kOrange);
    if (bPlotRatio){
      m_rp1->SetFillColor(kGray);
      m_rp2->SetFillColor(kGray);
    }
  }

  m_pad1->Draw();
  m_pad2->Draw();

  m_rp1_top->Draw(); 
  m_rp2_top->Draw();
  
  if (bPlotRatio){
    m_rp1->Draw();
    m_rp2->Draw();
  }

  return;

}

void SPlotter::SetupCanvasForEPS()
{
  // set up a canvas for single EPS files
  // optimised plots for including in theses or publications and documents
  // different possibilities 
  // ratio/no ratio plots

  Int_t CanWidth;
  Int_t CanHeight;
  CanWidth = 400;
  CanHeight = 400;

  // set up the canvas
  m_can = new TCanvas("canvas","Control Plots", CanWidth, CanHeight);

  Float_t yplot = 0.65;
  Float_t yratio = 0.34;

                                                //  coordinates:
  // set up the coordinates of the two pads:    //  			 
  Float_t y1, y2, y3;                           //  y3 +-------------+	
  y3 = 0.99;                                    //     |             |	
  y2 = y3-yplot;                                //     |     pad1    |	
  y1 = y2-yratio;                               //  y2 |-------------|	
  Float_t x1, x2;                               //     |     rp1     |	
  x1 = 0.01;                                    //  y1 +-------------+	
  x2 = 0.99;                                    //     x1            x2	
                                                // 			
                                                // No Pad 2!            
                                                                                                                                             
      
  m_rp1_top = new TPad("pad1", "Control Plots 2", x1, y2, x2, y3);
  m_rp1 = new TPad("rp1", "Ratio2", x1, y1, x2, y2);
  m_pad1 = new TPad("pad1", "Control Plots 2", x1, y1, x2, y3);
 
  m_rp2_top = new TPad("pad1", "Control Plots 2", x1, y2, x2, y3);
  m_rp2 = new TPad("rp1", "Ratio2", x1, y1, x2, y2);
  m_pad2 = new TPad("pad1", "Control Plots 2", x1, y1, x2, y3);

      
  m_pad1->SetTopMargin(0.05); m_pad1->SetBottomMargin(0.16);  m_pad1->SetLeftMargin(0.19); m_pad1->SetRightMargin(0.05);
  m_pad2->SetTopMargin(0.05); m_pad2->SetBottomMargin(0.16);  m_pad2->SetLeftMargin(0.19); m_pad2->SetRightMargin(0.05);
  
  m_rp1_top->SetTopMargin(0.065); m_rp1_top->SetBottomMargin(0.01);  m_rp1_top->SetLeftMargin(0.19); m_rp1_top->SetRightMargin(0.05);
  m_rp2_top->SetTopMargin(0.065); m_rp2_top->SetBottomMargin(0.01);  m_rp2_top->SetLeftMargin(0.19); m_rp2_top->SetRightMargin(0.05);
  m_rp1->SetTopMargin(0.01);    m_rp1->SetBottomMargin(0.35);  m_rp1->SetLeftMargin(0.19);  m_rp1->SetRightMargin(0.05);
  m_rp2->SetTopMargin(0.01);    m_rp2->SetBottomMargin(0.35);  m_rp2->SetLeftMargin(0.19);  m_rp2->SetRightMargin(0.05);    
  
  if (debug){
    m_rp1_top->SetFillColor(kYellow);
    m_rp2_top->SetFillColor(kOrange);
    if (bPlotRatio){
      m_rp1->SetFillColor(kGray);
      m_rp2->SetFillColor(kGray);
    }
  }

  m_pad1->Draw();
  m_pad2->Draw();

  m_rp1_top->Draw(); 
  m_rp2_top->Draw();
  
  if (bPlotRatio){
    m_rp1->Draw();
    m_rp2->Draw();
  }

  return;

}

void SPlotter::OpenPostscript(TString dir, TString hname)
{
  // create a new ps file with the directory in the name
  // optional: for EPS files add the name of the histogram

  TString filename(m_ps_name);
  filename.ReplaceAll(".ps","");
  filename.Append("_");
  filename.Append(dir);
  filename.Append(".ps");

  if (bSingleEPS){
    filename.ReplaceAll(".ps","");
    filename.Append("_");
    filename.Append(hname);
    filename.Append(".eps");
    
  } else {
    
    TString text(dir);
    text.Prepend("Plotting all histograms in directory ");
    cout << "\n+-------------------------- SFrame Plotter ---------------------------+" << endl;
    cout <<   "| " << setw(60)<< text                                    << "        |" << endl;
    cout <<   "+---------------------------------------------------------------------+" << endl;
    m_page = 0;
  }
      
  m_ps = NULL;
  if (bSingleEPS){
    m_ps = new TPostScript(filename, 113); // eps output
  } else {
    if (bPortrait){
      m_ps = new TPostScript(filename, 111);  // ps output
      m_ps->Range(20.0, 30.0);
    } else {
      m_ps = new TPostScript(filename, 112);  // ps output
      m_ps->Range(27.0, 18.0);
    }
  }

}

void SPlotter::ClosePostscript()
{
  // close the ps file and set page number to 0  
  if (m_ps){
    m_ps->Close();
    delete m_ps;
    m_ps = NULL;
  }
  m_page = 0;
}

void SPlotter::ProcessAndPlot(std::vector<TObjArray*> histarr)
{
  // loop over all arrays in the input array and plot them 

  if (histarr.size()<1){
    cerr << "SPlotter::ProcessAndPlot: No arrays of histograms given. Abort." << endl;
    exit(EXIT_FAILURE);
  }

  if (histarr[0]->GetEntries()<1){
    cerr << "SPlotter::ProcessAndPlot: No histograms given. Abort." << endl;
    exit(EXIT_FAILURE);
  }
  
  if (bPlotRatio && histarr.size()==1){
    cerr << "SPlotter::ProcessAndPlot: Only one process given, can not plot " 
	 << " ratio. Steering correct?" << endl;
    exit(EXIT_FAILURE);
  }

  SetupGlobalStyle();

  TString psname = m_ps_name;
  TString current_dir = "";

  TString namebase = psname;
  namebase.ReplaceAll( ".ps", "" );

  // loop over all histograms and plot them!
  int iplot = 1;
  bool bleg = true;
  for (int i=0; i<histarr[0]->GetEntries(); ++i){

    // get the histograms for the different processes
    vector<SHist*> hists = GetPlottableHists(histarr, i);    

    // no plottable hists found at position i
    if (debug) cout << "Number of plottable hists at index " << i << " = " << hists.size() << endl;
    if (hists.size()==0) continue;

    if (bShapeNorm) ShapeNormalise(hists);

    int ipad = GetCurrentPad(iplot);
    
    if (debug) cout << "Plotting histograms " << hists[0]->GetName() 
		    << " iplot = " << iplot << " ipad = " << ipad << endl;
    
    // new directory? create new ps file for ps-book!
    if (!bSingleEPS){
      TString dir = hists[0]->GetDir();
      if (dir.CompareTo(current_dir)!=0){
	      if (iplot!=1) DrawPageNum();
        	Cleanup();
        	SetupCanvas();
        	OpenPostscript(dir);
        	current_dir = dir;
        	iplot = 1;
        	bleg = true;
      }

      // new page every second plot
      if (iplot%2==1){
      	if (debug) cout << "Creating new page with number " << m_page << endl;
        DrawPageNum();
      	if (need_update) m_can->Update();
      	m_ps->NewPage();
      	++m_page;
      }

    // new file for each plot in single EPS mode
    } else {
      TString dir = hists[0]->GetDir();
      TString hname = hists[0]->GetName();
      Cleanup();
      SetupCanvasForEPS();
      if (debug) cout << "Creating new eps file with name " << dir << "_" << hname << endl;
      OpenPostscript(dir, hname);
      current_dir = dir;
      iplot = 1;
      bleg = true;
    }

    // cosmetics
    DoCosmetics(hists);

    // ---------- do what we set out to do: plot! ----------------

    if (hists[0]->IsYieldPlot()){  // special treatment for lumi yield plot
    
      PlotLumiYield(hists[0], ipad);

    } else { // usual plots

      PlotHists(hists, ipad);
      // draw a legend     
      if (bleg){
      	DrawLegend(GetHistsAtIndex(histarr, i));
      	if (!bDrawLegend) bleg = false;
      }
      // draw lumi information
      TString hname = hists[0]->GetName();
      if (hname.Contains("mjhtt")){
      	if (bDrawLumi) DrawLumi(18.3);
      } else {
      	if (bDrawLumi) DrawLumi();
      }
      // finally: redraw axes
      gPad->RedrawAxis();

      // draw the ratio
      if (bPlotRatio){
      	if (bZScoreInRatio) PlotZScore(hists, ipad);
      	else PlotRatios(hists, ipad);
      }

      m_can->SaveAs( namebase + "_" + hname + ".root");
    }

    ++iplot;
  }

  // done!
  if (!bSingleEPS) DrawPageNum();

  if (need_update) m_can->Update();

  
  Cleanup(); 
  
}

void SPlotter::PlotLumiYield(SHist* hist, int ipad)
{
  // plot the lumi yield histogram

  if (ipad==1) m_pad1->cd();
  if (ipad==2) m_pad2->cd();

  hist->Draw();

  // calculate the average
  TH1* h = hist->GetHist();
  double sum=0;
  int bins=0;
  for (int i=1; i<h->GetNbinsX()+1; ++i){
    if (h->GetBinContent(i)>0){
      sum += h->GetBinContent(i);
      bins++;
    }
  }
  double av = sum / bins;

  // calculate average with outlier-rejection (4sigma)
  sum=0;
  bins=0;
  for (int i=1; i<h->GetNbinsX()+1; ++i){
    if (h->GetBinContent(i)>0){
      double dev = TMath::Abs( (h->GetBinContent(i) - av)/h->GetBinError(i) );
      if (dev<4){
      	sum += h->GetBinContent(i);
      	bins++;
      } else {
      	cout << "Lumi yield: outlier in bin " << i << " with content " << h->GetBinContent(i) << " average = " << av << endl;
      } 
    }
  }
  av = sum / bins;

  // calculate error on mean and chi2
  double dev = 0;
  double chi2 = 0;
  for (int i=1; i<h->GetNbinsX()+1; ++i){
    if (h->GetBinContent(i)>0){
      double pull = (h->GetBinContent(i) - av)/h->GetBinError(i);
      if (TMath::Abs(pull)<4){
      	dev += TMath::Power(h->GetBinContent(i)-av, 2);
      	chi2 += pull*pull;
      }
    }
  }
  double err = TMath::Sqrt(dev/bins);

  // highlight points with deviations of more than 3, 4 and 5 sigma 
  double xr = h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin();
  double wi = gPad->GetAbsWNDC() * (1 - gPad->GetLeftMargin() - gPad->GetRightMargin());
  double he = gPad->GetAbsHNDC() * (1 - gPad->GetTopMargin() - gPad->GetBottomMargin());
  double ar = wi/he;
  double fudge = 1.;
  if (bSingleEPS) fudge = 1.2;
  double r1 = 0.02*xr*fudge;
  double yr = h->GetMaximum()-h->GetMinimum();
  double r2 = 0.016*yr*ar*fudge;
  for (int i=1; i<h->GetNbinsX()+1; ++i){
    if (h->GetBinContent(i)>0){
      double pull = (h->GetBinContent(i) - av)/h->GetBinError(i);
      if (TMath::Abs(pull)>5){
      	TEllipse* circ = new TEllipse(h->GetXaxis()->GetBinCenter(i), h->GetBinContent(i), r1, r2);
      	circ->SetFillColor(kWhite);
      	circ->SetLineColor(kRed);
      	circ->Draw();
      } else if (TMath::Abs(pull)>4){
      	TEllipse* circ = new TEllipse(h->GetXaxis()->GetBinCenter(i), h->GetBinContent(i), r1, r2);
      	circ->SetFillColor(kWhite);
      	circ->SetLineColor(kOrange);
      	circ->Draw();
      } else if (TMath::Abs(pull)>3){
      	TEllipse* circ = new TEllipse(h->GetXaxis()->GetBinCenter(i), h->GetBinContent(i), r1, r2);
      	circ->SetFillColor(kWhite);
      	circ->SetLineColor(kSpring);
      	circ->Draw();
      }
    }
  }

  // draw the average
  TF1* f = new TF1("average", "[0]", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  f->SetLineColor(kAzure+1);
  f->SetLineWidth(1);
  f->SetParameter(0, av);
  f->Draw("same");


  TF1* fup = new TF1("up", "[0]", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  TF1* fdown = new TF1("down", "[0]", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  fup->SetParameter(0, av+err);
  fdown->SetParameter(0, av-err);
  fup->SetLineColor(kAzure+1);
  fdown->SetLineColor(kAzure+1);
  fup->SetLineWidth(1);
  fdown->SetLineWidth(1);
  fup->SetLineStyle(kDashed);
  fdown->SetLineStyle(kDashed);
  fup->Draw("same");
  fdown->Draw("same");

  TLatex* text = new TLatex();
  text->SetTextFont(42);
  text->SetNDC();
  text->SetTextColor(kBlack);
  text->SetTextSize(0.05);
  if (bSingleEPS)  text->SetTextSize(0.04);
  text->SetTextAlign(11);
  TString info = TString::Format("#chi^{2} / ndf");
  text->DrawLatex(0.5, 0.30, info.Data());
  info = TString::Format("%3.1f / %d", chi2, bins-1);
  text->DrawLatex(0.65, 0.30, info.Data());
  info = TString::Format("average");
  text->DrawLatex(0.5, 0.23, info.Data());
  info = TString::Format("%4.1f #pm %4.1f", av, err);
  text->DrawLatex(0.65, 0.23, info.Data());

  hist->Draw("same");

  return;

}


void SPlotter::PlotHists(vector<SHist*> hists, int ipad)
{
  // plot all histograms in the array

  if (ipad==1){
    if (bPlotRatio) m_rp1_top->cd();
    else m_pad1->cd();
  }
  if (ipad==2){
    if (bPlotRatio) m_rp2_top->cd();
    else m_pad2->cd();
  }

  bool isok = SetMinMax(hists);
  if (isok) SetLogAxes(hists);

  // first get some basic histograms
  SHist* sstack = SelStack(hists);
  SHist* sdata  = SelData(hists);

  // first, draw data if it exists
  int ndrawn = 0;
  if (sdata){
    sdata->Draw();
    ++ndrawn;
  }

  if (debug){
    if (sdata){
      cout << "\nHist name = " << sdata->GetName() << " process = " << sdata->GetProcessName() << endl;
      cout << "hists, entries = " << hists.size() << endl;
      cout << "Data entries = " << sdata->GetHist()->Integral() << endl;
    }
    if (sstack){
      double stack_entries = 0;
      TObjArray* arr = sstack->GetStack()->GetStack();
      TH1* h = (TH1*)arr->At(arr->GetEntries()-1);
      stack_entries = h->Integral();
      cout << "Stack entries = " << stack_entries << endl;
      TList* hists = sstack->GetStack()->GetHists();
      // calculate individual area
      for (int i=0; i<hists->GetSize(); ++i){
      	TH1* h = (TH1*) hists->At(i);
      	int iend = h->GetNbinsX();
      	double area = h->Integral(1,iend);
      	cout << "  entries of histogram " << i << " in stack = " << area << endl;
      }
    }
  }

  if (m_printout){
    if (sdata) cout << "\nDATA: N = " << sdata->GetHist()->Integral() << endl;
  }
  m_written = 0;

  // first round
  int nh = hists.size();

  for (int i=0; i<nh; ++i){
    SHist* sh = hists[i];
    if (sh->IsStack()) continue;
    if (sh==sdata) continue;
    if (ndrawn==0) sh->Draw();
    else sh->Draw("same");
    ++ndrawn;
  }
 
  // now draw the stack
  if (sstack){
    if (ndrawn==0){
      sstack->Draw();
      need_update = false;    
    } else {
      sstack->Draw("same");
    }
  }
 
  // second round
  for (int i=0; i<nh; ++i){
    SHist* sh = hists[i];
    if (sh->IsStack()) continue;
    if (sh==sdata) continue;
    sh->Draw("same");
  }

  // draw normalisation error if it is given
  if (sstack){
    DrawSysError(sstack);
  }

  // draw data on top
  if (sdata){ 
    sdata->Draw("same");
    //cout << "name = " << sdata->GetName() << endl;
    // for data, set to draw poissonian coverage if required
    if (bPubStyleErrors){
      int lastnonzero = 99999;
      if (sstack){
        for (int i=1; i<sstack->GetStack()->GetHistogram()->GetNbinsX()+1; ++i){
          TH1D* h = (TH1D*) sstack->GetStack()->GetStack()->At(sstack->GetStack()->GetStack()->GetLast());
          //cout << "bin " << i << " stack entries = " << h->GetBinContent(i) << endl;
          if (h->GetBinContent(i)>0) lastnonzero = i;
        }
      }
      DrawPoissonCoverage(sdata, lastnonzero);
    }
  }

  gPad->RedrawAxis();
  
}

void SPlotter::DrawPoissonCoverage(SHist* data, int lastbin)
{

  double alpha = 1 - 0.6827;

  TGraphAsymmErrors* ge = new TGraphAsymmErrors();

  TH1D* dh = (TH1D*) data->GetHist();
  int np = 0;
  bool first = false;
  TString name = dh->GetName();
  for (int i=1; i<dh->GetNbinsX()+1; ++i){
    int N = dh->GetBinContent(i);
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
    if (N>0) first = true; 
    if (name=="muo_T1B0_L1chi2lo_mtt" && i==5) first = true;
    if (first && N<5 && i<=lastbin){
      ge->SetPoint(np,dh->GetXaxis()->GetBinCenter(i), N);
      ge->SetPointEYlow(np, N-L);
      ge->SetPointEYhigh(np, U-N);
      ++np;
    }
  }
  
  ge->SetLineColor(kBlack);
  ge->SetLineWidth(2);
  ge->Draw("Z0 same");

}

void SPlotter::DrawPoissonCoverageInRatio(vector<SHist*> hists)
{

  // get the data and the stack
  SHist* sstack = SelStack(hists);
  SHist* sdata  = SelData(hists);
  TH1D* mc = (TH1D*) sstack->GetStack()->GetStack()->At(sstack->GetStack()->GetStack()->GetLast());

  double alpha = 1 - 0.6827;

  TGraphAsymmErrors* ge = new TGraphAsymmErrors();

  TH1D* dh = (TH1D*) sdata->GetHist();
  TString name = dh->GetName();
  int np = 0;
  bool first = false;
  bool print = false;
  for (int i=1; i<dh->GetNbinsX()+1; ++i){
    int N = dh->GetBinContent(i);
    if (N>0) first = true;
    if (name=="muo_T1B0_L1chi2lo_mtt" && i==5) first = true; 
    if (name=="muo_T1B0_L1chi2lo_mtt" && i==7) print = true; 
    if (N>=5) continue;
    if (!first) continue;

    double pred = mc->GetBinContent(i);
    if (pred==0) continue;

    double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1);
    double rU = U / mc->GetBinContent(i);    
    double L = 0.; 
    double rL = 0.;
    if (N>0){ 
      L = ROOT::Math::gamma_quantile(alpha/2,N,1.);
      rL = L / mc->GetBinContent(i);
    }

    if (print){
      cout << "bin = " << i << " " << dh->GetXaxis()->GetBinLowEdge(i) << " < mtt < " << dh->GetXaxis()->GetBinUpEdge(i) << endl;
      cout << "N = " << N << endl;
      cout << "U = " << U << endl;
      cout << "L = " << L << endl;
      cout << "mc bin content = " << mc->GetBinContent(i) << endl;
      cout << "rU = " << rU << endl;
      cout << "rL = " << rL << endl;
    }

    if (first){
      double r = N / mc->GetBinContent(i);
      ge->SetPoint(np,dh->GetXaxis()->GetBinCenter(i), r);
      if (print){
        cout << "setting point to " << dh->GetXaxis()->GetBinCenter(i) << " and " << r << endl;
        cout << "setting point error to " << rL << " and " << rU << endl;
        print = false;
      }
      ge->SetPointEYlow(np, r-rL);
      ge->SetPointEYhigh(np, rU-r);
      ++np;
    }
  }
  
  ge->SetLineColor(kBlack);
  ge->SetLineWidth(2);
  ge->Draw("Z0 same");


}

void SPlotter::DrawSysError(SHist* stack)
{
  // plot an error band corresponding to the overall
  // systematic error
  // if a theta file is used as input with shape variations, 
  // also the error from those is drawn
  TH1* h = (TH1*) stack->GetStack()->GetStack()->At( stack->GetStack()->GetStack()->GetEntries()-1 );
  //TH1* e = (TH1*) h->Clone();
  TGraphAsymmErrors* eAsym = new TGraphAsymmErrors();

  for (Int_t i=1; i<h->GetNbinsX()+1; ++i){
    Double_t sys = 0; 
    if (m_syserr>0) sys = m_syserr*h->GetBinContent(i);
    Double_t stat = h->GetBinError(i);
    Double_t norm_err = CalcNormErrorForBin(stack, i);
    Double_t sys_err_plus = CalcShapeSysErrorForBinFromTheta(stack, i, "up");
    Double_t sys_err_minus = CalcShapeSysErrorForBinFromTheta(stack, i, "down");
    Double_t ey_low = TMath::Sqrt(norm_err*norm_err + sys*sys + stat*stat + sys_err_minus*sys_err_minus);
    Double_t ey_up = TMath::Sqrt(norm_err*norm_err + sys*sys + stat*stat + sys_err_plus*sys_err_plus);
    Double_t ex_low = (h->GetXaxis()->GetBinCenter(i)) - (h->GetXaxis()->GetBinLowEdge(i));
    Double_t ex_up =  (h->GetXaxis()->GetBinUpEdge(i))-h->GetXaxis()->GetBinCenter(i);
    eAsym -> SetPoint(i, h->GetXaxis()->GetBinCenter(i), h->GetBinContent(i)); 
    eAsym -> SetPointError(i, ex_low, ex_up, ey_low, ey_up); 

    if (m_printout){
      cout << "Bin " << i << ": total background " << stack->GetName() << " N = " << h->GetBinContent(i) << " +- " << (ey_low+ey_up)/2. << endl;
      cout << "( stat = " << stat << " xs_rate = " << norm_err << " norm = " << sys << " sys = " << (sys_err_plus + sys_err_minus)/2. << " ) " << endl;
    }
   
  }
   
  static Int_t LightGray     = TColor::GetColor( "#aaaaaa" );
  //h->SetFillColor(kGray+2);
  eAsym->SetFillColor(LightGray);
  eAsym->SetLineWidth(1);
  eAsym->SetFillStyle(3245);
  eAsym->Draw("E2 same");


}

double SPlotter::CalcNormErrorForBin(SHist* stack, int ibin)
{
  // calculate the normalisation uncertainty of a single bin in the stack
  // due to normalisation uncertainty on different processes
  
  double err = 0;
  for (int i=0; i<stack->GetStack()->GetStack()->GetEntries(); ++i){
    TH1* h = (TH1*) stack->GetStack()->GetHists()->At(i);
    err += h->GetBinContent(ibin)*stack->GetUnc(i);
  }
  return err;
}



double SPlotter::CalcShapeSysErrorForBinFromTheta(SHist* stack, int ibin, TString sign)
{
  double absoluteerr = 0;
  double squarederr = 0;
  double err = 0;
  if (m_shapesys_arr.size()==0)//no systamtics given in theta file
    return err;

  if (sign!="up" && sign!="down"){
    cout << "error in call to CalcShapeSysErrorForBinFromTheta: sign can only be 'up' or 'down', no systematic error will be calculated" << endl;
    return err;
  }
    
  // loop over all background samples to find the process
  for (int i=0; i<stack->GetStack()->GetStack()->GetEntries(); ++i){  
    TH1* h = (TH1*) stack->GetStack()->GetHists()->At(i);
    TString histname = h->GetName(); //e.g. HT__QCD
    histname.ReplaceAll("__", "#");
    TObjArray* histnamePieces = histname.Tokenize("#");
    TString variableName =  ((TObjString*)histnamePieces->At(0))-> GetString(); //this is the sample name e.g. HT
    TString sampleName =  ((TObjString*)histnamePieces->At(1))-> GetString(); //this is the sample name e.g. QCD or TTbar
    
    Double_t sample_err2 = 0;

    // loop over all systematic error samples
    for (unsigned int syst = 0; syst < m_shapesys_arr.size(); ++syst){
      TObjArray* arr = m_shapesys_arr[syst];
      for (int nplots = 0; nplots < arr->GetEntries(); ++nplots){
	SHist* hSys = (SHist*)arr->At(nplots);
	TH1* hSyst = hSys->GetHist();
	TString systFullName = hSys -> GetProcessName();//e.g. QCD__uncert__plus

	systFullName.ReplaceAll("__","#");
	TObjArray* systFullNamePieces = systFullName.Tokenize("#");
	TString systVariableName = hSys -> GetName();       

	// continue if the channel of the systematic sample has the same name as the channel of the background process
	if (variableName == systVariableName){

	  // check if the the sign is the same as requested
	  TString syssign = ((TObjString*) systFullNamePieces->At(2))->GetString();
	  if (syssign == sign){
	 
	    // check if systematic uncertainty comes from the same sample as the background (e.g. ttbar)
	    if (systFullNamePieces->Contains(sampleName)){
	      double fac = 1.0;
	      if (systFullName.Contains("ttbar")) fac = 0.94;      
        if (systFullName.Contains("wjet")) fac = 0.98;              
        if (systFullName.Contains("sitop")) fac = 0.94;              
	      //cout << "warning! factor of 0.95 for systematics!!! (line 956)" << endl;
	      absoluteerr = (hSyst->GetBinContent(ibin)*fac)-(h->GetBinContent(ibin));

	      // the second one contains the name of the uncertainty: check if the error should be reduced
	      TString sysname = ((TObjString*) systFullNamePieces->At(1))->GetString();
	      
	      // loop over systematics that should be reduced, find the right factor
	      for (Int_t j=0; j<m_ScaleSysUncName->GetEntries(); ++j){
      		TString sysname_to_red = ((TObjString*) m_ScaleSysUncName->At(j))->GetString();
      		if (sysname == sysname_to_red){		
      		  absoluteerr *= m_sysweight.At(j);
      		}
	      }
	      
	      // got it: add to the total error in quadrature
	      squarederr += absoluteerr*absoluteerr;	   
	      sample_err2 += absoluteerr*absoluteerr;	   
	    }
	  }	  
	}
      }
    }
    // -------------- only for output
    if (m_written==0){
      if (m_printout){
	cout << fixed << setprecision(0);
	double stat2 = TMath::Power(h->GetBinError(ibin),2);
	double xs_norm_err = h->GetBinContent(ibin)*stack->GetUnc(i);
	double norm_err = h->GetBinContent(ibin)*m_syserr;
	double tot_err = TMath::Sqrt(sample_err2 + stat2 + norm_err*norm_err + xs_norm_err*xs_norm_err);
	cout << "sample " << sampleName << " N = " << h->GetBinContent(ibin) << " +- " << tot_err 
	     << " (stat = " << TMath::Sqrt(stat2) << " sys = " << TMath::Sqrt(sample_err2) <<  " norm (lumi) = " 
	     << norm_err << " norm (xs) = " << xs_norm_err << " ) " << endl;
      }
    }
  }
  m_written = 1;
  err = TMath::Sqrt(squarederr);
  return err;
  
}

void SPlotter::PlotRatios(vector<SHist*> hists, int ipad)
{
  // plot all histograms in the array

  if (ipad==1) m_rp1->cd();
  if (ipad==2) m_rp2->cd();

  // calculate ratios
  vector<SHist*> ratios = CalcRatios(hists);

  gPad->SetLogx(0);
  gPad->SetLogy(0);

  int ndrawn = 0;
  int nh = ratios.size();
  double xmin = 0;
  double xmax = 0;
  for (int i=0; i<nh; ++i){
    SHist* rh = ratios[i];
    rh->DrawNoErrorX(false);
    TString name = rh->GetName();
    if (name.Contains("_lx")) gPad->SetLogx(1);
    if (ndrawn==0) rh->Draw();
    else rh->Draw("same");
    ++ndrawn;
    int first = rh->GetHist()->GetXaxis()->GetFirst();
    xmin = rh->GetHist()->GetXaxis()->GetBinLowEdge(first);
    int last = rh->GetHist()->GetXaxis()->GetLast();
    xmax = rh->GetHist()->GetXaxis()->GetBinUpEdge(last);
  }

  // draw line at 1:
  TLine* unity = new TLine(xmin, 1., xmax, 1.);
  unity->SetLineColor(kBlack);
  unity->SetLineStyle(kDashed);
  unity->Draw();

  // draw coverage for empty bins
  DrawPoissonCoverageInRatio(hists);

  //  cout << "Drawing reweighting function with errors..." << endl;
  bool draw_top_pt_rew = false;
  if (draw_top_pt_rew){
    TF1* fl = new TF1("rew_lo", "TMath::Exp(0.156-0.00137*x)",0, 400);
    TF1* fu = new TF1("rew_up", "TMath::Exp(0.156-0.00137*400)", 400, 2000);
    fl->SetLineWidth(1);
    fu->SetLineWidth(1);
    fl->SetLineStyle(kSolid);
    fu->SetLineStyle(kSolid);
    fl->SetLineColor(kRed);
    fu->SetLineColor(kRed);
    
    TF1* fl2 = new TF1("rew_lo", "TMath::Power(TMath::Exp(0.156-0.00137*x),2)",0, 400);
    TF1* fu2 = new TF1("rew_up", "TMath::Power(TMath::Exp(0.156-0.00137*400),2)", 400, 2000);
    fl2->SetLineWidth(1);
    fu2->SetLineWidth(1);
    fl2->SetLineStyle(kDotted);
    fu2->SetLineStyle(kDotted);
    fl2->SetLineColor(kRed);
    fu2->SetLineColor(kRed);
    
    TF1* f0 = new TF1("rew_lo", "1",0, 2000);
    f0->SetLineWidth(1);
    f0->SetLineStyle(kDotted);
    f0->SetLineColor(kRed);

    fl->Draw("same");
    fu->Draw("same");
    fl2->Draw("same");
    fu2->Draw("same");
    f0->Draw("same");
  }
  

  // draw theory uncertainty due to scale variations
  bool draw_scale_unc = false;
  bool electron = false;
  bool muon = false;
  if (draw_scale_unc){
  
    // electron
    if (electron){
      TString name = hists[0]->GetName();
      if (name.Contains("lep")){
	TF1* fdn1 = new TF1("fdn", "pol4", 120., 600);  
	TF1* fup1 = new TF1("fup", "[0]+[1]*(x-600.)", 600., 1400);
	fdn1->SetParameters(0.818924, 0.00322385, -9.86616e-06, 1.08823e-08, -3.46706e-12);
	fup1->SetParameters(1.10266e+00, 2.81655e-03);
	fdn1->SetLineWidth(1);
	fup1->SetLineWidth(1);
	fdn1->SetLineColor(kBlack);
	fup1->SetLineColor(kBlack);
	fdn1->Draw("same");
	fup1->Draw("same");
	
	TF1* fdn2 = new TF1("fdn", "pol4", 120., 600);  
	TF1* fup2 = new TF1("fup", "[0]+[1]*(x-600.)", 600., 1400);
	fdn2->SetParameters(1.18108, -0.00322385, 9.86615e-06, -1.08822e-08, 3.46703e-12);
	fup2->SetParameters(8.97344e-01, -2.81655e-03);
	fdn2->SetLineWidth(1);
	fup2->SetLineWidth(1);
	fdn2->SetLineColor(kBlack);
	fup2->SetLineColor(kBlack);
	fdn2->Draw("same");
	fup2->Draw("same");
      }
      
      if (name.Contains("had")){
	TF1* fdn1 = new TF1("fdn", "pol4", 0., 600);  
	TF1* fup1 = new TF1("fup", "[0]+[1]*(x-600.)", 600., 1400);
	fdn1->SetParameters(1.08361, 0.00180409, -8.69107e-06, 1.30316e-08, -4.14883e-12);
	fup1->SetParameters(1.31441e+00, 1.56129e-03);
	fdn1->SetLineWidth(1);
	fup1->SetLineWidth(1);
	fdn1->SetLineColor(kBlack);
	fup1->SetLineWidth(kBlack);
	fdn1->Draw("same");
	fup1->Draw("same");
	
	TF1* fdn2 = new TF1("fdn", "pol4", 0., 600);  
	TF1* fup2 = new TF1("fup", "[0]+[1]*(x-600.)", 600., 1400);
	fdn2->SetParameters(0.916392, -0.00180409, 8.69107e-06, -1.30316e-08, 4.14883e-12);
	fup2->SetParameters(6.85593e-01, -1.56129e-03);
	fdn2->SetLineWidth(1);
	fup2->SetLineWidth(1);
	fdn2->SetLineColor(kBlack);
  fup2->SetLineColor(kBlack);
	fdn2->Draw("same");
	fup2->Draw("same");
      }
    }

    // muon
    if (muon){
      TString name = hists[0]->GetName();
      if (name.Contains("lep")){
	TF1* fdn1 = new TF1("fdn", "pol4", 0., 600);  
	TF1* fup1 = new TF1("fup", "[0]+[1]*(x-600.)", 600., 1400);
	fdn1->SetParameters(1.96862, -0.0116395, 5.31737e-05, -1.04026e-07, 7.53312e-11);
	fup1->SetParameters(1.42066e+00, 1.68181e-03);
	fdn1->SetLineWidth(1);
	fup1->SetLineWidth(1);
	fdn1->SetLineColor(kBlack);
	fup1->SetLineColor(kBlack);
	fdn1->Draw("same");
	fup1->Draw("same");
	
	TF1* fdn2 = new TF1("fdn", "pol4", 0., 600);  
	TF1* fup2 = new TF1("fup", "[0]+[1]*(x-600.)", 600., 1400);
	fdn2->SetParameters(0.0313793, 0.0116395, -5.31737e-05, 1.04026e-07, -7.53312e-11);
	fup2->SetParameters(5.79338e-01, -1.68181e-03);
	fdn2->SetLineWidth(1);
	fup2->SetLineWidth(1);
	fdn2->SetLineColor(kBlack);
	fup2->SetLineColor(kBlack);
	fdn2->Draw("same");
	fup2->Draw("same");
      }
      
      if (name.Contains("had")){
	TF1* fdn1 = new TF1("fdn", "pol4", 0., 600);  
	TF1* fup1 = new TF1("fup", "[0]+[1]*(x-600.)", 600., 1400);
	fdn1->SetParameters(1.03754, 0.00132089, -8.34035e-06, 1.88731e-08, -1.16808e-11);
	fup1->SetParameters(1.39030e+00, 1.82452e-03);
	fdn1->SetLineWidth(1);
	fup1->SetLineWidth(1);
	fdn1->SetLineColor(kBlack);
	fup1->SetLineColor(kBlack);
	fdn1->Draw("same");
	fup1->Draw("same");
	
	TF1* fdn2 = new TF1("fdn", "pol4", 0., 600);  
	TF1* fup2 = new TF1("fup", "[0]+[1]*(x-600.)", 600., 1400);
	fdn2->SetParameters(0.962457, -0.00132089, 8.34034e-06, -1.8873e-08, 1.16808e-11);
	fup2->SetParameters(6.09700e-01, -1.82452e-03);
	fdn2->SetLineWidth(1);
	fup2->SetLineWidth(1);
	fdn2->SetLineColor(kBlack);
	fup2->SetLineColor(kBlack);
	fdn2->Draw("same");
	fup2->Draw("same");
      }

    }

  }

 
  gPad->RedrawAxis();
  
}


vector<SHist*> SPlotter::CalcRatios(vector<SHist*> hists)
{
  // build ratios from the array 'hists'
  // by default it is checked if a data histogram exists,
  // which is then divided by the stack
  // steerable: which histograms should be calculated for the ratio

  // first get the basic histograms
  SHist* sstack = SelStack(hists);
  SHist* sdata  = SelData(hists);
  
  vector<SHist*> ratios;

  // TODO: ratio if neither stack nor data exist
  if (!sstack || !sdata){    
    return ratios;
  }
  
  SHist* rd = (SHist*) sdata->Duplicate();
  TH1D*  rdhist = (TH1D*) rd->GetHist();

  // get the denominator: the last element in the stack is the sum of all
  TObjArray* arr = sstack->GetStack()->GetStack();
  TH1D* denom = (TH1D*) arr->At(arr->GetEntries()-1);

  rdhist->Divide(denom);

  // set the error to display only the error on the data
  for (Int_t ibin=1;ibin<denom->GetNbinsX()+1; ++ibin){
    Double_t val = sdata->GetHist()->GetBinContent(ibin);
    Double_t err = sdata->GetHist()->GetBinError(ibin);
    Double_t rel_err = err / val;
    rdhist->SetBinError(ibin, rel_err * rdhist->GetBinContent(ibin) );
  }
  //rdhist->GetYaxis()->SetTitle(rd->GetProcessName() + " / BG");
  rdhist->GetYaxis()->SetTitle("Data / Bkg");
  if (bSingleEPS){
    SingleEPSRatioCosmetics(rdhist);
  } else {
    RatioCosmetics(rdhist);
  }

  // one histogram for the MC statistical error and one for the total error
  SHist* mcerr = new SHist(rdhist);
  mcerr->GetHist()->SetName("MCstat");
  mcerr->SetProcessName("MCstat");
  TH1D* MCstat = (TH1D*)mcerr->GetHist();

  SHist* mctot = new SHist(rdhist);
  mctot->GetHist()->SetName("MCtot");
  mctot->SetProcessName("MCtot");
  TH1D* MCtot = (TH1D*)mctot->GetHist();
  TGraphAsymmErrors* eAsym = new TGraphAsymmErrors();

  for (Int_t ibin=1;ibin<denom->GetNbinsX()+1; ++ibin){

    Double_t val = denom->GetBinContent(ibin);
    Double_t err = denom->GetBinError(ibin);
    MCstat->SetBinContent(ibin,  1.0);
    MCstat->SetBinError(ibin,  err/val);    

    Double_t sys = 0;
    if (m_syserr>0) sys = m_syserr;
    Double_t norm_err = CalcNormErrorForBin(sstack, ibin)/val;

    Double_t tot = TMath::Sqrt(norm_err*norm_err + sys*sys + err/val*err/val);
    MCtot->SetBinContent(ibin, 1.0);
    MCtot->SetBinError(ibin, tot);

    Double_t sys_err_plus = CalcShapeSysErrorForBinFromTheta(sstack, ibin, "up");
    Double_t sys_err_minus = CalcShapeSysErrorForBinFromTheta(sstack, ibin, "down");
    if (sys_err_plus < sys_err_minus){
      Double_t temp = sys_err_plus;
      sys_err_plus = sys_err_minus;
      sys_err_minus = temp;
    }

    Double_t ey_low = TMath::Sqrt(norm_err*norm_err + sys*sys + err/val*err/val + sys_err_minus/val*sys_err_minus/val);
    Double_t ey_up = TMath::Sqrt(norm_err*norm_err + sys*sys +err/val*err/val + sys_err_plus/val*sys_err_plus/val);
    Double_t ex_low = (denom->GetXaxis()->GetBinCenter(ibin)) - (denom->GetXaxis()->GetBinLowEdge(ibin));
    Double_t ex_up =  (denom->GetXaxis()->GetBinUpEdge(ibin))-denom->GetXaxis()->GetBinCenter(ibin);
    eAsym -> SetPoint(ibin, denom->GetXaxis()->GetBinCenter(ibin), 1.); 
    eAsym -> SetPointError(ibin, ex_low, ex_up, ey_low, ey_up); 

    // set error to 0 for empty bins
//    if (bIgnoreEmptyBins && val<0.05 && ibin<15){
    if (bIgnoreEmptyBins && val<0.005){
      //cout << "no MC in bin " << ibin << " lower = " << denom->GetXaxis()->GetBinLowEdge(ibin) << " upper = " << denom->GetXaxis()->GetBinUpEdge(ibin) << endl;
      MCstat->SetBinError(ibin, 0.);
      MCtot->SetBinError(ibin, 0.);
      MCstat->SetBinContent(ibin, 0.);
      MCtot->SetBinContent(ibin, 0.);
      eAsym -> SetPointError(ibin, ex_low, ex_up, 0, 0); 
      eAsym -> SetPoint(ibin, denom->GetXaxis()->GetBinCenter(ibin), 0.); 
    }
   
  }

  //static Int_t VLightGray    = TColor::GetColor( "#eeeeee" );
  static Int_t MLightGray    = TColor::GetColor( "#dddddd" );
  static Int_t LightGray     = TColor::GetColor( "#aaaaaa" );
  //static Int_t Gray          = TColor::GetColor( "#888888" );

  MCstat->SetMarkerStyle(0);
  MCstat->SetMarkerSize(0);
  MCstat->SetLineColor(MLightGray);
  if (bIgnoreEmptyBins){
    MCstat->SetLineColor(kBlack);
    MCstat->SetLineStyle(kDashed);
  }
  MCstat->SetFillColor(MLightGray);

  MCtot->SetMarkerStyle(0);
  MCtot->SetMarkerSize(0);
  MCtot->SetLineColor(LightGray);
  if (bIgnoreEmptyBins){
    MCtot->SetLineColor(kBlack);
    MCtot->SetLineStyle(kDashed);
  } 
  MCtot->SetFillColor(LightGray);

  eAsym->SetMarkerStyle(0);
  eAsym->SetMarkerSize(0);
  eAsym->SetLineColor(LightGray);
  eAsym->SetFillColor(LightGray);

  if (m_shapesys_arr.size()!=0){
    mctot->SetAsymmErrors(eAsym);
  }

  ratios.push_back(mctot);
  ratios.push_back(mcerr);
  ratios.push_back(rd);	  
 
  return ratios;

}

void SPlotter::PlotZScore(vector<SHist*> hists, int ipad)
{
  // plot all histograms in the array

  if (ipad==1) m_rp1->cd();
  if (ipad==2) m_rp2->cd();

  // calculate ratios
  vector<SHist*> ratios = CalcZScore(hists);

  gPad->SetLogx(0);
  gPad->SetLogy(0);

  int ndrawn = 0;
  int nh = ratios.size();
  for (int i=0; i<nh; ++i){
    SHist* rh = ratios[i];
    TString name = rh->GetName();
    if (name.Contains("_lx")) gPad->SetLogx(1);
    if (ndrawn==0) rh->Draw();      
    else rh->Draw("same");
    
    ++ndrawn;
  }

  // some lines to guide the eye
  TH1* h = ratios[0]->GetHist();
  //double xmin = h->GetXaxis()->GetXmin();
  //double xmax = h->GetXaxis()->GetXmax();
  double xmin = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst());
  double xmax = h->GetXaxis()->GetBinUpEdge(h->GetXaxis()->GetLast());
  TLine* zerol = new TLine(xmin, 0, xmax, 0);
  TLine* upl = new TLine(xmin, 3, xmax, 3);
  TLine* downl = new TLine(xmin, -3, xmax, -3);

  zerol->SetLineColor(kBlack);
  upl->SetLineColor(kGray+1);
  downl->SetLineColor(kGray+1);

  upl->SetLineStyle(kDotted);
  downl->SetLineStyle(kDotted);
  
  zerol->Draw();
  upl->Draw();
  downl->Draw();
 
  gPad->RedrawAxis();
  
}

vector<SHist*> SPlotter::CalcZScore(vector<SHist*> hists)
{
  // build ratios from the array 'hists'
  // by default it is checked if a data histogram exists,
  // which is then divided by the stack
  // steerable: which histograms should be calculated for the ratio

  // first get the basic histograms
  SHist* sstack = SelStack(hists);
  SHist* sdata  = SelData(hists);
  
  vector<SHist*> scores;

  // TODO: score if neither stack nor data exist
  if (!sstack || !sdata){    
    return scores;
  }
  
  // observed
  TH1D*  obshist = (TH1D*) sdata->GetHist();

  // expected
  TObjArray* arr = sstack->GetStack()->GetStack();
  TH1D* exphist = (TH1D*) arr->At(arr->GetEntries()-1);

  // error on expected
  SHist experr(exphist);
  experr.GetHist()->SetName("ExpErr");
  experr.SetProcessName("ExpErr");
  TH1D* experrhist = (TH1D*)experr.GetHist();

  // the result
  SHist* zscore = sdata->Duplicate();
  TH1D* zscorehist = (TH1D*)zscore->GetHist();
  zscorehist->GetYaxis()->SetTitle("Z-Score");
  zscore->DrawNoErrorX(true);
  ZScoreCosmetics(zscorehist);

  // to guide the eye
  SHist* mcerr = zscore->Duplicate();
  mcerr->GetHist()->SetName("MCstat");
  mcerr->SetProcessName("MCstat");
  TH1D* MCstat = (TH1D*)mcerr->GetHist();

  SHist* mctot = zscore->Duplicate();
  mctot->GetHist()->SetName("MCtot");
  mctot->SetProcessName("MCtot");
  TH1D* MCtot = (TH1D*)mctot->GetHist();

  for (Int_t ibin=1;ibin<exphist->GetNbinsX()+1; ++ibin){
    Double_t val = exphist->GetBinContent(ibin);
    Double_t staterr = exphist->GetBinError(ibin);

    Double_t sys = 0;
    if (m_syserr>0) sys = m_syserr*val; // take absolute error
    Double_t norm_err = CalcNormErrorForBin(sstack, ibin);

    Double_t tot_err2 = norm_err*norm_err + sys*sys + staterr*staterr;

    Double_t sys_err_plus = CalcShapeSysErrorForBinFromTheta(sstack, ibin, "up");
    Double_t sys_err_minus = CalcShapeSysErrorForBinFromTheta(sstack, ibin, "down");
    // symmetrize
    double sys_err_tot = ( fabs(sys_err_plus) + fabs(sys_err_minus) ) / 2.;

    tot_err2 += sys_err_tot*sys_err_tot;
    
    experrhist->SetBinError(ibin, TMath::Sqrt(tot_err2));   
  }

  // calculate Z-score for each bin
  int ntoys = 20000;
  for (Int_t ibin = 1; ibin<exphist->GetNbinsX()+1; ++ibin){
    Double_t vobs = obshist->GetBinContent(ibin);
    Double_t vexp = exphist->GetBinContent(ibin);
    Double_t verr_exp = experrhist->GetBinError(ibin);
    Int_t nbins = max(int(vexp)*3,50);
    TH1D expdist("expdist","expdist", nbins, -0.5, nbins-0.5); // creates a histogram for the Z-score computation, based on the Expected no of events.
    TRandom3 rand(0);
    for (int i=0; i<ntoys; ++i){
      double x = rand.Gaus(vexp, verr_exp);
      if (x<0) continue;
      expdist.Fill(rand.Poisson(x), TMath::Gaus(x,vexp,verr_exp,true)); // fill a histogram with a poisson dist. weighted with gaussian shape
    }
    expdist.Scale(1./expdist.Integral());

    double integXtoInf = expdist.Integral(expdist.FindFixBin(vobs)+1,expdist.GetNbinsX());
    integXtoInf += expdist.GetBinContent(expdist.FindFixBin(vobs))*0.5; // add half of the bin content because of binning 
    double pvalue = 1-integXtoInf;
    double vzscore = TMath::NormQuantile(pvalue);
    //double zscore = (vobs-vexp)/verr_exp; // naive definition, would also work (and we would not need to dice)   

    // save the histograms for debugging
    TString fname = TString::Format("temp_%i.root", ibin);
    TFile* f = new TFile(fname, "RECREATE");
    expdist.Write();
    f->Write();
    f->Close();

    if (vobs<1){
      vzscore = -999;
    }

    zscorehist->SetBinContent(ibin, vzscore);
    zscorehist->SetBinError(ibin, 0);

    // just to guide the eye: 1sigma and 2sigma bands
    MCstat->SetBinContent(ibin,  0.0);
    MCstat->SetBinError(ibin,  1.0);

    MCtot->SetBinContent(ibin, 0.0);
    MCtot->SetBinError(ibin, 2.0);


  }

  static Int_t MLightGray    = TColor::GetColor( "#dddddd" );
  static Int_t LightGray     = TColor::GetColor( "#aaaaaa" );

  MCstat->SetMarkerStyle(0);
  MCstat->SetMarkerSize(0);
  MCstat->SetLineColor(MLightGray);
  MCstat->SetFillColor(MLightGray);

  MCtot->SetMarkerStyle(0);
  MCtot->SetMarkerSize(0);
  MCtot->SetLineColor(LightGray);
  MCtot->SetFillColor(LightGray);

  scores.push_back(mctot);
  scores.push_back(mcerr);
  scores.push_back(zscore);
 
  return scores;
}


void SPlotter::DrawLegend(vector<SHist*> hists)
{
  // draw a legend  

  int narr = hists.size();
  float yfrac = 0.06;
  if (!bPlotRatio) yfrac = 0.05;

  float top = 0.89;
  if (!bPlotRatio && bDrawLumi) top = 0.86;
  if (bSingleEPS){
    top = 0.88;
    if (bPlotRatio) top = 0.88;
    //if (bDrawLumi) top = 0.88;
  }
  float ysize = yfrac*narr;
  float xleft = 0.7;
  if (bSingleEPS) xleft = 0.55;
  float xright = 0.88;
  if (!bPortrait){
    top = 0.99;
    ysize = 0.07*narr;
    xleft = 0.72;
    xright = 0.96;
  }
  
  TString hname = hists[0]->GetName();
  if (hname.Contains("tau32")){
    xleft = 0.23;
    xright = 0.46;
    top = 0.78;
  }

  ysize = 0.06*4;
  TLegend *leg = new TLegend(xleft,top-ysize,xright,top, NULL,"brNDC");
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  if (bSingleEPS) leg->SetTextSize(0.05);

  // do the ordering by hand
  cout << "warning: plotting legend - check if ordering is ok (by hand, line 1530 in SPlotter!" << endl;
  //Int_t j[] = {0, 4, 1, 2, 3, 5, 6}; // dilepton
  //Int_t j[] = {0, 2, 1, 3, 4}; // CMSTT
  Int_t j[] = {0, 7, 1, 2, 3, 4, 5, 6, 8, 9}; // l+jets case

  //Int_t j[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}; // general case


  for (Int_t i=0; i<narr; ++i){

    //SHist* sh = hists[j[i]];
    SHist* sh = hists[i];
    if (sh->IsStack()) continue;

    TString legname = TString::Format("leg_entry_%i",i);
    TString legtitle = sh->GetLegName();
    TLegendEntry* entry = NULL;
    int marker = sh->GetHist()->GetMarkerStyle();
    int lstyle = sh->GetHist()->GetLineStyle();

    if (legtitle=="W+light") continue; //legtitle = "W(#rightarrow l #nu)+jets";
    if (legtitle=="W+c") continue;
    if (legtitle=="W+b") continue;
    if (legtitle=="single-top") legtitle = "Other";
    if (legtitle=="Z+jets") continue;
    if (legtitle=="diboson") continue;

    legtitle.ReplaceAll("TeV 1%", "TeV, 1%");
    legtitle.Prepend(" ");

    if (marker>0){
      entry = leg->AddEntry(legname, legtitle, "pe");
      entry->SetLineWidth(2);
      entry->SetLineColor(sh->GetHist()->GetLineColor());
      entry->SetMarkerColor(sh->GetHist()->GetLineColor());
      entry->SetMarkerStyle(marker);
      entry->SetMarkerSize(1.0);
      if (bSingleEPS) entry->SetMarkerSize(0.8);
    } else {
      
      if (sh->IsUsedInStack()){
      	entry = leg->AddEntry(legname, legtitle, "f");
      	entry->SetLineWidth(1);
      	entry->SetLineColor(sh->GetHist()->GetLineColor());
      	entry->SetFillColor(sh->GetHist()->GetLineColor());
      	entry->SetFillStyle(1001);

      } else {
      	entry = leg->AddEntry(legname, legtitle, "l");
      	entry->SetLineColor(sh->GetHist()->GetLineColor());
      	entry->SetMarkerStyle(0);
      	entry->SetMarkerSize(0);
      	entry->SetMarkerColor(sh->GetHist()->GetLineColor());
      	entry->SetLineWidth(2);
      	entry->SetLineStyle(lstyle);
      }
      entry->SetTextAlign(12);
      //entry->SetTextColor(fSampleColors.At(i));
    }
  }
  leg->Draw();


  // auxiliary text
  TString name = hists[0]->GetName();
  cout << "name = " << name << endl;
  TString infotext;
  if (name.Contains("ele_t0b0", TString::kIgnoreCase)) infotext = "e+jets, 0 t tag, 0 b tag";
  if (name.Contains("ele_t0b1", TString::kIgnoreCase)) infotext = "e+jets, 0 t tag, 1 b tag";
  if (name.Contains("ele_t1b0", TString::kIgnoreCase)) infotext = "e+jets, 1 t tag";

  if (name.Contains("muo_t0b0", TString::kIgnoreCase)) infotext = "#mu+jets, 0 t tag, 0 b tag";
  if (name.Contains("muo_t0b1", TString::kIgnoreCase)) infotext = "#mu+jets, 0 t tag, 1 b tag";
  if (name.Contains("muo_t1b0", TString::kIgnoreCase)) infotext = "#mu+jets, 1 t tag";

  if (name.Contains("lepton_0top0btag")) infotext = "e/#mu+jets, 0 t tag, 0 b tag";
  if (name.Contains("lepton_0top1btag")) infotext = "e/#mu+jets, 0 t tag, 1 b tag";
  if (name.Contains("lepton_1top")) infotext = "e/#mu+jets, 1 t tag";

  if (name.Contains("ee")) infotext = "ee";
  if (name.Contains("mumu")) infotext = "#mu#mu";
  if (name.Contains("emu")) infotext = "e#mu";

  if (name == "btag0") infotext = "|#Deltay| < 1.0; 0 b tag";
  if (name == "btag1") infotext = "|#Deltay| < 1.0; 1 b tag";
  if (name == "btag2") infotext = "|#Deltay| < 1.0; 2 b tag";
  if (name == "btag3") infotext = "|#Deltay| > 1.0; 0 b tag";
  if (name == "btag4") infotext = "|#Deltay| > 1.0; 1 b tag";
  if (name == "btag5") infotext = "|#Deltay| > 1.0; 2 b tag";

  if (name == "httbtag0") infotext = "H_{T} > 800 GeV, 0 b tag (low-mass)";
  if (name == "httbtag1") infotext = "H_{T} > 800 GeV, 1 b tag (low-mass)";
  if (name == "httbtag2") infotext = "H_{T} > 800 GeV, 2 b tag (low-mass)";
  if (name == "mjhttbtag0") infotext = "H_{T} < 800 GeV, 0 b tag (low-mass)";
  if (name == "mjhttbtag1") infotext = "H_{T} < 800 GeV, 1 b tag (low-mass)";
  if (name == "mjhttbtag2") infotext = "H_{T} < 800 GeV, 2 b tag (low-mass)";

  TLatex *text1 = new TLatex(3.5, 24, infotext);
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetX(0.20);
  text1->SetTextFont(42);
  text1->SetTextSize(0.06);
  text1->SetY(0.995);

  if (name.BeginsWith("btag") || name.BeginsWith("htt") || name.BeginsWith("mjhtt")){
    text1->SetX(0.19);
    text1->SetTextSize(0.053);
    text1->SetY(0.99);
  }
  if (name.BeginsWith("el_") || name.BeginsWith("mu_") || name.BeginsWith("lepton_")){
    text1->SetY(0.995);
  }
  text1->Draw();
  
}


void SPlotter::DrawLumi(double lumi)
{

  TString infotext;
  if (lumi < 0)
    infotext = TString::Format("%3.1f fb^{-1} (13 TeV)", m_lumi);
  else 
    infotext = TString::Format("%3.1f fb^{-1} (13 TeV)", lumi);

  TLatex *text1 = new TLatex(3.5, 24, infotext);
  text1->SetNDC();
  text1->SetTextAlign(33);
  text1->SetX(0.95);
  text1->SetTextFont(42);
  if (bPlotRatio){ 
    text1->SetTextSize(0.06);
    text1->SetY(1.);
  } else {
    text1->SetTextSize(0.045);
    text1->SetY(1.);
  }
  text1->Draw();

  if (bForPublication || bForPrelim){
    TString cmstext = "CMS";
    TLatex *text2 = new TLatex(3.5, 24, cmstext);
    text2->SetNDC();
    text2->SetTextAlign(13);
    text2->SetX(0.24);
    text2->SetTextFont(62);
    if (bPlotRatio){ 
      text2->SetTextSize(0.08);
      text2->SetY(0.87);
    } else {
      text2->SetTextSize(0.05);
      text2->SetY(0.87);
    }
    text2->Draw();
  }

  if (bForPrelim){
    TString preltext = "Preliminary";
    TLatex *text3 = new TLatex(3.5, 24, preltext);
    text3->SetNDC();
    text3->SetTextAlign(13);
    text3->SetX(0.24);
    text3->SetTextFont(52);
    if (bPlotRatio){ 
      text3->SetTextSize(0.06);
      text3->SetY(0.78);
    } else {
      text3->SetTextSize(0.035);
      text3->SetY(0.78);
    }
    text3->Draw();
  }
  
}

void SPlotter::DoCosmetics(vector<SHist*> hists)
{

  // loop over all histograms and make them pretty
  int nh = hists.size();
  for (int i=0; i<nh; ++i){
    SHist* sh = hists[i];
    if (sh->IsStack()) continue;
    GeneralCosmetics(sh->GetHist());
    if (bPortrait) PortraitCosmetics(sh->GetHist());
    if (!bPortrait) LandscapeCosmetics(sh->GetHist());
    if (bSingleEPS) SingleEPSCosmetics(sh->GetHist());
    if (sh->IsYieldPlot()) YieldCosmetics(sh->GetHist());
  }

}

SHist* SPlotter::SelStack(vector<SHist*> hists)
{
  // select the stack histogram from the array
  int narr = hists.size();
  SHist* h = NULL;
  for (int i=0; i<narr; ++i){
    if (hists[i]->IsStack()){
      h=hists[i];
      break;
    }
  }
  return h;
}

SHist* SPlotter::SelData(vector<SHist*> hists)
{
  // select the data histogram from the array
  int narr = hists.size();
  SHist* h = NULL;
  TString process;
  for (int i=0; i<narr; ++i){
    process = hists[i]->GetProcessName();
    if (process.Contains("data", TString::kIgnoreCase)){
      h = hists[i];
      break;
    }
  }
  return h;
}

bool SPlotter::SetMinMax(vector<SHist*> hists)
{
  // set minimum and maximum of all histograms
  int narr = hists.size();
  TString name = hists[0]->GetName();
  double max = 0;
  double min = FLT_MAX;
  for (int i=0; i<narr; ++i){
    if (max<hists[i]->GetMaximum()) max = hists[i]->GetMaximum();
    if (min>hists[i]->GetMinimum(1e-6)) min = hists[i]->GetMinimum(1e-6);
    double imin = hists[i]->GetMinimum(1e-10);
    if (min>imin){
      if (imin>1e-10){
	      min = imin;
      }
    }
  }

  bool isok = true;
  if (max<1e-6){
    isok = false;
    return isok;
  }

  bool islog = false;
  double uscale = 1.2;
  if (name.Contains("_lxy") || name.Contains("_ly") || bPlotLogy ){
    islog = true;
    uscale = 12.;
  }

  if (name == "el_1top_mttbar"){
    uscale = 20.;
  }

  if (name == "mu_1top_mttbar"){
    uscale = 20.;
  }

  if (name.Contains("lepton")){
    uscale = 40;
  }

  if (name == "lepton_1top_mttbar"){
    uscale = 50;
    min = 0.4;
  }

  if (name.Contains("mjhtt")){
    uscale = 50;
    min = 0.5;
  } else if (name.Contains("htt")){
    uscale = 30;
    min = 0.5;
  }

  for (int i=0; i<narr; ++i){
    SHist* h = hists[i];
    if (h->IsStack()){ 
      if (!islog){
	h->GetStack()->SetMinimum(0.0011);
      } else { 
	if (min>1e-10){
	  if (min<0.1){
	    h->GetStack()->SetMinimum(0.04);
	  } else {
	    h->GetStack()->SetMinimum(min);
	  }
	}
      }
      h->GetStack()->SetMaximum(uscale*max);
    } else {
      if (!islog){ 
	h->GetHist()->SetMinimum(0.0011);
      } else { 
	if (min>1e-10){
	  if (min<0.1){
	    h->GetHist()->SetMinimum(0.04);
	  } else {
	    h->GetHist()->SetMinimum(min);
	  }
	}
      }
      h->GetHist()->SetMaximum(uscale*max);
    }
  }  

  return isok;
}

void SPlotter::SetLogAxes(vector<SHist*> hists)
{
  // set log axes 
  TString name = hists[0]->GetName();
  gPad->SetLogx(0);
  gPad->SetLogy(0);
  if (name.Contains("_lxy")){
    gPad->SetLogx(1);
    gPad->SetLogy(1);
  } else if (name.Contains("_lx")){
    gPad->SetLogx(1);
  } else if (name.Contains("_ly")){
    gPad->SetLogy(1);
  } else {
    // do nothing, all fine
  }

  // override if requested
  if (bPlotLogy){
    gPad->SetLogy(1);
  }

  return;
}

int SPlotter::GetCurrentPad(int np)
{
  // get the current pad, depending on the number of 
  // already processed plots
  int ipad = 1;
  int rest = np%2;
  if (rest==0) ipad=2;      
  return ipad;
}

vector<SHist*> SPlotter::GetHistsAtIndex(std::vector<TObjArray*> histarr, int index)
{
  // fill an array with histograms at position index. 

  vector<SHist*> hists;
  int narr = histarr.size();
  for (int i=0; i<narr; ++i){
    SHist* hist = (SHist*)histarr[i]->At(index);
    hists.push_back(hist);
  }
  
  return hists;

}

vector<SHist*> SPlotter::GetPlottableHists(std::vector<TObjArray*> histarr, int index)
{
  // fill an array with plottable histograms at position index. 
  // if the first histogram in the array should be plotted (DoPlot flag), 
  // then take at first position in the array
  // otherwise look for the stack and plot it first
  // only then all other histograms are added

  if (debug) cout << "\nSPlotter: Collecting plottable hists for index " << index << endl;
  vector<SHist*> hists;
  bool gotstack = false;
  SHist* hist = (SHist*)histarr[0]->At(index);

  // check if the histogram is a 2D or 3D histogram, 
  // plotting not supported yet, to come
  if (hist->GetHist()->InheritsFrom(TH2::Class())){
    if (debug) cout << "Hist inherits from TH2, return without adding any to the array " << endl;
    return hists;
  }

  TString name = hist->GetName();
  TString process = hist->GetProcessName();
  if (process.Contains("data",TString::kIgnoreCase) 
      && name.Contains("_perlumibin", TString::kIgnoreCase)){
    hist->SetIsYieldPlot(true);
    hists.push_back(hist);
    return hists;
  }

  if (hist->DoDraw()){ // take first hist
    hists.push_back(hist);
    gotstack = false;
    if (debug) cout << "Adding hist " << hist->GetHist()->GetName()
		    << " from process " << hist->GetProcessName() 
		    << " and directory " << hist->GetDir() << " to array." << endl;
  } else { // try if stack exists
    TObjArray* stacks = GetStacks(histarr);
    if (stacks){
      hist = (SHist*)stacks->At(index);
      hists.push_back(hist);
      gotstack = true;
      if (debug) cout << "Adding stack " << hist->GetStack()->GetName()
		      << " from process " << hist->GetProcessName() 
		      << " and directory " << hist->GetDir() << " to array." << endl;
    }
  }
  
  // loop over the rest and add them to the array
  int narr = histarr.size();
  for (int i=1; i<narr; ++i){

    SHist* hist = (SHist*)histarr[i]->At(index);

    if (hist->DoDraw()){ // take it if it should be drawn

      if (hist->IsStack()){
	if (!gotstack){ // take the stack only if not already added
	  hists.push_back(hist);
	  if (debug) cout << "Adding stack " << hist->GetStack()->GetName()
			  << " from process " << hist->GetProcessName() 
			  << " and directory " << hist->GetDir() << " to array." << endl;
	}
      } else { // take the histogram if it's not the stack hist
	hists.push_back(hist);
	if (debug) cout << "Adding hist " << hist->GetHist()->GetName()
			<< " from process " << hist->GetProcessName() 
			<< " and directory " << hist->GetDir() << " to array." << endl;
      }
    }
  }

  if (debug) cout << "SPlotter: Done with collecting plottable hists for index " 
		  << index << ", got " << hists.size() << " histograms" << endl;
  
  return hists;

}

void SPlotter::DrawPageNum()
{
  
  m_can->cd();
  TPaveText* text; 
  TString s;
  s.Form("%i",m_page);
  if (bPortrait){
    text = new TPaveText(0.93, 0.00, 0.97, 0.03, "NDC");
  } else {
    text = new TPaveText(0.03,0.00, 0.06, 0.03, "NDC");
  }
  text->SetBorderSize(0);
  text->SetFillColor(0);
  text->AddText(s.Data());
  text->Draw("same");
  
}

void SPlotter::GeneralCosmetics(TH1* hist)
{
  // set Y-axis title
  TString ytitle = "Events";
  Double_t w = hist->GetXaxis()->GetBinWidth(1);
  if (w==99) w = 100;
  TString title = hist->GetTitle();
  if (title.Contains("GeV")){
    ytitle = TString::Format("Events / %i GeV", (Int_t) w);
  }

  hist->GetYaxis()->SetTitle(ytitle);  
  
  // set X-axis title
  hist->GetXaxis()->SetTitle(hist->GetTitle()); 

  hist->SetTitle("");

  if (bShapeNorm) {
    hist->GetYaxis()->SetTitle("#DeltaN/N");
  }

  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetLabelFont(42);

  hist->SetLineWidth(2);

}

void SPlotter::PortraitCosmetics(TH1* hist)
{
  // top histogram of the ratio plot
  if (bPlotRatio){
    
    // x-axis
    hist->GetXaxis()->SetTickLength(0.05);

    // y-axis
    hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetLabelSize(0.062);
    hist->GetYaxis()->SetLabelOffset(0.01);   
    hist->GetYaxis()->SetTitleOffset(0.8);
    hist->GetYaxis()->SetTickLength(0.02);
  
  // only this histogram
  } else {

    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetLabelOffset(0.008);
    hist->GetXaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.2);
    
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.045);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelOffset(0.011);

  }
  
}

void SPlotter::SingleEPSCosmetics(TH1* hist)
{

  hist->SetMarkerSize(1.0);

  // top histogram of the ratio plot
  if (bPlotRatio){

    // x-axis
    hist->GetXaxis()->SetTickLength(0.05);
    hist->GetXaxis()->SetLabelSize(0.0);
    hist->GetXaxis()->SetTitleSize(0.0);
    hist->GetXaxis()->SetTitleOffset(1.0);

    // y-axis
    hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetLabelSize(0.062);
    hist->GetYaxis()->SetLabelOffset(0.01);   
    hist->GetYaxis()->SetTitleOffset(1.15);
    hist->GetYaxis()->SetTickLength(0.02);
  
  // only this histogram
  } else {

    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetLabelOffset(0.008);
    hist->GetXaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.2);
    
    hist->GetYaxis()->SetTitleOffset(1.8);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.045);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelOffset(0.011);

  }
  
}

void SPlotter::SingleEPSRatioCosmetics(TH1* hist)
{

  hist->GetYaxis()->SetRangeUser(0.3, 1.7);
  //hist->GetYaxis()->SetRangeUser(0.05, 1.95);
  hist->SetMarkerSize(0.7);

  // cosmetics for portrait mode 
  if (bPortrait){
    hist->SetTitle("");
    
    // x-axis
    hist->GetXaxis()->SetLabelSize(0.12);
    hist->GetXaxis()->SetTickLength(0.08);
    hist->GetXaxis()->SetTitleSize(0.14);
    hist->GetXaxis()->SetTitleOffset(1.15);
    hist->GetXaxis()->SetNdivisions(1005);
	  
    // y-axis
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleSize(0.14);
    hist->GetYaxis()->SetTitleOffset(0.57);
    hist->GetYaxis()->SetLabelSize(0.11);
    //hist->GetYaxis()->SetNdivisions(210);
    hist->GetYaxis()->SetNdivisions(505);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelOffset(0.011);

    // cosmetics for landscape mode 
  } else {
    
    hist->SetTitle("");
    hist->SetTitleOffset(1.1, "X");
    hist->SetTitleOffset(0.5, "Y");
    hist->SetLabelOffset(0.02, "X");
    hist->SetLabelOffset(0.01, "Y");
	  
    hist->GetXaxis()->SetLabelSize(0.14);
    hist->GetXaxis()->SetTickLength(0.07);
    hist->GetXaxis()->SetTitleSize(0.15);
	  
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleSize(0.11);
    hist->GetYaxis()->SetLabelSize(0.12);
    hist->GetYaxis()->SetNdivisions(505);
    hist->GetYaxis()->SetTickLength(0.03);
    
  }

}


void SPlotter::YieldCosmetics(TH1* hist)
{
  // cosmetics for the lumi yield histogram
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetLabelOffset(0.008);
    hist->GetXaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.2);
    
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.045);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelOffset(0.011);

    if (bSingleEPS){
      hist->GetYaxis()->SetTitleOffset(1.8);
      hist->GetYaxis()->SetTitleSize(0.055);
      hist->GetYaxis()->SetLabelSize(0.05);
      hist->GetYaxis()->SetTickLength(0.02);
      hist->GetYaxis()->SetLabelOffset(0.011);
    }

    hist->GetXaxis()->SetTitle("integrated luminosity [fb^{-1}]");
    double dlum = hist->GetXaxis()->GetBinWidth(1);
    TString xtit = TString::Format("events per %3.1f fb^{-1}", dlum);
    hist->GetYaxis()->SetTitle(xtit);
}



void SPlotter::LandscapeCosmetics(TH1* hist)
{

  // FIXME: need to define sensible style
  hist->GetXaxis()->SetLabelSize(0.12);
  hist->GetXaxis()->SetTickLength(0.08);
  hist->GetXaxis()->SetTitleSize(0.12);
  hist->GetXaxis()->SetTitleOffset(1.25);
	  

}

void SPlotter::RatioCosmetics(TH1* hist)
{

  hist->GetYaxis()->SetRangeUser(0.3, 1.7);
  //hist->GetYaxis()->SetRangeUser(0.05, 1.95);
  hist->SetMarkerSize(0.7);

  // cosmetics for portrait mode 
  if (bPortrait){
    hist->SetTitle("");
    
    // x-axis
    hist->GetXaxis()->SetLabelSize(0.12);
    hist->GetXaxis()->SetTickLength(0.08);
    hist->GetXaxis()->SetTitleSize(0.12);
    hist->GetXaxis()->SetTitleOffset(1.25);
	  
    // y-axis
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleSize(0.12);
    hist->GetYaxis()->SetTitleOffset(0.46);
    hist->GetYaxis()->SetLabelSize(0.11);
    //hist->GetYaxis()->SetNdivisions(210);
    hist->GetYaxis()->SetNdivisions(505);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelOffset(0.011);

    // cosmetics for landscape mode 
  } else {
    
    hist->SetTitle("");
    hist->SetTitleOffset(1.1, "X");
    hist->SetTitleOffset(0.5, "Y");
    hist->SetLabelOffset(0.02, "X");
    hist->SetLabelOffset(0.01, "Y");
	  
    hist->GetXaxis()->SetLabelSize(0.14);
    hist->GetXaxis()->SetTickLength(0.07);
    hist->GetXaxis()->SetTitleSize(0.15);
	  
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleSize(0.11);
    hist->GetYaxis()->SetLabelSize(0.12);
    hist->GetYaxis()->SetNdivisions(505);
    hist->GetYaxis()->SetTickLength(0.03);
    
  }

}

void SPlotter::ZScoreCosmetics(TH1* hist)
{

  // use the usual ratio cosmetics first
  if (bSingleEPS){
    SingleEPSRatioCosmetics(hist);
  } else {
    RatioCosmetics(hist);
  }

  hist->GetYaxis()->SetRangeUser(-4.5, 4.5);

  return;

}

void SPlotter::ShapeNormalise(std::vector<SHist*> hists)
{
  for (unsigned int i=0; i<hists.size(); ++i){
    hists[i]->NormaliseToArea();
  }
}

