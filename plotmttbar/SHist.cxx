#include "SHist.h"
#include <TList.h>
#include <TExec.h>
#include <iostream>
#include <stdlib.h>

using namespace std;

ClassImp(SHist)

SHist::SHist(TH1* hist)
{
  // constructor, create a new histogram such that the 
  // SHist object is the owner of it

  if (!hist){
    cerr << "Called constructor of SHist (hist) with NULL pointer, abort." << endl;
    exit(EXIT_FAILURE);
  }

  m_hist = (TH1*) hist->Clone();
  m_stack = NULL;
  m_asymme = NULL;
  m_weight = 1.;
  m_is_stack = false;
  m_is_used_in_stack = false;
  m_draw_marker = true;
  m_draw = true;
  m_is_yield_plot = false;
  m_draw_noxerr = true; 
  
}

SHist::SHist(THStack* stack)
{
  // constructor with a THStack, creates a new stack such that
  // the SHist object is the owner of it

  if (!stack){
    cerr << "Called constructor of SHist (stack) with NULL pointer, abort." << endl;
    exit(EXIT_FAILURE);
  }

  m_stack = (THStack*)stack->Clone();
  m_hist = NULL;
  m_asymme = NULL;
  m_weight = 1.;
  m_is_stack = true;
  m_draw_marker = false;
  m_draw = true;
  m_is_yield_plot = false;
  m_draw_noxerr = true;
  
}

SHist::~SHist()
{
  if (m_hist){
    m_hist->Delete();
    m_hist = NULL;
  }
  if (m_stack){
    m_stack->Delete();
    m_stack = NULL;
  }
  if (m_asymme){
    delete m_asymme;
  }
}

SHist* SHist::Duplicate()
{
  SHist* s = NULL;
  if (m_hist) s = new SHist(m_hist);
  if (m_stack) s = new SHist(m_stack);
  if (m_asymme){
    TGraphAsymmErrors* t = new TGraphAsymmErrors(*m_asymme);
    s->SetAsymmErrors(t);
  }
  s->SetDir(m_dir);
  s->SetProcessName(m_process);
  s->SetLegName(m_leg_name);
  s->SetWeight(m_weight);
  s->SetDoDraw(m_draw);
  s->SetDrawMarker(m_draw_marker);
  s->SetIsUsedInStack(m_is_used_in_stack);
  return s;
}

const char* SHist::GetName() const
{
  if (m_is_stack) return m_stack->GetName();
  else return m_hist->GetName();
}

void SHist::SetName(TString name)
{
  if (m_is_stack) m_stack->SetName(name);
  else m_hist->SetName(name);
}

TH1* SHist::GetHist()
{
  return m_hist;
}

THStack* SHist::GetStack()
{
  return m_stack;
}

void SHist::SetProcessName(TString name)
{
  m_process = name;
}

TString SHist::GetProcessName()
{
  return m_process;
}

void SHist::SetLegName(TString name)
{
  m_leg_name = name;
  m_leg_name.ReplaceAll("SPACE", " ");
  m_leg_name.ReplaceAll("_", " ");
  m_leg_name.ReplaceAll("ttbar", "t#bar{t}");
  m_leg_name.ReplaceAll("~", ",");
  m_leg_name.ReplaceAll("[", "{");
  m_leg_name.ReplaceAll("]", "}");
  m_leg_name.ReplaceAll("_", " ");
  m_leg_name.ReplaceAll("TTbar", "t#bar{t}");
}

TString SHist::GetLegName()
{
  return m_leg_name;
}

void SHist::SetDir(TString name)
{
  m_dir = name;
}

TString SHist::GetDir()
{
  return m_dir;
}


void SHist::SetIsStack(bool flag)
{
  m_is_stack = flag;
}

bool SHist::IsStack()
{
  return m_is_stack;
}

void SHist::SetIsUsedInStack(bool flag)
{
  m_is_used_in_stack = flag;
}

bool SHist::IsUsedInStack()
{
  return m_is_used_in_stack;
}

void SHist::SetIsYieldPlot(bool flag)
{
  m_is_yield_plot = flag;
}

bool SHist::IsYieldPlot()
{
  return m_is_yield_plot;
}

void SHist::SetWeight(double weight)
{
  m_weight = weight;
}

double SHist::GetWeight()
{
  return m_weight;
}

void SHist::SetUnc(double unc, int i)
{
  // set the uncertainty on the normalisation of the sample
  // usually just one number (array at position 0), 
  // but can be an array in case the SHist is a stack
  if (i<m_unc_arr.GetSize()+1) m_unc_arr.Set(m_unc_arr.GetSize()+1);
  m_unc_arr.AddAt(unc, i);
}

double SHist::GetUnc(int i)
{
  if (m_unc_arr.GetSize()==0) return 0.;
  return m_unc_arr.At(i);
}

void SHist::SetDrawMarker(bool flag)
{
  m_draw_marker = flag;
}

bool SHist::DrawMarker()
{
  if (m_draw_marker) return true;
  else return false;
}

bool SHist::DrawLine()
{
  if (!m_draw_marker) return true;
  else return false;
}

void SHist::SetDoDraw(bool flag)
{
  m_draw = flag;
}

bool SHist::DoDraw()
{
  return m_draw;
}

void SHist::DrawNoErrorX(bool flag)
{
  m_draw_noxerr = flag;
}

bool SHist::GetNoErrorX()
{
  return m_draw_noxerr;
}


void SHist::Draw(Option_t *option)
{
  // here is all the fun: Draw the histogram

  TString dopt(option);

  if (m_process.Contains("MCstat") || m_process.Contains("MCtot")){
    //TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");
    //setex2->Draw();
    m_hist->DrawCopy("E2 " + dopt);
    if (m_asymme){
      m_asymme->Draw("E2 " + dopt);
      return;
    }
    // this was previously used to draw a line at 1, do this rather in SPlotter::PlotRatios
    //m_hist->SetFillColor(0);
    //m_hist->DrawCopy("HIST same");
    return;
  }

   //if (m_asymme){
    //  m_asymme->Draw("E2 " + dopt);
    //  return;
    //} else {
      //m_hist->DrawCopy("E2 " + dopt);
      //m_hist->SetFillColor(0);
      //m_hist->DrawCopy("HIST same");
      //return;
      //}
  
  if (IsStack()){
    m_stack->Draw("HIST " + dopt);
  } else {
    // skip 2d hists for now
    if (m_hist->InheritsFrom("TH2")) return;

    if (m_hist->GetMarkerStyle()>0){

      if (m_draw_noxerr){
      	m_hist->Draw("E0 X0 " + dopt);
      } else {
      	for (Int_t i=1; i<m_hist->GetNbinsX()+1; ++i){
      	  if (m_hist->GetBinContent(i)==0) m_hist->SetBinError(i,0);
      	}
	      m_hist->Draw("E0 X0 " + dopt);
      }
    } else {
      m_hist->Draw("HIST " + dopt);
    }
  }

}

double SHist::GetMinimum(double minval)
{
  // get minimum
  if (IsStack()){
    return m_stack->GetMinimum();
  } else {
    return m_hist->GetMinimum(minval);
  }

}

double SHist::GetMaximum()
{
  // get maximum
  if (IsStack()){
    return m_stack->GetMaximum();
  } else {
    return m_hist->GetMaximum();
  }

}

void SHist::SetAsymmErrors(TGraphAsymmErrors* as)
{
  // give the SHist an asymmetric error
  m_asymme = as;  
}

TGraphAsymmErrors* SHist::GetAsymmErrors()
{
  // return the asymmetric errors
  return m_asymme;
}

void SHist::NormaliseToArea()
{
  // noramlise the histogram according to its area
  // the resulting histogram has area = 1

  double area = 0.;
  int ibeg = 1; // integrate from bin 1 
  int iend = 1;

  if (IsStack()){
    TList* hists = m_stack->GetHists();
    // calculate total area
    for (int i=0; i<hists->GetSize(); ++i){
      TH1* h = (TH1*) hists->At(i);
      iend = h->GetNbinsX();
      area += h->Integral(ibeg,iend);
    }
    // scale each histogram
    if (area>0){
      for (int i=0; i<hists->GetSize(); ++i){
	TH1* h = (TH1*) hists->At(i);
	h->Scale(1./area);
      }
      m_stack->Modified();
    }

  } else {

    iend = m_hist->GetNbinsX();
    area = m_hist->Integral(ibeg,iend);
    if (area>0){
      m_hist->Scale(1./area);
    }
  }

}
