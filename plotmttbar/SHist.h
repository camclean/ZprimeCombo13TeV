#ifndef SHIST_H
#define SHIST_H

#include <TH1.h>
#include <THStack.h>
#include <TGraphAsymmErrors.h>
#include <TH2.h>
#include <TString.h>
#include <TObject.h>

class SHist : public TObject
{

 public:
  
  SHist(TH1* hist);
  SHist(THStack* stack);
  ~SHist();

  const char* GetName() const;
  void SetName(TString name);

  SHist* Duplicate();

  void SetProcessName(TString);
  TString GetProcessName();
  
  void SetLegName(TString);
  TString GetLegName();

  void SetDir(TString);
  TString GetDir();
  
  void SetIsUsedInStack(bool);
  bool IsUsedInStack();

  void SetIsStack(bool);
  bool IsStack();

  void SetIsYieldPlot(bool);
  bool IsYieldPlot();

  void SetWeight(double);
  double GetWeight();

  void SetUnc(double unc, int i=0);
  double GetUnc(int i=0);

  void SetDoDraw(bool);
  bool DoDraw();

  void SetDrawMarker(bool);
  bool DrawMarker();
  bool DrawLine();

  TH1* GetHist();
  THStack* GetStack();

  double GetMinimum(double minval);
  double GetMaximum(); 

  void NormaliseToArea();

  void DrawNoErrorX(bool);
  bool GetNoErrorX();

  virtual void Draw(Option_t *option="");

  void SetAsymmErrors(TGraphAsymmErrors* as);
  TGraphAsymmErrors* GetAsymmErrors();

 private:
  TH1* m_hist;
  THStack* m_stack;
  TGraphAsymmErrors* m_asymme;
  double m_weight;
  TArrayD m_unc_arr;
  TString m_process;
  TString m_leg_name;
  TString m_dir;
  bool m_is_stack;
  bool m_is_used_in_stack;
  bool m_draw_marker;
  bool m_draw;
  bool m_is_yield_plot;
  bool m_draw_noxerr;

  ClassDef(SHist,0)  // SFrame histograms

};

#endif
