#include <iostream>
#include <stdio.h>
#include <sys/stat.h>
#include <TSystem.h>
#include <TObjString.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TClass.h>
#include <TMath.h>

#include "FileParser.h"
#include "SHist.h"

using namespace std;

FileParser::FileParser()
{
  m_file = NULL;
  m_hists = NULL;
  m_shapeSys = NULL;
  debug = false;
  m_do_cumulative = true;
}

FileParser::~FileParser()
{
  if (m_file){
    CloseFile();
  }
}

void FileParser::CloseFile()
{
  if (m_file){
    m_file->Close();
    delete m_file;
    m_file = NULL;
  }
}


void FileParser::Clear()
{
  // clear all arrays, but do not close file
  if (m_hists) m_hists->Clear();
  if (m_shapeSys) m_shapeSys->Clear();
}

bool FileParser::FileExists(TString filename)
{
  struct stat buf;
  if (stat(filename.Data(), &buf) != -1)
    {
      return true;
    }
  return false;
}


void FileParser::OpenFile(TString fname, TString cyclename)
{
  // open the root files with names given in the TObjArray

  if (m_file != NULL){
    cerr << "FileParser::OpenFile: Can not open new file, since file " 
	 << m_file->GetName() << " is still in memory. Abort." << endl;
    exit(EXIT_FAILURE);
  }

  if (cyclename.Sizeof()!=0){
    TString Prefix(cyclename);
    Prefix.Append(".");
    fname.Prepend(Prefix);
  }  

  // check if name consists of a wildcard, if so use hadd to combine histograms
  if (fname.Contains("*")){
    TString target(fname);
    target.ReplaceAll("*","");

    // check if target exists, delete if yes
    if (FileExists(target)){
      if (debug) cout << "Target exists, removing file " << target << endl;
      remove(target);
    }
      
    TString command = "hadd " + target + " " + fname;
    int res = gSystem->Exec(command);
    if(res != 0){
        cerr << "hadd command '" << command << "' failed with error code " << res << ", aborting." << endl;
        exit(EXIT_FAILURE);
    }
    fname = target;
  }

  // check if name consists of a wildcard, if so use hadd to combine histograms
  if (fname.Contains("?")){
    TString target(fname);
    target.ReplaceAll("?","");

    // check if target exists, delete if yes
    if (FileExists(target)){
      if (debug) cout << "Target exists, removing file " << target << endl;
      remove(target);
    }
      
    TString command = "hadd " + target + " " + fname;
    int res = gSystem->Exec(command);
    if(res != 0){
        cerr << "hadd command '" << command << "' failed with error code " << res << ", aborting." << endl;
        exit(EXIT_FAILURE);
    }
    fname = target;
  }

  if (debug) cout << "Opening file with name " << fname << "..." << endl;
  m_file = new TFile(fname, "READ");
  if (debug){
    cout << "... success! pointer = " << m_file << endl;
    cout << "name = " << m_file << endl;
    cout << " is open? " << m_file->IsOpen() << endl;
    m_file->ls();
  }
    
  if (!m_file->IsOpen()) {
    cout << endl << "FileParser: File " << fname << " does not exist!!!" << endl;
    exit(EXIT_FAILURE);
  } else { // success!
    cout << "FileParser: Successfully opened file " << fname << endl;
  }

  StoreProcessName(fname);

  // create a new TObjArray to store all histograms
  m_hists = new TObjArray();

  return;
}


void FileParser::OpenThetaFile(TString cyclename)
{
  // open the root files with names given in the TObjArray

  if (m_file != NULL){
    cerr << "FileParser::OpenFile: Can not open new file, since file " 
	 << m_file->GetName() << " is still in memory. Abort." << endl;
    exit(EXIT_FAILURE);
  }
  TString fname = "" ;
  if (cyclename.Sizeof()!=0){
    fname = cyclename;
      }  

  // fname = target;


  if (debug) cout << "Opening file with name " << fname << "..." << endl;
  m_file = new TFile(fname, "READ");
  if (debug){
    cout << "... success! pointer = " << m_file << endl;
    cout << "name = " << m_file << endl;
    cout << " is open? " << m_file->IsOpen() << endl;
    m_file->ls();
  }
    
  if (!m_file->IsOpen()) {
    cout << endl << "FileParser: File " << fname << " does not exist!!!" << endl;
    exit(EXIT_FAILURE);
  } else { // success!
    cout << "FileParser: Successfully opened file " << fname << endl;
  }

  StoreProcessName(fname);

  // create a new TObjArray to store all histograms
  m_hists = new TObjArray();
  m_shapeSys = new TObjArray();
  return;
}


void FileParser::StoreProcessName(TString name)
{
  
  TObjArray* pieces = name.Tokenize(".");
  for (int i=0; i<pieces->GetEntries(); ++i){
    TString piece = ((TObjString*)pieces->At(i))->GetString();
    if (piece.CompareTo("root")==0){
      m_process = ((TObjString*)pieces->At(i-1))->GetString();
      if (debug) cout << "Process in file = " << m_process << endl;
    }
  }
}


TObjArray* FileParser::FindSubdirs()
{
  // find all subdirectories (former histogram collections) in the open file
  // returns a TObjArray with the names of the subdirectories 

  m_file->cd();
  TObjArray* dirnames = new TObjArray();
  TString dirname(""); // empty directory, to stay in home dir first
  dirnames->Add(new TObjString(dirname));

  TKey *key;
  TIter nextkey( gDirectory->GetListOfKeys() );
  while ( (key = (TKey*)nextkey())) {
    TObject *obj = key->ReadObj();
    if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {    // found a subdirectory! 
      TString dirname(((TDirectory*) obj)->GetName());
      dirnames->Add(new TObjString(dirname));
      if (debug) cout << "Found directory " << dirname << endl;
    }
  }
  return dirnames;

}

void FileParser::BrowseFile()
{

  if (!m_file){
    cerr << "FileParser::BrowseFile: No file open. Abort." << endl;
    exit(EXIT_FAILURE);
  }
  TObjArray* dirs = FindSubdirs();

  // loop over all directories and get the histograms
  for (Int_t i=0; i<dirs->GetEntries(); ++i){

    TString dirname = ((TObjString*)dirs->At(i))->GetString();
    if (debug) cout << "Getting all histograms from directory " << dirname << endl;

    m_file->cd();
    gDirectory->Cd(dirname);

    // loop over all histograms in the directory and store them
    TKey *key;
    TIter nextkey( gDirectory->GetListOfKeys() );
    while ( (key = (TKey*)nextkey())) {

      TObject *obj = key->ReadObj();

      if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {

	// histogram found
	TH1* thist = (TH1*) obj;

	TString name = thist->GetName();
	//if (name.Contains("chi2hi")) continue;
	//if (!name.Contains("top")) continue;

	if (m_do_cumulative) MakeCumulativeHist(thist);
	TH1* rebinned = Rebin(thist, dirname);
	SHist* shist = NULL;
	if (rebinned){
	  shist = new SHist(rebinned);
	} else {
	  shist = new SHist(thist);
	}
	shist->SetProcessName(m_process);
	if (dirname==""){
	  shist->SetDir("Main");
	} else {
	  shist->SetDir(dirname);
	}
	if (debug) cout << "Adding hist " << shist->GetHist()->GetName() 
			<< " (process = " << m_process << ")" << endl;
	m_hists->Add(shist);	
      }
      
      delete obj;

    }

  }

  return;

}


void FileParser::BrowseThetaFile(TString sample)
{

  if (!m_file){
    cerr << "FileParser::BrowseFile: No file open. Abort." << endl;
    exit(EXIT_FAILURE);
  }
   
  m_file->cd();
  TKey *key;
  TIter nextkey( m_file->GetListOfKeys() );

  while ( (key = (TKey*)nextkey())) {
    
    TString histName = key->GetName();

    histName.ReplaceAll("__", "#");
    TObjArray* pieces = histName.Tokenize("#");

    if (((TObjString*)pieces->At(1))->GetString() == sample){

      TObject *obj = key->ReadObj();
      
      if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {

	// histogram found
	TH1* thist = (TH1*) obj;

  thist->SetName(histName);
  TString name = thist->GetName();

  if (name.Contains("chi2hi")) continue;
  if (name.Contains("mll")) continue;

  // debugging
  //if (name.Contains("ele_t0b0", TString::kIgnoreCase)) continue; 
  //if (name.Contains("ele_t0b1", TString::kIgnoreCase)) continue; 
  //if (name.Contains("ele_t1b0", TString::kIgnoreCase)) continue; 

  //if (name.Contains("muo_t0b0", TString::kIgnoreCase)) continue;
  //if (name.Contains("muo_t0b1", TString::kIgnoreCase)) continue;
  //if (name.Contains("muo_t1b0", TString::kIgnoreCase)) continue;

  // ele HLT and IDs nuisance parameters
  if (name.Contains("ele", TString::kIgnoreCase) && (!name.Contains("DATA"))){
    thist->Scale(0.94);
  } 
  // ttag light nuisance parameters
  if (name.Contains("T1", TString::kIgnoreCase) && (!name.Contains("DATA")) 
      && (!name.Contains("ttbar", TString::kIgnoreCase))) {
    thist->Scale(0.85);
  } 

	if (m_do_cumulative) MakeCumulativeHist(thist);
	TH1* rebinned = Rebin(thist, "");


	//TH1* rebinned = NULL;
	SHist* shist = NULL;
	if (rebinned){
	  shist = new SHist(rebinned);
	}
	else {
	  shist = new SHist(thist);
	}
	TString proc_name = ((TObjString*)pieces->At(1))->GetString();
	shist->SetProcessName(proc_name);
	SetProcessName(proc_name);
	TString hname = ((TObjString*)pieces->At(0))->GetString();
	shist->SetName(hname);
	  
	shist->SetDir("Main");
	
	if (pieces->GetEntries()>2){
	  m_shapeSys->Add(shist);
	  proc_name = ((TObjString*)pieces->At(1))->GetString() + "__" + ((TObjString*)pieces->At(2))->GetString() + "__" + ((TObjString*)pieces->At(3))->GetString();
	  shist->SetProcessName(proc_name);	  
	  SetProcessName(proc_name);
	  if (debug) cout << "Adding hist to systematic sample: " << shist->GetHist()->GetName() 
	  		  << " (process = " << m_process << ")" << endl;
	} else {
	  m_hists->Add(shist);
	  if (debug) cout << "Adding hist " << shist->GetHist()->GetName() 
			  << " (process = " << m_process << ")" << endl;
	}
      }
      
      delete obj;
    }
  }

  return;

}

void FileParser::MakeCumulativeHist(TH1* hist)
{
  for (Int_t i=1; i<hist->GetNbinsX()+1; ++i){
    Double_t sum = 0;
    Double_t sumw2 = 0;
    for (int j=i; j<hist->GetNbinsX()+1; ++j){
      sum += hist->GetBinContent(j);
      sumw2 += hist->GetSumw2()->At(j);
    }
    hist->SetBinContent(i, sum);
    hist->SetBinError(i, TMath::Sqrt(sumw2));
  }

}


TH1* FileParser::Rebin(TH1* hist, TString dirname)
{						

  TString name(hist->GetName());
  TString title(hist->GetTitle());

  if (name.Contains("mtt")) {
    //  TH1* rebinned = hist->Rebin(2);
    TH1* rebinned = hist->Rebin(2);
    rebinned->GetXaxis()->SetRangeUser(0,3500);
    rebinned->SetTitle("M_{t#bart} [GeV]");

    if (name.BeginsWith("ele_T1B0_L1chi2lo_mtt")){
      //cout << "e+jets, 1top, Name = " << name << endl;
      //rebinned->SetBinContent(32,0);
      //rebinned->SetBinContent(33,0);
      //rebinned->SetBinContent(34,0);
      //rebinned->SetBinContent(35,0);
    }
    return rebinned;
  }


  if (name == "Pt_toplep_rec") {
    TH1* rebinned = hist->Rebin(2);
    rebinned->GetXaxis()->SetRangeUser(0,1000);
    rebinned->SetTitle("P_{T, top}^{lept} [GeV]");
    return rebinned;
  }

  if (name == "Pt_tophad_rec") {
    TH1* rebinned = hist->Rebin(2);
    rebinned->GetXaxis()->SetRangeUser(0,1000);
    rebinned->SetTitle("P_{T, top}^{had} [GeV]");
    return rebinned;
  }


  if (name.CompareTo("toptags")==0){// && dirname.Contains("cutflow6") && title.Contains("electron")){
   
    Double_t binsx[] = {0, 960, 1020, 1080, 1140, 1200, 1260, 1320, 1380, 1440, 1500, 1560, 1620, 1680, 1740, 1800, 1860, 1920, 1980, 2040, 2100, 2400, 3000};
    name.Append("_rebin_lx");
    TH1* rebinned = hist->Rebin(22, name, binsx);
    rebinned->SetTitle("HT [GeV]");
    return rebinned;

  } else if (name.BeginsWith("mu_0top0btag_mttbar")) {
    
    TH1* rebinned = hist->Rebin(2);
    if (name=="mu_0top0btag_mttbar__wlight"){
      rebinned->SetBinContent(35,0.1589);
      rebinned->SetBinError(35,0.3867);
    }
    rebinned->GetXaxis()->SetRangeUser(0,3500);
    rebinned->SetTitle("M_{t#bart} [GeV]");
    return rebinned;

  } else if (name.BeginsWith("mu_0top1btag_mttbar")) {
    
    TH1* rebinned = hist->Rebin(2);

    rebinned->SetBinContent(32,0);
    rebinned->SetBinContent(33,0);
    rebinned->SetBinContent(34,0);
    rebinned->SetBinContent(35,0);
    rebinned->SetBinContent(37,0);
    rebinned->SetBinContent(40,0);
    rebinned->SetBinContent(43,0);

    rebinned->GetXaxis()->SetRangeUser(0,3500);
    rebinned->SetTitle("M_{t#bart} [GeV]");
    return rebinned;

  } else if (name.BeginsWith("mu_1top_mttbar")) {
    
    //Double_t binsx[] = {0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3500, 5000};
    //Double_t binsx[] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 
    //1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 
    //3000, 3500, 5000};
    //name.Append("_rebin_lx");
    //TH1* rebinned = hist->Rebin(30, name, binsx);
    TH1* rebinned = hist->Rebin(2);

    rebinned->SetBinContent(29,0);
    rebinned->SetBinContent(30,0);
    rebinned->SetBinContent(31,0);
    rebinned->SetBinContent(32,0);
    rebinned->SetBinContent(33,0);
   
    rebinned->GetXaxis()->SetRangeUser(0,3500);
    rebinned->SetTitle("M_{t#bart} [GeV]");
    return rebinned;

  } else if (name.BeginsWith("el_0top0btag_mttbar")) {
    cout << "name = " << name << endl;
    TH1* rebinned = hist->Rebin(2);
    if (name=="el_0top0btag_mttbar__wlight")
    {
      rebinned->SetBinContent(32, 0.31089);
      rebinned->SetBinError(32, 0.642);
      rebinned->SetBinContent(33, 0.3073);
      rebinned->SetBinError(33, 0.9042);
    }
    rebinned->GetXaxis()->SetRangeUser(0,3500);
    rebinned->SetTitle("M_{t#bart} [GeV]");
    return rebinned;

  } else if (name.BeginsWith("el_0top1btag_mttbar")) {
    
    TH1* rebinned = hist->Rebin(2);
    if (name=="el_0top1btag_mttbar__wlight")
    {
      rebinned->SetBinContent(33, 0.1038);
      rebinned->SetBinError(33, 0.338);
      rebinned->SetBinContent(35, 0.095);
      rebinned->SetBinError(35, 0.2975);
    }
    rebinned->GetXaxis()->SetRangeUser(0,3500);
    rebinned->SetTitle("M_{t#bart} [GeV]");
    return rebinned;

  } else if (name.BeginsWith("el_1top_mttbar")) {

    //Double_t binsx[] = {0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3500, 5000};
    //Double_t binsx[] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 
    //		        1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 
    //			3000, 3100, 3500, 5000};
    //name.Append("_rebin_lx");
    //TH1* rebinned = hist->Rebin(30, name, binsx);   
    TH1* rebinned = hist->Rebin(2);

    cout << "rebinned e+jets mttbar hist! Nbins = " << rebinned->GetNbinsX() << endl;
    rebinned->SetBinContent(30, 0);
    rebinned->SetBinContent(31, 0);
    rebinned->SetBinContent(32, 0);
    rebinned->SetBinContent(34, 0);

    rebinned->GetXaxis()->SetRangeUser(0,3500);
    rebinned->SetTitle("M_{t#bart} [GeV]");
    return rebinned;

  } else if (name.BeginsWith("lepton_0top0btag_mttbar")) {
    
    TH1* rebinned = hist->Rebin(2);
    rebinned->GetXaxis()->SetRangeUser(0,3500);
    rebinned->SetTitle("M_{t#bart} [GeV]");
    return rebinned;

  } else if (name.BeginsWith("lepton_0top1btag_mttbar")) {
    
    TH1* rebinned = hist->Rebin(2);
    rebinned->GetXaxis()->SetRangeUser(0,3500);
    rebinned->SetTitle("M_{t#bart} [GeV]");
    return rebinned;

  } else if (name.BeginsWith("lepton_1top_mttbar")) {
    
    //Double_t binsx[] = {0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 5000};
    //name.Append("_rebin_lx");
    //TH1* rebinned = hist->Rebin(24, name, binsx);
    TH1* rebinned = hist->Rebin(2);
    rebinned->GetXaxis()->SetRangeUser(0,3500);
    rebinned->SetTitle("M_{t#bart} [GeV]");
    return rebinned;

  } else if (name.BeginsWith("btag")) {

    TH1* rebinned = hist->Rebin(1);

    if (name.BeginsWith("btag2")){
      rebinned->GetXaxis()->SetRangeUser(500,2500);
    } else if (name.BeginsWith("btag3") || name.BeginsWith("btag4")){
      rebinned->GetXaxis()->SetRangeUser(500,4500);
    } else {
      rebinned->GetXaxis()->SetRangeUser(500,3500);      
    }
    rebinned->SetTitle("M_{t#bart} [GeV]");

    for (int i=1; i<rebinned->GetNbinsX(); ++i){

      if (rebinned->GetBinContent(i)<0.){
        rebinned->SetBinContent(i, fabs(rebinned->GetBinContent(i-1)-0.05));
        rebinned->SetBinError(i, rebinned->GetBinError(i-1)+0.05);        
      }      
      if (rebinned->GetBinContent(i)<0.1){
        rebinned->SetBinContent(i, 0.);
        rebinned->SetBinError(i, 0.);        
      }
    }    
    if (name.Contains("btag3") && name.Contains("qcd")){
      rebinned->SetBinContent(11, 0.);
      rebinned->SetBinError(11, 0.);
    }   
    if (name.Contains("btag4") && name.Contains("qcd")){
      rebinned->SetBinContent(43, 0.101);
      rebinned->SetBinError(43, 0.101);
      rebinned->SetBinContent(12, 0.);
      rebinned->SetBinError(12, 0.);      
    }
    if (name.Contains("btag5") && name.Contains("qcd")){
      rebinned->SetBinContent(28, 0.1205);
      rebinned->SetBinError(28, 0.1044);
    }   

    return rebinned;     

    return rebinned;

/*
  } else if (name.BeginsWith("btag")) {

    //Double_t binsx[] = {0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3500, 5000};
    name.Append("_rebin_lx");
    //TH1* rebinned = hist->Rebin(19, name, binsx);   
    //TH1* rebinned = hist->Rebin(1);
    TH1* rebinned = new TH1D(name, hist->GetTitle(), 50, 0, 5000);
    for (int i=1; i<51; ++i){
      if ((i<20 && hist->GetBinContent(i)<0.4) || ((i>25 && hist->GetBinContent(i)<0.15))) continue;
      rebinned->SetBinContent(i, hist->GetBinContent(i));
      rebinned->SetBinError(i, hist->GetBinError(i));
    }

    // do some filling up of empty bins to display errors correctly (at all)
    if (name.BeginsWith("btag1__")){
      rebinned->SetBinContent(8, 0);
      rebinned->SetBinError(8, 0);
    }
    cout << "name = " << name << endl;
    if (name == "btag2__ttbar_rebin_lx"){
      rebinned->SetBinContent(20, 0.30363);
    }
    if (name == "btag1__qcd_rebin_lx"){
      rebinned->SetBinContent(26, 0.31576);
      rebinned->SetBinError(26, 0.4583);
    }
    if (name == "btag4__qcd_rebin_lx"){
      rebinned->SetBinContent(34, 0.3128);
      rebinned->SetBinError(34, 0.4379);
    }
    if (name == "btag5__qcd_rebin_lx"){
      rebinned->SetBinContent(36, 0);
      rebinned->SetBinError(36, 0);
      rebinned->SetBinContent(26, 0.3113);
      rebinned->SetBinError(26, 0.3620);
      rebinned->SetBinContent(27, 0.381);
      rebinned->SetBinError(27, 0.3823);
    }

    rebinned->GetXaxis()->SetRangeUser(0,3500);
    rebinned->SetTitle("M_{t#bart} [GeV]");
    return rebinned;
*/

  } else if (name.Contains("htt")) {

    //Double_t binsx[] = {0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3500, 5000};
    name.Append("_rebin_lx");
    //TH1* rebinned = hist->Rebin(19, name, binsx);   
    TH1* rebinned = hist->Rebin(2);
    for (int i=1; i<rebinned->GetNbinsX()+1; ++i){
      if ((i<10 && hist->GetBinContent(i)<0.1) || ((i>25 && hist->GetBinContent(i)<0.07))){
        rebinned->SetBinContent(i, 0);
        rebinned->SetBinError(i, 0);
      }
    }

    // some corrections by hand
    if (name.Contains("httbtag2")){
      rebinned->SetBinContent(27, 0.);
      rebinned->SetBinError(27, 0.);
      rebinned->SetBinContent(31, 0.);
      rebinned->SetBinError(31, 0.);
    }
    if (name=="mjhttbtag0__qcd_rebin_lx"){
      rebinned->SetBinContent(33, 0.126);
      rebinned->SetBinError(33, 0.2591);
    }
    if (name.Contains("mjhttbtag1")){
      rebinned->SetBinContent(31, 0.);
      rebinned->SetBinError(31, 0.);
    }
    if (name.Contains("mjhttbtag2")){
      rebinned->SetBinContent(29, 0.);
      rebinned->SetBinError(29, 0.);
    }

    rebinned->GetXaxis()->SetRangeUser(0,3500);
    rebinned->SetTitle("M_{t#bart} [GeV]");
    return rebinned;

  } else if (name.Contains("mumu") || name.Contains("ee") || name.Contains("emu")) {

    double binsx[]={0,100,200,300,400,500,600,700,800,900,1000,
		    1100,1200,1300,
		    1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600, 2700, 2800, 2900, 3000,3500}; 
    name.Append("_rebin_lx");
    TH1* rebinned = hist->Rebin(30, name, binsx);  
    //TH1* rebinned = hist->Rebin(3);

    if (name.Contains("ee")){
      rebinned->SetBinContent(30, 0);
      rebinned->SetBinError(30, 0);
    }
    if (name.Contains("mumu")){
      rebinned->SetBinContent(3, 0);
      rebinned->SetBinError(3, 0);
      rebinned->SetBinContent(29, 0);
      rebinned->SetBinError(29, 0);
    }

    rebinned->GetXaxis()->SetRangeUser(0,3000);
    rebinned->SetTitle("M_{t#bart} [GeV]");
    return rebinned;

  } else if (name.Contains("dilepton")) {

    double binsx[]={0,100,200,300,400,500,600,700,800,900,1000,
		    1100,1200,1300,
		    1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600, 2700, 2800, 2900, 3000,3500}; 
    name.Append("_rebin_lx");
    TH1* rebinned = hist->Rebin(30, name, binsx);   
    //TH1* rebinned = hist->Rebin(3);
    rebinned->GetXaxis()->SetRangeUser(0,3000);
    rebinned->SetTitle("M_{t#bart} [GeV]");
    return rebinned;



  } else {
    return NULL;
  }

}

void FileParser::SetInfo(TString legname, double weight, int colour, int marker, float unc)
{
  
  for (int i=0; i<m_hists->GetEntries(); ++i){
    SHist* sh = (SHist*)m_hists->At(i);
    sh->SetLegName(legname);
    sh->SetWeight(weight);
    sh->SetUnc(unc);
    if (weight>0) sh->GetHist()->Scale(weight);
    sh->GetHist()->SetMarkerColor(colour);
    sh->GetHist()->SetLineColor(colour);

    if (marker > 1 ){
      sh->SetDrawMarker(true);
      sh->GetHist()->SetMarkerStyle(marker);
    } else {
      sh->SetDrawMarker(false);
      sh->GetHist()->SetMarkerStyle(0);
      sh->GetHist()->SetMarkerSize(0);
      sh->GetHist()->SetLineWidth(2);   
    }

    // histogram is transparent if marker < 0  
    if (marker < 0 ){
      // change line style
      if (marker==-1) sh->GetHist()->SetLineStyle(kDashed);
      if (marker==-2) sh->GetHist()->SetLineStyle(kDotted);
      if (marker==-3) sh->GetHist()->SetLineStyle(kDashDotted);
      if (marker==-4) sh->GetHist()->SetLineStyle(kDashDotted);    
    }
  }
}
