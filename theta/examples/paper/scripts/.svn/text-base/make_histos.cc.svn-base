#include "TFile.h"
#include "TH1.h"

#include <cmath>

using namespace std;

int main(){
    TFile f("results/templates.root", "recreate");
    TH1 * signal = new TH1D("signal", "signal", 100, 0, 500);
    TH1 * signal_plus = new TH1D("signal_plus", "signal", 100, 0, 500);
    TH1 * signal_minus = new TH1D("signal_minus", "signal", 100, 0, 500);
    TH1 * bkg = new TH1D("bkg", "bkg", 100, 0, 500);
    for(int i=1; i<=100; ++i){
        double x = signal->GetBinCenter(i);
        signal->SetBinContent(i, exp(-pow(x - 250, 2) / (2*50*50)));
        signal_plus->SetBinContent(i, exp(-pow(x - 250 * 1.1, 2) / (2*52*52)));
        signal_minus->SetBinContent(i, exp(-pow(x - 250 * 0.9, 2) / (2*48*48)));
        bkg->SetBinContent(i, exp(-0.0012*x));
    }
    f.Write();
}



