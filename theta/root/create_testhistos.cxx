#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"

int main(){
   TFile file("testhistos.root", "recreate");
   TH1D * histo1d = new TH1D("histo1d", "histo1d", 24, -4, 20);
   for(int i=0; i<=25; ++i){
      histo1d->SetBinContent(i, i + 12.0); // content is lower bin edge + 17.
   }
   TH2D* histo2d = new TH2D("histo2d", "histo2d", 10, 0, 10, 11, 0, 11);
   for(int i=0; i<=11; ++i){
      for(int j=0; j<=12; ++j){
           histo2d->SetBinContent(i, j, (i + 0.78) * (j + 3.02));
      }
   }
   TH3D * histo3d = new TH3D("histo3d", "histo3d", 10, 0, 10, 11, 0, 11, 12, 0, 12);
   for(int i=0; i<=11; ++i){
      for(int j=0; j<=12; ++j){
         for(int k=0; k<=13; ++k){
           histo3d->SetBinContent(i, j, k, (i + 0.12) * (j + 1.34) * (k + 5.67));
         }
      }
   }
   file.Write();
   file.Close();
}


