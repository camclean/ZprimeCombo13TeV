#include "TFile.h"
#include "TTree.h"

int main(){
   TFile file("testtree.root", "recreate");
   TTree * tree = new TTree("tree", "tree");
   float o1, o2;
   //assume o1 and o2 might be observables of two analyses named 'analysis 1' and 'analysis 2'. Both
   // observables have range [0, 1].
   //Events 0--199 are only selected by analysis1, events 200--399 only by analsysis2, so the other
   // observable is set to -1 for these events to indicate that it has not been selected.
   tree->Branch("o1", &o1, "data/F");
   tree->Branch("o2", &o2, "data/F");
   for(size_t i=0; i<1000; ++i){
      o1 = o2 = 0.5;
      if(i<=200){
         o2 = -1;
      }
      else if(i<=400){
         o1 = -1;
      }
      tree->Fill();
   }
   file.Write();
   file.Close();
}

