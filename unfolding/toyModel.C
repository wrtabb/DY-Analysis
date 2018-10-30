#include "TRandom.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TTimeStamp.h"

void counter(Long64_t i, Long64_t N);
const int nBins = 43;
const double massbins[] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 
                           91,96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 
                           185,200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 
                           1500,3000};

const double pi = TMath::Pi();
const TString toyModelName = "toyData.root";
const int nEvents = 10e6;

void toyModel()
{
  //Data to store
  float massTrue,massMeasured;

  gRandom->SetSeed(0);
  //Setting up root file and tree for storing data
  TFile fToyData(toyModelName,"recreate");
  TTree*tree = new TTree("toyData","Toy Data");
  tree->Branch("massTrue",&massTrue,"massTrue/f");
  tree->Branch("massMeasured",&massMeasured,"massMeasured/f");

  //Defining mass distribution model
  //and resolution model
  TF1*fMassModel = new TF1("fMassModel","([0]/x)+gaus(1)",15,3000);
  fMassModel->SetParameters(1000,100,90,7);
  TF1*fResolutionModel = new TF1("fPhiResolutionModel","gaus(0)",-2,2);
  fResolutionModel->SetParameters(1,0,3);
  
  //Filling variables
  Long64_t nentries = 0;
  for(Long64_t i=0;i<nEvents;i++){
    counter(i,nEvents);
    massTrue = fMassModel->GetRandom();
    massMeasured = massTrue+fResolutionModel->GetRandom();
    if(massMeasured<15 || massMeasured>3000) continue;
    nentries++;
    tree->Fill();
  }
  //fMassModel->Draw(); 
  tree->Write(); 
  fToyData.Close();
  cout << "Events processed: " << nentries << endl;
}

//Counter for tracking program progress
void counter(Long64_t i, Long64_t N)
{
  int P = 100*(i)/(N);
  TTimeStamp eventTimeStamp;
  if(i%(N/100)==0) {
    cout << "toyModel.C " << "[Time: " << eventTimeStamp.AsString("s") << "] " << P
    << "%" << endl;
  }
return;
}


