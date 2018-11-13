#include "TRandom.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TTimeStamp.h"

void counter(Long64_t i, Long64_t N);
const double pi = TMath::Pi();
const TString mcDist = "/home/hep/wrtabb/git/DY-Analysis/plots/plotsDY.root";
const TString toyModelName = "toyData.root";
const int nEvents = 1e7;
enum ModelType{
  IN_MASS_RANGE,//events not allowed to migrate into/out of mass range (15,3000)
  OUT_MASS_RANGE//events are allowed to migrate into/out of mass range (15,3000)
};

void toyModel()
{
  //Data to store
  float massTrue,massMeasured,smear;
  int modelType = IN_MASS_RANGE;//define if bin migration can move into/out of mass range
  
  //Setting up root file and tree for storing data
  TFile fToyData(toyModelName,"recreate");
  TTree*tree = new TTree("toyData","");
  tree->Branch("massTrue",&massTrue,"massTrue/f");
  tree->Branch("massMeasured",&massMeasured,"massMeasured/f");
  
  TFile*file = new TFile(mcDist);
  //Defining mass distribution model
  //and resolution model
  TF1*fToyModel = new TF1("fToyModel","([0]/x)+gaus(1)",10,200);
  fToyModel->SetParameters(10,1,90,7);
  TF1*fResolutionModel = new TF1("fPhiResolutionModel","gaus(0)",-20,20);
  fResolutionModel->SetParameters(1,0,3);
  TH1D*hMassDist = (TH1D*)file->Get("hGenAllInvMass"); 
  //Filling variables
  Long64_t nentries = 0;
  for(Long64_t i=0;i<nEvents;i++){
    counter(i,nEvents);
    massTrue = fToyModel->GetRandom();
    smear = fResolutionModel->GetRandom();
    massMeasured = massTrue+smear;
    if(modelType == IN_MASS_RANGE){
      if(massMeasured<200&&massMeasured>10){
        tree->Fill();
        nentries++;
      }
    }//end modelType IN_MASS_RANGE
    else{
      tree->Fill();
      nentries++;
    }//end modelType OUT_MASS_RANGE
  }//end for loop
  fToyData.Close();
  cout << "Number of events processed: " << nentries << endl;
}//end toyModel()

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


