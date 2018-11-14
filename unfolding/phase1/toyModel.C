#include "TRandom.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TTimeStamp.h"

const int nLogBins = 43;
const int nLogBins2 = 2*nLogBins;
const float  massbins[] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81,
  86, 91,96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,200, 220, 243, 273,
  320, 380, 440, 510, 600, 700, 830, 1000, 1500,3000};
const float massbins2[] = {15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,
                            57.5,60,62,64,66,68,70,72,74,76,78.5,81,83.5,86,88.5,91,93.5,96,
                            98.5,101,103.5,106,108,110,112.5,115,117.5,120,123,126,129.5,133,
                            137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,210,220,
                            231.5,243,258,273,296.5,320,350,380,410,440,475,510,555,600,650,
                            700,765,830,915,1000,1250,1500,2250,3000};

void counter(Long64_t i, Long64_t N);
const double pi = TMath::Pi();
const TString mcDist = "/home/hep/wrtabb/git/DY-Analysis/plots/plotsDY.root";
const TString toyModelName = "toyData.root";
const TString histSaveName = "toyUnfold.root";
const int nBins = 40;
const int nBins2 = 80;
const double binLow = 20;
const double binHigh = 200;
const int nEvents = 1e7;
enum ModelType{
  IN_MASS_RANGE,//events not allowed to migrate into/out of mass range 
  OUT_MASS_RANGE//events are allowed to migrate into/out of mass range 
};

void toyModel()
{
  TH1::SetDefaultSumw2();
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0); 
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

  TH1D*hTrue = new TH1D("hTrue","",nLogBins,massbins);
  TH1D*hReco = new TH1D("hReco","",nLogBins2,massbins2);
  TH2D*hMatrix = new TH2D("hMatrix","",nLogBins,massbins,nLogBins2,massbins2);
  //Filling variables
  Long64_t nentries = 0;
  for(Long64_t i=0;i<nEvents;i++){
    counter(i,nEvents);
    massTrue = hMassDist->GetRandom();
    smear = fResolutionModel->GetRandom();
    massMeasured = massTrue+smear;
    if(modelType == IN_MASS_RANGE){
      if(massMeasured<3000&&massMeasured>15){
        tree->Fill();
        nentries++;
        hTrue->Fill(massTrue);
        hReco->Fill(massMeasured);
        hMatrix->Fill(massTrue,massMeasured);
      }
    }//end modelType IN_MASS_RANGE
    else{
      tree->Fill();
      nentries++;
    }//end modelType OUT_MASS_RANGE
  }//end for loop
  fToyData.Close();
  TCanvas*canvas=new TCanvas("canvas","",10,10,1000,1000);
  canvas->SetLogy();
  canvas->SetLogx();
  hMatrix->GetXaxis()->SetMoreLogLabels();
  hMatrix->GetXaxis()->SetNoExponent();
  hMatrix->GetYaxis()->SetMoreLogLabels();
  hMatrix->GetYaxis()->SetNoExponent();
  hMatrix->GetXaxis()->SetTitle("true mass [GeV]");
  hMatrix->GetYaxis()->SetTitle("reco mass [GeV]");
  hMatrix->Draw("colz");
  canvas->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase1Plots/step1MigrationMatrix_RecoInMassRange.png");
  TFile*file2 = new TFile(histSaveName,"recreate");
  file2->cd();
  hMatrix->Write();
  hReco->Write();
  hTrue->Write();
  file2->Write();
  file2->Close();

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


