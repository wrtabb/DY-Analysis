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
const TString mcSF =  "/home/hep/wrtabb/git/DY-Analysis/plots/dataVsMC.root";
const TString toyModelName = "toyData.root";
const TString histSaveName = "toyUnfold.root";
const int nBins = 40;
const int nBins2 = 80;
const double binLow = 20;
const double binHigh = 200;
const double massMax = 3000;
const double massMin = 15;
const int nEvents = 1e7;

const bool exactClosure = true;//set exact closure
const bool effInc = false; //include efficiency

void toyModel()
{
  TH1::SetDefaultSumw2();
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0); 
  //gROOT->SetBatch(kTRUE);
  //Data to store
  float massTrue,massMeasured,smear;
  float massMeasuredData,massMeasuredMC; 
  
  ofstream parameterFile;
  parameterFile.open("parameters.txt");
  parameterFile << exactClosure << " " << effInc << endl;
  parameterFile.close();
  //Setting up root file and tree for storing data
  TFile fToyData(toyModelName,"recreate");
  TFile*file = new TFile(mcDist);

  TFile*fileSF = new TFile(mcSF);
  TProfile*profileSF = (TProfile*)fileSF->Get("hSFvsInvMassAll_pfx");

  //Defining mass distribution model
  //and resolution model
  TF1*fToyModel = new TF1("fToyModel","([0]/x)+gaus(1)",10,200);
  fToyModel->SetParameters(10,1,90,7);
  TF1*fResolutionModel = new TF1("fPhiResolutionModel","gaus(0)",-20,20);
  fResolutionModel->SetParameters(1,0,3);

  //Efficiency Plot
  TEfficiency*efficiency = (TEfficiency*)file->Get("Efficiency");
  //MC Mass Distribution for filling toy model
  TH1D*hMassDist = (TH1D*)file->Get("hGenInvMass"); 
  hMassDist->SetName("hMassDist");

  double effData,effMC,rand,rho;
  TRandom3*random = new TRandom3();

  const int nHists = 3;
  TH1D*hReco[nHists];
  TH1D*hTrue[nHists];
  TH2D*hMatrix[nHists];
  const TString recoName[] = {"hReco","hMCReco","hAltReco"};
  const TString trueName[] = {"hTrue","hMCTrue","hAltTrue"};
  const TString matrixName[] = {"hMatrixDontUse","hMatrix","hAltMatrix"};
  Long64_t N = 0;
  for(int j=0;j<nHists;j++){
    if(!exactClosure){
      gRandom->SetSeed(j+1);
      random->SetSeed(j+1);
    }
    else{
      gRandom->SetSeed(1);
      random->SetSeed(1);
    }
    hReco[j] = new TH1D(recoName[j],"",nLogBins2,massbins2);
    hTrue[j] = new TH1D(trueName[j],"",nLogBins,massbins);
    hMatrix[j] = new TH2D(matrixName[j],"",nLogBins,massbins,nLogBins2,massbins2);
    for(Long64_t i=0;i<nEvents;i++){
      counter(N,nHists*nEvents);
      N++;
      massTrue = hMassDist->GetRandom();
      smear = fResolutionModel->GetRandom();
      massMeasured = massTrue+smear;
      bool seenInData = kTRUE;
      bool seenInMC = kTRUE;

      if(effInc){
        rho = profileSF->GetBinContent(profileSF->FindBin(massTrue));
        effMC = efficiency->GetEfficiency(hMassDist->FindBin(massMeasured));
        effData = effMC*rho;
        rand = random->Rndm(); 
        if(rand>effMC) seenInMC = false;
        if(rand>effData) seenInData = false; 
      }
      else rho = 1.0;
      massMeasuredData = massMeasuredMC = massMeasured;
      if(!seenInData) massMeasuredData = 0;
      if(!seenInMC) massMeasuredMC = 0;
      hReco[j]->Fill(massMeasuredData);
      hTrue[j]->Fill(massTrue);
      if(seenInMC){
        hMatrix[j]->Fill(massTrue,massMeasuredMC, rho);
        hMatrix[j]->Fill(massTrue,0.0,1-rho);
      }
      else hMatrix[j]->Fill(massTrue,massMeasuredMC,1.0);
    }
  }//end j loop
  
  fToyData.Close();
  TString saveName;
  saveName = 
    "/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase1Plots/step1MigrationMatrix";
  if(effInc) saveName += "_EffInc";
  else saveName += "_NoEff";
  if(exactClosure) saveName += "_ClosureTest";

  TCanvas*canvas=new TCanvas("canvas","",10,10,1000,1000);
  canvas->SetLogy();
  canvas->SetLogx();
  hMatrix[1]->GetXaxis()->SetMoreLogLabels();
  hMatrix[1]->GetXaxis()->SetNoExponent();
  hMatrix[1]->GetYaxis()->SetMoreLogLabels();
  hMatrix[1]->GetYaxis()->SetNoExponent();
  hMatrix[1]->GetXaxis()->SetTitle("true mass [GeV]");
  hMatrix[1]->GetYaxis()->SetTitle("reco mass [GeV]");
  hMatrix[1]->GetYaxis()->SetTitleOffset(1.5);
  hMatrix[1]->Draw("colz");
  
  saveName += ".png";
  canvas->SaveAs(saveName);
  TFile*file2 = new TFile(histSaveName,"recreate");
  file2->cd();
  for(int i=0;i<nHists;i++){
    hMatrix[i]->Write();
    hReco[i]->Write();
    hTrue[i]->Write();
  }
  hMassDist->Write();
  file2->Write();
  file2->Close();

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


