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
const TString backName[] = {"hBack1","hBack2","hBack3"};
const double massMax = 3000;
const double massMin = 15;
const int nEvents = 1e7;
const bool backInc = true;//include background
const bool exactClosure = true;//set exact closure
const bool effInc = false; //include efficiency

void backToyModel()
{
  TH1::SetDefaultSumw2();
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  ofstream parameterFile;
  parameterFile.open("parameters.txt");
  parameterFile << exactClosure << " " << effInc << " " << backInc << endl;
  parameterFile.close();
  TFile fToyData(toyModelName,"recreate");
  TFile*file = new TFile(mcDist);
  TFile*fileSF = new TFile(mcSF);
  TProfile*profileSF = (TProfile*)fileSF->Get("hSFvsInvMassAll_pfx");
  
  TF1*fResolutionModel = new TF1("fPhiResolutionModel","gaus(0)",-20,20);
  fResolutionModel->SetParameters(1,0,3);
  TEfficiency*efficiency = (TEfficiency*)file->Get("Efficiency");

  TH1D*hMassDist = (TH1D*)file->Get("hGenInvMass");
  hMassDist->SetName("hMassDist");
  double backDistInt,massDistInt;
  double massTrue,massReco,massTrueSig,massRecoSig,massBack;
  double weight,backWeight,signalWeight;
  TH1D*hBackgroundDist;
  TH1D*hMassWithBackground;
  if(backInc){
    TH1D*hFakes = (TH1D*)fileSF->Get("hFakesInvMass");
    TH1D*hEW = (TH1D*)fileSF->Get("hEWInvMass");
    TH1D*hTops = (TH1D*)fileSF->Get("hTopsInvMass");
    for(int i=0;i<87;i++){
      if(hFakes->GetBinContent(i)<0.0) hFakes->SetBinContent(i,0.0);
      if(hEW->GetBinContent(i)<0.0) hEW->SetBinContent(i,0.0);
      if(hTops->GetBinContent(i)<0.0) hTops->SetBinContent(i,0.0);
    }
    hBackgroundDist = (TH1D*)hFakes->Clone("hBackgroundDist");
     hBackgroundDist->Add(hEW);
     hBackgroundDist->Add(hTops);
     hBackgroundDist->SetMarkerStyle(20);
     hBackgroundDist->SetMarkerColor(kBlack);
     hBackgroundDist->SetLineColor(kBlack);
    hMassWithBackground = (TH1D*)hMassDist->Clone();
     hMassWithBackground->Add(hBackgroundDist);
     hMassWithBackground->SetMarkerStyle(20);
     hMassWithBackground->SetMarkerColor(kBlack);
     hMassWithBackground->SetLineColor(kBlack);
  weight = hMassWithBackground->Integral()/nEvents; 
  backWeight = weight*hBackgroundDist->Integral()/hMassWithBackground->Integral();
  signalWeight = hMassDist->Integral()/nEvents;
  }//end if backInc

 TH1D*hReco = new TH1D("hReco","",nLogBins2,massbins2);
 TH1D*hTrue = new TH1D("hTrue","",nLogBins,massbins);
  hTrue->SetMarkerStyle(20);
  hTrue->SetMarkerColor(kBlack);
 TH1D*hBack = new TH1D("hBack","",nLogBins2,massbins2);
  hBack->SetFillColor(kOrange+1);
  hBack->SetLineColor(kOrange+1);
 TH2D*hMatrix = new TH2D("hMatrix","",nLogBins,massbins,nLogBins2,massbins2);
 TH2D*hMatrixSig = new TH2D("hMatrixSig","",nLogBins,massbins,nLogBins2,massbins2);
 TH1D*hRecoSig = new TH1D("hRecoSig","",nLogBins2,massbins2);
 TH1D*hTrueSig = new TH1D("hTrueSig","",nLogBins,massbins);
  hTrueSig->SetFillColor(kRed+2);
  hTrueSig->SetLineColor(kRed+2);

 TRandom3*random = new TRandom3();
 Long64_t N=0;
 double smear;
 for(int i=0;i<nEvents;i++){
   counter(N,nEvents);
   N++;
   smear = fResolutionModel->GetRandom();
   massTrue = hMassWithBackground->GetRandom();//signal + background
   massReco = massTrue+smear;   
   massBack = hBackgroundDist->GetRandom();//Background only
   massTrueSig = hMassDist->GetRandom();//Signal only
   massRecoSig = massTrueSig+smear;

   hBack->Fill(massBack,backWeight);
   hTrue->Fill(massTrue,weight);
   hReco->Fill(massReco,weight);
   hTrueSig->Fill(massTrueSig,signalWeight);
   hRecoSig->Fill(massRecoSig,signalWeight);
   hMatrix->Fill(massTrue,massReco,weight);
 }

  hBack->Rebin(2);
  THStack*hStack = new THStack("hStack","");
   hStack->Add(hBack);
   hStack->Add(hTrueSig);  
  TCanvas*canvas = new TCanvas("canvas","",10,10,1000,1400);
   canvas->SetLogy();
   canvas->SetLogx();
   canvas->SetGrid();
  hStack->Draw("hist");
  hBackgroundDist->Rebin(2);
  hMassWithBackground->Rebin(2);
  hMassWithBackground->Draw("PE,same"); 
  hBackgroundDist->Draw("PE,same");

  TFile*file2 = new TFile(histSaveName,"recreate");
  file2->cd();
  hBackgroundDist->Write();
  hMassDist->Write();
  hBack->Write();
  hTrue->Write();
  hReco->Write();
  hMatrix->Write();
  file2->Write();
  file2->Close();

}


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

