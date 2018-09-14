#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TEfficiency.h"

//File Locations
const TString dataFileName = "./plots/unfoldingMatrices.root";
const TString effFileName = "./plots/plotsDY.root";

//Parameters
Long64_t luminosity = 35867;//luminosity of data
const int nMassBins = 43;
const double massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86,
                             91,96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171,
                             185,200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000,
                             1500,3000};

void xSecCalc()
{
  gStyle->SetOptStat(0);
  //Opening files
  TFile*fData = new TFile(dataFileName);
  TFile*fEff = new TFile(effFileName);
  
  //Defining Histograms From Files
  TH1F*hData = (TH1F*)fData->Get("hMassUnfolded");//Unfolded data from unfolded.C
  TEfficiency*hEff = (TEfficiency*)fEff->Get("Efficiency");//Total efficiency
  TEfficiency*hAcceptance = (TEfficiency*)fEff->Get("Acceptance");//Acceptance

  //Cross Section Histogram
  TH1F*hXsec = new TH1F("hXsec","",nMassBins,massbins);
  hXsec->Sumw2();
  hXsec->SetTitle("Drell-Yan Cross Section");  
  hXsec->GetXaxis()->SetTitle("m [GeV]");
  hXsec->GetYaxis()->SetTitle("d#sigma/dm [pb/GeV]");
  hXsec->SetLineColor(kBlack);
  hXsec->SetMarkerColor(kBlack);
  hXsec->SetMarkerStyle(20);
  hXsec->GetXaxis()->SetMoreLogLabels();
  hXsec->GetXaxis()->SetNoExponent();
  hXsec->SetMinimum(0.0000001);
  hXsec->SetMaximum(1000);
  double acc, eff, sig, xSec, binWidth;
  for(int i=1;i<=nMassBins;i++){
    binWidth = hData->GetBinWidth(i);
    acc = hAcceptance->GetEfficiency(i);
    eff = hEff->GetEfficiency(i);
    sig = hData->GetBinContent(i);
    xSec = sig/(acc*eff*luminosity*binWidth);
    hXsec->SetBinContent(i,xSec);    
  }
  TCanvas*canvas = new TCanvas("canvas","",10,10,1200,1200);
  canvas->SetLogx();
  canvas->SetLogy();
  canvas->SetGrid();
  hXsec->Draw("PE");
  canvas->SaveAs("./plots/xSecInvMass.png");
}
