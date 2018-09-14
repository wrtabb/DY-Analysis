#include "TUnfold.h"
#include "TUnfoldSys.h"
#include "TUnfoldDensity.h"
#include "TUnfoldBinning.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
using namespace std;

//File Names
const TString testSampleFileName = "unfoldMCSample.root";
const TString testMatrixFileName = "migrationMatrix.root";

void testTUnfold()
{
  //Load the files
  TFile*testSampleFile = new TFile(testSampleFileName);
  TFile*testMatrixFile = new TFile(testMatrixFileName);

  //Define hisograms
  TH1F*hObserved = (TH1F*)testSampleFile->Get("hMCInvMass0");
  hObserved->GetXaxis()->SetNoExponent();
  hObserved->GetXaxis()->SetMoreLogLabels();

  TH2F*hMigrationMatrix = (TH2F*)testMatrixFile->Get("migMatrixGENFSvsReco");
  double bkg = 1.0;//amount to place in under/overflow bins
  int overBin = 44;//n=43, so this is the overflow bin
  int underBin = 0;//underflow bin
  int underFlowLow = hMigrationMatrix->GetBin(underBin,underBin);//bins for the four corners
  int underFlowHigh = hMigrationMatrix->GetBin(underBin,overBin);//bins for the four corners
  int overFlowLow = hMigrationMatrix->GetBin(overBin,underBin);//bins for the four corners
  int overFlowHigh = hMigrationMatrix->GetBin(overBin,overBin);//bins for the four corners
  hMigrationMatrix->SetBinContent(underFlowLow,bkg);//setting bin content in over/underflow
  hMigrationMatrix->SetBinContent(underFlowHigh,bkg);//setting bin content in over/underflow
  hMigrationMatrix->SetBinContent(overFlowLow,bkg);//setting bin content in over/underflow
  hMigrationMatrix->SetBinContent(overFlowHigh,bkg);//setting bin content in over/underflow

  hMigrationMatrix->GetXaxis()->SetNoExponent();
  hMigrationMatrix->GetXaxis()->SetMoreLogLabels();
  
  cout << endl;
  cout << "Migration Matrix Information:" << endl;
  cout << "-----------------------------" << endl;
  cout << "nx Bins = " << hMigrationMatrix->GetNbinsX() << endl;
  cout <<  "ny Bins = " << hMigrationMatrix->GetNbinsY() << endl;
  cout << endl;
  cout << "Observed Matrix Information:" << endl;
  cout << "----------------------------" << endl;
  cout << "nbins = " << hObserved->GetNbinsX() << endl;
  cout << endl;

  //Begin TUnfold
  TUnfold unfold(hMigrationMatrix,TUnfold::kHistMapOutputVert);
  Double_t tau=1.E-4;
  Double_t biasScale=0.0;
  unfold.DoUnfold(tau,hObserved,biasScale);
  //The Unfolded Distribution
  TH1*hUnfolded = unfold.GetOutput("Unfolded");
  //A distribution showing the overall correlations
  //TH2*hCorrelations = unfold.GetRhoIJ("Correlations");

  //Define canvases
  //TCanvas*canvas = new TCanvas("canvas","",10,10,1200,1200);
  //canvas->SetLogy();
  //canvas->SetLogx();
 
  //Draw Plots
  //hObserved->Draw("hist");
  //hUnfolded->Draw("PE,same");
  //hMigrationMatrix->Draw("colz");
}
