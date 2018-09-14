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

void testTUnfoldDensity()
{
  //Load the files
  TFile*testSampleFile = new TFile(testSampleFileName);
  TFile*testMatrixFile = new TFile(testMatrixFileName);

  //Define hisograms
  TH1F*hObserved = (TH1F*)testSampleFile->Get("hMCInvMass0");
  hObserved->GetXaxis()->SetNoExponent();
  hObserved->GetXaxis()->SetMoreLogLabels();

  TH2F*hMigrationMatrix = (TH2F*)testMatrixFile->Get("migMatrixGENFSvsReco");
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
  TUnfold::ERegMode regMode =
    TUnfold::kRegModeCurvature;
  TUnfold::EConstraint constraintMode=
    TUnfold::kEConstraintArea;
  TUnfoldDensity::EDensityMode densityFlags=
    TUnfoldDensity::kDensityModeBinWidth;
  TUnfoldDensity unfold(hMigrationMatrix,TUnfold::kHistMapOutputVert,
    regMode,constraintMode,densityFlags);
  Double_t tau=1.E-4;
  Double_t biasScale=0.0;
  unfold.DoUnfold(tau,hObserved,biasScale);
  
  //The Unfolded Distribution
  TH1*hUnfolded = unfold.GetOutput("Unfolded");
  //A distribution showing the overall correlations
 // TH2*hCorrelations = unfold.GetRhoIJ("Correlations");

  //Define canvases
  //TCanvas*canvas = new TCanvas("canvas","",10,10,1200,1200);
  //canvas->SetLogy();
  //canvas->SetLogx();
 
  //Draw Plots
  //hObserved->Draw("hist");
  //hUnfolded->Draw("hist,same");
  //hMigrationMatrix->Draw("colz");
}
