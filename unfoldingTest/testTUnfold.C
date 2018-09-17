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
const double massbins[44] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,
                             110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,
                             380,440,510,600,700,830,1000,1500,3000};
const int nBins = 43;
//File Names
const TString fileName= "unfoldingClosureTest.root";

void testTUnfold()
{
  //Load the files
  TFile*file= new TFile(fileName);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  //Define hisograms
  TH1F*hReco = (TH1F*)file->Get("migMatrixGENFSvsReco_py");
  TH1F*hGen = (TH1F*)file->Get("migMatrixGENFSvsReco_px");
  TH2F*hMatrix = (TH2F*)file->Get("migMatrixGENFSvsReco");
  hReco->GetXaxis()->SetNoExponent();
  hReco->GetXaxis()->SetMoreLogLabels();
  hReco->SetMarkerStyle(20);
  hReco->SetMarkerColor(kBlack);
  hReco->SetLineColor(kBlack);
  hGen->GetXaxis()->SetNoExponent();
  hGen->GetXaxis()->SetMoreLogLabels();
  hGen->SetFillColor(kRed+2);
  hGen->SetLineColor(kRed+2);
  hMatrix->GetXaxis()->SetNoExponent();
  hMatrix->GetXaxis()->SetMoreLogLabels();
  hMatrix->GetYaxis()->SetNoExponent();
  hMatrix->GetYaxis()->SetMoreLogLabels();

  cout << endl;
  cout << "Migration Matrix Information:" << endl;
  cout << "-----------------------------" << endl;
  cout << "nx Bins = " << hMatrix->GetNbinsX() << endl;
  cout << "ny Bins = " << hMatrix->GetNbinsY() << endl;
  cout << endl;
  cout << "Observed Matrix Information:" << endl;
  cout << "----------------------------" << endl;
  cout << "nbins = " << hReco->GetNbinsX() << endl;
  cout << endl;

  double binContent;
  int globalBin,globalBinFlipped;
  TH2F*hMatrixFlipped=new TH2F("hMatrixFlipped","",nBins,massbins,nBins,massbins);
  hMatrixFlipped->SetTitle("X and Y flipped");
  hMatrixFlipped->GetXaxis()->SetTitle("Reco");
  hMatrixFlipped->GetYaxis()->SetTitle("Gen");
  for(int i=0;i<nBins;i++){
    for(int j=0;j<nBins;j++){
    globalBin = hMatrix->GetBin(i+1,j+1);
    binContent = hMatrix->GetBinContent(globalBin);
    globalBinFlipped = hMatrix->GetBin(j+1,i+1);
    hMatrixFlipped->SetBinContent(globalBinFlipped,binContent);
    }
  }
  //Begin TUnfold

  TUnfoldDensity unfold(hMatrixFlipped,TUnfold::kHistMapOutputHoriz);
  Double_t tau=0.0;
  Double_t biasScale=0.0;
  unfold.DoUnfold(tau,hReco,biasScale);
  //The Unfolded Distribution
  TH1*hUnfolded = unfold.GetOutput("Unfolded");
  hUnfolded->SetMarkerStyle(22);
  hUnfolded->SetMarkerColor(kBlue+2);

  //Define canvases
  TCanvas*canvas1 = new TCanvas("canvas1","",10,10,1200,1200);
  canvas1->SetLogy();
  canvas1->SetLogx();

  TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
  legend->SetTextSize(0.02);
  legend->AddEntry(hGen,"True Distribution");
  legend->AddEntry(hReco,"Measured Distribution");
  legend->AddEntry(hUnfolded,"Unfolded Distribution");
  hGen->GetXaxis()->SetTitle("Invariant mass [GeV]");
  hGen->Draw("hist");
  hReco->Draw("PE,same");
  hUnfolded->Draw("PE,same");
  legend->Draw("same");

  canvas1->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/testUnfolding.jpg");
}
