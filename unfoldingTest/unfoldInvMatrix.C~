#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

// Constants

const TString dataFileName = "sampleMockData.root";
const TString simFileName = "sampleSimulation.root";
const TString mockDataTreeName = "mockData";
const TString simTreeName = "sim";

const int nMassBins = 15;
const int massMin = 20;
const int massMax = 200;

// Declarations
void setHistAttributes(TH1F *hist, Color_t color);

// Main function

void unfoldInvMatrix(){

  // Variables. Note that the variables are used both for two different
  // input trees. This is ok because this is not done concurrently.
  float massTrue, massMeasured;
  TBranch *b_massTrue;
  TBranch *b_massMeasured;

  // Set up files and trees.
  TFile *fMockData = new TFile(dataFileName);
  TTree *mockDataTree = (TTree*)fMockData->Get(mockDataTreeName);
  mockDataTree->SetBranchAddress("massTrue"    , &massTrue,     &b_massTrue);
  mockDataTree->SetBranchAddress("massMeasured", &massMeasured, &b_massMeasured);

  TFile *fSim = new TFile(simFileName);
  TTree *simTree = (TTree*)fSim->Get(simTreeName);
  simTree->SetBranchAddress("massTrue"    , &massTrue,     &b_massTrue);
  simTree->SetBranchAddress("massMeasured", &massMeasured, &b_massMeasured);

  // 
  // Loop over mock data and fill histograms
  //
  TH1F *hMassTrue     = new TH1F("hMassTrue"    ,"",nMassBins, massMin, massMax);
  TH1F *hMassMeasured = new TH1F("hMassMeasured","",nMassBins, massMin, massMax);
  for(int i=0; i < mockDataTree->GetEntries(); i++){
    mockDataTree->GetEntry(i);
    hMassTrue->Fill(massTrue);
    hMassMeasured->Fill(massMeasured);
  }
  // Also pack data into vectors
  TVectorD yieldsMeasured(nMassBins);
  TVectorD yieldsTrue(nMassBins);
  for(int i=1; i<=nMassBins; i++){
    yieldsMeasured(i-1) = hMassMeasured->GetBinContent(i);
    yieldsTrue    (i-1) = hMassTrue  ->GetBinContent(i);
  }

  // 
  // Prepare for unfolding, compose all the matrices
  //
  // Loop over simulation and fill the migration matrix histogram
  TH2F *hMigrationMatrix = new TH2F("hMigrationMatrix","migration matrix",
				    nMassBins, massMin, massMax,
				    nMassBins, massMin, massMax);
  for(int i=0; i < simTree->GetEntries(); i++){
    simTree->GetEntry(i);
    hMigrationMatrix->Fill(massTrue, massMeasured);
  }
  // Switch to working with matrices
  TMatrixD migrationMatrix(nMassBins,nMassBins);
  for(int i=1; i<=nMassBins; i++){
    for(int j=1; j<=nMassBins; j++){
      migrationMatrix(i-1,j-1) = hMigrationMatrix->GetBinContent(i,j);
    }
  }
  // Create the response matrix: normalize in columns (true mass) 
  TMatrixD responseMatrix(nMassBins,nMassBins);
  for(int i=0; i<nMassBins; i++){
    // Compute the sum of all elements for all observed masses for this true mass
    float sum=0;
    for(int j=0; j<nMassBins; j++){
      sum += migrationMatrix(i,j);
    }
    // Fill the response matrix
    for(int j=0; j<nMassBins; j++){
      if( sum != 0){
	responseMatrix(i,j) = migrationMatrix(i,j)/sum;
      }else{
	printf("can't properly normalize, sum=0\n");
	responseMatrix(i,j) = 0;
      }
    }
  }
  // "Unfolding" matrix: inverted response matrix
  TMatrixD unfoldingMatrix = responseMatrix;
  Double_t det;
  unfoldingMatrix.Invert(&det);
  // Repack the response and unfolding matrices into 2D histograms for plotting
  TH2F *hResponseMatrix = (TH2F*)hMigrationMatrix->Clone("hResponseMatrix");
  hResponseMatrix->SetTitle("response matrix");
  TH2F *hUnfoldingMatrix = (TH2F*)hMigrationMatrix->Clone("hUnfoldingMatrix");
  hUnfoldingMatrix->SetTitle("unfolding matrix");
  for(int i=0; i<nMassBins; i++){
    for(int j=0; j<nMassBins; j++){
      hResponseMatrix->SetBinContent(i+1, j+1, responseMatrix(i,j));
      hUnfoldingMatrix->SetBinContent(i+1, j+1, unfoldingMatrix(i,j));
    }
  }

  // 
  // Perform unfolding with the matrix inversion method
  //
  TVectorD yieldsUnfolded = (unfoldingMatrix.T()) * yieldsMeasured;
  // Save into histogram
  TH1F *hMassUnfolded = (TH1F*)hMassTrue->Clone("hMassUnfolded");
  hMassUnfolded->Reset();
  for(int i=0; i<nMassBins; i++){
    hMassUnfolded->SetBinContent(i+1, yieldsUnfolded(i) );
  }

  //
  // Plot the results
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  c1->cd();
  gStyle->SetOptStat(0);
  float max1 = hMassTrue->GetMaximum();
  hMassTrue->GetYaxis()->SetRangeUser(0, max1*1.3);
  setHistAttributes(hMassTrue, kBlack);
  setHistAttributes(hMassMeasured, kBlue);
  setHistAttributes(hMassUnfolded, kRed);
  hMassTrue->Draw("hist");
  hMassMeasured->Draw("same,PE");
  hMassUnfolded->Draw("same,hist");

  TCanvas *c2 = new TCanvas("c2","c2",100,10,600,600);
  c2->cd();
  gStyle->SetPalette(1);
  hMigrationMatrix->Draw("COLZ");

  TCanvas *c3 = new TCanvas("c3","c3",200,10,600,600);
  c3->cd();
  gStyle->SetPalette(1);
  hResponseMatrix->Draw("COLZ");

  TCanvas *c4 = new TCanvas("c4","c4",300,10,600,600);
  c4->cd();
  gStyle->SetPalette(1);
  hUnfoldingMatrix->Draw("COLZ");

  

}

void setHistAttributes(TH1F *hist, Color_t color){

  hist->SetLineWidth(2);
  hist->SetLineColor(color);
}

