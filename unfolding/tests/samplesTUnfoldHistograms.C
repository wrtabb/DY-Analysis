
#include <iostream>
#include <map>
#include <cmath>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include "TUnfoldBinning.h"

using namespace std;

const TString dataFileName = "sampleMockData.root";
const TString simFileName = "sampleSimulation.root";
const TString simFileName2 = "sampleSimulation2.root";

const TString dataTreeName = "mockData";
const TString simTreeName = "sim";

const int nMassBins =15;
const int massMin = 20;
const int massMax = 200;

void samplesTUnfoldHistograms()
{
  TH1::SetDefaultSumw2();
  float massTrue, massMeasured;  
  TBranch*b_massTrue;
  TBranch*b_massMeasured;

  TFile*dataFile = new TFile(dataFileName);
  TTree*dataTree = (TTree*)dataFile->Get(dataTreeName);
  dataTree->SetBranchAddress("massTrue",&massTrue,&b_massTrue);
  dataTree->SetBranchAddress("massMeasured", &massMeasured, &b_massMeasured);

  TFile*simFile = new TFile(simFileName);
  TTree*simTree = (TTree*)simFile->Get(simTreeName);
  simTree->SetBranchAddress("massTrue",&massTrue,&b_massTrue);
  simTree->SetBranchAddress("massMeasured",&massMeasured,&b_massMeasured);

  TFile*simFile2 = new TFile(simFileName2);
  TTree*simTree2 = (TTree*)simFile2->Get(simTreeName);
  simTree2->SetBranchAddress("massTrue",&massTrue,&b_massTrue);
  simTree2->SetBranchAddress("massMeasured",&massMeasured,&b_massMeasured);

  //Data distributions
  TH1F*hDataTrue = new TH1F("hDataTrue","",nMassBins,massMin,massMax);
  TH1F*hDataMeasured = new TH1F("hDataMeasured","",nMassBins,massMin,massMax);
  Long64_t nentries,mentries;
  nentries = dataTree->GetEntries();
  for(int i=0;i<nentries;i++){
    dataTree->GetEntry(i);
    hDataTrue->Fill(massTrue);
    hDataMeasured->Fill(massMeasured);
  }
  
  //MC distributions
  TH2F*hMigrationMatrix = new TH2F("hMigrationMatrix","",nMassBins,massMin,massMax,nMassBins,
    massMin,massMax); 
  TH1F*hMCTrue = new TH1F("hMCTrue","",nMassBins,massMin,massMax);
  TH1F*hMCMeasured = new TH1F("hMCMeasured","",nMassBins,massMin,massMax);
 
  TH2F*hMigrationMatrix2 = new TH2F("hMigrationMatrix2","",nMassBins,massMin,massMax,nMassBins,
    massMin,massMax); 
  TH1F*hMCTrue2 = new TH1F("hMCTrue2","",nMassBins,massMin,massMax);
  TH1F*hMCMeasured2 = new TH1F("hMCMeasured2","",nMassBins,massMin,massMax);

  mentries = simTree->GetEntries();
  for(int i=0;i<nentries;i++){
    simTree->GetEntry(i);
    hMCTrue->Fill(massTrue);
    hMCMeasured->Fill(massMeasured);
    hMigrationMatrix->Fill(massMeasured,massTrue);

    simTree2->GetEntry(i);
    hMCTrue2->Fill(massTrue);
    hMCMeasured2->Fill(massMeasured);
    hMigrationMatrix2->Fill(massMeasured,massTrue);
  }
  TFile*histFile = new TFile("sampleHistos.root","recreate");
  hMCTrue->Write();
  hMCMeasured->Write();
  hMigrationMatrix->Write();
  hMCTrue2->Write();
  hMCMeasured2->Write();
  hMigrationMatrix2->Write();
  histFile->Write();
  histFile->Close();
 
}
