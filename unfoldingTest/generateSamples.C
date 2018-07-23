#include "TRandom.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

// Declarations
enum Model {
  MODEL1
};

// Constants
const int nDataEvents = 10000; // events mimicking real data
const int nSimEvents = 100000; // events mimicking MC used for unfolding procedure

const Model model = MODEL1;
const TString outDataFileName = "sampleMockData.root";
const TString outSimFileName = "sampleSimulation.root";

// Main function

void generateSamples(){

  //
  // Models
  //
  TF1 *funcTruthDataShape = new TF1("funcTruthDataShape","([0]/x)+gaus(1)",10,200);
  funcTruthDataShape->SetParameters(10,1,90,7);

  TF1 *funcResolutionModel = new TF1("funcResolutionModel","gaus(0)", -20, 20);
  funcResolutionModel->SetParameters(1,0,3);

  TF1 *funcSimModelShape = 0;
  TF1 *funcSimResolutionModel = 0;

  if( model == MODEL1 ){
    // Example: in "simulation" the true shape and resolution match real life
    funcSimModelShape = (TF1*)funcTruthDataShape->Clone("funcSimModelShape");
    funcSimResolutionModel = (TF1*)funcResolutionModel->Clone("funcSimResolutionModel");
    printf("Setting up simulation model: exact match to \"real\" data\n");
  }else{
    printf("unknown model requested\n");
    return;
  }
  
  // Variables we will store
  float massTrue, massMeasured;
  
  // 
  // Generate and save mock data sample
  //
  printf("Generating mock data\n");
  TFile fMockData(outDataFileName,"recreate");
  TTree *treeData = new TTree("mockData", "Mock Data");
  treeData->Branch("massTrue", &massTrue, "massTrue/F");
  treeData->Branch("massMeasured", &massMeasured, "massMeasured/F");
  for(int i=0; i<nDataEvents; i++){
    massTrue = funcTruthDataShape->GetRandom();
    // the "true mass" is smeared according to the resolutino function
    // that simulates mismeasurement with a real detector
    massMeasured = massTrue + funcResolutionModel->GetRandom();
    treeData->Fill();
  }
  treeData->Write();
  fMockData.Close();


  // 
  // Generate and save simulated sample
  //
  printf("Generating simulation\n");
  TFile fSim(outSimFileName,"recreate");
  TTree *treeSim = new TTree("sim", "Simulated sample");
  treeSim->Branch("massTrue", &massTrue, "massTrue/F");
  treeSim->Branch("massMeasured", &massMeasured, "massMeasured/F");
  for(int i=0; i<nSimEvents; i++){
    massTrue = funcSimModelShape->GetRandom();
    // the "true mass" is smeared according to the resolutino function
    // that simulates mismeasurement with a real detector
    massMeasured = massTrue + funcSimResolutionModel->GetRandom();
    treeSim->Fill();
  }
  treeSim->Write();
  fSim.Close();

  printf("Done\n");
}
