#include "TRandom.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TLorentzVector.h"

const double pi = TMath::Pi();
double eMass = 0.000511;
const double phiRange = pi;
const double etaRange = 2.4;
const double ptHigh = 28;
const double ptLow =17;
const TString fileName = "ntuple_skim_638.root";
const TString treeName = "recoTree/DYtree";
const TString toyModelName = "toyData.root";
const int nEvents = 10000;
const int nBinsPhi = 50;
const float avgLeadpT = 50;
const float avgSubpT = 35;

void toyModel()
{
  TRandom*random = new TRandom3(time(0));
  //Data to store
  float phiTrue,etaTrue,ptTrue,massTrue;
  float phiMeasured,etaMeasured,ptMeasured,massMeasured;
  float sign;
  TFile fToyData(toyModelName,"recreate");
  TTree*tree = new TTree("toyData","Toy Data");
  tree->Branch("phiTrue",&phiTrue,"phiTrue/F");
  tree->Branch("etaTrue",&etaTrue,"etaTrue/F");
  tree->Branch("ptTrue",&ptTrue,"ptTrue/F");
  tree->Branch("massTrue",&massTrue,"massTrue/f");
  tree->Branch("phiMeasured",&phiMeasured,"phiMeasured/F");
  tree->Branch("etaMeasured",&etaMeasured,"etaMeasured/F");
  tree->Branch("ptMeasured",&ptMeasured,"ptMeasured/F");
  tree->Branch("massMeasured",&massMeasured,"massMeasured/f");

  TF1*ptShapeLead = new TF1("ptShape","[0]/x+gaus(50)",10,500);
  TF1*ptShapeSub = new TF1("ptShape","[0]/x+gaus(35)",10,500);
  TF1*fResolutionModel = new TF1("fPhiResolutionModel","gaus(0)",-20,20);
  ptShapeLead->SetParameter(0,1.0);
  ptShapeSub->SetParameter(0,1.0);
  fResolutionModel->SetParameters(1,0,3);
  TLorentzVector ele1;
  TLorentzVector ele2;
  TLorentzVector ele1Measured;
  TLorentzVector ele2Measured;

  for(int i=0;i<nEvents;i++){
    //Electron 1
    if(random->Rndm()<=0.5) sign = 1;
    else sign = -1;
    phiTrue = sign*phiRange*random->Rndm(); 
    if(random->Rndm()<=0.5) sign = 1;
    else sign = -1;
    etaTrue = sign*etaRange*random->Rndm();
    ptTrue = ptShapeLead->GetRandom();
    ele1.SetPtEtaPhiM(ptTrue,etaTrue,phiTrue,eMass);

    phiMeasured = phiTrue+fResolutionModel->GetRandom(-0.5,0.5);
    etaMeasured = etaTrue+fResolutionModel->GetRandom(-0.5,0.5);
    ptMeasured = ptTrue+fResolutionModel->GetRandom(-5,5);
    ele1Measured.SetPtEtaPhiM(ptMeasured,etaMeasured,phiMeasured,eMass);
    
    //Electron 2
    if(random->Rndm()<=0.5) sign = 1;
    else sign = -1;
    phiTrue = sign*phiRange*random->Rndm(); 
    if(random->Rndm()<=0.5) sign = 1;
    else sign = -1;
    etaTrue = sign*etaRange*random->Rndm();
    ptTrue = ptShapeSub->GetRandom();
    ele2.SetPtEtaPhiM(ptTrue,etaTrue,phiTrue,eMass);

    phiMeasured = phiTrue+fResolutionModel->GetRandom(-0.5,0.5);
    etaMeasured = etaTrue+fResolutionModel->GetRandom(-0.5,0.5);
    ptMeasured = ptTrue+fResolutionModel->GetRandom(-5,5);
    ele2Measured.SetPtEtaPhiM(ptMeasured,etaMeasured,phiMeasured,eMass);

    massTrue = (ele1+ele2).M();
    massMeasured = (ele1Measured+ele2Measured).M();
  tree->Write();
    if(massTrue<15 || massTrue>3000) continue;
    if(!((ele1.Pt()>ptLow&&ele2.Pt()>ptHigh)||(ele1.Pt()>ptHigh&&ele2.Pt()>ptLow)))
      continue;
  }
  fToyData.Close();
}
