#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"
void skimFile(TString fileName,TString suffix);
double _prefiringweight;

void skimNtuples()
{
 ifstream fileList("skimFile.txt");
 TString suffix = "_skim";
 TString rootFile;
 while(!fileList.eof()){
  fileList >> rootFile;
  TString fileName = rootFile;
  skimFile(fileName,suffix);
 }
}

void skimFile(TString fileName, TString suffix){
 TFile*file = new TFile(fileName);
 TTree*tree = (TTree*)file->Get(treeName);
 TChain*chain = new TChain(treeName);
 chain->Add(fileName.Data());
 chain->SetBranchStatus("*",0);

 TBranch*b_Nelectrons;
 TBranch*b_Electron_pT;
 TBranch*b_Electron_eta;
 TBranch*b_Electron_phi;
 TBranch*b_Electron_passMediumID;
 TBranch*b_HLT_ntrig;
 TBranch*b_HLT_trigType;
 TBranch*b_HLT_trigFired;
 TBranch*b_GENEvt_weight;
 TBranch*b_nVertices;
 TBranch*b_nPileUp;
 TBranch*b_GENnPair;
 TBranch*b_GENLepton_eta;
 TBranch*b_GENLepton_phi;
 TBranch*b_GENLepton_pT;
 TBranch*b_GENLepton_ID;
 TBranch*b_GENLepton_isHardProcess;
 TBranch*b_GENLepton_fromHardProcessFinalState;
 TBranch*b_prefiringweight;

 chain->SetBranchStatus("Nelectrons",1);
 chain->SetBranchStatus("nVertices",1);
 chain->SetBranchStatus("nPileUp",1);
 chain->SetBranchStatus("Electron_pT",1);
 chain->SetBranchStatus("Electron_eta",1);
 chain->SetBranchStatus("Electron_phi",1);
 chain->SetBranchStatus("Electron_passMediumID",1);
 chain->SetBranchStatus("HLT_ntrig",1);
 chain->SetBranchStatus("HLT_trigType",1);
 chain->SetBranchStatus("HLT_trigFired",1);
 chain->SetBranchStatus("HLT_trigName",1);
 chain->SetBranchStatus("GENEvt_weight",1);
 chain->SetBranchStatus("GENnPair",1);
 chain->SetBranchStatus("GENLepton_eta",1);
 chain->SetBranchStatus("GENLepton_phi",1);
 chain->SetBranchStatus("GENLepton_pT",1);
 chain->SetBranchStatus("GENLepton_ID",1);
 chain->SetBranchStatus("GENLepton_isHardProcess",1);
 chain->SetBranchStatus("GENLepton_fromHardProcessFinalState",1);
 chain->SetBranchStatus("_prefiringweight",1);

 TString newFileName = "./skim1/";
 newFileName += fileName.ReplaceAll(".root","")+suffix+TString(".root");
 TFile*newFile = new TFile(newFileName,"recreate");
 TTree*newTree = chain->CloneTree(0);

 newTree->Write();
 newFile->Close();
}

