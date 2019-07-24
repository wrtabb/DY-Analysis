//This script opens the v2.6 ntuples from hadoop
//Then it takes only the branches I need to keep
//Saves them to a new tree
//Saves the tree to my working area
//
//After saving all root files in a particular directory
//It then copies them to hadoop
//And deletes the files from the working directory
//It is done this way to avoid copying huge numbers
//Of files, taking up lots of hard drive space,
//To the working directory

#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"
#include "header.h"

void skimFile(TString fileName,TString fileNameNoDir,TString subdirectory,TString suffix);
void skim(TString directory,TString subdirectory);

double _prefiringweight;

void skimNtuples()
{
 TString directory[nDir];
 for(int i=0;i<2;i++){
  directory[i] = baseDir+subdirectory[i];
  //-----Open root files, create skimmed versions, save to working dorectory-----//
  skim(directory[i],subdirectory[i]);
  //-----Copy newly created files to hadoop-----//
  gSystem->Exec("cp -r ./skimmed_files/"+subdirectory[i]+"skim1 "+directory[i]);
  //-----Delete newly created files from working directory-----//
  gSystem->Exec("rm -r ./skimmed_files/"+subdirectory[i]+"skim1/*.root");
 }
}

void skim(TString directory,TString subdirectory){
 ifstream fileList(directory+"skim_file.txt");
 ifstream fileListNoDir(directory+"temp_file.txt");
 TString suffix = "_skim";
 TString fileName;
 TString fileNameNoDir;
 while(true){
  fileList >> fileName;
  fileListNoDir >> fileNameNoDir;
  if(fileList.eof()) break;
  skimFile(fileName,fileNameNoDir,subdirectory,suffix);
 }
}

void skimFile(TString fileName,TString fileNameNoDir,TString subdirectory,TString suffix){
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

 TString newFileName = "./skimmed_files/"+subdirectory+"skim1/";
 newFileName += fileNameNoDir.ReplaceAll(".root","")+suffix+TString(".root");
 cout << newFileName << endl;
 TFile*newFile = new TFile(newFileName,"recreate");
 TTree*newTree = chain->CloneTree(0);

 newTree->Write();
 newFile->Close();
}

