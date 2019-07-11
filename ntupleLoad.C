#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"
#include "/home/hep/wrtabb/git/DY-Analysis/headers/NtuplesV2P6Location.h"
#include "/home/hep/wrtabb/git/DY-Analysis/headers/DrellYanCuts.h"
#include "/home/hep/wrtabb/git/DY-Analysis/headers/Functions.h"

void ntupleLoad()
{
 TH1::SetDefaultSumw2();
 TTimeStamp ts_start;
 cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
 TStopwatch totaltime;
 totaltime.Start();
 bool isMC;//is Monte Carlo
 gStyle->SetOptStat(0);

 ofstream output("ntupleEntries.txt");

 //-----Defining branches-----//
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

 cout << "Loading ntuples" << endl;
 //The names of every directory being loaded
 TString dirNames[numChains] = {LLM10to50,LLM50toInf,LLM100to200,LLM200to400,LLM400to500,
  LLM500to700,LLM700to800,LLM800to1000,LLM1000to1500,LLM1500to2000,LLM2000to3000};

 TChain*chains[numChains];
 vector <TString> *subFiles[numChains];
 for(int iChain=0;iChain<numChains;iChain++){
  subFiles[iChain] = new vector<TString>;
  if(iChain==MC10to50){
   subFiles[iChain]->push_back(dirNames[iChain]+"v1/");
   subFiles[iChain]->push_back(dirNames[iChain]+"v2/");
   subFiles[iChain]->push_back(dirNames[iChain]+"ext1v1/");
  }
  else if(iChain==MC50to100){
   subFiles[iChain]->push_back(dirNames[iChain]+"base/");
   //subFiles[iChain]->push_back(dirNames[iChain]+"madgraph/");
   //subFiles[iChain]->push_back(dirNames[iChain]+"madgraph_ext/");
  }
  else subFiles[iChain]->push_back(dirNames[iChain]);
 }//end loop over chains 

 TString files;
 Long64_t subDirectorySize;
 Long64_t totalentries = 0;
 for(int iChain=0;iChain<numChains;iChain++){
  chains[iChain] = new TChain(treeName);
  subDirectorySize = subFiles[iChain]->size();
  for(int k=0;k<subDirectorySize;k++){
   TFileCollection filecoll("dum");//Object for creating a list of files in a directory
   files=subFiles[iChain]->at(k);
   files+="*.root";
   filecoll.Add(files);
   chains[iChain]->AddFileInfoList(filecoll.GetList());
   cout << files << endl;
   cout << chains[iChain]->GetEntries() << " events loaded" << endl;
   output << files << ", " << chains[iChain]->GetEntries() << endl;
   if(chains[iChain]->GetEntries()==0){
    cout << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "ERROR: Broken files or files not found in: " << endl;
    cout << files << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << endl;
    return;
   }
  }//end loop over files

  totalentries=totalentries+chains[iChain]->GetEntries();

  //-----Setting addresses for branches-----//
  chains[iChain]->SetBranchAddress("Nelectrons", &Nelectrons, &b_Nelectrons);
  chains[iChain]->SetBranchAddress("nVertices", &nVertices, &b_nVertices);
  chains[iChain]->SetBranchAddress("nPileUp", &nPileUp, &b_nPileUp);
  chains[iChain]->SetBranchAddress("Electron_pT", &Electron_pT, &b_Electron_pT);
  chains[iChain]->SetBranchAddress("Electron_eta",&Electron_eta, &b_Electron_eta);
  chains[iChain]->SetBranchAddress("Electron_phi",&Electron_phi, &b_Electron_phi);
  chains[iChain]->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,
    &b_Electron_passMediumID);
  chains[iChain]->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
  chains[iChain]->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
  chains[iChain]->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
  chains[iChain]->SetBranchAddress("HLT_trigName",&pHLT_trigName);   
  chains[iChain]->SetBranchAddress("GENEvt_weight",&GENEvt_weight,&b_GENEvt_weight);
  chains[iChain]->SetBranchAddress("GENnPair", &GENnPair, &b_GENnPair);
  chains[iChain]->SetBranchAddress("GENLepton_eta", &GENLepton_eta, &b_GENLepton_eta);
  chains[iChain]->SetBranchAddress("GENLepton_phi",&GENLepton_phi, &b_GENLepton_phi);
  chains[iChain]->SetBranchAddress("GENLepton_pT",&GENLepton_pT, &b_GENLepton_pT);
  chains[iChain]->SetBranchAddress("GENLepton_ID",&GENLepton_ID, &b_GENLepton_ID);
  chains[iChain]->SetBranchAddress("GENLepton_isHardProcess",&GENLepton_isHardProcess,
   &b_GENLepton_isHardProcess);
  chains[iChain]->SetBranchAddress
   ("GENLepton_fromHardProcessFinalState",&GENLepton_fromHardProcessFinalState,
   &b_GENLepton_fromHardProcessFinalState);

 }//end iChain loop
  
 cout << "Total Events Loaded: " << totalentries << endl;
 cout << endl;
 output.close(); 
 totaltime.Stop();
 Double_t TotalCPURunTime = totaltime.CpuTime();
 Double_t TotalRunTime = totaltime.RealTime();
 TTimeStamp ts_end;
 cout << endl;
 cout << "**************************************************************************" << endl;
 cout << "Total CPU RunTime: " << TotalCPURunTime/60 << " minutes" << endl;
 cout << "Total Real RunTime: " << TotalRunTime/60 << " minutes" << endl;
 cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;   
 cout << "**************************************************************************" << endl;
 cout << endl;
}

