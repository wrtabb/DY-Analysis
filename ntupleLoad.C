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

 cout << "Loading ntuples" << endl;
 //The names of every directory being loaded
 TString dirNames[numChains] = {DYLL_M10to50,DYLL_M50toInf,DYLL_M100to200,DYLL_M200to400,
  DYLL_M400to500,DYLL_M500to700,DYLL_M700to800,DYLL_M800to1000,DYLL_M1000to1500,
  DYLL_M1500to2000,DYLL_M2000to3000};

 TChain*chains[numChains];
 vector <TString> *subFiles[numChains];
 for(int iChain=0;iChain<numChains;iChain++){
  subFiles[iChain] = new vector<TString>;
  if(iChain==MC10to50){
   subFiles[iChain]->push_back(dirNames[iChain]+"/v1/");
   subFiles[iChain]->push_back(dirNames[iChain]+"/v2/");
   subFiles[iChain]->push_back(dirNames[iChain]+"/ext1v1/");
  }
 else subFiles[iChain]->push_back(dirNames[iChain]);
 }

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

