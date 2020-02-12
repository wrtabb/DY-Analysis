#include "NtuplesV2P6Location.h"
#include "VariableList.h"

Long64_t LoadTrees(std::vector<TString>dirNames,SampleType sampleType,LepType lepType);
void InitializeBranches(TChain*chain,bool isMC,LepType lepType);
//And which lepton is being analyzed
//And which samples to load
//For samples, choose one of these:
// LL		: Leptons
// EW		: electroweak (background)
// TT		: tops (background)
// DATA		: data
void GetHists(SampleType sampleType,LepType lepType)
{
 std::vector<TString> dirNames;
 if(sampleType==LL) dirNames = dirNamesLL;
 else if(sampleType==EW) dirNames = dirNamesEW;
 else if(sampleType==TT) dirNames = dirNamesTT;
 else if(sampleType==DATA) dirNames = dirNamesData;

 int dirSize = dirNames.size();
 if(sampleType==LL){
  if(lepType==ELE){
   for(int i=0;i<dirSize;i++){
    dirNames.at(i)+= "/EE";
   }
  }
  else if(lepType==MUON){
   for(int i=0;i<dirSize;i++){
    dirNames.at(i)+= "/MuMu";
   }
  }
  else if(lepType==TAU){
   for(int i=0;i<dirSize;i++){
    dirNames.at(i)+= "/TauTau";
   }
  }
 }//end if sampletype
 LoadTrees(dirNames,sampleType,lepType);
}
 
Long64_t LoadTrees(std::vector<TString>dirNames,SampleType sampleType,LepType lepType)
{
 TTimeStamp ts_start;
 cout << "Begin loading trees:" << endl;
 cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
 TStopwatch totaltime;
 totaltime.Start();

 bool isMC = true;
 if(sampleType==DATA){
  isMC = false;
 }

 const int numChains = dirNames.size();
 TString files;
 Long64_t subDirectorySize;
 Long64_t totalentries = -1;

 TString fileNames;
 fileNames = "/skims_0002/*.root";
 vector <TString> *subFiles[numChains];
 for(int iChain=0;iChain<numChains;iChain++){
  subFiles[iChain] = new vector<TString>;
  if(sampleType==LL && iChain==M10to50){
   subFiles[iChain]->push_back(dirNames.at(iChain)+"/ext1v1");
   subFiles[iChain]->push_back(dirNames.at(iChain)+"/v1");
   subFiles[iChain]->push_back(dirNames.at(iChain)+"/v2");
  }
  else if(sampleType==LL && iChain==M50to100){
   subFiles[iChain]->push_back(dirNames.at(iChain)+"/base");
  }
  else if(sampleType==EW && iChain==W_PLUS_JETS) {
   subFiles[iChain]->push_back(dirNames.at(iChain));
   subFiles[iChain]->push_back(dirNames.at(iChain)+"_ext");
   subFiles[iChain]->push_back(dirNames.at(iChain)+"_ext2v5");
  }
  else if(sampleType==DATA && iChain==RUN_H) {
   subFiles[iChain]->push_back(dirNames.at(iChain)+"ver2");
   subFiles[iChain]->push_back(dirNames.at(iChain)+"ver3");
  }
  else if(sampleType==TT && iChain==TT0to700) {
   subFiles[iChain]->push_back(dirNames.at(iChain));
   subFiles[iChain]->push_back(dirNames.at(iChain)+"Backup");
  }
  else subFiles[iChain]->push_back(dirNames.at(iChain));

 }//end loop over iChains 

 totalentries = 0;
 TChain*chains[numChains];
 for(int iChain=0;iChain<numChains;iChain++){
  chains[iChain] = new TChain(treeName);
  subDirectorySize = subFiles[iChain]->size();
  for(int k=0;k<subDirectorySize;k++){
   files=subFiles[iChain]->at(k);
   files+=fileNames;
   chains[iChain]->Add(files);
   cout << files << endl;
   cout << chains[iChain]->GetEntries() << " events loaded" << endl;
   if(chains[iChain]->GetEntries()==0){
    cout << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "ERROR: Broken files or files not found in: " << endl;
    cout << files << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << endl;
    return 0;
   }
  }//end loop over files
 totalentries=totalentries+chains[iChain]->GetEntries();
 InitializeBranches(chains[iChain],isMC,lepType);
 }

 cout << "Total Events Loaded: " << totalentries << endl;
 cout << endl;

 totaltime.Stop();
 Double_t TotalCPURunTime = totaltime.CpuTime();
 Double_t TotalRunTime = totaltime.RealTime();
 TTimeStamp ts_end;
 cout << endl;
 cout << "End loading trees:" << endl;
 cout << "**************************************************************************" << endl;
 cout << "Total CPU RunTime: " << TotalCPURunTime/60 << " minutes" << endl;
 cout << "Total Real RunTime: " << TotalRunTime/60 << " minutes" << endl;
 cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
 cout << "**************************************************************************" << endl;
 cout << endl;

 
 
 return totalentries;
}


void InitializeBranches(TChain*chain,bool isMC,LepType lepType)
{
 //-----HLT Branches-----//
 chain->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
 chain->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
 chain->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
 chain->SetBranchAddress("HLT_trigName",&pHLT_trigName);

 //-----Reco-level branches-----//
 chain->SetBranchAddress("nVertices", &nVertices, &b_nVertices);
 chain->SetBranchAddress("nPileUp", &nPileUp, &b_nPileUp);
 if(lepType==ELE){
  chain->SetBranchAddress("Nelectrons", &Nelectrons, &b_Nelectrons);
  chain->SetBranchAddress("Electron_pT", &Electron_pT, &b_Electron_pT);
  chain->SetBranchAddress("Electron_eta",&Electron_eta, &b_Electron_eta);
  chain->SetBranchAddress("Electron_phi",&Electron_phi, &b_Electron_phi);
  chain->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,
                                    &b_Electron_passMediumID);
 }
 else if(lepType==MUON){
  cout << "The muon branches for reco-muons do not exist in current ntuples" << endl;
  cout << "This will be fixed as soon as possible." << endl;
  return ;
 }

 //-----Gen-level branches-----//
 if(isMC){
  chain->SetBranchAddress("GENEvt_weight",&GENEvt_weight,&b_GENEvt_weight);
  chain->SetBranchAddress("GENnPair", &GENnPair, &b_GENnPair);
  chain->SetBranchAddress("GENLepton_eta", &GENLepton_eta, &b_GENLepton_eta);
  chain->SetBranchAddress("GENLepton_phi",&GENLepton_phi, &b_GENLepton_phi);
  chain->SetBranchAddress("GENLepton_pT",&GENLepton_pT, &b_GENLepton_pT);
  chain->SetBranchAddress("GENLepton_ID",&GENLepton_ID, &b_GENLepton_ID);
  chain->SetBranchAddress("GENLepton_isHardProcess",&GENLepton_isHardProcess,
                                   &b_GENLepton_isHardProcess);
  chain->SetBranchAddress("GENLepton_fromHardProcessFinalState",
                                   &GENLepton_fromHardProcessFinalState,
                                   &b_GENLepton_fromHardProcessFinalState);
 }//end isMC
 return;
}//end Initialize Branches
