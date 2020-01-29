#include "NtuplesV2P6Location.h"
#include "VariableList.h"

Long64_t LoadTrees(std::vector<TString>dirNames,SampleType sampleType,LepType lepType);
//Want to set the ntuple version being used
//And which lepton is being analyzed
//And which samples to load
//For samples, choose one of these:
// EE		: electrons
// EE_RECO	: electrons, reco only
// MUMU		: muons (!!!!!not included yet!!!!!!)
// TAUTAU	: taus (used for background)
// EW		: electroweak (background)
// TT		: tops (background)
// DATA		: data
void tempLoadNtuples(SampleType sampleType,LepType lepType)
{
 std::vector<TString> dirNames;
 if(sampleType==LL) dirNames = dirNamesLL;
 else if(sampleType==EW) dirNames = dirNamesEW;
 else if(sampleType==TT) dirNames = dirNamesTT;
 else if(sampleType==DATA) dirNames = dirNamesData;

 int dirSize = dirNames.size();
 if(lepType==ELE){
  for(int i=0;i<dirSize;i++){
   dirNames.at(i)+= "EE/";
  }
 }
 else if(lepType==MUON){
  for(int i=0;i<dirSize;i++){
   dirNames.at(i)+= "MuMu/";
  }
 }
 else if(lepType==TAU){
  for(int i=0;i<dirSize;i++){
   dirNames.at(i)+= "TauTau/";
  }
 }
 
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
 fileNames = "skims_0002/*.root";
 vector <TString> *subFiles[numChains];
 for(int iChain=0;iChain<numChains;iChain++){
  subFiles[iChain] = new vector<TString>;
  if(sampleType==LL && iChain==M10to50){
   subFiles[iChain]->push_back(dirNames.at(iChain)+"ext1v1/");
   subFiles[iChain]->push_back(dirNames.at(iChain)+"v1/");
   subFiles[iChain]->push_back(dirNames.at(iChain)+"v2/");
  }
  else if(sampleType==LL && iChain==M50to100){
   subFiles[iChain]->push_back(dirNames.at(iChain)+"base/");
  }
  else if(sampleType==EW && iChain==W_PLUS_JETS) {
   subFiles[iChain]->push_back(dirNames.at(iChain));
   subFiles[iChain]->push_back(dirNames.at(iChain)+"_ext/");
   subFiles[iChain]->push_back(dirNames.at(iChain)+"_ext2v5/");
  }
  else if(sampleType==DATA && iChain==RUN_H) {
   subFiles[iChain]->push_back(dirNames.at(iChain)+"ver2/");
   subFiles[iChain]->push_back(dirNames.at(iChain)+"ver3/");
  }
  else if(sampleType==TT && iChain==TT0to700) {
   subFiles[iChain]->push_back(dirNames.at(iChain));
   subFiles[iChain]->push_back(dirNames.at(iChain)+"Backup/");
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
