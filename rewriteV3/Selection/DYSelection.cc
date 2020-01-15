#include "DYSelection.hh"
#include <TBranch.h>
#include <iostream>
#include <TStopwatch.h>
#include <TTimeStamp.h>

//Constructor
// Want to set the ntuple version being used
// And which lepton is being analyzed
// And which samples to load
// For samples, choose one of these:
// LL           : Leptons
// LL_RECO      : Leptons, reco only
// MUMU         : muons (!!!!!not included yet!!!!!!)
// TAUTAU       : taus (used for background)
// EW           : electroweak (background)
// TT           : tops (background)
// DATA         : data
DYSelection::DYSelection(NtupleVersion ntup, LepType lepType, SampleType sampleType)
{
 std::vector<TString> dirNames;
 if(sampleType==LL) dirNames = dirNamesEE;
 else if(sampleType==LL_RECO) dirNames = dirNamesEEReco;
 else if(sampleType==TAUTAU) dirNames = dirNamesTAUTAU;
 else if(sampleType==EW) dirNames = dirNamesEW;
 else if(sampleType==TT) dirNames = dirNamesTT;
 else if(sampleType==DATA) dirNames = dirNamesData;

 LoadTrees(ntup,dirNames,sampleType,lepType);
}

Long64_t DYSelection::LoadTrees(NtupleVersion ntup,std::vector<TString>dirNames,SampleType sampleType,LepType lepType)
{
 TTimeStamp ts_start;
 cout << "Begin loading trees:" << endl;
 cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
 TStopwatch totaltime;
 totaltime.Start();

 bool isMC = true;
 bool isReco = true;//this refers to the reco only samples
 if(sampleType==DATA){
  isMC = false;
  isReco = false;
 }
 if(sampleType==LL) isReco = false;
 cout << "isReco = " << isReco << endl;

 const int numChains = dirNames.size();
 TString files;
 Long64_t subDirectorySize;
 Long64_t totalentries = -1;
 if(ntup==V2P6){
  cout << "V2P6 is not yet implemented" << endl;
  return 0;
 }

 else if(ntup==V2P3){
  TString fileNames;
  if(!isReco && isMC) fileNames = "/*.root";
  else fileNames = "/skims_0002/*.root";
  vector <TString> *subFiles[numChains];
  for(int iChain=0;iChain<numChains;iChain++){
   subFiles[iChain] = new vector<TString>;
   if(((sampleType==LL || sampleType==LL_RECO) && iChain==EE10to50) ||
      (sampleType==TAUTAU && iChain==TAUTAU10to50)) {
    subFiles[iChain]->push_back(dirNames.at(iChain)+"_ext1v1");
    subFiles[iChain]->push_back(dirNames.at(iChain)+"_v1");
    subFiles[iChain]->push_back(dirNames.at(iChain)+"_v2");
   }
   else if((sampleType==LL && iChain==EE100to200) ||
           (sampleType==TAUTAU && iChain==TAUTAU100to200) ||
           (sampleType==EW && iChain==WJETS)) {
    subFiles[iChain]->push_back(dirNames.at(iChain));
    subFiles[iChain]->push_back(dirNames.at(iChain)+"_ext");
   }
   else if(sampleType==DATA && iChain==RUNH) {
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
 }//end if ntup==V2P3

 else{
  cout << "ERROR: Ntuple version not correctly chosen!" << endl;
  return 0;
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

 //-----Initialize branches-----//
 InitBranches(numChains,isMC,isReco,lepType);

 //-----Open all needed files and load historams-----//
 LoadHistograms();
 
 return totalentries;
}//end LoadTrees()
   
void DYSelection::InitBranches(int numChains,bool isMC,bool isReco,LepType lepType)
{
 cout << "numChains = " << numChains << endl;
 for(int iChain=0;iChain<numChains;iChain++){

  //-----HLT Branches-----//
  chains[iChain]->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
  chains[iChain]->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
  chains[iChain]->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
  chains[iChain]->SetBranchAddress("HLT_trigName",&pHLT_trigName);

  //-----Reco-level branches-----//
  chains[iChain]->SetBranchAddress("Nelectrons", &Nelectrons, &b_Nelectrons);
  chains[iChain]->SetBranchAddress("nVertices", &nVertices, &b_nVertices);
  chains[iChain]->SetBranchAddress("nPileUp", &nPileUp, &b_nPileUp);
  chains[iChain]->SetBranchAddress("Electron_pT", &Electron_pT, &b_Electron_pT);
  chains[iChain]->SetBranchAddress("Electron_eta",&Electron_eta, &b_Electron_eta);
  chains[iChain]->SetBranchAddress("Electron_phi",&Electron_phi, &b_Electron_phi);
  chains[iChain]->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,
                                   &b_Electron_passMediumID);

  //only MC samples get gen weights
  if(isMC){
   chains[iChain]->SetBranchAddress("GENEvt_weight",&GENEvt_weight,&b_GENEvt_weight);
   //samples that are reco only do not get gen level branches
   if(!isReco){
    //-----Gen-level branches-----//
    chains[iChain]->SetBranchAddress("GENnPair", &GENnPair, &b_GENnPair);
    chains[iChain]->SetBranchAddress("GENLepton_eta", &GENLepton_eta, &b_GENLepton_eta);
    chains[iChain]->SetBranchAddress("GENLepton_phi",&GENLepton_phi, &b_GENLepton_phi);
    chains[iChain]->SetBranchAddress("GENLepton_pT",&GENLepton_pT, &b_GENLepton_pT);
    chains[iChain]->SetBranchAddress("GENLepton_ID",&GENLepton_ID, &b_GENLepton_ID);
    chains[iChain]->SetBranchAddress("GENLepton_isHardProcess",&GENLepton_isHardProcess,
                                     &b_GENLepton_isHardProcess);
    chains[iChain]->SetBranchAddress("GENLepton_fromHardProcessFinalState",
                                     &GENLepton_fromHardProcessFinalState,
                                     &b_GENLepton_fromHardProcessFinalState);
   }//end if !isReco   
  }//end isMC
 }//end loop over chains
}//end InitBranches()

//-----Load files and histograms-----//
void DYSelection::LoadHistograms()
{
 //-----For pileup weights-----//
 pileupRatioFile  = new TFile(pileupRatioName);
 hPileupRatio = (TH1F*)pileupRatioFile->Get("hPileupRatio");
 //-----For Scale Factors-----//
 fileLeg2SF  = new TFile(leg2SFName);
 fileMedIDSF = new TFile(medIDSFName);
 fileRecoSF  = new TFile(recoSFName);
 hLeg2SF  = (TH2F*) fileLeg2SF->Get("EGamma_SF2D");
 hMedIDSF = (TH2F*)fileMedIDSF->Get("EGamma_SF2D");
 hRecoSF  = (TH2F*) fileRecoSF->Get("EGamma_SF2D");
}

