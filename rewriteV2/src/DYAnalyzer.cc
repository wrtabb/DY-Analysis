#include "../include/DYAnalyzer.hh"
#include <TBranch.h>
#include <iostream>
#include <TStopwatch.h>
#include <TTimeStamp.h>

using namespace std;

std::vector<TString> dirNamesEE = {
 //Electrons
 EEM10to50,
 EEM50to100,
 EEM100to200,
 EEM200to400,
 EEM400to500,
 EEM500to700,
 EEM700to800,
 EEM800to1000,
 EEM1000to1500,
 EEM1500to2000,
 EEM2000to3000
};

std::vector<TString> dirNamesEEReco = {
 //Electrons: Reco only
 EEM10to50Reco,
 EEM50to100Reco,
 EEM100to200Reco,
 EEM200to400Reco,
 EEM400to500Reco,
 EEM500to700Reco,
 EEM700to800Reco,
 EEM800to1000Reco,
 EEM1000to1500Reco,
 EEM1500to2000Reco,
 EEM2000to3000Reco
};

std::vector<TString> dirNamesTAUTAU = {
 //TauTau
 Taus10to50Reco,
 Taus50to100Reco,
 Taus100to200Reco,
 Taus200to400Reco,
 Taus400to500Reco,
 Taus500to700Reco,
 Taus700to800Reco,
 Taus800to1000Reco,
 Taus1000to1500Reco,
 Taus1500to2000Reco,
 Taus2000to3000Reco
};
std::vector<TString> dirNamesEW = {
 //EW
 WJetsReco,
 WWReco,
 WWtoLNuReco,
 ZZtoLReco,
 WZReco,
 WZtoLNuReco
};

std::vector<TString> dirNamesTT = {
 //Tops
 TT0to700Reco,
 TT700to1000Reco,
 TT1000andUpReco,
 tWReco,
 tWantiReco
};

std::vector<TString> dirNamesData = {
 //Data
 runB,
 runC,
 runD,
 runE,
 runF,
 runG,
 runH
};

//Constructor
//Want to set the ntuple version being used
//And which lepton is being analyzed
//And which samples to load
//For samples, choose one of these:
//    EE 	: Electrons
//    EE_RECO	: Electrson reco only
//    TAUTAU	: Taus
//    EW	: Electroweak
//    TT	: Tops
//    DATA	: Data
//
//    NOTE:	EE is for trees with gen and reco level
//    		EE_RECO are skimmed to only contain reco events
DYAnalyzer::DYAnalyzer(NtupleVersion ntup, LepType lepType, SampleType sampleType)
{
 std::vector<TString> dirNames;
 if(sampleType==EE) dirNames = dirNamesEE;
 else if(sampleType==EE_RECO) dirNames = dirNamesEEReco;
 else if(sampleType==TAUTAU) dirNames = dirNamesTAUTAU;
 else if(sampleType==EW) dirNames = dirNamesEW;
 else if(sampleType==TT) dirNames = dirNamesTT;
 else if (sampleType==DATA) dirNames = dirNamesData;

 //Calling an instance of this class automatically begins loading trees which can take
 //up to several minutes
 LoadTrees(ntup,dirNames,sampleType);
 //GetVars(sampleType);
} 

//Load all trees
Long64_t DYAnalyzer::LoadTrees(NtupleVersion ntup,std::vector<TString>dirNames,SampleType sampleType)
{
 TTimeStamp ts_start;
 cout << "Begin loading trees:" << endl;
 cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
 TStopwatch totaltime;
 totaltime.Start();

 bool isMC = true;
 bool isReco = true;//this refers to the reco only samples
 if(sampleType==DATA) isMC = false;
 if(sampleType==EE || sampleType==DATA) isReco = false;
 const int numChains = dirNames.size();
 TString files;
 Long64_t subDirectorySize;
 Long64_t totalentries = -1;
 //The loading process is different for different ntuples versions because there are slightly
 //different samples included as well as different directory structure
 if(ntup==V2P6){
  cout << "V2P6 is not yet implemented" << endl;
  return 0;
 }

 else if(ntup==V2P3){
  TString fileNames;
  if(isMC && !isReco) fileNames = "/*.root";
  else fileNames = "/skims_0002/*.root";
  vector <TString> *subFiles[numChains];
  for(int iChain=0;iChain<numChains;iChain++){
   subFiles[iChain] = new vector<TString>;
   if(((sampleType==EE || sampleType==EE_RECO) && iChain==EE10to50) || 
      (sampleType==TAUTAU && iChain==TAUTAU10to50)) {
    subFiles[iChain]->push_back(dirNames.at(iChain)+"_ext1v1");
    subFiles[iChain]->push_back(dirNames.at(iChain)+"_v1");
    subFiles[iChain]->push_back(dirNames.at(iChain)+"_v2");
   }
   else if((sampleType==EE && iChain==EE100to200) ||
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

  }//end loop over chains 

  totalentries = 0;
  for(int iChain=0;iChain<numChains;iChain++){
   cout << "Chain: " << iChain << endl;
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

 //--Initialize branches-----//
 InitBranches(isMC,isReco);
 //-----Open all needed files and load histograms-----//
 LoadHistograms();

 return totalentries;
}//end LoadTrees()

void DYAnalyzer::InitBranches(bool isMC,bool isReco)
{
 for(int iChain=0;iChain<numChains;iChain++){

  //-----HLT Branches-----//
  chains[iChain]->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
  chains[iChain]->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
  chains[iChain]->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
  chains[iChain]->SetBranchAddress("HLT_trigName",&pHLT_trigName);

  chains[iChain]->SetBranchAddress("GENEvt_weight",&GENEvt_weight,&b_GENEvt_weight);

  //-----Reco-level branches-----//
  chains[iChain]->SetBranchAddress("Nelectrons", &Nelectrons, &b_Nelectrons);
  chains[iChain]->SetBranchAddress("nVertices", &nVertices, &b_nVertices);
  chains[iChain]->SetBranchAddress("nPileUp", &nPileUp, &b_nPileUp);
  chains[iChain]->SetBranchAddress("Electron_pT", &Electron_pT, &b_Electron_pT);
  chains[iChain]->SetBranchAddress("Electron_eta",&Electron_eta, &b_Electron_eta);
  chains[iChain]->SetBranchAddress("Electron_phi",&Electron_phi, &b_Electron_phi);
  chains[iChain]->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,
    &b_Electron_passMediumID);
 
  if(!isReco && isMC){
   //-----Gen-level branches-----//
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
  }
 }
}//end InitBranches()

//-----Lepton selection-----//
//-----Get exactly two leptons from hard process and FSR-----//
int DYAnalyzer::GetGenLeptons(LepType lepType,int &idxHardLep1,int &idxHardLep2,
                              int &idxFSRLep1,int &idxFSRLep2)
{
 //PDG ID codes
 int lepID = 0;
 if      (lepType==ELE)  lepID = 11;
 else if (lepType==MUON) lepID = 13;
 else {
  cout << "ERROR: Appropriate lepton not selected." << endl;
  cout << "NOTE: This analysis does not handle taus." << endl;
  return 0;
 }
 int nDileptons = 0;

 for(int iLep=0;iLep<GENnPair;iLep++){
  for(int jLep=iLep+1;jLep<GENnPair;jLep++){
   //lepton ID selection
   if(!(abs(GENLepton_ID[iLep])==lepID && abs(GENLepton_ID[jLep])==lepID)) continue;
   //Require opposite sign for electrons
   if(GENLepton_ID[iLep]*GENLepton_ID[jLep]>0 && lepID==11) continue;
   if(GENLepton_isHardProcess[iLep]==1 & GENLepton_isHardProcess[jLep]==1){
    idxHardLep1 = iLep;
    idxHardLep2 = jLep;
    nDileptons++;
   }//end if hard process
   if(GENLepton_fromHardProcessFinalState[iLep]==1 && 
    GENLepton_fromHardProcessFinalState[jLep]==1){
    idxFSRLep1 = iLep;
    idxFSRLep2 = jLep;
   }//end if FSR
  }//end jLep loop
 }//end iLep loop
 return nDileptons;
}

int DYAnalyzer::GetRecoElectrons(int &leadEle,int &subEle)
{
 int numDileptons = 0;

 for(int iEle = 0; iEle < Nelectrons; iEle++) {
  if(!Electron_passMediumID[iEle]) continue;
  for(int jEle = iEle+1; jEle < Nelectrons; jEle++) {
   if(!Electron_passMediumID[jEle]) continue;
   if(AcceptanceCut(Electron_pT[iEle],Electron_pT[jEle],Electron_eta[iEle],
    Electron_eta[jEle])){
    numDileptons++;
    if(Electron_pT[iEle]>Electron_pT[jEle]){
     leadEle = iEle;
     subEle = jEle;
    }
    else {
     leadEle = jEle;
     subEle = iEle;
    }
   }
  }//end jEle loop
 }//end iEle loop
 return numDileptons;
}

bool DYAnalyzer::AcceptanceCut(double pt1,double pt2,double eta1,double eta2)
{
  if(abs(eta1)>etaGapLow && abs(eta1)<etaGapHigh) return false;
  if(abs(eta2)>etaGapLow && abs(eta2)<etaGapHigh) return false;
  if(abs(eta1)>etaHigh||abs(eta2)>etaHigh) return false;
  if(!((pt1>ptLow && pt2>ptHigh)||(pt1>ptHigh && pt2>ptLow))) return false;
  return true;
}

bool DYAnalyzer::GenToRecoMatchCut(int genIndex,int &recoIndex)
{
 double dR,deta,dphi;
 float dRMin = 100000;
 recoIndex=-1;
 for(int iEle=0;iEle<Nelectrons;iEle++){
  deta=Electron_eta[iEle]-GENLepton_eta[genIndex];
  dphi=abs(Electron_phi[iEle]-GENLepton_phi[genIndex]);
  if(dphi>pi) dphi=2*pi-dphi;
  dR=sqrt(deta*deta+dphi*dphi);
  if(dR<dRMin){
   recoIndex=iEle;
   dRMin=dR;
  }
 }//end for loop

 bool matchFound = true;
 if(dRMin>=dRMinCut){
  recoIndex=-1;
  matchFound=false;
 }
  
 return matchFound;
}

//-----Medium ID Cuts-----//
bool DYAnalyzer::MediumIDCut(bool passID1,bool passID2)
{
 if(passID1 && passID2) return true;
 else return false;
}

//-----HLT Cuts-----//
bool DYAnalyzer::HLTCut()
{
 //pHLT_trigName is a pointer to the vector HLT_trigName contained in the tree
 int trigNameSize = pHLT_trigName->size();
 bool passHLT = false;
 for(int iHLT=0;iHLT<trigNameSize;iHLT++) {
  trigName = pHLT_trigName->at(iHLT);
  if(trigName.CompareTo(triggerUsed)==0 && HLT_trigFired[iHLT]==1){
   passHLT = true;
   break;
  }
 }
 return passHLT;
}

//-----Load files and histograms-----//
void DYAnalyzer::LoadHistograms()
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


//-----Adding weights to events-----//
double DYAnalyzer::GetTotalWeight(bool isReco,int iChain,double genWeight,double xSecWeight,
                                  double eta1,double eta2,double pt1,double pt2)
{

 //Taking care of edge cases
 if(pt1<ptBinLow)  pt1 = ptBinLow;
 if(pt2<ptBinLow)  pt2 = ptBinLow;
 if(pt1>ptBinHigh) pt1 = ptBinHigh;
 if(pt2>ptBinHigh) pt2 = ptBinHigh;
 
 //Initialize weights to 1.0
 sfWeight = 1.0;
 pileupWeight = 1.0;
 totalWeight = 1.0;

 pileupWeight = hPileupRatio->GetBinContent(hPileupRatio->FindBin(nPileUp));

 sfReco1=hRecoSF->GetBinContent(hRecoSF->FindBin(eta1,pt1));
 sfReco2=hRecoSF->GetBinContent(hRecoSF->FindBin(eta2,pt2));
 sfID1=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eta1,pt1));
 sfID2=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eta2,pt2));
 sfHLT=(hLeg2SF->GetBinContent(hLeg2SF ->FindBin(eta1,pt1)))*
              (hLeg2SF->GetBinContent(hLeg2SF ->FindBin(eta2,pt2)));
 if(isReco) sfWeight = sfReco1*sfReco2*sfID1*sfID2*sfHLT;

 totalWeight = genWeight*xSecWeight*pileupWeight*sfWeight;
 return totalWeight;
}

double DYAnalyzer::GetGenWeightSum(int iChain)
{
 TTimeStamp ts_start;
 cout << "Begin getting gen weights:" << endl;
 cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
 TStopwatch totaltime;
 totaltime.Start();

 sumGenWeight = 0.0;
 sumRawGenWeight = 0.0;
 varGenWeight = 0.0;
 Long64_t localEntry;
 for(Long64_t i=0;i<GetDYEntries(iChain);i++){
  localEntry = chains[iChain]->LoadTree(i);
  b_GENEvt_weight->GetEntry(localEntry);
  genWeight = GENEvt_weight/fabs(GENEvt_weight);//normalized genweight
  sumGenWeight += genWeight;
  varGenWeight += GENEvt_weight*GENEvt_weight;//variance of genweights
  sumRawGenWeight += GENEvt_weight; 
 }  

 totaltime.Stop();
 Double_t TotalCPURunTime = totaltime.CpuTime();
 Double_t TotalRunTime = totaltime.RealTime();
 TTimeStamp ts_end;
 cout << endl;
 cout << "End Getting Gen Weights:" << endl;
 cout << "**************************************************************************" << endl;
 cout << "Total CPU RunTime: " << TotalCPURunTime/60 << " minutes" << endl;
 cout << "Total Real RunTime: " << TotalRunTime/60 << " minutes" << endl;
 cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
 cout << "**************************************************************************" << endl;
 cout << endl;

 return sumGenWeight;
}//end GetGenWeight

double DYAnalyzer::GetXsecWeight(int iChain,bool useGenWeight)
{
 if(useGenWeight) return dataLuminosity*(xSec[iChain]/1.0);
 else return dataLuminosity*(xSec[iChain]/GetDYEntries(iChain));
}//end GetXsecWeight

//-----Other functions-----//
void DYAnalyzer::Counter(Long64_t i,Long64_t N,TString name)
{
 int P = 100*(i)/(N);
 TTimeStamp eventTimeStamp;
 if(i%(N/100)==0) {
  cout << name << ": [Time: " << eventTimeStamp.AsString("s") << "] " << P << "%" << endl;
 }
 return;
}

double DYAnalyzer::CalcInvMass(double pt1,double eta1,double phi1,double m1,double pt2,
                             double eta2,double phi2,double m2)
{
 TLorentzVector v1;
 TLorentzVector v2;
 v1.SetPtEtaPhiM(pt1,eta1,phi1,m1);
 v2.SetPtEtaPhiM(pt2,eta2,phi2,m2);
 return (v1+v2).M();
}

void DYAnalyzer::GetEfficiencies(TH1*hist0,TH1*hist1,TString name)
{
 //Having counts in the underflow bin causes issues with TEfficiency
 //So remove them first
 TString saveFileName = "data/efficiencies.root";
 hist0->SetBinContent(0,0);
 hist1->SetBinContent(0,0);

 TEfficiency*eff = new TEfficiency((*hist1),(*hist0));
 eff->SetTitle(name);
 eff->SetMarkerStyle(20);
 eff->SetMarkerSize(0.5);
 eff->SetName(name);  
 eff->SetStatisticOption(TEfficiency::kFNormal);
 
 TFile*saveFile = new TFile(saveFileName,"update");
 eff->Write();
 saveFile->Close();
}

//Define invariant mass histograms
//Default linear binning of 598 corresponds to 5 GeV per bin
TH1D*DYAnalyzer::DefineMassHist(BinType type,TString histName,int nBins = 598)
{
 float lowBin = 10;
 float highBin = 3000;
 TH1D*hist;
 if(type==LOG){ 
  hist = new TH1D(histName,"",nLogBins,massbins);
 }
 else if(type==LINEAR){
  hist = new TH1D(histName,"",nBins,lowBin,highBin);
 }
 else{
  hist = new TH1D("INVALID","",0,0,0);
  cout << "ERROR: Histogram binning not defined!!!!!!!!!!!" << endl;
 }
 return hist;
}
/*
vector<vector<double>> GetVars(LepType lepType,TChain chains)
{
 int iHard1 = -1;
 int iHard2 = -1;
 int iFSR1 = -1; 
 int iFSR2 = -1;

 double lepMass;
 double pt1,pt2,eta1,eta2,phi1,phi2,invMass;

 std::vector<std::vector<double>> vars;
 if(lepType==ELE) lepMass = eMass;
 else if(lepType==MUON) lepMass = muMass;
 for(int iChain=0;iChain<numchains;iChain++){
  for(Long64_t iEvent=0;iEvent<chains[iChain]->GetEntries();iEvent++){
   if(isMC && !isReco){  
    chains[iChain]->GetEntry(iEvent);
    pt1  = GENLepton_pT[iHard1];
    pt2  = GENLepton_pT[iHard2];
    eta1 = GENLepton_eta[iHard1];
    eta2 = GENLepton_eta[iHard2];
    phi1 = GENLepton_phi[iHard1];
    phi2 = GENLepton_phi[iHard2];

    invMass = -1;
    GetGenLeptons(lepType,iHard1,iHard2,iFSR1,iFSR2);
    invMass = CalcInvMass(pt1,eta1,phi1,lepMass,pt2,eta2,phi2,lepMass);

    vars.at(iEvent) = {pt1,pt2,eta1,eta2,phi1,phi2,invMass};
   }//end if isMC && !isReco
  }//end event loop
 }//end chain loop
}
