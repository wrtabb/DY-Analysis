#include "VariableList.h"

Long64_t LoadTrees(std::vector<TString>dirNames,SampleType sampleType,LepType lepType);
void InitializeBranches(TChain*chain,bool isMC,LepType lepType);
void EventLoop(TChain*chain,SampleType sampleType,LepType lepType,double genWeightSum,
               double xSecWeight);
int GetGenLeptons(LepType lepType,int &idxHardLep1,int &idxHardLep2,int &idxFSRLep1,
                  int &idxFSRLep2);
int GetRecoLeptons();
double CalcInvMass(double pt1,double eta1,double phi1,double m1,double pt2,double eta2,
                   double phi2,double m2);
bool PassAcceptance(double pt1,double pt2,double eta1,double eta2);
bool GenToRecoMatch(int genIndex,int &recoIndex);
bool PassMediumID(bool passID1,bool passID2);
bool HLTCut();
double GetGenWeightSum(TChain*chain);
double GetTotalWeight(bool isReco,double genWeight,double xSecWeight,
                      double eta1,double eta2,double pt1,double pt2);
double GetXsecWeight(int iChain,bool useGenWeight);
TH1D*DefineMassHist(BinType type,TString histName,int nBins);
TH2D*DefineMatrixHist(BinType type,TString histName,int nBins);

//And which lepton is being analyzed
//And which samples to load
//For samples, choose one of these:
// LL		: Leptons
// EW		: electroweak (background)
// TT		: tops (background)
// FAKES 	: Fakes, W+Jets (background)
// DATA		: data
void GetHists(SampleType sampleType,LepType lepType)
{
 std::vector<TString> dirNames;
 if(sampleType==LL){
  dirNames = dirNamesLL;
  xSec = xSecLL;
 }
 else if(sampleType==EW){
  dirNames = dirNamesEW;
  xSec = xSecEW;
 }
 else if(sampleType==TT){
  dirNames = dirNamesTT;
  xSec = xSecTT;
 }
 else if(sampleType==FAKES){
  dirNames = dirNamesFakes;
  xSec = xSecFakes;
 }
 else if(sampleType==DATA){
  dirNames = dirNamesData;
  xSec = xSecData;
 }
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
 //Temporarily stop loading skims and instead load full trees
 //The skims are missing important GEN-level branches and need to be fixed
 //fileNames = "/skims_0002/*.root";
 fileNames = "/*.root";
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
 double genWeightSum = GetGenWeightSum(chains[iChain]);
 bool useGenWeights = true;
 double xSecWeight = GetXsecWeight(iChain,useGenWeights); 
 EventLoop(chains[iChain],sampleType,lepType,genWeightSum,xSecWeight);
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
 cout << "---------------------" << endl;
 cout << "Initializing Branches" << endl;

 //-----HLT Branches-----//
 chain->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
 chain->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
 chain->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
 chain->SetBranchAddress("HLT_trigName",&pHLT_trigName);

 chain->SetBranchAddress("nVertices", &nVertices, &b_nVertices);
 chain->SetBranchAddress("nPileUp", &nPileUp, &b_nPileUp);

 //-----Reco-level branches-----//
 if(lepType==ELE){
  cout << "Initializing RECO branches" << endl;
  chain->SetBranchAddress("Nelectrons", &Nelectrons, &b_Nelectrons);
  chain->SetBranchAddress("Electron_pT", &Electron_pT, &b_Electron_pT);
  chain->SetBranchAddress("Electron_eta",&Electron_eta, &b_Electron_eta);
  chain->SetBranchAddress("Electron_phi",&Electron_phi, &b_Electron_phi);
  chain->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,
                                    &b_Electron_passMediumID);
 }
 else if(lepType==MUON){
  cout << "The muon branches for reco-muons do not exist in current ntuple skims" << endl;
  cout << "This will be fixed as soon as possible." << endl;
  return ;
 }

 //-----Gen-level branches-----//
 if(isMC){
  cout << "Initializing GEN branches" << endl;
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
 cout << "---------------------" << endl;
 return;
}//end Initialize Branches

void EventLoop(TChain*chain,SampleType sampleType,LepType lepType,double genWeightSum,
               double xSecWeight)
{
 Long64_t nEntries = chain->GetEntries();
 double lepMass;
 if     (lepType == ELE)  lepMass = eMass;
 else if(lepType == MUON) lepMass = muMass;
 else if (lepType == TAU) lepMass = tauMass; 
 else {
  lepMass = -1;
  cout << "Lepton type not properly selected!" << endl;
  cout << "Must be 'ELE', 'MUON', or 'TAU'" << endl;
  return;
 }

 int nHists = 5;
 TH1D*hMass[nHists];
 TH2D*hMatrix[nHists-1];
 for(int i=0;i<nHists;i++){
  TString histNameMass = "hMass";
  TString histNameMatrix = "hMatrix";
  histNameMass += i;
  histNameMatrix += i;
  hMass[i] = DefineMassHist(LOG,histNameMass,0);
  if(i<nHists-1) hMatrix[i] = DefineMatrixHist(LOG,histNameMatrix,0);
 }

 for(Long64_t i=0;i<nEntries;i++){
  Long64_t event = chain->GetEntry(i);
  double massHard = -1;
  double massFSR = -1;
  double massReco = -1;
  int iHard1 = -1;
  int iHard2 = -1;
  int iFSR1 = -1;
  int iFSR2 = -1;
  int leadLep = -1;
  int subLep = -1;
  int closestTrack1 = -1;
  int closestTrack2 = -1;
  int nGenDileptons;
  int nRecoDileptons;
 
  //Select leptons in each event
  nGenDileptons = GetGenLeptons(lepType,iHard1,iHard2,iFSR1,iFSR2);
  //nRecoDileptons = GetRecoLeptons(lepType,leadLep,subLep);
 
  //weights
  double weightReco = 1.0;
  double weightHard = 1.0;
  double genWeight = (GENEvt_weight/fabs(GENEvt_weight))/genWeightSum;
  double xSecWeight = 1.0;
  double pt1 = Electron_pT[leadLep];
  double pt2 = Electron_pT[subLep];
  double eta1 = Electron_eta[leadLep];
  double eta2 = Electron_eta[subLep];

  weightReco = GetTotalWeight(true,genWeight,xSecWeight,eta1,eta2,pt1,pt2);
  weightHard = GetTotalWeight(false,genWeight,xSecWeight,eta1,eta2,pt1,pt2);
  //Determine which leptons pass cuts
  bool passAcceptance = PassAcceptance(GENLepton_pT[iHard1],GENLepton_pT[iHard2],
                                       GENLepton_eta[iHard1],GENLepton_eta[iHard2]);
  bool passRecoMatch1 = GenToRecoMatch(iFSR1,closestTrack1);
  bool passRecoMatch2 = GenToRecoMatch(iFSR2,closestTrack2);
  bool passRecoMatch = passRecoMatch1 && passRecoMatch2;
  bool passMediumID = PassMediumID(Electron_passMediumID[closestTrack1],
                                  Electron_passMediumID[closestTrack2]);
  bool passHLT = HLTCut();

  if(iHard1>=0 && iHard2>=0){
   massHard = CalcInvMass(GENLepton_pT[iHard1],GENLepton_eta[iHard1],
                          GENLepton_phi[iHard1],lepMass,GENLepton_pT[iHard2],
                          GENLepton_eta[iHard2],GENLepton_phi[iHard2],lepMass);
  }
  if(iFSR1>=0 && iFSR2>=0){
   massFSR = CalcInvMass(GENLepton_pT[iFSR1],GENLepton_eta[iFSR1],
                         GENLepton_phi[iFSR1],lepMass,GENLepton_pT[iFSR2],
                         GENLepton_eta[iFSR2],GENLepton_phi[iFSR2],lepMass);
  }
  if(leadLep>=0 && subLep>=0){
   massReco = CalcInvMass(Electron_pT[leadLep],Electron_eta[leadLep],
                          Electron_phi[leadLep],lepMass,Electron_pT[subLep],
                          Electron_eta[subLep],Electron_phi[subLep],lepMass);
  }
  
  hMass[0]->Fill(massHard,weightHard);
  if(!passAcceptance) massReco = 0;
  hMass[1]->Fill(massReco,weightReco);
  hMatrix[0]->Fill(massReco,massHard);
  if(!passRecoMatch) massReco = 0;
  hMass[2]->Fill(massReco,weightReco);
  hMatrix[1]->Fill(massReco,massHard);
  if(!passMediumID) massReco = 0;
  hMass[3]->Fill(massReco,weightReco);
  hMatrix[2]->Fill(massReco,massHard);
  if(!passHLT) massReco = 0;
  hMass[4]->Fill(massReco,weightReco);
  hMatrix[3]->Fill(massReco,massHard);

  TFile*histSave = new TFile("unfoldingHists.root","recreate");
  for(int i=0;i<nHists;i++){
   hMass[i]->Write();
   if(i<nHists-1) hMatrix[i]->Write();
  }
  histSave->Write();
  histSave->Close();
 }//end event loop
}//end EventLoop()

int GetGenLeptons(LepType lepType,int &idxHardLep1,int &idxHardLep2,int &idxFSRLep1,
                  int &idxFSRLep2)
{
 int lepID;
 if      (lepType==ELE)  lepID = 11;
 else if (lepType==MUON) lepID = 13;
 else if (lepType==TAU)  lepID = 15;
 else {
  cout << "ERROR: Appropriate lepton not selected" << endl;
  return 0;
 }
 int nDileptons = 0;

 for(int iLep=0;iLep<GENnPair;iLep++){
  for(int jLep=iLep+1;jLep<GENnPair;jLep++){
   if(!(abs(GENLepton_ID[iLep])==lepID && abs(GENLepton_ID[jLep])==lepID)) continue;
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

double CalcInvMass(double pt1,double eta1,double phi1,double m1,double pt2,double eta2,
                   double phi2,double m2)
{
 TLorentzVector v1;
 TLorentzVector v2;
 v1.SetPtEtaPhiM(pt1,eta1,phi1,m1);
 v2.SetPtEtaPhiM(pt2,eta2,phi2,m2);
 return (v1+v2).M();
}

bool PassAcceptance(double pt1,double pt2,double eta1,double eta2)
{
 if(abs(eta1)>etaGapLow && abs(eta1)<etaGapHigh) return false;
 if(abs(eta2)>etaGapLow && abs(eta2)<etaGapHigh) return false;
 if(abs(eta1)>etaHigh||abs(eta2)>etaHigh) return false;
 if(!((pt1>ptLow && pt2>ptHigh)||(pt1>ptHigh && pt2>ptLow))) return false;
 return true;
} 

bool GenToRecoMatch(int genIndex,int &recoIndex)
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

bool PassMediumID(bool passID1,bool passID2)
{
 if(passID1 && passID2) return true;
 else return false;
}

bool HLTCut()
{
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

double GetGenWeightSum(TChain*chain)
{
 TTimeStamp ts_start;
 cout << "Begin getting gen weights:" << endl;
 cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
 TStopwatch totaltime;
 totaltime.Start();

 double genWeight;
 double sumGenWeight = 0.0;
 double sumRawGenWeight = 0.0;
 double varGenWeight = 0.0;
 Long64_t localEntry;
 for(Long64_t i=0;i<chain->GetEntries();i++){
  localEntry = chain->LoadTree(i);
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

double GetTotalWeight(bool isReco,double genWeight,double xSecWeight,
                      double eta1,double eta2,double pt1,double pt2)
{
 if(pt1<ptBinLow)  pt1 = ptBinLow;
 if(pt2<ptBinLow)  pt2 = ptBinLow;
 if(pt1>ptBinHigh) pt1 = ptBinHigh;
 if(pt2>ptBinHigh) pt2 = ptBinHigh;

 TFile*pileupRatioFile  = new TFile(pileupRatioName);
 TH1F*hPileupRatio = (TH1F*)pileupRatioFile->Get("hPileupRatio");
 TFile*fileLeg2SF  = new TFile(leg2SFName);
 TFile*fileMedIDSF = new TFile(medIDSFName);
 TFile*fileRecoSF  = new TFile(recoSFName);
 TH2F*hLeg2SF  = (TH2F*) fileLeg2SF->Get("EGamma_SF2D");
 TH2F*hMedIDSF = (TH2F*)fileMedIDSF->Get("EGamma_SF2D");
 TH2F*hRecoSF  = (TH2F*) fileRecoSF->Get("EGamma_SF2D");
 double sfWeight = 1.0;
 double pileupWeight = 1.0;
 double totalWeight = 1.0;

 pileupWeight = hPileupRatio->GetBinContent(hPileupRatio->FindBin(nPileUp));

 double sfReco1=hRecoSF->GetBinContent(hRecoSF->FindBin(eta1,pt1));
 double sfReco2=hRecoSF->GetBinContent(hRecoSF->FindBin(eta2,pt2));
 double sfID1=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eta1,pt1));
 double sfID2=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eta2,pt2));
 double sfHLT=(hLeg2SF->GetBinContent(hLeg2SF ->FindBin(eta1,pt1)))*
              (hLeg2SF->GetBinContent(hLeg2SF ->FindBin(eta2,pt2)));
 if(isReco) sfWeight = sfReco1*sfReco2*sfID1*sfID2*sfHLT;

 totalWeight = genWeight*xSecWeight*pileupWeight*sfWeight;
 return totalWeight;
}

double GetXsecWeight(int iChain,TChain*chain,bool useGenWeight)
{
 if(useGenWeight) return dataLuminosity*(xSec.at(iChain)/1.0);
 else return dataLuminosity*(xSec.at(iChain)/chain->GetEntries());
}//end GetXsecWeight

TH1D*DefineMassHist(BinType type,TString histName,int nBins)
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

TH2D*DefineMatrixHist(BinType type,TString histName,int nBins)
{
 float lowBin = 10;
 float highBin = 3000;
 TH2D*hist;
 if(type==LOG){
  hist = new TH2D(histName,"",nLogBins,massbins,nLogBins2,massbins2);
 }
 else if(type==LINEAR){
  hist = new TH2D(histName,"",nBins,lowBin,highBin,nBins,lowBin,highBin);
 }
 else{
  hist = new TH2D("INVALID","",0,0,0,0,0,0);
  cout << "ERROR: Histogram binning not defined!!!!!!!!!!!" << endl;
 }
 return hist;
}
