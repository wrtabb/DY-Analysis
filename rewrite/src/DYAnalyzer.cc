#include "../include/DYAnalyzer.hh"
#include <TBranch.h>
#include <iostream>
#include <TStopwatch.h>
#include <TTimeStamp.h>

using namespace std;

TString dirNames[numChains] = {
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

Long64_t DYAnalyzer::LoadTrees()
{
 TTimeStamp ts_start;
 cout << "Begin loading trees:" << endl;
 cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
 TStopwatch totaltime;
 totaltime.Start();

 vector <TString> *subFiles[numChains];
 for(int iChain=0;iChain<numChains;iChain++){
  subFiles[iChain] = new vector<TString>;
  if(iChain==MC10to50){
   subFiles[iChain]->push_back(dirNames[iChain]+"ext1v1");
   subFiles[iChain]->push_back(dirNames[iChain]+"v1");
   subFiles[iChain]->push_back(dirNames[iChain]+"v2");
  }
  else if(iChain==MC100to200){
   subFiles[iChain]->push_back(dirNames[iChain]);
   subFiles[iChain]->push_back(dirNames[iChain]+"_ext");
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
   files=subFiles[iChain]->at(k);
   files+="/*.root";
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

 //--Initialize branches-----//
 InitBranches();
 //-----Open all needed files and load histograms-----//
 LoadHistograms();

 return totalentries;
}//end LoadTrees()

void DYAnalyzer::InitBranches()
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
}//end InitBranches()

//-----Return simple values-----//
Long64_t DYAnalyzer::GetDYEntries(int iChain)
{
 return chains[iChain]->GetEntries();
}

Long64_t DYAnalyzer::GetDYEntry(int iChain,Long64_t iEntry)
{
 return chains[iChain]->GetEntry(iEntry);
}

//-----Lepton selection-----//
//-----Get exactly two leptons from hard process and FSR-----//
int DYAnalyzer::GetGenLeptons(LepType lepType,int &idxHardLep1,int &idxHardLep2,
                              int &idxFSRLep1,int &idxFSRLep2)
{
 //PDG ID codes
 int lepID = 0;
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
