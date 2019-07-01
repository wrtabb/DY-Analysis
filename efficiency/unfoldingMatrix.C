#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"
#include "/home/hep/wrtabb/git/DY-Analysis/headers/ntupleSkimLocation.h"

const int ptBinHigh = 499;
const int ptBinLow = 26;
int nVertices;
void counter(Long64_t i, Long64_t N);
bool findGenToRecoMatch(int genIndex,int &recoIndex);

TLorentzVector getDielectronP4(double pt1,double eta1,double phi1,double m1,double pt2,
 double eta2,double phi2,double m2);
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2);
enum chainNum{
 MC10to50,
 MC50to100,
 MC100to200,
 MC200to400,
 MC400to500,
 MC500to700,
 MC700to800,
 MC800to1000,
 MC1000to1500,
 MC1500to2000,
 MC2000to3000
};

void unfoldingMatrix()
{
 TTimeStamp ts_start;
 cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
 TStopwatch totaltime;
 totaltime.Start();
 bool isMC; //is Monte Carlo
 gStyle->SetOptStat(0);
 //Defining branches
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
 
 //Loading ntuples
 cout << "Loading ntuples" << endl;
 //The names of every directory being loaded
 TString dirNames[numChains] = {EEM10to50,EEM50to100,EEM100to200,EEM200to400,EEM400to500,
  EEM500to700,EEM700to800,EEM800to1000,EEM1000to1500,EEM1500to2000,EEM2000to3000};


 TChain*chains[numChains];
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
   TFileCollection filecoll("dum");//Object for creating a list of files in a directory
   files=subFiles[iChain]->at(k);
   files+="/*.root";
   filecoll.Add(files);
   chains[iChain]->AddFileInfoList(filecoll.GetList());
   cout << files << endl;
   cout << chains[iChain]->GetEntries() << " events loaded" << endl;
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

  //Setting addresses for branches
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
 TH1D*hMC = new TH1D("hMC","",nLogBins2,massbins2);
 TH1D*hTrue = new TH1D("hTrue","",nLogBins,massbins);
 TH2D*hMatrix = new TH2D("hMatrix","",nLogBins,massbins,nLogBins2,massbins2);

 TFile*pileupRatioFile  = new TFile(pileupRatioName);
 TH1F*hPileupRatio = (TH1F*)pileupRatioFile->Get("hPileupRatio");
 TFile*fileLeg2SF = new TFile(leg2SFName);
 TH2F*hLeg2SF = (TH2F*)fileLeg2SF->Get("EGamma_SF2D");
 TFile*fileMedIDSF = new TFile(medIDSFName);
 TH2F*hMedIDSF = (TH2F*)fileMedIDSF->Get("EGamma_SF2D");
 TFile*fileRecoSF = new TFile(recoSFName);
 TH2F*hRecoSF = (TH2F*)fileRecoSF->Get("EGamma_SF2D");
 
 cout << "Starting Event Loop" << endl;
 double varGenWeight,lumiEffective,nEffective,localEntry,sumGenWeight,sumRawGenWeight, 
  totalWeight,sfWeight,weightNoPileup,xSecWeight,genWeight,pileupWeight,xSecWeightAlone;
 Long64_t nentries;
 Long64_t count = 0;
 TString compareHLT = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
 TString trigName;
 int trigNameSize;
 double lumi = dataLuminosity;//luminosity for xsec weighting
 double sfReco1,sfReco2,sfID1,sfID2,sfHLT;//efficiency scale factors
 double eEta1, eEta2, ePt1, ePt2;

 double binx,binWidth;   
 Long64_t nData = 0;
 //Loop over samples
 for(int iChain=0;iChain<numChains;iChain++) {

  cout << endl;
  cout << "Processing chain: " << dirNames[iChain] << endl;
  cout << endl;
  
  nentries = chains[iChain]->GetEntries();
  cout << "Number of Entries: " << nentries << endl; 
  //Finding normalized genWeights,sums,variances
  sumGenWeight = 0.0;
  sumRawGenWeight = 0.0;
  varGenWeight = 0.0;
  for(Long64_t i=0;i<nentries;i++){
   localEntry = chains[iChain]->LoadTree(i);
   b_GENEvt_weight->GetEntry(localEntry);
   genWeight = GENEvt_weight/fabs(GENEvt_weight);	//normalized genweight
   sumGenWeight += genWeight;
   varGenWeight += GENEvt_weight*GENEvt_weight; //variance of genweights
   sumRawGenWeight += GENEvt_weight; 
  }          
  nEffective = (sumRawGenWeight*sumRawGenWeight)/varGenWeight;
  lumiEffective = nEffective/xSec[iChain];
   
  //Event loop
  for(Long64_t i=0;i<nentries;i++) {      
   counter(count,totalentries);
   count = count+1; 
   chains[iChain]->GetEntry(i);
    
   TLorentzVector hardP4;
   TLorentzVector fsrP4;
   int idxGenEle1 = -1;
   int idxGenEle2 = -1;
   int idxGenEleFS1 = -1;
   int idxGenEleFS2 = -1;
   int nGenDielectrons = 0;
   double invMassHard = 0;
   double invMassFS = 0;
  
   //Gen loop
   for(int kLep=0;kLep<GENnPair;kLep++){
    for(int lLep=kLep+1;lLep<GENnPair;lLep++){
     if(!(abs(GENLepton_ID[kLep])==11 && abs(GENLepton_ID[lLep])==11))
     continue;
     if(GENLepton_ID[kLep]*GENLepton_ID[lLep]>0) continue;
     if(GENLepton_isHardProcess[kLep]==1 && GENLepton_isHardProcess[lLep]==1){
      idxGenEle1 = kLep;
      idxGenEle2 = lLep;
      nGenDielectrons++;
     }//end if hard process
     if(GENLepton_fromHardProcessFinalState[kLep]==1 && 
      GENLepton_fromHardProcessFinalState[lLep]==1){
      idxGenEleFS1 = kLep;
      idxGenEleFS2 = lLep;
     }
    }//end lLep loop
   }//end kLep loop
   if(nGenDielectrons!=1){
    cout << "!!!!!!!!!!WARNING!!!!!!!!!!" << endl;
    cout << "Gen level produces too many or too few electron pairs" << endl;
    continue;
   }
   bool passAcceptance = true;
   //Make sure all events are within acceptance region
   if(!passDileptonKinematics(GENLepton_pT[idxGenEle1],GENLepton_pT[idxGenEle2],
    GENLepton_eta[idxGenEle1],GENLepton_eta[idxGenEle2])) passAcceptance = false;

   //Get 4-momentum for gen-level electron pairs
   hardP4 = getDielectronP4(GENLepton_pT[idxGenEle1],GENLepton_eta[idxGenEle1],
    GENLepton_phi[idxGenEle1],eMass,GENLepton_pT[idxGenEle2],GENLepton_eta[idxGenEle2],            GENLepton_phi[idxGenEle2],eMass);
   fsrP4 = getDielectronP4(GENLepton_pT[idxGenEleFS1],GENLepton_eta[idxGenEleFS1],
    GENLepton_phi[idxGenEleFS1],eMass,GENLepton_pT[idxGenEleFS2],GENLepton_eta[idxGenEleFS2],
    GENLepton_phi[idxGenEleFS2],eMass);

   //calculate gen-level invariant masses
   if(passAcceptance) invMassHard = hardP4.M();
   invMassFS = fsrP4.M();

   //HLT cut
   trigNameSize = pHLT_trigName->size();
   bool passHLT = kFALSE;	  
   for(int iHLT=0;iHLT<trigNameSize;iHLT++) {
    trigName = pHLT_trigName->at(iHLT);	  
    if(trigName.CompareTo(compareHLT)==0) {
     if(HLT_trigFired[iHLT]==1) passHLT = kTRUE;
     else passHLT = kFALSE;		     
     break; 
    }
   }//end loop over triggers 
  
   int numDielectrons = 0;
   int subEle = -1;
   int leadEle = -1;
   double invMass = 0;
   TLorentzVector dielectronP4;

   //Reco Electron loop
   for(int iEle = 0; iEle < Nelectrons; iEle++) {
    if(!Electron_passMediumID[iEle]) continue;
    for(int jEle = iEle+1; jEle < Nelectrons; jEle++) {
     if(!Electron_passMediumID[jEle]) continue;
     if(passDileptonKinematics(Electron_pT[iEle],Electron_pT[jEle],Electron_eta[iEle],
      Electron_eta[jEle])){
      numDielectrons++;
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
    
   dielectronP4 = getDielectronP4(Electron_pT[leadEle],Electron_eta[leadEle],
    Electron_phi[leadEle],eMass,Electron_pT[subEle],Electron_eta[subEle],
    Electron_phi[subEle],eMass);

   //Calculate reco invariant mass
   invMass = dielectronP4.M();
   int closestTrackLep1, closestTrackLep2;
   closestTrackLep1 = closestTrackLep2 = -1;
   bool genToRecoMatchedLep1 = findGenToRecoMatch(idxGenEleFS1,closestTrackLep1);
   bool genToRecoMatchedLep2 = findGenToRecoMatch(idxGenEleFS2,closestTrackLep2);

   //All cuts
   if(!(genToRecoMatchedLep1 && genToRecoMatchedLep2)) invMass=0;
   if(!Electron_passMediumID[closestTrackLep1]) invMass=0;
   if(!Electron_passMediumID[closestTrackLep2]) invMass=0;
   if(numDielectrons!=1) invMass = 0;
   if(leadEle<0||subEle<0) invMass = 0;
   if(!passHLT) invMass = 0;

   eEta1 = Electron_eta[leadEle];
   eEta2 = Electron_eta[subEle];
   ePt1 = Electron_pT[leadEle];
   ePt2 = Electron_pT[subEle];
  
   if(ePt1<ptBinLow) ePt1 = ptBinLow;
   if(ePt2<ptBinLow) ePt2 = ptBinLow;
   if(ePt1>ptBinHigh) ePt1 = ptBinHigh;
   if(ePt2>ptBinHigh) ePt2 = ptBinHigh;

   //Determining weighting factors
   xSecWeight=lumi*(xSec[iChain]/1.0);//xSecWeight when used with genWeight 
   xSecWeightAlone = lumi*(xSec[iChain]/nentries);//xSecWeight when used without genWeight
   pileupWeight = hPileupRatio->GetBinContent(hPileupRatio->FindBin(nPileUp));
   genWeight = (GENEvt_weight/fabs(GENEvt_weight))/sumGenWeight;

   sfReco1=hRecoSF->GetBinContent(hRecoSF->FindBin(eEta1,ePt1));
   sfReco2=hRecoSF->GetBinContent(hRecoSF->FindBin(eEta2,ePt2));
   sfID1=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eEta1,ePt1));
   sfID2=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eEta2,ePt2));
   sfHLT=(hLeg2SF->GetBinContent(hLeg2SF->FindBin(eEta1,ePt1)))*
     (hLeg2SF->GetBinContent(hLeg2SF->FindBin(eEta2,ePt2)));
   sfWeight = sfReco1*sfReco2*sfID1*sfID2*sfHLT;
   totalWeight = genWeight*xSecWeight*pileupWeight;

   //Fill the histograms
   hMC->Fill(invMass,totalWeight*sfWeight);
   hTrue->Fill(invMassHard,totalWeight);
   hMatrix->Fill(invMassHard,invMass,totalWeight*sfWeight);
   hMatrix->Fill(invMassHard,0.0,totalWeight*(1-sfWeight));
  }//end event loop   
 }//end chain loop 
  
  TFile *rootFile = new TFile("/home/hep/wrtabb/git/DY-Analysis/data/unfoldIn.root","RECREATE");
  rootFile->cd();
  hMC->Write();
  hTrue->Write();
  hMatrix->Write();
  rootFile->Write();
  rootFile->Close();
  
  totaltime.Stop();
  Double_t TotalCPURunTime = totaltime.CpuTime();
  Double_t TotalRunTime = totaltime.RealTime();
  TTimeStamp ts_end;
  cout << endl;
  cout << "**************************************************************************" << endl;
  cout << "Total CPU RunTime: " << TotalCPURunTime/60 << " minutes" << endl;
  cout << "Total Real RunTime: " << TotalRunTime/60 << " minutes" << endl;
  cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;   
  cout << "Number of Events Processed: " << count << endl;
  cout << "Number of data events Processed: " << nData << endl;
  cout << "**************************************************************************" << endl;
  cout << endl;
  
}//end main function

//Counter for tracking program progress
void counter(Long64_t i, Long64_t N)
{
  int P = 100*(i)/(N);  
  TTimeStamp eventTimeStamp;
  if(i%(N/100)==0) {
    cout << "dataVsMC.C " << "[Time: " << eventTimeStamp.AsString("s") << "] " << P 
      << "%" << endl;
  }
  return;
}

//Kinematic cuts
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2)
{
  if(abs(eta1)>etaGapLow && abs(eta1)<etaGapHigh) return kFALSE;
  if(abs(eta2)>etaGapLow && abs(eta2)<etaGapHigh) return kFALSE; 
  if(abs(eta1)>etaHigh||abs(eta2)>etaHigh) return kFALSE; 
  if(!((pt1>ptLow && pt2>ptHigh)||(pt1>ptHigh && pt2>ptLow))) return kFALSE;
  return kTRUE;
}

TLorentzVector getDielectronP4(double pt1,double eta1,double phi1,double m1,double pt2,
  double eta2,double phi2,double m2)
{
  TLorentzVector vElectron1;
  TLorentzVector vElectron2;
  vElectron1.SetPtEtaPhiM(pt1,eta1,phi1,m1);
  vElectron2.SetPtEtaPhiM(pt2,eta2,phi2,m2);
  return vElectron1+vElectron2;
}

bool findGenToRecoMatch(int genIndex, int &recoIndex)
{
  double dR,deta,dphi;
  float dRMin = 100000;
  recoIndex=-1;
  for(int iEle=0;iEle<Nelectrons;iEle++)
    {
      deta=Electron_eta[iEle]-GENLepton_eta[genIndex];
      dphi=abs(Electron_phi[iEle]-GENLepton_phi[genIndex]);
      if(dphi>pi) dphi=2*pi-dphi;
      dR=sqrt(deta*deta+dphi*dphi);

      if(dR<dRMin)
        { 
          recoIndex=iEle;
          dRMin=dR;
        }
    }
  bool matchFound = kTRUE;
  if(dRMin>=dRMinCut)
    {
      recoIndex=-1;
      matchFound=kFALSE;
    }
  return matchFound;
}
