#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"

void counter(Long64_t i, Long64_t N);
double calcInvMass(double pt1,double eta1,double phi1,double m1,double pt2,double eta2,double phi2,double m2);
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2);
bool passPromptGenElectron(int ID, int fromfinalstate);
bool passHardProcess(int ID, int hardProces);
bool findGenToRecoMatch(int genIndex,int &recoIndex);

enum chainNum 
  {        
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
enum InvMassHist 
{
  KINEMATIC_CUTS,
  RECO_MATCHED,
  ID_CUTS,
  ALL_ELE,
  HLT_CUTS,
  RECO_ELE
};

////////////////////////////////////////////////////////////////////////
//  Here we have number of MC and Data events from the ntuples        //
//  which pass all cuts. The ratio of these is used to fill           //
//  distribution with a realistic number of events to test unfolding  //
////////////////////////////////////////////////////////////////////////
const double numMC = 21344924;
const double numData = 12549278;
const double dataMCRatio = numData/numMC;//ratio of data to MC events 

const TString treeName = "recoTree/DYTree";
const TString pileupRatioName = "/home/hep/wrtabb/git/DY-Analysis/plots/pileup.root"; 
const TString leg2SFName = "/home/hep/wrtabb/git/DY-Analysis/SFs/Leg2_SF.root"; 
const TString medIDSFName = "/home/hep/wrtabb/git/DY-Analysis/SFs/MediumID_SF.root";
const TString recoSFName = "/home/hep/wrtabb/git/DY-Analysis/SFs/Reco_SF.root";
const TString effName = "/home/hep/wrtabb/git/DY-Analysis/plots/efficiencies.root";

const int nSubSamples10to50 = 3;
const int nSubSamples100to200 = 2;
const int ptBinHigh = 499;
const int ptBinLow = 26;
const int nMatrixHistos = 3;
const int nInvMassHistos = 6;
const TString histMatrixNames[nMatrixHistos] = {"hGenHardvsGenFS","hGenFSvsReco",
 "hGenHardvsReco"};
const TString matrixYaxisTitles[nMatrixHistos] = {"Gen-Level Final State Invariant Mass [GeV]",
 "Reco Invariant Mass [GeV]","Reco Invariant Mass [GeV]"};
const TString matrixXaxisTitles[nMatrixHistos] = {"Gen-Level Hard Process Invariant Mass [GeV]" ,"Gen-Level Final State Invariant Mass [GeV]","Gen-Level Hard Process Invariant Mass [GeV]"};
const TString histInvMassNames[nInvMassHistos] = {"hGenInvMass","hGenMatchedInvMass",
 "hGenPassIDInvMass","hGenAllInvMass","hHLTGenInvMass","hRecoInvMass"};
const TString histInvMassTitles[nInvMassHistos] = {"Only Kinematic Cuts","Reco-Gen Matched",
 "Medium ID Cuts","Final State: No cuts","HLT Cut","Reconstructed"};

void makeUnfDists()
{
 TTimeStamp ts_start;
 cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
 TStopwatch totaltime;
 totaltime.Start();
 gROOT->SetBatch(kTRUE);
 gStyle->SetOptStat(0);
 gStyle->SetPalette(1);

 //Defining branches
 TBranch*b_GENnPair;
 TBranch*b_GENLepton_eta;
 TBranch*b_GENLepton_phi;
 TBranch*b_GENLepton_pT;
 TBranch*b_GENLepton_ID;
 TBranch*b_GENLepton_isHardProcess;
 TBranch*b_GENLepton_fromHardProcessFinalState;
 TBranch*b_GENEvt_weight;
 TBranch*b_Nelectrons;
 TBranch*b_Electron_pT;
 TBranch*b_Electron_eta;
 TBranch*b_Electron_phi;
 TBranch*b_Electron_passMediumID;
 TBranch*b_HLT_ntrig;
 TBranch*b_HLT_trigType;
 TBranch*b_HLT_trigFired;
 TBranch*b_nPileUp;

 //Loading ntuples
 cout << "Loading ntuples" << endl;

 TString dirNames[numChains] = {      
  "/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE/crab_DYLL_M10to50_",
  "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_truncated_M50To100/EE",
  "/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE/crab_DYLL_M100to200",
  "/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
  "/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
  "/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
  "/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
  "/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
  "/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
  "/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
  "/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE"
 }; 

 TString baseDirectory = 
   "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/ikrav/DrellYan_13TeV_2016/v2p3/DYJetsToLL_allMasses_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"; 
 TChain*chains[numChains];
 vector <TString> *subFiles[numChains]; 
 for(int iChain=0;iChain<numChains;iChain++) {
  subFiles[iChain] = new vector<TString>;
  if(iChain==MC10to50) {
   subFiles[iChain]->push_back(dirNames[iChain]+"ext1v1");
   subFiles[iChain]->push_back(dirNames[iChain]+"v1");
   subFiles[iChain]->push_back(dirNames[iChain]+"v2");
  }
  else if(iChain==MC100to200) {
   subFiles[iChain]->push_back(dirNames[iChain]);
   subFiles[iChain]->push_back(dirNames[iChain]+"_ext");
  }
  else subFiles[iChain]->push_back(dirNames[iChain]);      
 }//end iChain loop 

 TString files;  
 Long64_t subDirectorySize;
 Long64_t totalentries = 0;
 for(int iChain=0;iChain<numChains;iChain++) {
  chains[iChain] = new TChain(treeName);
  subDirectorySize = subFiles[iChain]->size();
  for(int k=0;k<subDirectorySize;k++) {	  	      
   TFileCollection filecoll("dum");//Object for creating a list of files in a directory
   files = baseDirectory;
   files+=subFiles[iChain]->at(k);
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
   }//end if
 }
  
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
  chains[iChain]->SetBranchAddress("Nelectrons", &Nelectrons, &b_Nelectrons);
  chains[iChain]->SetBranchAddress("nPileUp", &nPileUp, &b_nPileUp);
  chains[iChain]->SetBranchAddress("Electron_pT", &Electron_pT, &b_Electron_pT);
  chains[iChain]->SetBranchAddress("Electron_eta",&Electron_eta, &b_Electron_eta);
  chains[iChain]->SetBranchAddress("Electron_phi",&Electron_phi, &b_Electron_phi);
  chains[iChain]->SetBranchAddress
   ("Electron_passMediumID",&Electron_passMediumID,&b_Electron_passMediumID);
  chains[iChain]->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
  chains[iChain]->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
  chains[iChain]->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
  chains[iChain]->SetBranchAddress("HLT_trigName",&pHLT_trigName);
  
  totalentries=totalentries+chains[iChain]->GetEntries();      
 }//end loading ntuples

 cout << endl;
 cout << "Total Events Loaded: " << totalentries << endl;

 //Histograms for unfolding test
 TH1D*hReco = new TH1D("hReco","",nLogBins2,massbins2);
 TH1D*hData = new TH1D("hData","",nLogBins2,massbins2);
 TH1D*hTrue = new TH1D("hTrue","",nLogBins,massbins);
 TH2D*hMatrix = new TH2D("hMatrix","",nLogBins2,massbins2,nLogBins,massbins);
 
 //Event Loop
 cout << "Starting Event Loop" << endl;
 double invMassFSR, rapidity, dileptonPt, dileptonEta, xSecWeight, weightNoPileup, genWeight, 
  varGenWeight, totalWeight, lumiEffective, nEffective, localEntry, sumGenWeight, 
  sumRawGenWeight, pileupWeight, sfReco1, sfReco2, sfID1, sfID2, sfHLT;
 double invMassHardProcess,sfWeight;
 TRandom3*random = new TRandom3();
 double rand;
 Long64_t nentries;
 Long64_t count = 0;
 double lumi = dataLuminosity;
 TString HLTname = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
 TString trigName;
 int trigNameSize;
 long int nTooManyDielectrons = 0;
 long int nTooManyDielectronsFS = 0; 
 double eEta1, eEta2, ePt1, ePt2;
 
 TFile*pileupRatioFile  = new TFile(pileupRatioName);
 TH1F*hPileupRatio = (TH1F*)pileupRatioFile->Get("hPileupRatio");
 TFile*fileLeg2SF = new TFile(leg2SFName);
 TH2F*hLeg2SF = (TH2F*)fileLeg2SF->Get("EGamma_SF2D");
 TFile*fileMedIDSF = new TFile(medIDSFName);
 TH2F*hMedIDSF = (TH2F*)fileMedIDSF->Get("EGamma_SF2D");
 TFile*fileRecoSF = new TFile(recoSFName);
 TH2F*hRecoSF = (TH2F*)fileRecoSF->Get("EGamma_SF2D");
 TH2D*hEffSFvsMass = new TH2D("hEffSFvsMass","",nLogBins,massbins,100,0,1);
 Long64_t nDataEvents = 0;
 Long64_t nMCEvents = 0;

 for(int iChain=0;iChain<numChains;iChain++)  {
  cout << endl;
  cout << "Processing chain: " << dirNames[iChain] << endl;
  cout << endl;

  nentries = chains[iChain]->GetEntries();
  xSecWeight=lumi*(xSec[iChain]/1.0);      
  sumGenWeight = 0;
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
  pileupWeight = hPileupRatio->GetBinContent(hPileupRatio->FindBin(nPileUp));
  genWeight = GENEvt_weight/fabs(GENEvt_weight);
  genWeight = genWeight/sumGenWeight;

  for(Long64_t i=0;i<nentries;i++){      
   chains[iChain]->GetEntry(i);
   counter(count,totalentries);
   count = count+1;

   int idxGenEle1, idxGenEle2, idxGenEleFS1, idxGenEleFS2;
   idxGenEle1 = idxGenEle2 = idxGenEleFS1 = idxGenEleFS2 = -1;
   int nGenDielectrons = 0;//number of hard process dielectrons
   int nGenDielectronsFS = 0;//number of final stat dielectrons
   double invMassReco = 0;//reconstructed invariant mass
   double invMassTrue = 0;//true invariant mass

   // Loop over gen leptons and find the electron pair at the isHardProcess
   // and isHardProcessFinalState level.
   for(int kLep=0;kLep<GENnPair;kLep++){
    for(int lLep=kLep+1;lLep<GENnPair;lLep++){
     //keep only electron pairs   
     if(!(abs(GENLepton_ID[kLep])==11 && abs(GENLepton_ID[lLep])==11))
      continue;
     //require opposite sign
     if(GENLepton_ID[kLep]*GENLepton_ID[lLep]>0) continue;
     //Found dielectron pair from hard process 
     if(GENLepton_isHardProcess[kLep]==1 && GENLepton_isHardProcess[lLep]==1){
      idxGenEle1 = kLep;
      idxGenEle2 = lLep;
      nGenDielectrons++;
     }
    } // end inner loop over gen leptons
   } // end outer loop over gen leptons

   //Calculate true invariant mass from hard process
   invMassTrue = calcInvMass(GENLepton_pT[idxGenEle1],GENLepton_eta[idxGenEle1],
    GENLepton_phi[idxGenEle1],eMass,GENLepton_pT[idxGenEle2],GENLepton_eta[idxGenEle2],
    GENLepton_phi[idxGenEle2],eMass);		  
	  	  
   //Reco loop
   int idxRecoEle1,idxRecoEle2;
   idxRecoEle1=idxRecoEle2=-1;
   if(Nelectrons<2) continue;
   int nEle = 0;
   for(int iEle = 0; iEle < Nelectrons; iEle++){
    for(int jEle = iEle+1; jEle < Nelectrons; jEle++){	  
     //Find dielectron which passes kinemtatic constraints
     if(passDileptonKinematics(Electron_pT[iEle],Electron_pT[jEle],
     Electron_eta[iEle],Electron_eta[jEle])) nEle++;
     //Reco electrons which passed cuts
     if(nEle==1){
      //keeping only pairs of electrons per event
      idxRecoEle1 = iEle;
      idxRecoEle2 = jEle;		    
     }
    }//end inner reco loop	   
   }//end reco loop

   if(nEle!=1) continue;//Must have two electrons only
   if(idxRecoEle1>=0&&idxRecoEle2>=0)
    //Reconstructed invariant mass if dielectron pair found
    invMassReco=calcInvMass(Electron_pT[idxRecoEle1],Electron_eta[idxRecoEle1],
     Electron_phi[idxRecoEle1],eMass,Electron_pT[idxRecoEle2],
     Electron_eta[idxRecoEle2],Electron_phi[idxRecoEle2],eMass);	  
   if(nGenDielectrons==0) continue; // must be DY->mumu or tautau event, skip it
   if(nGenDielectrons>=2){
    // Strange, there should be only two electrons from hard process
    printf("More than two hard process dielectrons found\n");
    // skip event, but count the number of cases
    nTooManyDielectrons++;
    continue;
   }
  
  //Medium ID cuts
  if(!Electron_passMediumID[idxRecoEle1]) invMassReco = 0;//iLep electron ID cut
  if(!Electron_passMediumID[idxRecoEle2]) invMassReco = 0;//jLep electron ID cut
  //Kinematic cuts
  if(!passDileptonKinematics(Electron_pT[idxGenEle1],Electron_pT[idxGenEle2],
   Electron_eta[idxGenEle1],Electron_eta[idxGenEle2])) invMassTrue = 0; 	
   
  //Weights
  pileupWeight = hPileupRatio->GetBinContent(hPileupRatio->FindBin(nPileUp));
  genWeight = GENEvt_weight/fabs(GENEvt_weight);
  genWeight = genWeight/sumGenWeight;
  eEta1 = Electron_eta[idxRecoEle1];
  eEta2 = Electron_eta[idxRecoEle2];
  ePt1 = Electron_pT[idxRecoEle1];
  ePt2 = Electron_pT[idxRecoEle2];

  if(ePt1<ptBinLow) ePt1 = ptBinLow;//pull this information from the histograms
  if(ePt2<ptBinLow) ePt2 = ptBinLow;//raise bin
  if(ePt1>ptBinHigh) ePt1 = ptBinHigh;//lower bin
  if(ePt2>ptBinHigh) ePt2 = ptBinHigh;//

  sfReco1=hRecoSF->GetBinContent(hRecoSF->FindBin(eEta1,ePt1));
  sfReco2=hRecoSF->GetBinContent(hRecoSF->FindBin(eEta2,ePt2));
  sfID1=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eEta1,ePt1));
  sfID2=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eEta2,ePt2));
  sfHLT=(hLeg2SF->GetBinContent(hLeg2SF->FindBin(eEta1,ePt1)))*
    (hLeg2SF->GetBinContent(hLeg2SF->FindBin(eEta2,ePt2)));
  sfWeight = sfReco1*sfReco2*sfID1*sfID2*sfHLT;
  totalWeight = genWeight*xSecWeight*pileupWeight;

  const Long64_t N = nentries;

  //HLT cut
  trigNameSize = pHLT_trigName->size();
  bool passHLT = kFALSE;	  
  for(int iHLT=0;iHLT<trigNameSize;iHLT++) {
   trigName = pHLT_trigName->at(iHLT);
   if(trigName.CompareTo(HLTname)==0){
    if(HLT_trigFired[iHLT]==1){
     passHLT = true;	
    }
    else passHLT = false;
    break; 
   } 
  }// end loop over triggers

  ///////////////////////////////////////////////////////////
  // Apply matching to reconstructed electrons requirement //
  ///////////////////////////////////////////////////////////
  //int closestLep1, closestLep2;
  //closestLep1 = closestLep2 = -1;
  //bool genToRecoMatchedLep1 = findGenToRecoMatch(idxGenEleFS1,closestLep1);	      
  //if(closestLep1!=idxRecoEle1&&closestLep1!=idxRecoEle2) invMassReco = 0; 
  //bool genToRecoMatchedLep2 = findGenToRecoMatch(idxGenEleFS2,closestLep2);	      
  //if(closestLep2!=idxRecoEle1&&closestLep2!=idxRecoEle2) invMassReco = 0; 
  //if(!(genToRecoMatchedLep1 && genToRecoMatchedLep2)){
  // invMassReco = 0;
  //}

  // Apply HLT requirement
  if(!passHLT){ 
   invMassReco = 0;
  }
  rand = random->Rndm();
  if(invMassReco != 0) nMCEvents++; 
  if(rand>dataMCRatio){
   hReco->Fill(invMassReco,totalWeight*sfWeight);
   hMatrix->Fill(invMassReco,invMassTrue,totalWeight*sfWeight);
   hMatrix->Fill(0.0,invMassTrue,totalWeight*(1.0-sfWeight));
  }
  else{
   hData->Fill(invMassReco,totalWeight*sfWeight);
   hTrue->Fill(invMassTrue,totalWeight);
   if(invMassReco != 0) nDataEvents++;
  }
  }//end event loop   
 }//end chain loop 
  
 //Writing data to file
 TFile*fileUnf = new TFile("unfOut.root","recreate");
 hReco->Write();
 hData->Write();
 hTrue->Write();
 hMatrix->Write();
 fileUnf->Write();
 fileUnf->Close();
 
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

 ofstream eventFile;
 eventFile.open("nEventsNew.txt");
 eventFile << "MC Events: " << nMCEvents << ", Data Events: " << nDataEvents << endl; 
 eventFile.close();
  
}//end main function

//Counter for tracking program progress
void counter(Long64_t i, Long64_t N)
{
  int P = 100*(i)/(N);  
  TTimeStamp eventTimeStamp;
  if(i%(N/100)==0)
    {
      cout << "makeUnfDists.C " << "[Time: " << eventTimeStamp.AsString("s") << "] " << P << "%" << endl;
    }
  return;
}

//Finding correspondence between 
//reconstructed and Generated leptons
bool findGenToRecoMatch(int genIndex, int &recoIndex)
{
  double dR,deta,dphi;
  float dRMin = 100000;
  recoIndex=-1;
  //Matching gen level electrons with reconstructed electrons
  for(int iEle=0;iEle<Nelectrons;iEle++) {
   deta=Electron_eta[iEle]-GENLepton_eta[genIndex];
   dphi=abs(Electron_phi[iEle]-GENLepton_phi[genIndex]);
   if(dphi>pi) dphi=2*pi-dphi;
   dR=sqrt(deta*deta+dphi*dphi);

   if(dR<dRMin){
    dRMin=dR;
    recoIndex = iEle;
    }
  }//end of Loop 
  bool matchFound = true;
  if(dRMin>=dRMinCut) {
      matchFound=false;
      recoIndex = -1;
  }
  return matchFound;
}//end findGetToRecoMatch

//Kinematic cuts
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2)
{
  if(abs(eta1)>etaGapLow && abs(eta1)<etaGapHigh) return kFALSE;//eta cut
  if(abs(eta2)>etaGapLow && abs(eta2)<etaGapHigh) return kFALSE; //eta cut
  if(abs(eta1)>etaHigh||abs(eta2)>etaHigh) return kFALSE; //eta cut
  if(!((pt1>ptLow && pt2>ptHigh)||(pt1>ptHigh && pt2>ptLow))) return kFALSE;
  return kTRUE;
}

//Lepton ID and Final state cuts
bool passPromptGenElectron(int ID, int fromfinalstate)
{
  if(abs(ID)!=11) return kFALSE;		
  if(fromfinalstate!=1) return kFALSE;
  return kTRUE;
}

//Lepton ID and Hard Process Cuts
bool passHardProcess(int ID, int hardProcess)
{
  if(abs(ID)!=11) return kFALSE;		
  if(hardProcess!=1) return kFALSE;
  return kTRUE;
}

//Invariant mass calculator
double calcInvMass(double pt1,double eta1,double phi1,double m1,double pt2,double eta2,double phi2,double m2)
{  
  TLorentzVector vGenElectron1;
  TLorentzVector vGenElectron2;
  vGenElectron1.SetPtEtaPhiM(pt1,eta1,phi1,m1);
  vGenElectron2.SetPtEtaPhiM(pt2,eta2,phi2,m2);
  double invMass = (vGenElectron1+vGenElectron2).M();
  return invMass;
}//end calcInvMass
