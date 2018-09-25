#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "THStack.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TLegend.h"
#include <fstream>
#include <iostream>
#include <vector>
#include "TStyle.h"
#include "TEfficiency.h"
#include "TString.h"
#include "TLine.h"
#include "TTimeStamp.h"
#include "TFileCollection.h"
#include "THashList.h"
#include "TGraphAsymmErrors.h"

void counter(Long64_t i, Long64_t N);
double calcInvMass(double pt1,double eta1,double phi1,double m1,double pt2,double eta2,double phi2,double m2);
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2);
bool passPromptGenElectron(int ID, int fromfinalstate);
bool passHardProcess(int ID, int hardProces);
bool findGenToRecoMatch(int genIndex,int &recoIndex);

//Defining variables and arrays
const int MPSIZE = 2000;
int GENnPair, Nelectrons, HLT_ntrig, nPileUp;
double GENEvt_weight;
double GENLepton_phi[MPSIZE],GENLepton_eta[MPSIZE],GENLepton_pT[MPSIZE],GENLepton_Px[MPSIZE],GENLepton_Py[MPSIZE];
double GENLepton_Pz[MPSIZE],GENLepton_E[MPSIZE];
int GENLepton_ID[MPSIZE],GENLepton_isHardProcess[MPSIZE],GENLepton_fromHardProcessFinalState[MPSIZE];
double Electron_pT[MPSIZE], Electron_eta[MPSIZE], Electron_phi[MPSIZE];
double Electron_Energy[MPSIZE], Electron_Px[MPSIZE];
double Electron_Py[MPSIZE], Electron_Pz[MPSIZE], Electron_charge[MPSIZE];
bool Electron_passMediumID[MPSIZE];
int HLT_trigType[MPSIZE],HLT_trigFired[MPSIZE];
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
std::vector<std::string> HLT_trigName;
std::vector<std::string> *pHLT_trigName = &HLT_trigName;
const double massbins[44] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106, 
                             110,115,120,126,133,141,150,160,171,185,200,220,243,273,320, 
			     380,440,510,600,700,830,1000,1500,3000};
const double pi=TMath::Pi();
const int numChains = 11;
const int nLogBins = 43;
//Cross sections obtained from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SNUCMSYooDYntuple
const float xSec[numChains] = {6016.88,1873.52,76.2401,2.67606,0.139728,0.0792496,0.0123176,
                               0.01042,0.00552772,0.000741613,0.000178737};
const float etaHigh = 2.4;
const float etaGapHigh = 1.566; 
const float etaGapLow = 1.4442;
const float ptHigh = 28;
const float ptLow = 17;
const float eMass = 0.000511;
const TString treeName = "recoTree/DYTree";
const TString pileupRatioName = "./plots/pileup.root";
const TString leg2SFName = "./data/Leg2_SF.root";
const TString medIDSFName = "./data/MediumID_SF.root";
const TString recoSFName = "./data/Reco_SF.root";
const float dRMinCut = 0.3;
const int dataLuminosity = 35867; //Run2016B to Run2016H JSON. unit: /pb, Updated at 2017.07.30

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
const TString matrixXaxisTitles[nMatrixHistos] = {"Gen-Level Hard Process Invariant Mass [GeV]"  ,"Gen-Level Final State Invariant Mass [GeV]","Gen-Level Hard Process Invariant Mass [GeV]"};
const TString histInvMassNames[nInvMassHistos] = {"hGenInvMass","hGenMatchedInvMass",
  "hGenPassIDInvMass","hGenAllInvMass","hHLTGenInvMass","hRecoInvMass"};
const TString histInvMassTitles[nInvMassHistos] = {"Only Kinematic Cuts","Reco-Gen Matched",
  "Medium ID Cuts","Final State: No cuts","HLT Cut","Reconstructed"};

void efficiencies()
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

    TString dirNames[numChains] = 
      {      
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
 for(int iChain=0;iChain<numChains;iChain++)
   {
     subFiles[iChain] = new vector<TString>;
     if(iChain==MC10to50) 
       {
	 subFiles[iChain]->push_back(dirNames[iChain]+"ext1v1");
	 subFiles[iChain]->push_back(dirNames[iChain]+"v1");
	 subFiles[iChain]->push_back(dirNames[iChain]+"v2");
       }
     else if(iChain==MC100to200) 
       {
	 subFiles[iChain]->push_back(dirNames[iChain]);
	 subFiles[iChain]->push_back(dirNames[iChain]+"_ext");
       }
     else subFiles[iChain]->push_back(dirNames[iChain]);      
   } 
 TString files;  
 Long64_t subDirectorySize;
 Long64_t totalentries = 0;
 for(int iChain=0;iChain<numChains;iChain++)
   {
     chains[iChain] = new TChain(treeName);
     subDirectorySize = subFiles[iChain]->size();
      for(int k=0;k<subDirectorySize;k++)
	{	  	      
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
	  }
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

 TH1F*histInvMass[nInvMassHistos];
 for(int i=0;i<nInvMassHistos;i++){
   histInvMass[i]=new TH1F(histInvMassNames[i],"",nLogBins,massbins);
   histInvMass[i]->Sumw2();
   histInvMass[i]->GetXaxis()->SetTitle("invariant mass [GeV]");
   histInvMass[i]->GetXaxis()->SetMoreLogLabels();
   histInvMass[i]->GetXaxis()->SetNoExponent();
   histInvMass[i]->SetTitle(histInvMassTitles[i]);
 }
 
 TH2F*hpTvsMass = new TH2F("hpTvsMass","",43,massbins,598,10,3000);
 hpTvsMass->GetXaxis()->SetMoreLogLabels();
 hpTvsMass->GetXaxis()->SetNoExponent();
 hpTvsMass->GetYaxis()->SetTitle("p_{T} [GeV]"); 
 hpTvsMass->GetXaxis()->SetTitle("m_{ee} [GeV]"); 

 TH2F*hMatrix[nMatrixHistos];
 for(int i=0;i<nMatrixHistos;i++){
   hMatrix[i]=new TH2F(histMatrixNames[i],"",nLogBins,massbins,nLogBins,massbins);
   hMatrix[i]->GetYaxis()->SetTitle(matrixYaxisTitles[i]);
   hMatrix[i]->GetXaxis()->SetTitle(matrixXaxisTitles[i]);
   hMatrix[i]->GetYaxis()->SetNoExponent();
   hMatrix[i]->GetYaxis()->SetMoreLogLabels();
   hMatrix[i]->GetXaxis()->SetNoExponent();
   hMatrix[i]->GetXaxis()->SetMoreLogLabels();
 } 

 TH2F*migMatrixGENisHardvsGENFS = 
   new TH2F("migMatrixGENisHardvsGENFS","",nLogBins,massbins,nLogBins,massbins);
 migMatrixGENisHardvsGENFS->
   SetTitle("Migration Matrix: Gen-Level Final State vs. Gen-Level Hard Process");
 migMatrixGENisHardvsGENFS->GetYaxis()->
   SetTitle("Gen-Level Final State Dielectron Invariant Mass [GeV]");
 migMatrixGENisHardvsGENFS->GetXaxis()->
   SetTitle("Gen-Level Hard Process Dielecron Invariant mass [GeV]");
 migMatrixGENisHardvsGENFS->GetXaxis()->SetNoExponent();
 migMatrixGENisHardvsGENFS->GetXaxis()->SetMoreLogLabels();
 migMatrixGENisHardvsGENFS->GetYaxis()->SetNoExponent();
 migMatrixGENisHardvsGENFS->GetYaxis()->SetMoreLogLabels();
 TH2F*migMatrixGENFSvsReco = 
   new TH2F("migMatrixGENFSvsReco","",nLogBins,massbins,nLogBins,massbins);
 migMatrixGENFSvsReco->SetTitle("Migration Matrix: Reconstructed vs. Gen-Level Final State");
 migMatrixGENFSvsReco->GetXaxis()->
   SetTitle("Gen-Level Final State Dielectron Invariant Mass [GeV]");
 migMatrixGENFSvsReco->GetYaxis()->SetTitle("Reconstructed Dielecron Invariant mass [GeV]");
 migMatrixGENFSvsReco->GetXaxis()->SetNoExponent();
 migMatrixGENFSvsReco->GetXaxis()->SetMoreLogLabels();
 migMatrixGENFSvsReco->GetYaxis()->SetNoExponent();
 migMatrixGENFSvsReco->GetYaxis()->SetMoreLogLabels();
 TH2F*migMatrixGENisHardvsReco = 
   new TH2F("migMatrixGENisHardvsReco","",nLogBins,massbins,nLogBins,massbins);
 migMatrixGENisHardvsReco->
   SetTitle("Migration Matrix: Reconstructed vs. Gen-Level Hard Process");
 migMatrixGENisHardvsReco->GetXaxis()->
   SetTitle("Gen-Level Hard Process Dielectron Invariant Mass [GeV]");
 migMatrixGENisHardvsReco->GetYaxis()->
   SetTitle("Reconstructed Dielecron Invariant mass [GeV]");
 migMatrixGENisHardvsReco->GetXaxis()->SetNoExponent();
 migMatrixGENisHardvsReco->GetXaxis()->SetMoreLogLabels();
 migMatrixGENisHardvsReco->GetYaxis()->SetNoExponent();
 migMatrixGENisHardvsReco->GetYaxis()->SetMoreLogLabels();
 
 TH1F*hHardProcess[numChains];
 TString histbasename = "hHardProcess";
 TString histname;
 for(int jChain=0;jChain<numChains;jChain++)
   {
     histname = histbasename;
     histname+=jChain;
     hHardProcess[jChain] = new TH1F(histname,"",598,10,3000);
     hHardProcess[jChain]->SetFillColor(jChain+1);
     hHardProcess[jChain]->GetXaxis()->
       SetTitle("Gen-Level Dielectron Mass (isHardProcess) [GeV]"); 
     hHardProcess[jChain]->GetYaxis()->SetRangeUser(0.000001,1000000000);
     hHardProcess[jChain]->GetXaxis()->SetNoExponent();
     hHardProcess[jChain]->GetXaxis()->SetMoreLogLabels();
   }
 
 TFile *rootFile = new TFile("./plots/plotsDY.root","RECREATE");
 
 //Event Loop
 cout << "Starting Event Loop" << endl;
 double invMassFSR, rapidity, dileptonPt, dileptonEta, xSecWeight, weightNoPileup, genWeight, 
   varGenWeight, totalWeight, lumiEffective, nEffective, localEntry, sumGenWeight, 
   sumRawGenWeight, pileupWeight, sfReco1, sfReco2, sfID1, sfID2, sfHLT;
 double invMassHardProcess,sfWeight;
 Long64_t nentries;
 Long64_t count = 0;
 double nEvents = 250000;
 double lumi = dataLuminosity;
 //double lumi = nEvents/xSec[MC50to100];//luminosity of 50to100
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

 for(int iChain=0;iChain<numChains;iChain++)
   {
     cout << endl;
     cout << "Processing chain: " << dirNames[iChain] << endl;
     cout << endl;

      nentries = chains[iChain]->GetEntries();
      //nentries = nEvents;
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

      for(Long64_t i=0;i<nentries;i++)
	{      
	  chains[iChain]->GetEntry(i);
	  counter(count,totalentries);
	  count = count+1;
	  
	  pileupWeight = hPileupRatio->GetBinContent(hPileupRatio->FindBin(nPileUp));
	  genWeight = GENEvt_weight/fabs(GENEvt_weight);
	  genWeight = genWeight/sumGenWeight;
	  totalWeight = genWeight*xSecWeight*pileupWeight;

	  // Loop over gen leptons and find the electron pair at the isHardProcess
	  // and isHardProcessFinalState level.
	  int idxGenEle1, idxGenEle2, idxGenEleFS1, idxGenEleFS2;
	  idxGenEle1 = idxGenEle2 = idxGenEleFS1 = idxGenEleFS2 = -1;
	  int nGenDielectrons = 0;
	  int nGenDielectronsFS = 0;
	  for(int kLep=0;kLep<GENnPair;kLep++)
	    {
	    for(int lLep=kLep+1;lLep<GENnPair;lLep++)
	      {
	      // Require a dielectron
	      if(!(abs(GENLepton_ID[kLep])==11 && abs(GENLepton_ID[lLep])==11))
		continue;
	      // Require opposite signs
	      if(GENLepton_ID[kLep]*GENLepton_ID[lLep]>0) continue;
	      if(GENLepton_isHardProcess[kLep]==1 && GENLepton_isHardProcess[lLep]==1)
		{
		// Found a dielectron from hard process
		idxGenEle1 = kLep;
		idxGenEle2 = lLep;
		nGenDielectrons++;
		}
	      if(GENLepton_fromHardProcessFinalState[kLep]==1 && 
                 GENLepton_fromHardProcessFinalState[lLep]==1)
		{
		  // Found a dielectron from final state
		  idxGenEleFS1 = kLep;
		  idxGenEleFS2 = lLep;
		  nGenDielectronsFS++;
		}
	      } // end inner loop over gen leptons
	    } // end outer loop over gen leptons
	  	  
	  //Reco loop
	  double invMassReco;
	  invMassReco=0;
	  int idxRecoEle1,idxRecoEle2;
	  idxRecoEle1=idxRecoEle2=-1;
	  if(Nelectrons<2) continue;
	  int nEle = 0;
	  for(int iEle = 0; iEle < Nelectrons; iEle++)
	    {
	      for(int jEle = iEle+1; jEle < Nelectrons; jEle++)
		{	  
		  if(!passDileptonKinematics(Electron_pT[iEle],Electron_pT[jEle],
                    Electron_eta[iEle],Electron_eta[jEle])) continue; 	
		  if(!Electron_passMediumID[iEle]) continue;//iLep electron ID cut
		  if(!Electron_passMediumID[jEle]) continue;//jLep electron ID cut
		  nEle++;
		  //Reco electrons which passed cuts
		  if(nEle==1)//keeping only pairs of electrons per event
		    {
		      idxRecoEle1 = iEle;
		      idxRecoEle2 = jEle;		    
		    }
		}//end inner reco loop	   
	    }//end reco loop

	  if(idxRecoEle1>=0&&idxRecoEle2>=0)
	    invMassReco=calcInvMass(Electron_pT[idxRecoEle1],Electron_eta[idxRecoEle1],
              Electron_phi[idxRecoEle1],eMass,Electron_pT[idxRecoEle2],
              Electron_eta[idxRecoEle2],Electron_phi[idxRecoEle2],eMass);	  
	  if(nGenDielectrons==0) 
	    continue; // must be DY->mumu or tautau event, skip it
	  
	  if(nGenDielectrons>=2)
	    {
	      // Strange, there should be only two electrons from hard process
	      printf("More than two hard process dielectrons found\n");
	      // skip event, but count the number of cases
	      nTooManyDielectrons++;
	      continue;
	    }
	  
	  if(nGenDielectronsFS!=1)
	    {
	      // Odd, by now we should have only one pair 
	      printf("More than two hard process final state dielectrons found\n");
	      // skip event, but count the number of cases
	      nTooManyDielectronsFS++;
	      continue;
	    }	    
	  
	  invMassHardProcess = calcInvMass(GENLepton_pT[idxGenEle1],GENLepton_eta[idxGenEle1],
	    GENLepton_phi[idxGenEle1],eMass,GENLepton_pT[idxGenEle2],GENLepton_eta[idxGenEle2],            GENLepton_phi[idxGenEle2],eMass);		  
	  invMassFSR = calcInvMass(GENLepton_pT[idxGenEleFS1],GENLepton_eta[idxGenEleFS1],
            GENLepton_phi[idxGenEleFS1],eMass,GENLepton_pT[idxGenEleFS2],
            GENLepton_eta[idxGenEleFS2],GENLepton_phi[idxGenEleFS2],eMass);		  

	  hHardProcess[iChain]->Fill(invMassHardProcess,totalWeight);
	  
	  //HLT cut
	  trigNameSize = pHLT_trigName->size();
	  bool passHLT = kFALSE;	  
	  for(int iHLT=0;iHLT<trigNameSize;iHLT++)
            {
              trigName = pHLT_trigName->at(iHLT);
              if(trigName.CompareTo(HLTname)==0)
		{
                   if(HLT_trigFired[iHLT]==1) 
		     {
		       passHLT = kTRUE;	
		     }
		   else
		     {
		       passHLT = kFALSE;
		     }		     
                   break; 
		} 
	    }// end loop over triggers

	  // Fill histograms for acceptance and efficiency
	  // First, fill histogram for all dielectrons
	  histInvMass[ALL_ELE]->Fill(invMassFSR,totalWeight);
            
	  // Apply kinematic acceptance criteria
	  if(!passDileptonKinematics(GENLepton_pT[idxGenEleFS1],GENLepton_pT[idxGenEleFS2],
				     GENLepton_eta[idxGenEleFS1], GENLepton_eta[idxGenEleFS2])) 
	    continue;	      
	  // Both electrons are in kinematic acceptance at gen level
	  histInvMass[KINEMATIC_CUTS]->Fill(invMassFSR,totalWeight);
	  hpTvsMass->Fill(invMassFSR,GENLepton_pT[idxGenEleFS1],totalWeight);
	  hpTvsMass->Fill(invMassFSR,GENLepton_pT[idxGenEleFS2],totalWeight);

	  // Apply matching to reconstructed electrons requirement
	  int closestTrackLep1, closestTrackLep2;
	  closestTrackLep1 = closestTrackLep2 = -1;
	  bool genToRecoMatchedLep1 = findGenToRecoMatch(idxGenEleFS1,closestTrackLep1);	      
	  bool genToRecoMatchedLep2 = findGenToRecoMatch(idxGenEleFS2,closestTrackLep2);	      
	  if(!(genToRecoMatchedLep1 && genToRecoMatchedLep2)) continue;

	  // Both electrons are reconstructed
	  histInvMass[RECO_MATCHED]->Fill(invMassFSR,totalWeight);
	  
	  // Apply ID criteria:
	  // Dilepton pair at gen level matched to reco and passing ID at reco level
	  if(!Electron_passMediumID[closestTrackLep1]) continue;
	  if(!Electron_passMediumID[closestTrackLep2]) continue;
	  // Both electrons pass ID
	  histInvMass[ID_CUTS]->Fill(invMassFSR,totalWeight);
	  
	  // Apply HLT requirement
	  if(!passHLT) continue;
          
          eEta1 = Electron_eta[idxRecoEle1];
          eEta2 = Electron_eta[idxRecoEle2];
          ePt1 = Electron_pT[idxRecoEle1];
          ePt2 = Electron_pT[idxRecoEle2];

          if(ePt1<ptBinLow) ePt1 = ptBinLow;//pull this information from the histograms
          if(ePt2<ptBinLow) ePt2 = ptBinLow;//raise bin
          if(ePt1>ptBinHigh) ePt1 = ptBinHigh;//lower bin
          if(ePt2>ptBinHigh) ePt2 = ptBinHigh;//


	  // Event passed HLT cut
	  sfReco1=hRecoSF->GetBinContent(hRecoSF->FindBin(eEta1,ePt1));
          sfReco2=hRecoSF->GetBinContent(hRecoSF->FindBin(eEta2,ePt2));
          sfID1=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eEta1,ePt1));
          sfID2=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eEta2,ePt2));
          sfHLT=(hLeg2SF->GetBinContent(hLeg2SF->FindBin(eEta1,ePt1)))*
            (hLeg2SF->GetBinContent(hLeg2SF->FindBin(eEta2,ePt2)));
          sfWeight = sfReco1*sfReco2*sfID1*sfID2*sfHLT;

	  histInvMass[HLT_CUTS]->Fill(invMassFSR,totalWeight);
	  histInvMass[RECO_ELE]->Fill(invMassReco,totalWeight);
	  migMatrixGENisHardvsGENFS->Fill(invMassHardProcess,invMassFSR,totalWeight);
	  migMatrixGENFSvsReco->Fill(invMassFSR,invMassReco,totalWeight*sfWeight);
	  migMatrixGENisHardvsReco->Fill(invMassHardProcess,invMassReco,totalWeight);
	}//end event loop   
      
    }//end chain loop 
  
  TEfficiency*hMatchedEfficiency = 
    new TEfficiency((*histInvMass[RECO_MATCHED]),(*histInvMass[KINEMATIC_CUTS]));
  hMatchedEfficiency->SetTitle("Reconstruction Efficiency");
  hMatchedEfficiency->SetMarkerStyle(20);
  hMatchedEfficiency->SetMarkerSize(0.5);
  hMatchedEfficiency->SetName("RecoEfficiency");
  hMatchedEfficiency->SetStatisticOption(TEfficiency::kFNormal);
  TEfficiency* hIDEfficiency = 
    new TEfficiency((*histInvMass[ID_CUTS]),(*histInvMass[RECO_MATCHED]));
  hIDEfficiency->SetTitle("ID Efficiency");
  hIDEfficiency->SetMarkerStyle(20);
  hIDEfficiency->SetMarkerSize(0.5);
  hIDEfficiency->SetName("IDEfficiency");
  hIDEfficiency->SetStatisticOption(TEfficiency::kFNormal);
  TEfficiency* hAcceptance = 
    new TEfficiency((*histInvMass[KINEMATIC_CUTS]),(*histInvMass[ALL_ELE]));
  hAcceptance->SetTitle("Acceptance");
  hAcceptance->SetMarkerStyle(20);
  hAcceptance->SetMarkerSize(0.5);
  hAcceptance->SetName("Acceptance");  
  hAcceptance->SetStatisticOption(TEfficiency::kFNormal);
  TEfficiency* hHLTEfficiency = 
    new TEfficiency((*histInvMass[HLT_CUTS]),(*histInvMass[ID_CUTS]));
  hHLTEfficiency->SetTitle("HLT Efficiency");
  hHLTEfficiency->SetMarkerStyle(20);
  hHLTEfficiency->SetMarkerSize(0.5);
  hHLTEfficiency->SetName("HLTEfficiency");
  hHLTEfficiency->SetStatisticOption(TEfficiency::kFNormal);
  TEfficiency* hEfficiency = 
    new TEfficiency((*histInvMass[HLT_CUTS]),(*histInvMass[KINEMATIC_CUTS]));
  hEfficiency->SetTitle("Total Efficiency");
  hEfficiency->SetMarkerStyle(20);
  hEfficiency->SetMarkerSize(0.5);
  hEfficiency->SetName("Efficiency");
  hEfficiency->SetStatisticOption(TEfficiency::kFNormal);
  
  rootFile->cd();
  for(int i=0;i<nInvMassHistos;i++){
    histInvMass[i]->Write();
  } 
  hAcceptance->Write();
  hEfficiency->Write();
  hHLTEfficiency->Write();
  hpTvsMass->Write();
  migMatrixGENisHardvsGENFS->Write();
  migMatrixGENFSvsReco->Write();
  migMatrixGENisHardvsReco->Write();
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
  cout << "**************************************************************************" << endl;
  cout << endl;
  
}//end main function

//Counter for tracking program progress
void counter(Long64_t i, Long64_t N)
{
  int P = 100*(i)/(N);  
  TTimeStamp eventTimeStamp;
  if(i%(N/100)==0)
    {
      cout << "efficiencies.C " << "[Time: " << eventTimeStamp.AsString("s") << "] " << P << "%" << endl;
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
    }//end of Loop 
  bool matchFound = kTRUE;
  if(dRMin>=dRMinCut)
    {
      recoIndex=-1;
      matchFound=kFALSE;
    }
  return matchFound;
}//end findGetToRecoMatch

//Kinematic cuts
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2)
{
  //if(abs(eta1)>etaGapLow && abs(eta1)<etaGapHigh) return kFALSE;//eta cut
  //if(abs(eta2)>etaGapLow && abs(eta2)<etaGapHigh) return kFALSE; //eta cut
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
