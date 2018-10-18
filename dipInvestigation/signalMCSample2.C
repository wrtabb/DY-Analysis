#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "THStack.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
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
#include "TRatioPlot.h"
#include "THashList.h"
#include "TProfile.h"

void counter(Long64_t i, Long64_t N);
TLorentzVector getDielectronP4(double pt1,double eta1,double phi1,double m1,double pt2,
  double eta2,double phi2,double m2);
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2);

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

const int nHistos = 3;
enum histoTypes {
  GENFS,
  GENHARD,
  RECO
};
const TString histTitle[nHistos] = {"Gen-Level Final State","Gen-Level Hard Process",
  "Reco-Level"};

const TString histSave[nHistos] = {"_GENFS","_GENHARD","_RECO"};
//cutting parameters
const float etaHigh = 2.4;
const float etaGapHigh = 1.566; 
const float etaGapLow = 1.4442;
const float ptHigh = 28;
const float ptLow = 17;
const float dRMinCut = 0.3;
const int ptBinHigh = 499;
const int ptBinLow = 26;
const int invMassHigh = 120;
const int invMassLow = 60;
const float eMass = 0.000511;
const int dataLuminosity = 35867; //Run2016B to Run2016H JSON. unit: /pb, Updated at 2017.07.30
const TString treeName = "recoTree/DYTree";
const TString pileupRatioName = "/home/hep/wrtabb/git/DY-Analysis/plots/pileup.root";
const TString leg2SFName = "/home/hep/wrtabb/git/DY-Analysis/data/Leg2_SF.root";
const TString medIDSFName = "/home/hep/wrtabb/git/DY-Analysis/data/MediumID_SF.root";
const TString recoSFName = "/home/hep/wrtabb/git/DY-Analysis/data/Reco_SF.root";

const int numChains = 48; const double pi=TMath::Pi(); const float axisLow = 0.0001;

//InvMass
const int nBins = 2*46;
//const double massbins[] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 
//                             91,96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 
//                             185,200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 
//                             1500,3000};
//const double massbins[] = {100,105,107.5,108.5,109,109.5,110,110.5,111,112,115,120};
const double massbins[] = {15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,
                            57.5,60,62,64,66,68,70,72,74,76,78.5,81,83.5,86,88.5,91,93.5,96,
                            98.5,101,103.5,106,108.25,108.5,108.75,109,109.25,109.5,109.75,
                            110,112.5,115,117.5,120,123,126,129.5,133,
                            137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,210,220,
                            231.5,243,258,273,296.5,320,350,380,410,440,475,510,555,600,650,
                            700,765,830,915,1000,1250,1500,2250,3000};
const float binLowInvMass = 0;
const float binHighInvMass = 3000;
//InvMass Linear Plot
const int nBinsInvMassLinear = 60;
const float binLowInvMassLinear = 60;
const float binHighInvMassLinear = 120;
//nVertices
const int nBinsVert = 50;
const float binLowVert = 0;
const float binHighVert = 50;
//pT 
const int npTBins = 100;
const float binLowpT = 0;
const float binHighpT = 500;
//eta
const int nEtaBins = 50;
const float binLowEta = -2.5;
const float binHighEta = 2.5;
//rapidity
const int nYBins = 50;
const float binLowY = -2.5;
const float binHighY = 2.5;

const float xSec = 1921.8*3;

void signalMCSample2()
{
  TTimeStamp ts_start;
  cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
  TStopwatch totaltime;
  totaltime.Start();
  gStyle->SetOptStat(0);

  //Defining branches
  TBranch*b_GENnPair;
  TBranch*b_GENLepton_eta;
  TBranch*b_GENLepton_phi;
  TBranch*b_GENLepton_pT;
  TBranch*b_GENLepton_ID;
  TBranch*b_GENLepton_isHardProcess;
  TBranch*b_GENLepton_fromHardProcessFinalState;
  TBranch*b_Nelectrons;
  TBranch*b_Electron_pT;
  TBranch*b_Electron_eta;
  TBranch*b_Electron_phi;
  
  //Loading ntuples
  cout << "Loading ntuples" << endl;
  //The names of every directory being loaded
  TString dirName = 
    "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8";

  TString baseDir =
    "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/ikrav/DrellYan_13TeV_2016/v2p3";
  TString files;  
  Long64_t nentries = 0;
  TChain*chains = new TChain(treeName);
  int nFiles = 0;
  int maxFiles = 100;//maximum number of files to load
  cout << "****************************************************" << endl;
  cout << "Loading ntuples" << endl;
    for(int k=0;k<1000;k++) {	  	      
      if(nFiles>=maxFiles) break;
      files = baseDir;      
      files += dirName;
      //files += "/skims_0002/ntuple_skim_";
      files += "/ntuple_skim_";
      files += k;
      files += ".root";
      std::ifstream testFileStream(files);
      if(!(bool)testFileStream) continue;
      chains->Add(files);
      nFiles++;
      cout << files << endl;
      cout << chains->GetEntries() << " events loaded" << endl;	 
      
      if(chains->GetEntries()==0){//error message if no events loaded
	cout << endl;
	cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	cout << "ERROR: Broken files or files not found in: " << endl;
	cout << files << endl;
	cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	cout << endl;
	return;
      }
    }
    nentries = chains->GetEntries(); 

    //Setting addresses for branches
    chains->SetBranchAddress("GENnPair", &GENnPair, &b_GENnPair);
    chains->SetBranchAddress("GENLepton_eta", &GENLepton_eta, &b_GENLepton_eta);
    chains->SetBranchAddress("GENLepton_phi",&GENLepton_phi, &b_GENLepton_phi);
    chains->SetBranchAddress("GENLepton_pT",&GENLepton_pT, &b_GENLepton_pT);
    chains->SetBranchAddress("GENLepton_ID",&GENLepton_ID, &b_GENLepton_ID);
    chains->SetBranchAddress("GENLepton_isHardProcess",&GENLepton_isHardProcess,
                                     &b_GENLepton_isHardProcess);
    chains->SetBranchAddress
      ("GENLepton_fromHardProcessFinalState",&GENLepton_fromHardProcessFinalState,
                                     &b_GENLepton_fromHardProcessFinalState);
    chains->SetBranchAddress("Nelectrons", &Nelectrons, &b_Nelectrons);
    chains->SetBranchAddress("Electron_pT", &Electron_pT, &b_Electron_pT);
    chains->SetBranchAddress("Electron_eta",&Electron_eta, &b_Electron_eta);
    chains->SetBranchAddress("Electron_phi",&Electron_phi, &b_Electron_phi);
  
  cout << "Total Events Loaded: " << nentries << endl;
  cout << endl;
  TH1F*hSignalMC[nHistos];
  for(int i=0;i<nHistos;i++){
    TString histname = "hSignalMC";
    histname += i;
    hSignalMC[i] = new TH1F(histname,"",nBins,massbins);
    hSignalMC[i]->Sumw2();
    hSignalMC[i]->SetFillColor(kBlue+2);
    hSignalMC[i]->SetTitle(histTitle[i]);
    hSignalMC[i]->GetXaxis()->SetTitle("Dielectron invariant mass [GeV]");
    hSignalMC[i]->GetXaxis()->SetRangeUser(100,125);
  }
  //Event loop
  double invMass;
  int leadEle,subEle;
  cout << endl;
  cout << "Starting Event Loop" << endl;
  for(Long64_t i=0;i<nentries;i++) {      
    counter(i,nentries);
    chains->GetEntry(i);
    TLorentzVector DielectronP4;

    // Loop over gen leptons and find the electron pair at the isHardProcess
    // and isHardProcessFinalState level.
    int idxGenEle1, idxGenEle2, idxGenEleFS1, idxGenEleFS2;
    idxGenEle1 = idxGenEle2 = idxGenEleFS1 = idxGenEleFS2 = -1;
    int nGenDielectrons = 0;
    int nGenDielectronsFS = 0;
    for(int kLep=0;kLep<GENnPair;kLep++){
      for(int lLep=kLep+1;lLep<GENnPair;lLep++){
        // Require a dielectron
        if(!(abs(GENLepton_ID[kLep])==11 && abs(GENLepton_ID[lLep])==11)) continue;
        // Require opposite signs
        if(GENLepton_ID[kLep]*GENLepton_ID[lLep]>0) continue;

        if(GENLepton_isHardProcess[kLep]==1 && GENLepton_isHardProcess[lLep]==1){
          // Found a dielectron from hard process
          idxGenEle1 = kLep;
          idxGenEle2 = lLep;
          nGenDielectrons++;
        }
        if(GENLepton_fromHardProcessFinalState[kLep]==1 &&
          GENLepton_fromHardProcessFinalState[lLep]==1){
          // Found a dielectron from final state
          idxGenEleFS1 = kLep;
          idxGenEleFS2 = lLep;
          nGenDielectronsFS++;
        }
      }// end inner loop
    }// end outer loop

    // Calculate and fill Gen-Level Final State
    DielectronP4 = getDielectronP4(GENLepton_pT[idxGenEleFS1],GENLepton_eta[idxGenEleFS1],
    GENLepton_phi[idxGenEleFS1],eMass,GENLepton_pT[idxGenEleFS2],GENLepton_eta[idxGenEleFS2],
    GENLepton_phi[idxGenEleFS2],eMass);
    invMass=DielectronP4.M();
    hSignalMC[GENFS]->Fill(invMass);
    
    // Calculate and fill Gen-Level from hard process
    DielectronP4 = getDielectronP4(GENLepton_pT[idxGenEle1],GENLepton_eta[idxGenEle1],
    GENLepton_phi[idxGenEle1],eMass,GENLepton_pT[idxGenEle2],GENLepton_eta[idxGenEle2],
    GENLepton_phi[idxGenEle2],eMass);
    invMass=DielectronP4.M();
    hSignalMC[GENHARD]->Fill(invMass);

    //Reco loop
    if(Nelectrons<2) continue;   	  
    int numDielectrons = 0;
    int subEle = -1;
    int leadEle = -1;
    for(int iEle = 0; iEle < Nelectrons; iEle++) {
      for(int jEle = iEle+1; jEle < Nelectrons; jEle++) {
        if(Electron_pT[iEle]>Electron_pT[jEle]){
          leadEle = iEle; subEle = jEle;
        }
        else {
          leadEle = jEle; subEle = iEle;
        }	  
        
      }//end jEle loop
    }//end iEle loop
     
    if(leadEle<0||subEle<0) continue;

    // Calculate and fill Reco-level
    DielectronP4 = getDielectronP4(Electron_pT[leadEle],Electron_eta[leadEle],
      Electron_phi[leadEle],eMass,Electron_pT[subEle],Electron_eta[subEle],
      Electron_phi[subEle],eMass);
    invMass=DielectronP4.M();
    hSignalMC[RECO]->Fill(invMass);
  }//end event loop   

  // Draw to canvases and save
  TCanvas*canvas[nHistos];
  for(int i=0;i<nHistos;i++){
    canvas[i] = new TCanvas("canvas","",10,10,1000,1000);
    canvas[i]->SetGrid();
    canvas[i]->SetLogy();
    hSignalMC[i]->Draw("hist");
    TString saveName = "noSkim_noCuts_noWeights";
    saveName += histSave[i];
    saveName += ".png";
    canvas[i]->SaveAs(saveName);
  }

  totaltime.Stop();
  Double_t TotalCPURunTime = totaltime.CpuTime();
  Double_t TotalRunTime = totaltime.RealTime();
  TTimeStamp ts_end;
  cout << endl;
  cout << "**************************************************************************" << endl;
  cout << "Total CPU RunTime: " << TotalCPURunTime/60 << " minutes" << endl;
  cout << "Total Real RunTime: " << TotalRunTime/60 << " minutes" << endl;
  cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;   
  cout << "Number of Events Processed: " << nentries << endl;
  cout << "**************************************************************************" << endl;
  cout << endl;
  
}//end main function

//Counter for tracking program progress
void counter(Long64_t i, Long64_t N)
{
  int P = 100*(i)/(N);  
  TTimeStamp eventTimeStamp;
  if(i%(N/100)==0) {
    cout << "signalMCSample.C " << "[Time: " << eventTimeStamp.AsString("s") << "] " << P 
      << "%" << endl;
  }
  return;
}

//Kinematic cuts
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2)
{
  //if(abs(eta1)>etaGapLow && abs(eta1)<etaGapHigh) return kFALSE;
  //if(abs(eta2)>etaGapLow && abs(eta2)<etaGapHigh) return kFALSE; 
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
