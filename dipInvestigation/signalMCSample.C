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
int Nelectrons, HLT_ntrig, nVertices, nPileUp;
double GENEvt_weight;
double Electron_pT[MPSIZE], Electron_eta[MPSIZE], Electron_phi[MPSIZE];
double Electron_Energy[MPSIZE], Electron_Px[MPSIZE];
double Electron_Py[MPSIZE], Electron_Pz[MPSIZE], Electron_charge[MPSIZE];
bool Electron_passMediumID[MPSIZE];
int HLT_trigType[MPSIZE],HLT_trigFired[MPSIZE];
std::vector<std::string> HLT_trigName;
std::vector<std::string> *pHLT_trigName = &HLT_trigName;

const int nHistos = 5;
enum Histos {
  XSEC_WEIGHTS,
  SF_WEIGHTS,
  GEN_WEIGHTS,
  PU_WEIGHTS,
  TOTAL_WEIGHTS
};

TString histNames[]= {
  "hXsecWeight",
  "hSFweight",
  "hGenWeight",
  "hPUWeight",
  "hTotWeight"
};
const TString histTitles[] = {
  "Cross Section Weights Only",
  "Scale Factors Only",
  "Gen Weights Only",
  "Pileup Weights Only",
  "All Weights Applied"
};
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
                            98.5,101,103.5,106,108.25,108.5,108.75,109,109.25,109.5,109.75,110,112.5,115,117.5,120,123,126,129.5,133,
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

void signalMCSample()
{
  TTimeStamp ts_start;
  cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
  TStopwatch totaltime;
  totaltime.Start();
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
      files += "/skims_0002/ntuple_skim_";
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
    chains->SetBranchAddress("Nelectrons", &Nelectrons, &b_Nelectrons);
    chains->SetBranchAddress("nVertices", &nVertices, &b_nVertices);
    chains->SetBranchAddress("nPileUp", &nPileUp, &b_nPileUp);
    chains->SetBranchAddress("Electron_pT", &Electron_pT, &b_Electron_pT);
    chains->SetBranchAddress("Electron_eta",&Electron_eta, &b_Electron_eta);
    chains->SetBranchAddress("Electron_phi",&Electron_phi, &b_Electron_phi);
    chains->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,
      &b_Electron_passMediumID);
    chains->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
    chains->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
    chains->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
    chains->SetBranchAddress("HLT_trigName",&pHLT_trigName);   
    chains->SetBranchAddress("GENEvt_weight",&GENEvt_weight,&b_GENEvt_weight);
  
  cout << "Total Events Loaded: " << nentries << endl;
  cout << endl;

  TFile*pileupRatioFile  = new TFile(pileupRatioName);
  TH1F*hPileupRatio = (TH1F*)pileupRatioFile->Get("hPileupRatio");
  TFile*fileLeg2SF = new TFile(leg2SFName);
  TH2F*hLeg2SF = (TH2F*)fileLeg2SF->Get("EGamma_SF2D");
  TFile*fileMedIDSF = new TFile(medIDSFName);
  TH2F*hMedIDSF = (TH2F*)fileMedIDSF->Get("EGamma_SF2D");
  TFile*fileRecoSF = new TFile(recoSFName); 
  TH2F*hRecoSF = (TH2F*)fileRecoSF->Get("EGamma_SF2D");
 

  TH1F*hSignalMC[nHistos];
  for(int i=0;i<nHistos;i++){
    hSignalMC[i] = new TH1F(histNames[i],"",100,100,125); 
    hSignalMC[i]->Sumw2();
    hSignalMC[i]->SetFillColor(kBlue+2);
    hSignalMC[i]->GetXaxis()->SetNoExponent();
    hSignalMC[i]->GetXaxis()->SetMoreLogLabels();
    hSignalMC[i]->SetTitle(histTitles[i]);
    hSignalMC[i]->GetXaxis()->SetTitle("Dielectron invariant mass [GeV]");
  }

  double varGenWeight, lumiEffective, nEffective, localEntry, sumGenWeight, sumRawGenWeight, 
    totalWeight, sfWeight, weightNoPileup, xSecWeightGen, genWeight, pileupWeight, 
    xSecWeightAlone;
  TString compareHLT = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
  TString trigName;
  int trigNameSize;
  double lumi = dataLuminosity;//luminosity for xsec weighting
  double sfReco1,sfReco2,sfID1,sfID2,sfHLT;//efficiency scale factors
  double eEta1, eEta2, ePt1, ePt2;

  cout << endl;
  cout << "Starting to Normalize Gen Weights" << endl;
  //Finding normalized genWeights,sums,variances
  sumGenWeight = 0.0;
  sumRawGenWeight = 0.0;
  varGenWeight = 0.0;
  for(Long64_t i=0;i<nentries;i++){
    counter(i,nentries);
    localEntry = chains->LoadTree(i);
    b_GENEvt_weight->GetEntry(localEntry);
    genWeight = GENEvt_weight/fabs(GENEvt_weight);	//normalized genweight
    sumGenWeight += genWeight;
    varGenWeight += GENEvt_weight*GENEvt_weight; //variance of genweights
    sumRawGenWeight += GENEvt_weight; 
  }          
  nEffective = (sumRawGenWeight*sumRawGenWeight)/varGenWeight;
  lumiEffective = nEffective/xSec;
  
  //Event loop
  cout << endl;
  cout << "Starting Event Loop" << endl;
  for(Long64_t i=0;i<nentries;i++) {      
    counter(i,nentries);
    chains->GetEntry(i);
    if(Nelectrons<2) continue;   	  
    //HLT cut
    trigNameSize = pHLT_trigName->size();
    bool passHLT = kFALSE;	  

    for(int iHLT=0;iHLT<trigNameSize;iHLT++) {
      trigName = pHLT_trigName->at(iHLT);	  
      if(trigName.CompareTo(compareHLT)==0) {
        if(HLT_trigFired[iHLT]==1) {
          passHLT = kTRUE;	
        }
        else {
          passHLT = kFALSE;
        }		     
        break; 
      }
    } 
    if(!passHLT) continue;
   
    bool passNumEle = kFALSE;
    bool passKinematics = kFALSE;
    int numDielectrons = 0;
    int subEle = -1;
    int leadEle = -1;
    double invMassGen = -5000;
    double invMassReco = -5000;
    TLorentzVector recoDielectronP4;
    TLorentzVector genDielectronP4;

    //Electron loop
    for(int iEle = 0; iEle < Nelectrons; iEle++) {
      if(!Electron_passMediumID[iEle]) continue;
      for(int jEle = iEle+1; jEle < Nelectrons; jEle++) {
        if(!Electron_passMediumID[jEle]) continue;
        if(passDileptonKinematics(Electron_pT[iEle],Electron_pT[jEle],
          Electron_eta[iEle],Electron_eta[jEle])){
          numDielectrons++;
          if(Electron_pT[iEle]>Electron_pT[jEle]){
            leadEle = iEle; subEle = jEle;
          }
          else {
            leadEle = jEle; subEle = iEle;
          }	  
        }
      }//end jEle loop
    }//end iEle loop
     
   if(numDielectrons==1){
      passNumEle = kTRUE;
   }
   if(!passNumEle) continue;
   if(leadEle<0||subEle<0) continue;
   recoDielectronP4 = getDielectronP4(Electron_pT[leadEle],Electron_eta[leadEle],
     Electron_phi[leadEle],eMass,Electron_pT[subEle],Electron_eta[subEle],
     Electron_phi[subEle],eMass);
   genDielectronP4 = getDielectronP4(Electron_pT[leadEle],Electron_eta[leadEle],
     Electron_phi[leadEle],eMass,Electron_pT[subEle],Electron_eta[subEle],
     Electron_phi[subEle],eMass);
   
   invMassReco=recoDielectronP4.M();
   invMassGen=genDielectronP4.M();
   if(invMassReco<-1000) continue;
   
   eEta1 = Electron_eta[leadEle];
   eEta2 = Electron_eta[subEle];
   ePt1 = Electron_pT[leadEle];
   ePt2 = Electron_pT[subEle];
   
   if(ePt1<ptBinLow) ePt1 = ptBinLow;//pull this information from the histograms
   if(ePt2<ptBinLow) ePt2 = ptBinLow;//raise bin
   if(ePt1>ptBinHigh) ePt1 = ptBinHigh;//lower bin
   if(ePt2>ptBinHigh) ePt2 = ptBinHigh;//

   //Determining weighting factors
   xSecWeightAlone=lumi*(xSec/nentries);//xSecWeight when used alone
   xSecWeightGen = lumi*(xSec/1.0);//xSecWeight when used with Gen weight 
   pileupWeight = hPileupRatio->GetBinContent(hPileupRatio->FindBin(nPileUp));
   genWeight = (GENEvt_weight/fabs(GENEvt_weight))/sumGenWeight;

   sfReco1=hRecoSF->GetBinContent(hRecoSF->FindBin(eEta1,ePt1));
   sfReco2=hRecoSF->GetBinContent(hRecoSF->FindBin(eEta2,ePt2));
   sfID1=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eEta1,ePt1));
   sfID2=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eEta2,ePt2));
   sfHLT=(hLeg2SF->GetBinContent(hLeg2SF->FindBin(eEta1,ePt1)))*
     (hLeg2SF->GetBinContent(hLeg2SF->FindBin(eEta2,ePt2)));
   sfWeight = sfReco1*sfReco2*sfID1*sfID2*sfHLT;
   totalWeight = sfWeight*genWeight*xSecWeightGen*pileupWeight;

     hSignalMC[XSEC_WEIGHTS]->Fill(invMassReco,xSecWeightAlone);
     hSignalMC[GEN_WEIGHTS]->Fill(invMassReco,genWeight);
     hSignalMC[SF_WEIGHTS]->Fill(invMassReco,sfWeight);
     hSignalMC[PU_WEIGHTS]->Fill(invMassReco,pileupWeight);
     hSignalMC[TOTAL_WEIGHTS]->Fill(invMassReco,totalWeight);
   }//end event loop   

  TCanvas*canvas[nHistos];  
  TString canvasName;
  TString saveName;
  for(int i=0;i<nHistos;i++){
    canvasName = "canvas";
    canvasName += i;
    canvas[i] = new TCanvas(canvasName,"",10,10,1000,1000);
    canvas[i]->SetGrid();
    canvas[i]->SetLogy();
    //canvas[i]->SetLogx();
    hSignalMC[i]->GetXaxis()->SetRangeUser(100,125);
    hSignalMC[i]->Draw("hist");
    saveName = canvasName;
    saveName += "_ManyBins.png";
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
