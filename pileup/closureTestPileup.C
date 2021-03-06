

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "THStack.h"
#include "TH1F.h"
#include "TH1I.h"
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

void counter(Long64_t i, Long64_t N);
double calcInvMass(double pt1,double eta1,double phi1,double m1,double pt2,double eta2,double phi2,double m2);
double calcRapidity(double pt1,double eta1,double phi1,double m1,double pt2,double eta2,double phi2,double m2);
double calcDileptonPt(double pt1,double eta1,double phi1,double m1,double pt2,double eta2,double phi2,double m2);
double calcDileptonEta(double pt1,double eta1,double phi1,double m1,double pt2,double eta2,double phi2,double m2);
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2);

//Defining variables and arrays
const int MPSIZE = 2000;
const int nPileupBins = 100;
const float pileupBinLow = 0;
const float pileupBinHigh = nPileupBins;
int Nelectrons, HLT_ntrig, nVertices, nPileUp;
double GENEvt_weight;
double Electron_pT[MPSIZE], Electron_eta[MPSIZE], Electron_phi[MPSIZE];
double Electron_Energy[MPSIZE], Electron_Px[MPSIZE];
double Electron_Py[MPSIZE], Electron_Pz[MPSIZE], Electron_charge[MPSIZE];
bool Electron_passMediumID[MPSIZE];
int HLT_trigType[MPSIZE],HLT_trigFired[MPSIZE];
std::vector<std::string> HLT_trigName;
std::vector<std::string> *pHLT_trigName = &HLT_trigName;

//cutting parameters
const float etaHigh = 2.4;
const float etaGapHigh = 1.566; 
const float etaGapLow = 1.4442;
const float ptHigh = 28;
const float ptLow = 17;
const float dRMinCut = 0.3;

const float eMass = 0.000511;
const int dataLuminosity = 35867; //Run2016B to Run2016H JSON. unit: /pb, Updated at 2017.07.30
const TString treeName = "recoTree/DYTree";
const TString pileupRatioName = "./plots/pileup.root";
const int numChains = 48; const double pi=TMath::Pi(); const float axisLow = 0.0001;
enum ChainNum {
  QCD20to30,
  QCD30to50,
  QCD50to80,
  QCD80to120,
  QCD120to170,
  QCD170to300,
  QCD300toInf,  
  wJets,
  WW,    
  WWTo2L2Nu,
  ZZ,       
  ZZTo4L,
  WZ,
  WZTo3LNu,
  tt0to700,
  tt700to1000,
  tt1000toInf,
  tW,             
  tbarW,         
  EE10to50,
  EE50to100,
  EE100to200,
  EE200to400,
  EE400to500,
  EE500to700,
  EE700to800,
  EE800to1000,
  EE1000to1500,
  EE1500to2000,
  EE2000to3000,
  TAUTAU10to50,
  TAUTAU50to100,
  TAUTAU100to200,
  TAUTAU200to400,
  TAUTAU400to500,
  TAUTAU500to700,
  TAUTAU700to800,
  TAUTAU800to1000,
  TAUTAU1000to1500,
  TAUTAU1500to2000,
  TAUTAU2000to3000,
  DataRunB,
  DataRunC,
  DataRunD,
  DataRunE,
  DataRunF,
  DataRunG,
  DataRunH
}; 

//Histogram parameters
const int nHistoTypes = 11; 
const int nHistos = 5;
const TString histName[nHistos] = {
  "hFakes", "hEW", "hTops", "hMC", "hData"
};
const TString histTypeName[nHistoTypes] = {
  "InvMass", "InvMassLinear", "Vert", "VertWeighted", "pTLead", "pTSub", "pTDi", "EtaLead", "EtaSub", "EtaDi",
  "Rapidity"
};
const Color_t histFillColors[nHistos] = {
  kViolet+5, kRed+2, kBlue+2, kOrange-2, kWhite
};
const Color_t histLineColors[nHistos] = {
  kViolet+3, kRed+4, kBlue+3, kOrange+3, kBlack
};
const TString xAxisLabels[nHistoTypes] = {
  "Dielectron invariant mass [GeV]",
  "Dielectron invariant mass [GeV]",
  "Number of vertices",
  "Number of vertices"
  "p_{T} [GeV]",
  "p_{T} [GeV]",
  "p_{T} [GeV]",
  "#eta",
  "#eta",
  "#eta",
  "Y"
};
enum HistBins {
  FAKES,
  EW,
  TOPS,
  EE,
  DATA,
  UNDEF = -1
};
enum HistTypes {
  INV_MASS,
  INV_MASS_LINEAR,
  VERTICES,
  VERTICES_WEIGHTED,
  PT_LEAD,
  PT_SUB,
  PT_DI,
  ETA_LEAD,
  ETA_SUB,
  ETA_DI,
  RAPIDITY
};
//InvMass
const int nBinsInvMass = 43;
const float massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 
			    106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 
			    380, 440, 510, 600, 700, 830, 1000, 1500, 3000};
const float binLowInvMass = 0;
const float binHighInvMass = 3000;
//InvMass Linear Plot
const int nBinsInvMassLinear = 30;
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
const int nEtaBins = 100;
const float binLowEta = -2.5;
const float binHighEta = 2.5;
//rapidity
const int nYBins = 100;
const float binLowY = -2.5;
const float binHighY = 2.5;

const int nBins[nHistoTypes] = {nBinsInvMass,nBinsInvMassLinear,nBinsVert,nBinsVert,npTBins,npTBins,npTBins,
				nEtaBins,nEtaBins,nEtaBins,nYBins};
const float binLow[nHistoTypes] = {binLowInvMass,binLowInvMassLinear,binLowVert,binLowVert,binLowpT,binLowpT,
				   binLowpT,binLowEta,binLowEta,binLowEta,binLowY};
const float binHigh[nHistoTypes] = {binHighInvMass,binHighInvMassLinear,binHighVert,binHighVert,binHighpT,
				    binHighpT,binHighpT,binHighEta,binHighEta,binHighEta,binHighY};

//Cross sections calculated by Kyeongpil Lee
const float xSec[numChains] = {5352960,9928000,2890800,350000,62964,18810,1350,//QCD
			       61526.7,118.7,12.178,16.523,1.256,47.13,4.4297,//Bosons
			       734.577,76.605,20.578,35.85,35.85,//tops
			       6016.88,1873.52,76.2401,2.67606,0.139728,0.0792496,0.0123176,0.01042,//DYEE
			       0.00552772,0.000741613,0.000178737,//DYEE
			       6016.88,1873.52,76.2401,2.67606,0.139728,0.0792496,0.0123176,0.01042,//DYTAUTAU
			       0.00552772,0.000741613,0.000178737,//DYTAUTAU
			       1,1,1,1,1,1,1};//data (unweighted)

void closureTestPileup()
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
  
  //Loading ntuples
  cout << "Loading ntuples" << endl;
  //The names of every directory being loaded
  TString dirNames[numChains] = {
    "/QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8",
    "/QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt30to50",
    "/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt50to80",
    "/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt80to120",
    "/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCDEMEnriched_Pt120to170",
    "/QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8",
    "/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8",
    "/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo",
    "/WW_TuneCUETP8M1_13TeV-pythia8",
    "/WWTo2L2Nu_13TeV-powheg",
    "/ZZ_TuneCUETP8M1_13TeV-pythia8",
    "/ZZTo4L_13TeV_powheg_pythia8",
    "/WZ_TuneCUETP8M1_13TeV-pythia8",
    "/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8",
    "/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_truncated_M0To700/crab_ttbar",
    "/TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8",
    "/TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8",
    "/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
    "/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
    "/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE/crab_DYLL_M10to50",
    "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_truncated_M50To100/EE",
    "/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE/crab_DYLL_M100to200",
    "/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TauTau/crab_DYLL_M10to50",
    "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_truncated_M50To100/TauTau",
    "/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TauTau/crab_DYLL_M100to200",
    "/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TauTau",
    "/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TauTau",
    "/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TauTau",
    "/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TauTau",
    "/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TauTau",
    "/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TauTau",
    "/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TauTau",
    "/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TauTau",
    "/DoubleEG/crab_DoubleEG_RunB",
    "/DoubleEG/crab_DoubleEG_RunC",
    "/DoubleEG/crab_DoubleEG_RunD",
    "/DoubleEG/crab_DoubleEG_RunE",
    "/DoubleEG/crab_DoubleEG_RunF",
    "/DoubleEG/crab_DoubleEG_RunG",
    "/DoubleEG/crab_DoubleEG_RunH"
  };

  TString baseDirectory =
    "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/ikrav/DrellYan_13TeV_2016/v2p3";
  TString subDirectoryDY = "/DYJetsToLL_allMasses_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8";
  TChain*chains[numChains];

  vector <TString> *subFiles[numChains];  
  for(int iChain=0;iChain<numChains;iChain++) {
    if(iChain==WWTo2L2Nu||iChain==ZZTo4L||iChain==WZTo3LNu) continue;//Not using these in this analysis
    if(iChain==QCD20to30||iChain==QCD30to50||iChain==QCD50to80||iChain==QCD80to120||iChain==QCD120to170||
       iChain==QCD170to300||iChain==QCD300toInf) continue;//skipping QCD due to possible problems
    
    subFiles[iChain] = new vector<TString>;
    if(iChain==EE10to50||iChain==TAUTAU10to50) {
      subFiles[iChain]->push_back(dirNames[iChain]+"_ext1v1");
      subFiles[iChain]->push_back(dirNames[iChain]+"_v1");
      subFiles[iChain]->push_back(dirNames[iChain]+"_v2");
    }
    else if(iChain==EE100to200||iChain==TAUTAU100to200||iChain==wJets) {
      subFiles[iChain]->push_back(dirNames[iChain]);
      subFiles[iChain]->push_back(dirNames[iChain]+"_ext");
    }
    else if(iChain==DataRunH) {
      subFiles[iChain]->push_back(dirNames[iChain]+"ver2");
      subFiles[iChain]->push_back(dirNames[iChain]+"ver3");
    }
    else if(iChain==tt0to700) {
      subFiles[iChain]->push_back(dirNames[iChain]);
      subFiles[iChain]->push_back(dirNames[iChain]+"Backup");
    }
    else subFiles[iChain]->push_back(dirNames[iChain]);      
  } 
  int nFiles;
  const int maxFiles = 5;
  TString files;  
  Long64_t subDirectorySize;
  Long64_t totalentries = 0;
  for(int iChain=0;iChain<numChains;iChain++) {          
    if(iChain==WWTo2L2Nu||iChain==ZZTo4L||iChain==WZTo3LNu) continue;//not using these in this analysis
    if(iChain==QCD20to30||iChain==QCD30to50||iChain==QCD50to80||iChain==QCD80to120||iChain==QCD120to170||
       iChain==QCD170to300||iChain==QCD300toInf) continue;//skipping QCD due to possible problems
    chains[iChain] = new TChain(treeName);
    if(iChain==DataRunB||iChain==DataRunC||iChain==DataRunD||iChain==DataRunE||iChain==DataRunF||
       iChain==DataRunG||iChain==DataRunH) isMC = kFALSE;
    else isMC = kTRUE;
    if(!(iChain==EE50to100))continue;
    nFiles = 0;
    subDirectorySize = subFiles[iChain]->size();
    for(int k=0;k<subDirectorySize;k++){
      if(nFiles>=maxFiles) break;
      for(int j=0;j<10;j++) {	  	      
	if(nFiles>=maxFiles) break;
	files = baseDirectory;      
	if(iChain==EE10to50||iChain==EE50to100||iChain==EE100to200||iChain==EE200to400||
	   iChain==EE400to500||iChain==EE500to700||iChain==EE700to800||iChain==EE800to1000||
	   iChain==EE1000to1500||iChain==EE1500to2000||iChain==EE2000to3000||
	   iChain==TAUTAU10to50||iChain==TAUTAU50to100||iChain==TAUTAU100to200||iChain==TAUTAU200to400||
	   iChain==TAUTAU400to500||iChain==TAUTAU500to700||iChain==TAUTAU700to800||iChain==TAUTAU800to1000||
	   iChain==TAUTAU1000to1500||iChain==TAUTAU1500to2000||iChain==TAUTAU2000to3000) {
	  files+=subDirectoryDY;	 
      }      
	
      files+=subFiles[iChain]->at(k);
      files+="/skims_0002/";      
      files+="ntuple_skim_";
      files+=j;
      files+=".root";
      std::ifstream testFileStream(files);
      if(!(bool)testFileStream) continue;
      chains[iChain]->Add(files);
      nFiles++;
      cout << files << endl;
      cout << chains[iChain]->GetEntries() << " events loaded" << endl;	 
      }
      if(chains[iChain]->GetEntries()==0){//error message if no events loaded
	cout << endl;
	cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	cout << "ERROR: Broken files or files not found in: " << endl;
	cout << files << endl;
	cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	cout << endl;
	return;
      }
    }                
    totalentries+=chains[iChain]->GetEntries(); 

    //Setting addresses for branches
    chains[iChain]->SetBranchAddress("Nelectrons", &Nelectrons, &b_Nelectrons);
    chains[iChain]->SetBranchAddress("nVertices", &nVertices, &b_nVertices);
    chains[iChain]->SetBranchAddress("nPileUp", &nPileUp, &b_nPileUp);
    chains[iChain]->SetBranchAddress("Electron_pT", &Electron_pT, &b_Electron_pT);
    chains[iChain]->SetBranchAddress("Electron_eta",&Electron_eta, &b_Electron_eta);
    chains[iChain]->SetBranchAddress("Electron_phi",&Electron_phi, &b_Electron_phi);
    chains[iChain]->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,&b_Electron_passMediumID);
    chains[iChain]->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
    chains[iChain]->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
    chains[iChain]->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
    chains[iChain]->SetBranchAddress("HLT_trigName",&pHLT_trigName);   
    if(isMC) 
      chains[iChain]->SetBranchAddress("GENEvt_weight",&GENEvt_weight,&b_GENEvt_weight);
    else continue;
  }//end iChain loop
  
  cout << "Total Events Loaded: " << totalentries << endl;
  cout << endl;

  TFile*file = new TFile("./plots/pileup.root");
  TH1F*hDatanPileup = (TH1F*)file->Get("hPileupData");
  hDatanPileup->SetLineColor(kBlack);
  TH1F*hMCnPileup = new TH1F("hMCnPileup","",nPileupBins,pileupBinLow,pileupBinHigh);
  hMCnPileup->SetMarkerStyle(20);
  hMCnPileup->SetMarkerColor(kRed);
  hMCnPileup->SetLineColor(kRed);
  hMCnPileup->SetTitle("Pileup Unweighted");
  TH1F*hMCnPileupWeight = new TH1F("hMCnPileupWeight","",nPileupBins,pileupBinLow,pileupBinHigh);
  hMCnPileupWeight->SetMarkerStyle(20);
  hMCnPileupWeight->SetMarkerColor(kRed); 
  hMCnPileupWeight->SetLineColor(kRed);
  hMCnPileupWeight->SetTitle("Pileup with Weights");

  TFile*pileupRatioFile  = new TFile(pileupRatioName);
  TH1F*hPileupRatio = (TH1F*)pileupRatioFile->Get("hPileupRatio");
  
  cout << "Starting Event Loop" << endl;
  double invMass, rapidity, dileptonPt, dileptonEta, xSecWeight, weightNoPileup, genWeight, varGenWeight, 
    totalWeight, lumiEffective, nEffective, localEntry, sumGenWeight, sumRawGenWeight, pileupWeight;
  Long64_t nentries;
  Long64_t count = 0;
  int sampleCategory = -1;
  TString compareHLT = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
  TString trigName;
  int trigNameSize, subEle, leadEle;
  double lumi = dataLuminosity;//luminosity for xsec weighting

  //Loop over samples
  for(int iChain=0;iChain<numChains;iChain++) {
    if(iChain==WWTo2L2Nu||iChain==ZZTo4L||iChain==WZTo3LNu) continue;//not using these in this analysis
    if(iChain==QCD20to30||iChain==QCD30to50||iChain==QCD50to80||iChain==QCD80to120||iChain==QCD120to170||
       iChain==QCD170to300||iChain==QCD300toInf) continue;//skipping QCD due to possible problems
    if(iChain==DataRunB||iChain==DataRunC||iChain==DataRunD||iChain==DataRunE||iChain==DataRunF||
       iChain==DataRunG||iChain==DataRunH) isMC = kFALSE; //determine if chain is MC or Data
    else isMC = kTRUE;
     
    //Event loop
    nentries = chains[iChain]->GetEntries();   
    for(Long64_t i=0;i<nentries;i++) {      
      counter(count,totalentries);
      count = count+1; 
      chains[iChain]->GetEntry(i);	  	  

      pileupWeight = hPileupRatio->GetBinContent(hPileupRatio->FindBin(nPileUp));
      genWeight = GENEvt_weight/fabs(GENEvt_weight);
      genWeight = genWeight/sumGenWeight;      
      if(isMC) {//weights for MC 
	totalWeight = genWeight*xSecWeight*pileupWeight;      
	weightNoPileup = genWeight*xSecWeight;
      }
      else {//no weights for data
	totalWeight = 1.0;      
	weightNoPileup = 1.0;
      }
      if(isMC){
	hMCnPileup->Fill(nPileUp);
	hMCnPileupWeight->Fill(nPileUp,pileupWeight);
      }
    }//end event loop   
  }//end chain loop 

  double norm=1.0;
  double scaleData = norm/hDatanPileup->Integral();
  double scaleMC = norm/hMCnPileup->Integral();
  double scaleMCWeight = norm/hMCnPileupWeight->Integral();
  hDatanPileup->Scale(scaleData);
  hMCnPileup->Scale(scaleMC);
  hMCnPileupWeight->Scale(scaleMCWeight);
  for(int i=1;i<76;i++){
    double binN = i;
    double valueBefore = hMCnPileup->GetBinContent(i);
    double valueAfter = hMCnPileupWeight->GetBinContent(i);
    double valueData = hDatanPileup->GetBinContent(i);
  }

  hDatanPileup->GetXaxis()->SetRangeUser(0,75);
  hMCnPileup->GetXaxis()->SetRangeUser(0,75);
  hMCnPileupWeight->GetXaxis()->SetRangeUser(0,75);

  TCanvas*c=new TCanvas("c","",10,10,1500,700);
  hDatanPileup->Draw();
  TLegend*legend=new TLegend(0.65,0.9,0.9,0.7);
  legend->SetTextSize(0.02);
  legend->AddEntry(hDatanPileup,"Data");
  legend->AddEntry(hMCnPileup,"MC");
  c->Divide(2);
  c->cd(1);
  hDatanPileup->Draw("hist");
  hMCnPileup->Draw("P,same");
  legend->Draw("same");
  c->cd(2);
  hDatanPileup->Draw("hist");
  hMCnPileupWeight->Draw("P,same");
    
  c->SaveAs("./plots/pileupClosureTest.png");
  totaltime.Stop();
  Double_t TotalCPURunTime = totaltime.CpuTime();
  Double_t TotalRunTime = totaltime.RealTime();
  TTimeStamp ts_end;
  cout << endl;
  cout << "*****************************************************************************" << endl;
  cout << "Total CPU RunTime: " << TotalCPURunTime/60 << " minutes" << endl;
  cout << "Total Real RunTime: " << TotalRunTime/60 << " minutes" << endl;
  cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;   
  cout << "Number of Events Processed: " << count << endl;
  cout << "*****************************************************************************" << endl;
  cout << endl;
  
}//end main function

//Counter for tracking program progress
void counter(Long64_t i, Long64_t N)
{
  int P = 100*(i)/(N);  
  TTimeStamp eventTimeStamp;
  if(i%(N/100)==0) {
    cout << "dataVsMC.C " << "[Time: " << eventTimeStamp.AsString("s") << "] " << P << "%" << endl;
  }
  return;
}
