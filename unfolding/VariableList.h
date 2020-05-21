#ifndef VariableList_H
#define VariableList_H

#include "NtuplesV2P6Location.h"
#include <TLorentzVector.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

//-----File locations and tree name-----//
const TString treeName = "recoTree/DYTree";
const TString pileupRatioName = "/home/hep/wrtabb/DY-Analysis/data/pileup/pileup.root";
const TString leg2SFName = "/home/hep/wrtabb/DY-Analysis/data/SFs/Leg2_SF.root";
const TString medIDSFName = "/home/hep/wrtabb/DY-Analysis/data/SFs/MediumID_SF.root";
const TString recoSFName = "/home/hep/wrtabb/DY-Analysis/data/SFs/Reco_SF.root";
const TString pvzFileName = "/home/hep/wrtabb/DY-Analysis/data/PVz.root";

//-----Cut parameters-----//
const double etaGapLow = 1.4442;
const double etaGapHigh = 1.566;
const double etaHigh = 2.4;
const double ptLow = 17;
const double ptHigh = 28;
const float dRMinCut = 0.3;

//-----Branch variables-----//
const int MPSIZE = 2000;
int GENnPair;
int Nelectrons;
int HLT_ntrig;
int nPileUp;
double GENEvt_weight;
double GENLepton_phi[MPSIZE];
double GENLepton_eta[MPSIZE];
double GENLepton_pT[MPSIZE];
double GENLepton_Px[MPSIZE];
double GENLepton_Py[MPSIZE];
double GENLepton_Pz[MPSIZE];
double GENLepton_E[MPSIZE];
int GENLepton_ID[MPSIZE];
int GENLepton_isHardProcess[MPSIZE];
int GENLepton_fromHardProcessFinalState[MPSIZE];

double Electron_pT[MPSIZE]; 
double Electron_eta[MPSIZE];
double Electron_phi[MPSIZE];
double Electron_Energy[MPSIZE];
double Electron_Px[MPSIZE];
double Electron_Py[MPSIZE];
double Electron_Pz[MPSIZE];
double Electron_charge[MPSIZE];
bool Electron_passMediumID[MPSIZE];

int HLT_trigType[MPSIZE];
int HLT_trigFired[MPSIZE];
TString triggerUsed = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
TString trigName;
std::vector<std::string> HLT_trigName;
std::vector<std::string>*pHLT_trigName = &HLT_trigName;

const double pi=TMath::Pi();

const double eMass = 0.000511;
const double muMass = 0.105658;
const double tauMass = 1.7769;
const int dataLuminosity = 35867; //Run2016B to Run2016H JSON. unit: /pb, Updated at 2017.07.30
const int ptBinHigh = 499;
const int ptBinLow = 26;
int nVertices;

const int nYBins = 50;
const int nYBins2 = 2*nYBins;
const float binLowY = -2.4;
const float binHighY = 2.4;

const int nLogBinsMass = 43;
const int nLogBinsMass2 = 2*nLogBinsMass;
const int nLogBinsMass3 = nLogBinsMass+1;
const double massbins[] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,
 106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,
 1500,3000};
const double massbins2[] = {15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,
 50,52.5,55,57.5,60,62,64,66,68,70,72,74,76,78.5,81,83.5,86,88.5,91,93.5,96,98.5,101,103.5,
 106,108,110,112.5,115,117.5,120,123,126,129.5,133,137,141,145.5,150,155,160,165.5,171,178,
 185,192.5,200,210,220,231.5,243,258,273,296.5,320,350,380,410,440,475,510,555,600,650,700,
 765,830,915,1000,1250,1500,2250,3000};
const double massbins3[] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,
 106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,
 1500,2250,3000};

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/SNUCMSYooDYntuple (pb)
std::vector<double> xSec;
std::vector<double> xSecLL = {
 18610.0/3.0,//10to50 (NLO), NNLO value currently has huge uncertainty under investigation
 1923.26,//50to100 (NNLO)
 78.1258,//100to200 (NNLO)
 2.73309,//200to400 (NNLO)
 0.142945,//400to500 (NNLO)
 0.0809755,//500to700 (NNLO)
 0.0125589,//700to800 (NNLO)
 0.0105845,//800to1000 (NNLO)
 0.00556507,//1000to1500 (NNLO)
 0.000730495,//1500to2000 (NNLO)
 0.00016844//2000to3000 ((NNLO)
};

std::vector<double> xSecTT = {
 35.85,//ST_tW
 35.85,//ST_tbarW
 728.74,//ttbar_M0to700
 76.605,//ttbar_M700to1000
 20.578//ttbar_M1000toInf
};


std::vector<double> xSecEW = {
 118.7,//WW
 47.13,//WZ
 16.523,//ZZ
 18610.0/3.0,//10to50 (NLO), NNLO value currently has huge uncertainty under investigation
 1923.26,//50to100 (NNLO)
 78.1258,//100to200 (NNLO)
 2.73309,//200to400 (NNLO)
 0.142945,//400to500 (NNLO)
 0.0809755,//500to700 (NNLO)
 0.0125589,//700to800 (NNLO)
 0.0105845,//800to1000 (NNLO)
 0.00556507,//1000to1500 (NNLO)
 0.000730495,//1500to2000 (NNLO)
 0.00016844//2000to3000 ((NNLO)
};

std::vector<double> xSecData = {
 1,1,1,1,1,1,1,1
};

std::vector<double> xSecFakes = {
 61526.7//WJetsToLNu (NNLO)
};

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
TBranch*b_nVertices;

std::vector<TString> dirNamesLL = {
 DYLL_M10to50,   
 DYLL_M50to100,
 DYLL_M100to200,
 DYLL_M200to400,
 DYLL_M400to500,
 DYLL_M500to700,
 DYLL_M700to800,
 DYLL_M800to1000, 
 DYLL_M1000to1500,
 DYLL_M1500to2000,
 DYLL_M2000to3000
};

std::vector<TString> dirNamesEW = {
 WW_dir,
 WZ_dir,
 ZZ_dir,
 DYLL_M10to50,   
 DYLL_M50to100,
 DYLL_M100to200,
 DYLL_M200to400,
 DYLL_M400to500,
 DYLL_M500to700,
 DYLL_M700to800,
 DYLL_M800to1000, 
 DYLL_M1000to1500,
 DYLL_M1500to2000,
 DYLL_M2000to3000
};

std::vector<TString> dirNamesFakes = {
 WJets
};

std::vector<TString> dirNamesTT = {
 ST_tW,
 ST_tbarW,
 ttbar_M0to700,
 ttbar_M700to1000,
 ttbar_M1000toInf
};

std::vector<TString> dirNamesData = {
 DoubleEG_RunB,
 DoubleEG_RunC,
 DoubleEG_RunD,
 DoubleEG_RunE,
 DoubleEG_RunF,
 DoubleEG_RunG,
 DoubleEG_RunH,
};

std::vector<TString> dirNamesSM = {
 SM_2016B,
 SM_2016C,
 SM_2016D,
 SM_2016E,
 SM_2016F,
 SM_2016G,
 SM_2016H,
};

std::vector<TString> dirNamesZtoEE = {
 ZToEE_M50to120,
 ZToEE_M120to200,
 ZToEE_M200to400,
 ZToEE_M400to800,
 ZToEE_M800to1400,
 ZToEE_M1400to2300,
 ZToEE_M2300to3500,
 ZToEE_M3500to4500,
 ZToEE_M4500to6000,
 ZToEE_M6000toInf
};

enum SampleType{
 FAKES,
 EW,
 TT,
 LL,
 DATA
};

enum LepType{
 ELE,
 MUON,
 TAU
};

enum ChainLL{
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
 EE2000to3000
};

enum ChainEW{
 WW,
 WZ,
 ZZ,
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
 TAUTAU2000to3000
};
enum ChainFakes{
 W_PLUS_JETS
};

enum ChainTT{
 TW,
 T_BAR_W,
 TT0to700,
 TT700to1000,
 TT1000toInf
};

enum ChainData{
 RUN_B,
 RUN_C,
 RUN_D,
 RUN_E,
 RUN_F,
 RUN_G,
 RUN_H
};
enum BinType{
 LOG,
 LINEAR
};

#endif
