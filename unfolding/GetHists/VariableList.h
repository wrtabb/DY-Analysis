#include "NtuplesV2P6Location.h"
#include <TMath.h>

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
std::vector<std::string> HLT_trigName;
std::vector<std::string>*pHLT_trigName = &HLT_trigName;

const double pi=TMath::Pi();

const float eMass = 0.000511;
double muMass = 0.105658;
const int dataLuminosity = 35867; //Run2016B to Run2016H JSON. unit: /pb, Updated at 2017.07.30
const int ptBinHigh = 499;
const int ptBinLow = 26;
int nVertices;

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

const TString treeName = "recoTree/DYTree";
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
 WJets,
};

std::vector<TString> dirNamesTT = {
 ST_tW,
 ST_tbarW,
 ttbar,
 ttbarBackup_M0to700,
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
 LL,
 EW,
 TT,
 DATA
};

enum LepType{
 ELE,
 MUON,
 TAU
};

enum ChainLL{
 M10to50,   
 M50to100,
 M100to200,
 M200to400,
 M400to500,
 M500to700,
 M700to800,
 M800to1000, 
 M1000to1500,
 M1500to2000,
 M2000to3000
};

enum ChainEW{
 WW,
 WZ,
 ZZ,
 W_PLUS_JETS
};

enum ChainTT{
 TW,
 T_BAR_W,
 TTbar,
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
