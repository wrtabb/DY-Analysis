//-----Header file for the Drell-Yan analysis---------------//
//-----Contains variables that were used in many places-----//

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
#include "TRandom3.h"
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
#include "TF1.h"
#include "TUnfold.h"
#include "TUnfoldDensity.h"
#include "TSpline.h"
#include "TMatrixD.h"

//-----Histogram binning for mass plots
int nLogBins = 43;
const int nLogBins2 = 2*nLogBins;
const double massbins[] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106, 
 110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1500,
 3000};
const double massbins2[] = {15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,
 57.5,60,62,64,66,68,70,72,74,76,78.5,81,83.5,86,88.5,91,93.5,96,98.5,101,103.5,106,108,110,
 112.5,115,117.5,120,123,126,129.5,133,137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,
 210,220,231.5,243,258,273,296.5,320,350,380,410,440,475,510,555,600,650,700,765,830,915,1000,
 1250,1500,2250,3000};

//-----File locations and tree name
const TString treeName = "recoTree/DYTree";
const TString pileupRatioName = "/home/hep/wrtabb/git/DY-Analysis/data/pileup/pileup.root";
const TString leg2SFName = "/home/hep/wrtabb/git/DY-Analysis/data/SFs/Leg2_SF.root";
const TString medIDSFName = "/home/hep/wrtabb/git/DY-Analysis/data/SFs/MediumID_SF.root";
const TString recoSFName = "/home/hep/wrtabb/git/DY-Analysis/data/SFs/Reco_SF.root";
const TString pvzFileName = "/home/hep/wrtabb/git/DY-Analysis/data/PVz.root";

//-----Branch variables
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
std::vector<std::string> HLT_trigName;
std::vector<std::string> *pHLT_trigName = &HLT_trigName;

const double pi=TMath::Pi();
//Cross sections obtained from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SNUCMSYooDYntuple
const int numChains = 11;
const float xSec[numChains] = {6016.88,1873.52,76.2401,2.67606,0.139728,0.0792496,0.0123176,
 0.01042,0.00552772,0.000741613,0.000178737};
//Electron mass
const float eMass = 0.000511;
//Muon mass
double muMass = 0.105658;

const int dataLuminosity = 35867; //Run2016B to Run2016H JSON. unit: /pb, Updated at 2017.07.30

const int ptBinHigh = 499;
const int ptBinLow = 26;
int nVertices;

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

