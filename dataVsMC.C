////////////////////////////////////////////////////////////////////
//Data versus MC signal and backgrounds
//Robert Tabb
//wrtabb@huskers.unl.edu
////////////////////////////////////////////////////////////////////

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

void counter(Long64_t i, Long64_t N);
double calcInvMass(double pt1,double eta1,double phi1,double m1,double pt2,double eta2,double phi2,double m2);
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2);

//Defining variables and arrays
const int MPSIZE = 2000;
int Nelectrons, HLT_ntrig;
double Electron_pT[MPSIZE], Electron_eta[MPSIZE], Electron_phi[MPSIZE];
double Electron_Energy[MPSIZE], Electron_Px[MPSIZE];
double Electron_Py[MPSIZE], Electron_Pz[MPSIZE], Electron_charge[MPSIZE];
bool Electron_passMediumID[MPSIZE];
int HLT_trigType[MPSIZE],HLT_trigFired[MPSIZE];
std::vector<std::string> HLT_trigName;
std::vector<std::string> *pHLT_trigName = &HLT_trigName;
const double massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 
			     106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 
			     380, 440, 510, 600, 700, 830, 1000, 1500, 3000};
const double pi=TMath::Pi();
const int numChains = 7;
//Cross sections obtained from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SNUCMSYooDYntuple
//const float xSec[numChains] = {61526,5352960,9928000,2890800,350000,62964,18810,1350,118.7,16.523,47.13,
//			       831.76,35.85,35.85,18810.0,5705.9044344,226.6,7.77,0.4065,0.2334,0.03614,
//			       0.03047,0.01636,0.00218,0.0005156,1};
const float binLow = 60.001;
const float binHigh = 119.999;
const float axisLow = 0.0001;
const int nLinearBins = 30;
const int nLogBins = 43;

const float etaHigh = 2.4;
const float etaGapHigh = 1.566; 
const float etaGapLow = 1.4442;
const float ptHigh = 28;
const float ptLow = 17;
const float eMass = 0.000511;

const TString treeName = "recoTree/DYTree";
const float dRMinCut = 0.3;
const int maxFiles = 200;
/*
const TString whichChain[numChains] = 
    {"W + Jets","QCD pT 20 to 30","QCD pT 30 to 50","QCD pT 50 to 80","QCD pT 80 to 120","QCD pT 120 to 170",
     "QCD pT 170 to 300","QCD pT 300 to infinity","WW","ZZ","WZ","TT","tW","anti tW","MC 10 to 50","MC 50 to 100",
     "MC 100 to 200","MC 200 to 400","MC 400 to 500","MC 500 to 700","MC 700 to 800","MC 800 to 1000",
     "MC 1000 to 1500","MC 1500 to 2000","MC 2000 to 3000","Data"};
*/
void dataVsMC()
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

  ////////////////////////////////////////
  //iChain           Process            //
  //iChain==0        W+Jets             //
  //iChain==1        QCD20to30          //
  //iChain==2        QCD30to50          //
  //iChain==3        QCD50to80          //
  //iChain==4        QCD80to120         //
  //iChain==5        QCD120to170        //
  //iChain==6        QCD170to300        //
  //iChain==7        QCD300AndUp        //
  //iChain==8        WW                 //
  //iChain==9        ZZ                 //
  //iChain==10       WZ                 // 
  //iChain==11       tt                 //
  //iChain==12       tW top             //
  //iChain==13       tW antitop         //
  //iChain==14to24   MC Signal          //
  //iChain==25       Data               //
  ////////////////////////////////////////
      
  //Loading ntuples
  cout << "Loading ntuples" << endl;
  /*
  enum chainNum 
  {
    wJets = 0,
    QCD20to30 = 1,
    QCD30to50 = 2,
    QCD50to80 = 3,
    QCD80to120 = 4,        
    QCD120to170 = 5,       
    QCD170to300 = 6,       
    QCD300toInf = 7,       
    WW = 8,                 
    ZZ = 9,                
    WZ = 10,                 
    tt = 11,                
    tW = 12,             
    tbarW = 13,         
    MC10to50 = 14,
    MC50to100 = 15,
    MC100to200 = 16,
    MC200to400 = 17,
    MC400to500 = 18,
    MC500to700 = 19,
    MC700to800 = 20,
    MC800to1000 = 21,
    MC1000to1500 = 22,
    MC1500to2000 = 23,
    MC2000to3000 = 24,
    DataRunB = 25,
    DataRunC = 26,
    DataRunD = 27,
    DataRunE = 28,
    DataRunF = 29,
    DataRunG = 30,
    DataRunH = 31
  };*/

    enum chainNum 
  {
    DataRunB = 0,
    DataRunC = 1,
    DataRunD = 2,
    DataRunE = 3,
    DataRunF = 4,
    DataRunG = 5,
    DataRunH = 6
  };

  TString baseDirectory =  
    "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/ikrav/DrellYan_13TeV_2016/v2p3"; 
  TString ntupleFile = "/*.root";
  TChain*chains[numChains];
  vector <TString> *subFiles[numChains];  
  for(int iChain=0;iChain<numChains;iChain++)
    {
      subFiles[iChain] = new vector<TString>;

      if(iChain==DataRunB) subFiles[iChain]->push_back("/DoubleEG/crab_DoubleEG_RunB/skims_0001");
      if(iChain==DataRunC) subFiles[iChain]->push_back("/DoubleEG/crab_DoubleEG_RunC/skims_0001");
      if(iChain==DataRunD) subFiles[iChain]->push_back("/DoubleEG/crab_DoubleEG_RunD/skims_0001");
      if(iChain==DataRunE) subFiles[iChain]->push_back("/DoubleEG/crab_DoubleEG_RunE/skims_0001");
      if(iChain==DataRunF) subFiles[iChain]->push_back("/DoubleEG/crab_DoubleEG_RunF/skims_0001");
      if(iChain==DataRunG) subFiles[iChain]->push_back("/DoubleEG/crab_DoubleEG_RunG/skims_0001");
      if(iChain==DataRunH) 
	{
	  subFiles[iChain]->push_back("/DoubleEG/crab_DoubleEG_RunHver2/skims_0001");
	  subFiles[iChain]->push_back("/DoubleEG/crab_DoubleEG_RunHver3/skims_0001");
	}
    }
  
  
  TString files;  
  int nFiles;
  TFileCollection filecoll("dum");
  Long64_t subDirectorySize;
  Long64_t totalentries = 0;
  for(int iChain=0;iChain<numChains;iChain++)
    {
      chains[iChain] = new TChain(treeName);
      nFiles = 0;
      subDirectorySize = subFiles[iChain]->size();      
      

      //cout << subDirectorySize << endl;
      for(int k=0;k<subDirectorySize;k++)
	{	  	      
	  files = baseDirectory;
	  files+=subFiles[iChain]->at(k);
	  files+=ntupleFile;
	  filecoll.Add(files);
	  chains[iChain]->AddFileInfoList(filecoll.GetList());
	  cout << files << endl;
	  cout << chains[iChain]->GetEntries() << " events loaded" << endl;
	}                
	  
      //Setting addresses for all branches
      chains[iChain]->SetBranchAddress("Nelectrons", &Nelectrons, &b_Nelectrons);
      chains[iChain]->SetBranchAddress("Electron_pT", &Electron_pT, &b_Electron_pT);
      chains[iChain]->SetBranchAddress("Electron_eta",&Electron_eta, &b_Electron_eta);
      chains[iChain]->SetBranchAddress("Electron_phi",&Electron_phi, &b_Electron_phi);
      chains[iChain]->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,&b_Electron_passMediumID);
      chains[iChain]->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
      chains[iChain]->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
      chains[iChain]->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
      chains[iChain]->SetBranchAddress("HLT_trigName",&pHLT_trigName);
      
      totalentries=totalentries+chains[iChain]->GetEntries(); 
      
    }//end iChain loop
  
  cout << "Total Events Loaded: " << totalentries << endl;
  cout << endl;
  
  totaltime.Stop();
  Double_t TotalCPURunTime = totaltime.CpuTime();
  Double_t TotalRunTime = totaltime.RealTime();
  TTimeStamp ts_end;
  cout << endl;
  cout << "*****************************************************************************" << endl;
  cout << "Total CPU RunTime: " << TotalCPURunTime/60 << " minutes" << endl;
  cout << "Total Real RunTime: " << TotalRunTime/60 << " minutes" << endl;
  cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;   
  //cout << "Number of Events Processed: " << count << endl;
  cout << "*****************************************************************************" << endl;
  cout << endl;
}
