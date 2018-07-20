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
int GENnPair, Nelectrons, HLT_ntrig;
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
const double massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 
			     106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 
			     380, 440, 510, 600, 700, 830, 1000, 1500, 3000};
const double pi=TMath::Pi();
const int numChains = 26;
//Cross sections obtained from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SNUCMSYooDYntuple
const float xSec[numChains] = {61526,5352960,9928000,2890800,350000,62964,18810,1350,118.7,16.523,47.13,
			       831.76,35.85,35.85,18810.0,5705.9044344,226.6,7.77,0.4065,0.2334,0.03614,
			       0.03047,0.01636,0.00218,0.0005156,1};
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
const TString whichChain[numChains] = 
    {"W + Jets","QCD pT 20 to 30","QCD pT 30 to 50","QCD pT 50 to 80","QCD pT 80 to 120","QCD pT 120 to 170",
     "QCD pT 170 to 300","QCD pT 300 to infinity","WW","ZZ","WZ","TT","tW","anti tW","MC 10 to 50","MC 50 to 100",
     "MC 100 to 200","MC 200 to 400","MC 400 to 500","MC 500 to 700","MC 700 to 800","MC 800 to 1000",
     "MC 1000 to 1500","MC 1500 to 2000","MC 2000 to 3000","Data"};

void dataDYtoLL()
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
    WWTo2L2Nu = 9,
    ZZ = 10,       
    ZZTo4L = 11,
    WZ = 12,
    WZTo3LNu = 13,
    tt = 14,
    tt700to1000 = 15,
    tt1000toInf = 16,
    tW = 17,             
    tbarW = 18,         
    MC10to50 = 19,
    MC50to100 = 20,
    MC100to200 = 21,
    MC200to400 = 22,
    MC400to500 = 23,
    MC500to700 = 24,
    MC700to800 = 25,
    MC800to1000 = 26,
    MC1000to1500 = 27,
    MC1500to2000 = 28,
    MC2000to3000 = 29,
    DataRunB = 30,
    DataRunC = 31,
    DataRunD = 32,
    DataRunE = 33,
    DataRunF = 34,
    DataRunG = 35,
    DataRunH = 36
  };  
  
    TString subDirectory[numChains] = 
    {      
      "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_amcatnlo",
      "QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8",
      "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/",
      "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/",
      "QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/",
      "QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/",
      "QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/",
      "QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/",
      "WW_TuneCUETP8M1_13TeV-pythia8",
      "ZZ_TuneCUETP8M1_13TeV-pythia8",
      "WZ_TuneCUETP8M1_13TeV-pythia8",
      "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/",
      "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
      "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
      "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_",
      "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skims/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_skim.root",
      "DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200",
      "DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skims/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_skim.root",
      "DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skims/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_skim.root",
      "DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skims/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_skim.root",
      "DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skims/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_skim.root",
      "DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skims/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_skim.root",
      "DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skims/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_skim.root",
      "DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skims/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_skim.root",
      "DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skims/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_skim.root",
      "DoubleEG/crab_DoubleEG_RunB/ntuple_skim_"
    };

  TString sub10to50[3] =
    {
      "ext1v1/skims/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_crab_DYLL_M10to50_ext1v1_skim.root",
      "v1/skims/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_crab_DYLL_M10to50_v1_skim.root",
      "v2/skims/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_crab_DYLL_M10to50_v2_skim.root"
    };
  TString sub100to200[2] =
    {
      "/skims/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_crab_DYLL_M100to200_skim.root",
      "_ext/skims/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_crab_DYLL_M100to200_ext_skim.root"
    };
  TString subQCD30to50[2] =
    {
      "crab_QCDEMEnriched_Pt30to50",
      "crab_QCDEMEnriched_Pt30to50_ext1"
    };
  TString subQCD50to80[2] =
    {
      "crab_QCDEMEnriched_Pt50to80",
      "crab_QCDEMEnriched_Pt50to80_ext1"
    };
  TString subQCD80to120[2] =
    {
      "crab_QCDEMEnriched_Pt80to120",
      "crab_QCDEMEnriched_Pt80to120_ext1"
    };
  TString subQCD120to170[2] =
    {
      "crab_QCDEMEnriched_Pt120to170",
      "crab_QCDEMEnriched_Pt120to170_ext1"
    };
  TString subTT[2] =
    {
      "crab_ttbar",
      "crab_ttbarBackup"
    };
  TString inputBaseDirName =  
    "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/ikrav/DrellYan_13TeV_2016/v2p3/";  
  TChain*chains[numChains];
  TString files;  
  int nFiles;
  Long64_t totalentries = 0;
  for(int iChain=0;iChain<numChains;iChain++)
    {
      chains[iChain] = new TChain(treeName);
      nFiles = 0;
      for(int k=0;k<10000;k++)//Loop over individual files: ntuple_skim_k.root
	{	  
	  if(nFiles>=maxFiles) break;
	  if(iChain==wJets||iChain==QCD20to30||iChain==QCD170to300|iChain==QCD300toInf||iChain==ZZ||
	     iChain==WW||iChain==WZ||iChain==tW||iChain==tbarW)
	    {//wJets, QCD20to30, QCD170to300, QCD300toInf, ZZ, WW, WZ, tW, tbarW 
	      files = inputBaseDirName;
	      files+=subDirectory[iChain];
	      files+="/ntuple_skim_";
	      files+=k;
	      files+=".root";
	      std::ifstream testFileStream(files);
	      if(!(bool)testFileStream) continue;    	
	      chains[iChain]->Add(files); 	      
	      cout << "File: " << files << " loaded" << endl;	   
	      cout << "Entries: " << chains[iChain]->GetEntries() << endl;
	      nFiles++;	      
	      }  
	  if(iChain==tt)
	    {//tt 
	      for(int j=0;j<2;j++)
		{		  
		  if(nFiles>=maxFiles) break;
		  files = inputBaseDirName;
		  files+=subDirectory[iChain];
		  files+=subTT[j];
		  files+="/ntuple_skim_";
		  files+=k;
		  files+=".root";
		  std::ifstream testFileStream(files);
		  if(!(bool)testFileStream) continue;    	
		  chains[iChain]->Add(files); 		  
		  cout << "File: " << files << " loaded" << endl;
		  cout << "Entries: " << chains[iChain]->GetEntries() << endl;
		  nFiles++;		  
		}	  	    
	    }
	  if(iChain==QCD30to50)
	    {//QCD 30to50	    
	      for(int j=0;j<2;j++)
		{		  
		  files = inputBaseDirName;
		  files+=subDirectory[iChain];
		  files+=subQCD30to50[j];
		  files+="/ntuple_skim_";
		  files+=k;
		  files+=".root";
		  std::ifstream testFileStream(files);
		  if(!(bool)testFileStream) continue;    	
		  chains[iChain]->Add(files); 	       
		  cout << "File: " << files << " loaded" << endl;
		  cout << "Entries: " << chains[iChain]->GetEntries() << endl;
		  nFiles++;
		}	  	    
	    }
	  if(iChain==QCD50to80)
	    {//QCD 50to80	    
	      for(int j=0;j<2;j++)
		{		  
		  files = inputBaseDirName;
		  files+=subDirectory[iChain];
		  files+=subQCD50to80[j];
		  files+="/ntuple_skim_";
		  files+=k;
		  files+=".root";
		  std::ifstream testFileStream(files);
		  if(!(bool)testFileStream) continue;    	
		  chains[iChain]->Add(files); 		  
		  cout << "File: " << files << " loaded" << endl;
		  cout << "Entries: " << chains[iChain]->GetEntries() << endl;
		  nFiles++;
		}	  	    
	    }
	  if(iChain==QCD80to120)
	    {//QCD 80to120	    
	      for(int j=0;j<2;j++)
		{		  
		  files = inputBaseDirName;
		  files+=subDirectory[iChain];
		  files+=subQCD80to120[j];
		  files+="/ntuple_skim_";
		  files+=k;
		  files+=".root";
		  std::ifstream testFileStream(files);
		  if(!(bool)testFileStream) continue;    	
		  chains[iChain]->Add(files); 		  
		  cout << "File: " << files << " loaded" << endl;
		  cout << "Entries: " << chains[iChain]->GetEntries() << endl;
		  nFiles++;
		}	  	    
	    }
	  if(iChain==QCD120to170)
	    {//QCD 120to170	    
	      for(int j=0;j<2;j++)
		{		  
		  files = inputBaseDirName;
		  files+=subDirectory[iChain];
		  files+=subQCD120to170[j];
		  files+="/ntuple_skim_";
		  files+=k;
		  files+=".root";
		  std::ifstream testFileStream(files);
		  if(!(bool)testFileStream) continue;    	
		  chains[iChain]->Add(files); 		  
		  cout << "File: " << files << " loaded" << endl;
		  cout << "Entries: " << chains[iChain]->GetEntries() << endl;
		  nFiles++;
		}	  	    
	    }
	  if(iChain==Data)
	    {//Data sample
	      files = inputBaseDirName;
	      files+=subDirectory[iChain];
	      files += k;
	      files += ".root";
	      std::ifstream testFileStream(files);
	      if(!(bool)testFileStream) continue;    	      
	      chains[iChain]->Add(files);	      
	      cout << "File: " << files << " loaded" << endl;
	      cout << "Entries: " << chains[iChain]->GetEntries() << endl;
	      nFiles++;
	    }	
	}//end k loop over files
      	  if(iChain==MC10to50)
	    {//MC Sample 10to50
	      for(int j=0;j<3;j++)
		{		  
		  files = inputBaseDirName;
		  files+=subDirectory[iChain];
		  files+=sub10to50[j];
		  std::ifstream testFileStream(files);
		  if(!(bool)testFileStream) continue;    	
		  chains[iChain]->Add(files); 
		  cout << "File: " << files << " loaded" << endl;
		  cout << "Entries: " << chains[iChain]->GetEntries() << endl;
		}	  	      
	    }
	  if(iChain==MC100to200)
	    {//MC sample 100to200
	      for(int j=0;j<2;j++)
		{		  
		  files = inputBaseDirName;
		  files+=subDirectory[iChain];
		  files+=sub100to200[j];
		  std::ifstream testFileStream(files);
		  if(!(bool)testFileStream) continue;    	
		  chains[iChain]->Add(files); 
		  cout << "File: " << files << " loaded" << endl;
		  cout << "Entries: " << chains[iChain]->GetEntries() << endl;
		}  	      
	    }    
	  if(iChain==MC50to100||iChain==MC200to400||iChain==MC400to500||iChain==MC500to700||iChain==MC700to800||
	     iChain==MC800to1000||iChain==MC1000to1500||iChain==MC1500to2000||iChain==MC2000to3000)
	    {//MC samples other than 10to50 and 100to200
	      files = inputBaseDirName;
	      files += subDirectory[iChain];
	      std::ifstream testFileStream(files);
	      if(!(bool)testFileStream) continue;    	      
	      chains[iChain]->Add(files);
	      cout << "File: " << files << " loaded" << endl;
	      cout << "Entries: " << chains[iChain]->GetEntries() << endl;
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
  

  Long64_t dataEntries = chains[14]->GetEntries();
  cout << "Total Events Loaded: " << totalentries << endl;
  cout << endl;
    
  //Defining histograms
  TH1F*hMCInvMass = new TH1F("hMCInvMass","",nLogBins,massbins);
  hMCInvMass->Sumw2();
  hMCInvMass->SetFillColor(kOrange-2);
  hMCInvMass->SetLineColor(kOrange+3);
  hMCInvMass->GetXaxis()->SetTitle("m_{ee} [GeV]");
  hMCInvMass->GetXaxis()->SetMoreLogLabels();
  hMCInvMass->GetXaxis()->SetNoExponent();
  hMCInvMass->SetTitle("MC vs. Data Invariant Mass");
  hMCInvMass->SetMinimum(axisLow);
  TH1F*hMCInvMasslinear = new TH1F("hMCInvMasslinear","",nLinearBins,binLow,binHigh);
  hMCInvMasslinear->Sumw2();
  hMCInvMasslinear->SetFillColor(kOrange-2);
  hMCInvMasslinear->SetLineColor(kOrange+3);
  hMCInvMasslinear->GetXaxis()->SetTitle("m_{ee} [GeV]");
  hMCInvMasslinear->GetXaxis()->SetMoreLogLabels();
  hMCInvMasslinear->GetXaxis()->SetNoExponent();
  hMCInvMasslinear->SetTitle("MC vs. Data Invariant Mass");
  hMCInvMasslinear->SetMinimum(axisLow);
  TH1F*hDataInvMass = new TH1F("hDataInvMass","",nLogBins,massbins);
  hDataInvMass->Sumw2();
  hDataInvMass->SetLineColor(kBlack);
  hDataInvMass->SetMarkerColor(kBlack);
  hDataInvMass->SetMarkerSize(1);
  hDataInvMass->SetMarkerStyle(20);
  hDataInvMass->GetXaxis()->SetTitle("m_{ee} [GeV]");
  hDataInvMass->GetXaxis()->SetNoExponent();
  hDataInvMass->GetXaxis()->SetMoreLogLabels();
  hDataInvMass->SetMinimum(axisLow);
  TH1F*hDataInvMasslinear = new TH1F("hDataInvMasslinear","",nLinearBins,binLow,binHigh);
  hDataInvMasslinear->Sumw2();
  hDataInvMasslinear->SetLineColor(kBlack);
  hDataInvMasslinear->SetMarkerColor(kBlack);
  hDataInvMasslinear->SetMarkerSize(1);
  hDataInvMasslinear->SetMarkerStyle(20);
  hDataInvMasslinear->GetXaxis()->SetTitle("m_{ee} [GeV]");
  hDataInvMasslinear->GetXaxis()->SetNoExponent();
  hDataInvMasslinear->GetXaxis()->SetMoreLogLabels();
  hDataInvMasslinear->SetMinimum(axisLow);
  TH1F*hFakes = new TH1F("hFakes","",nLogBins,massbins);
  //W+Jets and QCD
  //iChain 0-7
  hFakes->Sumw2();
  hFakes->SetFillColor(kViolet+5);
  hFakes->SetLineColor(kViolet+3);
  hFakes->SetMinimum(axisLow);
  TH1F*hFakeslinear = new TH1F("hFakeslinear","",nLinearBins,binLow,binHigh);
  hFakeslinear->Sumw2();
  hFakeslinear->SetFillColor(kViolet+5);
  hFakeslinear->SetLineColor(kViolet+3);
  hFakeslinear->SetMinimum(axisLow);
  TH1F*hEW = new TH1F("hEW","",43,massbins);
  //WW, ZZ, WZ
  //iChain 8,9,10
  hEW->Sumw2();
  hEW->SetFillColor(kRed+2);
  hEW->SetLineColor(kRed+4);
  hEW->SetMinimum(axisLow);
  TH1F*hEWlinear = new TH1F("hEWlinear","",nLinearBins,binLow,binHigh);
  hEWlinear->Sumw2();
  hEWlinear->SetFillColor(kRed+2);
  hEWlinear->SetLineColor(kRed+4);
  hEWlinear->SetMinimum(axisLow);
  TH1F*hTops = new TH1F("hTops","",nLogBins,massbins);
  //tt, tW
  //iChain 11,12,13
  hTops->Sumw2();
  hTops->SetFillColor(kBlue+2);
  hTops->SetLineColor(kBlue+3);
  hTops->SetMinimum(axisLow);
  TH1F*hTopslinear = new TH1F("hTopslinear","",nLinearBins,binLow,binHigh);
  hTopslinear->Sumw2();
  hTopslinear->SetFillColor(kBlue+2);
  hTopslinear->SetLineColor(kBlue+3);
  hTopslinear->SetMinimum(axisLow);

  //Defining stacks
  THStack*hStack = new THStack("hStack","");
  hStack->Add(hFakes);
  hStack->Add(hEW);
  hStack->Add(hTops);
  hStack->Add(hMCInvMass);  
  THStack*hStacklinear = new THStack("hStacklinear","");
  hStacklinear->Add(hFakeslinear);
  hStacklinear->Add(hEWlinear);
  hStacklinear->Add(hTopslinear);
  hStacklinear->Add(hMCInvMasslinear); 

  TFile *rootFile = new TFile("./plots/dataVsMC.root","RECREATE");

  //Event Loop
  cout << "Starting Event Loop" << endl;
  double dpT, invMass, weight, dRMin;
  int dRMinIndex;  
  Long64_t nentries;
  Long64_t count = 0;
  TString compareHLT = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
  TString trigName;
  int trigNameSize;
  //nentries = 100000;
  double lumi = chains[MC50to100]->GetEntries()/xSec[MC50to100]; //luminosity of 50to100
  //double lumi = nentries/xSec[MC50to100];
  for(int iChain=0;iChain<numChains;iChain++)
    {
      //nentries = 100000;
      nentries = chains[iChain]->GetEntries();
      //if(chains[iChain]->GetEntries() < nentries) nentries = chains[iChain]->GetEntries();
      weight=lumi*(xSec[iChain]/nentries);
      cout << endl;
      cout << "Now processing files from " << "'" << whichChain[iChain] << "'" << endl;
      for(Long64_t i=0;i<nentries;i++)
	{      
	  chains[iChain]->GetEntry(i);
	  counter(count,totalentries);
	  //counter(count,26*nentries);
	  count = count+1;      
	    
	  //HLT cut
	  trigNameSize = pHLT_trigName->size();
	  bool passHLT = kFALSE;	  
	  for(int iHLT=0;iHLT<trigNameSize;iHLT++)
	    {
	      trigName = pHLT_trigName->at(iHLT);
	      
	      if(trigName.CompareTo(compareHLT)==0)
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
	    } 

	  //Electron loop
	  for(int iEle = 0; iEle < Nelectrons; iEle++)
	    {
	      for(int jEle = iEle+1; jEle < Nelectrons; jEle++)
		{	  
		  invMass=calcInvMass(Electron_pT[iEle],Electron_eta[iEle],Electron_phi[iEle],eMass,
				      Electron_pT[jEle],Electron_eta[jEle],Electron_phi[jEle],eMass);
		  
		  if(iChain==15 && invMass>100) continue;
		  if(!passDileptonKinematics(Electron_pT[iEle],Electron_pT[jEle],Electron_eta[iEle],
					     Electron_eta[jEle])) continue; 		  		 
		  if(!passHLT) continue;

		  if(!Electron_passMediumID[iEle]) continue;//iLep electron ID cut
		  if(!Electron_passMediumID[jEle]) continue;//jLep electron cut

		  if(iChain==Data)
		    {
		      hDataInvMass->Fill(invMass);
		      hDataInvMasslinear->Fill(invMass);
		    }
		  else if(iChain<WW) 
		    {
		      hFakes->Fill(invMass,weight);
		      hFakeslinear->Fill(invMass,weight);
		    }
		  else if(iChain>QCD300toInf&&iChain<tt) 
		    {
		      hEW->Fill(invMass,weight);
		      hEWlinear->Fill(invMass,weight);
		    }
		  else if(iChain>WZ&&iChain<MC10to50) 
		    {
		      hTops->Fill(invMass,weight);
		      hTopslinear->Fill(invMass,weight);
		    }		  
		  else if(iChain>tbarW&&iChain<Data)
		    {
		      hMCInvMass->Fill(invMass,weight);
		      hMCInvMasslinear->Fill(invMass,weight);
		    }
		}//end inner electron loop	   
	    }//end electron loop
	}//end event loop   
    }//end chain loop 

  double integralData, integralMC;
  
  integralData = 
    hDataInvMass->Integral(hDataInvMass->GetXaxis()->FindBin(binLow),hDataInvMass->GetXaxis()->FindBin(binHigh));
  integralMC = 
    hMCInvMass->Integral(hMCInvMass->GetXaxis()->FindBin(binLow),hMCInvMass->GetXaxis()->FindBin(binHigh))+
    hFakes->Integral(hFakes->GetXaxis()->FindBin(binLow),hFakes->GetXaxis()->FindBin(binHigh))+
    hEW->Integral(hEW->GetXaxis()->FindBin(binLow),hEW->GetXaxis()->FindBin(binHigh))+
    hTops->Integral(hTops->GetXaxis()->FindBin(binLow),hTops->GetXaxis()->FindBin(binHigh));

  double norm = integralData/integralMC;
  hMCInvMass->Scale(norm);
  hFakes->Scale(norm);
  hEW->Scale(norm);
  hTops->Scale(norm);

  hMCInvMasslinear->Scale(norm);
  hFakeslinear->Scale(norm);
  hEWlinear->Scale(norm);
  hTopslinear->Scale(norm);  

  TCanvas*canvas1 = new TCanvas("canvas1","",10,10,1000,1000);
  canvas1->SetLogx();
  canvas1->SetLogy();
  
  TLegend*legend = new TLegend(0.65,0.9,0.9,0.7);
  legend->SetTextSize(0.02);
  legend->AddEntry(hDataInvMass,"Data");
  legend->AddEntry(hMCInvMass,"#gamma^{*}/Z #rightarrow e^{-}e^{+}");
  legend->AddEntry(hTops,"t#bar{t}+tW+#bar{t}W");
  legend->AddEntry(hEW,"EW");
  legend->AddEntry(hFakes,"Fakes (Wjets+QCD");

  auto hDataMCRatio = new TRatioPlot(hStack,hDataInvMass);
  hDataMCRatio->GetXaxis()->SetTitle("m_{ee} [GeV]");  
  canvas1->cd();
  hDataMCRatio->Draw();
  hDataMCRatio->GetUpperPad()->cd();
  legend->Draw("same");
  canvas1->Update();
  canvas1->SaveAs("./plots/dataVsMClog.png");
  
  TCanvas*canvas2 = new TCanvas("canvas2","",10,10,1000,1000);
  canvas2->SetLogy();
  auto hDataMCRatiolinear = new TRatioPlot(hStacklinear,hDataInvMasslinear); 
  canvas2->cd();
  hDataMCRatiolinear->Draw();
  //TH1F*graph = (TH1F*)hDataMCRatiolinear->Clone();
  //graph->SetMinimum(0.001);
  //graph->Draw();
  hDataMCRatiolinear->GetUpperPad()->cd();
  legend->Draw("same");
 
  canvas2->SaveAs("./plots/dataVsMClinear.png");
  rootFile->cd();
  hStack->Write();
  hDataInvMass->Write();
  hMCInvMass->Write();
  hFakes->Write();
  hEW->Write();
  hTops->Write();
  hStacklinear->Write();
  hDataInvMasslinear->Write();
  hMCInvMasslinear->Write();
  hFakeslinear->Write();
  hEWlinear->Write();
  hTopslinear->Write();
  canvas1->Write();
  canvas2->Write();
  rootFile->Write();
  rootFile->Close();

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
  if(i%(N/100)==0)
    {
      cout << P << "%" <<  endl;
      cout << "[Time: " << eventTimeStamp.AsString("s") << "]" << endl;
    }
  return;
}

//Kinematic cuts
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2)
{
  //if(abs(eta1)>etaGapLow && abs(eta1)<etaGapHigh) return kFALSE;//eta cut
  //if(abs(eta2)>etaGapLow && abs(eta2)<etaGapHigh) return kFALSE; //eta cut
  if(abs(eta1)>etaHigh||abs(eta2)>etaHigh) return kFALSE; //eta cut
  if(!((pt1>ptLow && pt2>ptHigh)||(pt1>ptHigh && pt2>ptLow))) return kFALSE;
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

