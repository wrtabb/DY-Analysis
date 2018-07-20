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
#include "TFileCollection.h"
#include "TRatioPlot.h"

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
const int numChains = 37;
//Cross sections obtained from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SNUCMSYooDYntuple
const float xSec[numChains] = {5352960,9928000,2890800,350000,62964,18810,1350,//QCD
			       61526.7,118.7,12.178,16.523,1.256,47.13,4.4297,//Bosons
			       831.76,76.605,20.578,35.85,35.85,//tops
			       18810.0,5705.9044344,226.6,7.77,0.4065,0.2334,0.03614,//DY
			       0.03047,0.01636,0.00218,0.0005156,//DY
			       1,1,1,1,1,1,1};//data (unweighted)
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

  //Loading ntuples
  cout << "Loading ntuples" << endl;
  
  enum chainNum 
  {
    QCD20to30 = 0,
    QCD30to50 = 1,
    QCD50to80 = 2,
    QCD80to120 = 3,
    QCD120to170 = 4,
    QCD170to300 = 5,
    QCD300toInf = 6,  
    wJets = 7,
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

TString dirNames[numChains]=
  {//The names of every directory being loaded
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
    "/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttbar",
    "/TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8",
    "/TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8",
    "/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
    "/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
    "/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50",
    "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/skims_0001",
    "/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200",
    "/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "/DoubleEG/crab_DoubleEG_RunB",
    "/DoubleEG/crab_DoubleEG_RunC",
    "/DoubleEG/crab_DoubleEG_RunD",
    "/DoubleEG/crab_DoubleEG_RunE",
    "/DoubleEG/crab_DoubleEG_RunF",
    "/DoubleEG/crab_DoubleEG_RunG",
    "/DoubleEG/crab_DoubleEG_RunHver"
  };
 
 TString baseDirectory = 
   "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/ikrav/DrellYan_13TeV_2016/v2p3"; 
 TChain*chains[numChains];
 vector <TString> *subFiles[numChains];  
 for(int iChain=0;iChain<numChains;iChain++)
   {
     if(iChain==tt700to1000||iChain==tt1000toInf) continue;//skip these files for now
      if(iChain==WWTo2L2Nu||iChain==ZZTo4L||iChain==WZTo3Lnu) continue;//skip these files for now
     subFiles[iChain] = new vector<TString>;
     if(iChain==QCD30to50||iChain==QCD50to80||iChain==QCD80to120||iChain==QCD120to170) 
       {
	 subFiles[iChain]->push_back(dirNames[iChain]);
	 subFiles[iChain]->push_back(dirNames[iChain]+"_ext1");
       }
     else if(iChain==MC10to50) 
       {
	 subFiles[iChain]->push_back(dirNames[iChain]+"_ext1v1");
	 subFiles[iChain]->push_back(dirNames[iChain]+"_v1");
	 subFiles[iChain]->push_back(dirNames[iChain]+"_v2");
       }
     else if(iChain==MC100to200) 
       {
	 subFiles[iChain]->push_back(dirNames[iChain]);
	 subFiles[iChain]->push_back(dirNames[iChain]+"_ext");
       }
     else if(iChain==DataRunH) 
       {
	 subFiles[iChain]->push_back(dirNames[iChain]+"2");
	 subFiles[iChain]->push_back(dirNames[iChain]+"3");
       }
     else subFiles[iChain]->push_back(dirNames[iChain]);      
   }
 
 
 TString files;  
 TFileCollection filecoll("dum");//Object for creating a list of files in a directory
 Long64_t subDirectorySize;
 Long64_t totalentries = 0;
 for(int iChain=0;iChain<numChains;iChain++)
    {
      chains[iChain] = new TChain(treeName);
      subDirectorySize = subFiles[iChain]->size();    
      if(iChain==tt700to1000||iChain==tt1000toInf) continue;//skip these files for now
      if(iChain==WWTo2L2Nu||iChain==ZZTo4L||iChain==WZTo3Lnu) continue;//skip these files for now
      for(int k=0;k<subDirectorySize;k++)
	{	  	      
	  files = baseDirectory;
	  files+=subFiles[iChain]->at(k);
	  files+="/skims_0001/*.root";
	  filecoll.Add(files);
	  chains[iChain]->AddFileInfoList(filecoll.GetList());
	  cout << files << endl;
	  cout << chains[iChain]->GetEntries() << " events loaded" << endl;
	}                
	  
      //Setting addresses for branches
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
  double invMass, weight;
  Long64_t nentries;
  Long64_t count = 0;
  TString compareHLT = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
  TString trigName;
  int trigNameSize;
  //nentries = 10000;
  double lumi = chains[MC50to100]->GetEntries()/xSec[MC50to100]; //luminosity of 50to100
  //double lumi = nentries/xSec[MC50to100];
  for(int iChain=0;iChain<numChains;iChain++)
    {
      if(iChain==tt700to1000||iChain==tt1000toInf) continue;//skip these files for now
      if(iChain==WWTo2L2Nu||iChain==ZZTo4L||iChain==WZTo3Lnu) continue;//skip these files for now
      //nentries = 10000;
      nentries = chains[iChain]->GetEntries();
      //if(chains[iChain]->GetEntries() < nentries) nentries = chains[iChain]->GetEntries();
      weight=lumi*(xSec[iChain]/nentries);
      for(Long64_t i=0;i<nentries;i++)
	{      
	  chains[iChain]->GetEntry(i);
	  counter(count,totalentries);
	  //counter(count,7*nentries);
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
		  
		  if(iChain==MC50to100 && invMass>100) continue;
		  if(!passDileptonKinematics(Electron_pT[iEle],Electron_pT[jEle],Electron_eta[iEle],
					     Electron_eta[jEle])) continue; 		  		 
		  if(!passHLT) continue;

		  if(!Electron_passMediumID[iEle]) continue;//iLep electron ID cut
		  if(!Electron_passMediumID[jEle]) continue;//jLep electron cut

		  if(iChain==DataRunB||iChain==DataRunC||iChain==DataRunD||iChain==DataRunE||iChain==DataRunF||
		     iChain==DataRunG||iChain==DataRunH)
		    {
		      hDataInvMass->Fill(invMass);
		      hDataInvMasslinear->Fill(invMass);
		    }		  
		  else if(iChain==wJets||iChain==QCD20to30||iChain==QCD30to50||iChain==QCD50to80
			  ||iChain==QCD80to120||iChain==QCD120to170||iChain==QCD170to300||iChain==QCD300toInf) 
		    {
		      hFakes->Fill(invMass,weight);
		      hFakeslinear->Fill(invMass,weight);
		    }
		  else if(iChain==WW||iChain==ZZ||iChain==WZ||iChain==WWTo2L2Nu||iChain==ZZTo4L||iChain==WZTo3LNu) 
		    {
		      hEW->Fill(invMass,weight);
		      hEWlinear->Fill(invMass,weight);
		    }
		  else if(iChain==tt||iChain==tW||iChain==tbarW) 
		    {
		      hTops->Fill(invMass,weight);
		      hTopslinear->Fill(invMass,weight);
		    }		  
		  else if(iChain==MC10to50||iChain==MC50to100||iChain==MC100to200||iChain==MC200to400||
			  iChain==MC400to500||iChain==MC500to700||iChain==MC700to800||iChain==MC800to1000||
			  iChain==MC1000to1500||iChain==MC1500to2000||iChain==MC2000to3000)
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
