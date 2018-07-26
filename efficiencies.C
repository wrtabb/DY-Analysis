
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
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
#include "TStyle.h"
#include "TEfficiency.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TFileCollection.h"
#include "THashList.h"

void counter(Long64_t i, Long64_t N);
double calcInvMass(double pt1,double eta1,double phi1,double m1,double pt2,double eta2,double phi2,double m2);
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2);
bool passPromptGenElectron(int ID, int fromfinalstate);
bool passHardProcess(int ID, int hardProces);
bool findGenToRecoMatch(int genIndex,int &recoIndex);

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
const int numChains = 11;
const int nLogBins = 43;
//Cross sections obtained from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SNUCMSYooDYntuple
const float xSec[numChains] = {18810.0,5705.9044344,226.6,7.77,0.4065,0.2334,
			       0.03614,0.03047,0.01636,0.00218,0.0005156};
const float etaHigh = 2.5;
const float etaGapHigh = 1.566; 
const float etaGapLow = 1.4442;
const float ptHigh = 28;
const float ptLow = 17;
const float eMass = 0.000511;
const TString treeName = "recoTree/DYTree";
const float dRMinCut = 0.3;
const int nSubSamples10to50 = 3;
const int nSubSamples100to200 = 2;

void efficiencies()
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
  TBranch*b_Electron_passMediumID;
  TBranch*b_HLT_ntrig;
  TBranch*b_HLT_trigType;
  TBranch*b_HLT_trigFired;

  //Loading ntuples
  cout << "Loading ntuples" << endl;
  
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

    TString dirNames[numChains] = 
      {      
	"/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M10to50_",
	"/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_truncated_M50To100",
	"/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYLL_M100to200",
	"/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
	"/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
	"/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
	"/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
	"/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
	"/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
	"/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
	"/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"
      }; 
 TString baseDirectory = 
   "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/ikrav/DrellYan_13TeV_2016/v2p3"; 
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
	  files+="/skims/*.root";	  
	  filecoll.Add(files);
	  chains[iChain]->AddFileInfoList(filecoll.GetList());
	  cout << files << endl;
	  cout << chains[iChain]->GetEntries() << " events loaded" << endl;	  
	}
  
      //Setting addresses for all branches
      chains[iChain]->SetBranchAddress("GENnPair", &GENnPair, &b_GENnPair);
      chains[iChain]->SetBranchAddress("GENLepton_eta", &GENLepton_eta, &b_GENLepton_eta);
      chains[iChain]->SetBranchAddress("GENLepton_phi",&GENLepton_phi, &b_GENLepton_phi);
      chains[iChain]->SetBranchAddress("GENLepton_pT",&GENLepton_pT, &b_GENLepton_pT);
      chains[iChain]->SetBranchAddress("GENLepton_ID",&GENLepton_ID, &b_GENLepton_ID);
      chains[iChain]->SetBranchAddress("GENLepton_isHardProcess",&GENLepton_isHardProcess, 
				       &b_GENLepton_isHardProcess);
      chains[iChain]->SetBranchAddress("GENLepton_fromHardProcessFinalState",&GENLepton_fromHardProcessFinalState, 
				       &b_GENLepton_fromHardProcessFinalState);
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
    }//end loading ntuples
  
  cout << "Total Events Loaded: " << totalentries << endl;
  
  //Defining histograms
  TH1F*hGenDielectronInvMass = new TH1F("hDielectronInvMass","",nLogBins,massbins);
  hGenDielectronInvMass->Sumw2();
  hGenDielectronInvMass->SetLineColor(kRed);
  hGenDielectronInvMass->SetTitle("Only Kinematic Cut Dielectrons");
  hGenDielectronInvMass->GetXaxis()->SetTitle("m_{ee} (GeV)");
  hGenDielectronInvMass->GetXaxis()->SetMoreLogLabels();
  hGenDielectronInvMass->GetXaxis()->SetNoExponent();
  TH1F*hGenMatchedDielectronInvMass = new TH1F("hGenMatchedDielectronInvMass","",nLogBins,massbins);
  hGenMatchedDielectronInvMass->Sumw2();
  hGenMatchedDielectronInvMass->SetLineColor(kBlue);
  hGenMatchedDielectronInvMass->SetTitle("Reco Gen Matched Dielectrons");
  hGenMatchedDielectronInvMass->GetXaxis()->SetTitle("m_{ee} (GeV)");  
  hGenMatchedDielectronInvMass->GetXaxis()->SetMoreLogLabels();
  hGenMatchedDielectronInvMass->GetXaxis()->SetNoExponent();
  TH1F*hGenPassIDdielectronInvMass = new TH1F("hGenPasIDDielectronInvMass","",nLogBins,massbins);
  hGenPassIDdielectronInvMass->Sumw2();
  hGenPassIDdielectronInvMass->SetLineColor(kRed);
  hGenPassIDdielectronInvMass->SetTitle("Medium ID Cut Dielectrons");
  hGenPassIDdielectronInvMass->GetXaxis()->SetTitle("m_{ee} (GeV)");  
  hGenPassIDdielectronInvMass->GetXaxis()->SetMoreLogLabels();
  hGenPassIDdielectronInvMass->GetXaxis()->SetNoExponent();
  TH1F*hGenAllDielectronInvMass = new TH1F("hGenAllDielectronInvMass","",nLogBins,massbins);
  hGenAllDielectronInvMass->Sumw2();
  hGenAllDielectronInvMass->SetTitle("Final State Dielectrons: No Other cuts");
  hGenAllDielectronInvMass->SetLineColor(kBlack);
  hGenAllDielectronInvMass->GetXaxis()->SetTitle("m_{ee} (GeV)");  
  hGenAllDielectronInvMass->GetXaxis()->SetMoreLogLabels();
  hGenAllDielectronInvMass->GetXaxis()->SetNoExponent();
  TH1F*hHLTGenDielectronInvMass = new TH1F("hHLTGenDielectronInvMass","",nLogBins,massbins);
  hHLTGenDielectronInvMass->Sumw2();
  hHLTGenDielectronInvMass->SetTitle("HLT Cut");
  hHLTGenDielectronInvMass->SetLineColor(kBlack);
  hHLTGenDielectronInvMass->GetXaxis()->SetTitle("m_{ee} (GeV)"); 
  hHLTGenDielectronInvMass->GetXaxis()->SetMoreLogLabels();
  hHLTGenDielectronInvMass->GetXaxis()->SetNoExponent();
  TH2F*hpTvsMass = new TH2F("hpTvsMass","",43,massbins,598,10,3000);
  hpTvsMass->GetXaxis()->SetMoreLogLabels();
  hpTvsMass->GetXaxis()->SetNoExponent();
  hpTvsMass->GetYaxis()->SetTitle("p_{T} [GeV]"); 
  hpTvsMass->GetXaxis()->SetTitle("m_{ee} [GeV]"); 

  TH2F*migMatrixGENisHardvsGENFS = new TH2F("migMatrixGENisHardvsGENFS","",nLogBins,massbins,nLogBins,massbins);
  migMatrixGENisHardvsGENFS->SetTitle("Migration Matrix: Gen-Level Hard Process vs. Gen-Level Final State");
  migMatrixGENisHardvsGENFS->GetXaxis()->SetTitle("Gen-Level Final State Dielectron Invariant Mass [GeV]");
  migMatrixGENisHardvsGENFS->GetYaxis()->SetTitle("Gen-Level Hard Process Dielecron Invariant mass [GeV]");
  migMatrixGENisHardvsGENFS->GetXaxis()->SetNoExponent();
  migMatrixGENisHardvsGENFS->GetXaxis()->SetMoreLogLabels();
  migMatrixGENisHardvsGENFS->GetYaxis()->SetNoExponent();
  migMatrixGENisHardvsGENFS->GetYaxis()->SetMoreLogLabels();
  TH2F*migMatrixGENFSvsReco = new TH2F("migMatrixGENFSvsReco","",nLogBins,massbins,nLogBins,massbins);
  migMatrixGENFSvsReco->SetTitle("Migration Matrix: Gen-Level Final State vs. Reconstructed");
  migMatrixGENFSvsReco->GetYaxis()->SetTitle("Gen-Level Final State Dielectron Invariant Mass [GeV]");
  migMatrixGENFSvsReco->GetXaxis()->SetTitle("Reconstructed Dielecron Invariant mass [GeV]");
  migMatrixGENFSvsReco->GetXaxis()->SetNoExponent();
  migMatrixGENFSvsReco->GetXaxis()->SetMoreLogLabels();
  migMatrixGENFSvsReco->GetYaxis()->SetNoExponent();
  migMatrixGENFSvsReco->GetYaxis()->SetMoreLogLabels();
  TH2F*migMatrixGENisHardvsReco = new TH2F("migMatrixGENisHardvsReco","",nLogBins,massbins,nLogBins,massbins);
  migMatrixGENisHardvsReco->SetTitle("Migration Matrix: Gen-Level Hard Process vs. Reconstructed");
  migMatrixGENisHardvsReco->GetYaxis()->SetTitle("Gen-Level Hard Process Dielectron Invariant Mass [GeV]");
  migMatrixGENisHardvsReco->GetXaxis()->SetTitle("Reconstructed Dielecron Invariant mass [GeV]");
  migMatrixGENisHardvsReco->GetXaxis()->SetNoExponent();
  migMatrixGENisHardvsReco->GetXaxis()->SetMoreLogLabels();
  migMatrixGENisHardvsReco->GetYaxis()->SetNoExponent();
  migMatrixGENisHardvsReco->GetYaxis()->SetMoreLogLabels();

  TH1I*hnEle = new TH1I("hnEle","",10,0,10);
  hnEle->SetTitle("Number of Dielectrons Passing Cuts");

  TH1F*hHardProcess[numChains];
  TString histbasename = "hHardProcess";
  TString histname;
  for(int jChain=0;jChain<numChains;jChain++)
    {
      histname = histbasename;
      histname+=jChain;
      hHardProcess[jChain] = new TH1F(histname,"",598,10,3000);
      hHardProcess[jChain]->SetFillColor(jChain+1);
      hHardProcess[jChain]->GetXaxis()->SetTitle("Gen-Level Dielectron Mass (isHardProcess) [GeV]"); 
      hHardProcess[jChain]->GetYaxis()->SetRangeUser(0.000001,1000000000);
      hHardProcess[jChain]->GetXaxis()->SetNoExponent();
      hHardProcess[jChain]->GetXaxis()->SetMoreLogLabels();
    }

  TLegend*legend = new TLegend(0.70,0.9,0.95,0.80);
  legend->SetTextSize(0.03);
  TLegend*legend2 = new TLegend(0.65,0.9,0.9,0.7);
  legend2->SetTextSize(0.02);

  TFile *rootFile = new TFile("./plots/plotsDY.root","RECREATE");
  TCanvas*canvas1 = new TCanvas("cInvMassHardProcess","",10,10,1000,1000);
  canvas1->SetLogx();
  canvas1->SetLogy();
  canvas1->cd();

  //Event Loop
  cout << "Starting Event Loop" << endl;
  double dpT, invMass, invMassHardProcess, weight, dRMin;
  int dRMinIndex;  
  Long64_t nentries;
  Long64_t count = 0;
  double nEvents = 250000;
  double lumi = chains[MC50to100]->GetEntries()/xSec[MC50to100];//luminosity of 50to100
  //double lumi = nEvents/xSec[MC50to100];//luminosity of 50to100
  TString HLTname = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
  TString trigName;
  int trigNameSize;
  long int nTooManyDielectrons = 0;
  long int nTooManyDielectronsFS = 0;

  for(int iChain=0;iChain<numChains;iChain++)
    {
      nentries = chains[iChain]->GetEntries();
      //nentries = nEvents;
      weight=lumi*(xSec[iChain]/nentries);
      for(Long64_t i=0;i<nentries;i++)
	{      
	  chains[iChain]->GetEntry(i);
	  counter(count,totalentries);
	  //counter(count,11*nentries);
	  count = count+1;

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
	      if(GENLepton_fromHardProcessFinalState[kLep]==1 && GENLepton_fromHardProcessFinalState[lLep]==1)
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
		  if(!passDileptonKinematics(Electron_pT[iEle],Electron_pT[jEle],Electron_eta[iEle],
					     Electron_eta[jEle])) continue; 	
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
	  hnEle->Fill(nEle);

	  if(idxRecoEle1>=0&&idxRecoEle2>=0)
	    invMassReco=calcInvMass(Electron_pT[idxRecoEle1],Electron_eta[idxRecoEle1],Electron_phi[idxRecoEle1],
				    eMass,Electron_pT[idxRecoEle2],Electron_eta[idxRecoEle2],
				    Electron_phi[idxRecoEle2],eMass);	  
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
					   GENLepton_phi[idxGenEle1],eMass,GENLepton_pT[idxGenEle2],
					   GENLepton_eta[idxGenEle2],GENLepton_phi[idxGenEle2],eMass);		  
	  invMass = calcInvMass(GENLepton_pT[idxGenEleFS1],GENLepton_eta[idxGenEleFS1],
				GENLepton_phi[idxGenEleFS1],eMass,GENLepton_pT[idxGenEleFS2],
				GENLepton_eta[idxGenEleFS2],GENLepton_phi[idxGenEleFS2],eMass);		  

	  hHardProcess[iChain]->Fill(invMassHardProcess,weight);//Fill if hard process gen-electrons
	  
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
	  hGenAllDielectronInvMass->Fill(invMass,weight);

	  // Apply kinematic acceptance criteria
	  if(!passDileptonKinematics(GENLepton_pT[idxGenEleFS1],GENLepton_pT[idxGenEleFS2],
				     GENLepton_eta[idxGenEleFS1], GENLepton_eta[idxGenEleFS2])) 
	    continue;	      
	  // Both electrons are in kinematic acceptance at gen level
	  hGenDielectronInvMass->Fill(invMass,weight);
	  hpTvsMass->Fill(invMass,GENLepton_pT[idxGenEleFS1],weight);
	  hpTvsMass->Fill(invMass,GENLepton_pT[idxGenEleFS2],weight);

	  // Apply matching to reconstructed electrons requirement
	  int closestTrackLep1, closestTrackLep2;
	  closestTrackLep1 = closestTrackLep2 = -1;
	  bool genToRecoMatchedLep1 = findGenToRecoMatch(idxGenEleFS1,closestTrackLep1);	      
	  bool genToRecoMatchedLep2 = findGenToRecoMatch(idxGenEleFS2,closestTrackLep2);	      
	  if(!(genToRecoMatchedLep1 && genToRecoMatchedLep2)) continue;

	  // Both electrons are reconstructed
	  hGenMatchedDielectronInvMass->Fill(invMass,weight);
	  
	  // Apply ID criteria:
	  // Dilepton pair at gen level matched to reco and passing ID at reco level
	  if(!Electron_passMediumID[closestTrackLep1]) continue;
	  if(!Electron_passMediumID[closestTrackLep2]) continue;
	  // Both electrons pass ID
	  hGenPassIDdielectronInvMass->Fill(invMass,weight);
	  
	  // Apply HLT requirement
	  if(!passHLT) continue;
	  // Event passed HLT cut
	  hHLTGenDielectronInvMass->Fill(invMass,weight);
	  migMatrixGENisHardvsGENFS->Fill(invMassHardProcess,invMass,weight);
	  migMatrixGENFSvsReco->Fill(invMass,invMassReco,weight);
	  migMatrixGENisHardvsReco->Fill(invMassHardProcess,invMassReco,weight);
	}//end event loop   

      if(iChain==0)hHardProcess[iChain]->Draw("Bar");      
      else hHardProcess[iChain]->Draw("Barsame");	 	
    }//end chain loop 

  legend2->AddEntry(hHardProcess[0],"DYEE M10-50");
  legend2->AddEntry(hHardProcess[1],"DYEE M50-100");
  legend2->AddEntry(hHardProcess[2],"DYEE M100-200");
  legend2->AddEntry(hHardProcess[3],"DYEE M200-400");
  legend2->AddEntry(hHardProcess[4],"DYEE M400-500");
  legend2->AddEntry(hHardProcess[5],"DYEE M500-700");
  legend2->AddEntry(hHardProcess[6],"DYEE M700-800");
  legend2->AddEntry(hHardProcess[7],"DYEE M800-1000");
  legend2->AddEntry(hHardProcess[8],"DYEE M1000-1500");
  legend2->AddEntry(hHardProcess[9],"DYEE M1500-2000");
  legend2->AddEntry(hHardProcess[10],"DYEE M2000-3000");
  legend2->Draw("same");

  TEfficiency* hMatchedEfficiency = new TEfficiency((*hGenMatchedDielectronInvMass),(*hGenDielectronInvMass));
  hMatchedEfficiency->SetTitle("Reconstruction Efficiency");
  hMatchedEfficiency->SetMarkerStyle(20);
  hMatchedEfficiency->SetMarkerSize(0.5);
  hMatchedEfficiency->SetName("RecoEfficiency");
  hMatchedEfficiency->SetStatisticOption(TEfficiency::kFNormal);
  TEfficiency* hIDEfficiency = new TEfficiency((*hGenPassIDdielectronInvMass),(*hGenMatchedDielectronInvMass));
  hIDEfficiency->SetTitle("ID Efficiency");
  hIDEfficiency->SetMarkerStyle(20);
  hIDEfficiency->SetMarkerSize(0.5);
  hIDEfficiency->SetName("IDEfficiency");
  hIDEfficiency->SetStatisticOption(TEfficiency::kFNormal);
  TEfficiency* hAcceptance = new TEfficiency((*hGenDielectronInvMass),(*hGenAllDielectronInvMass));
  hAcceptance->SetTitle("Acceptance");
  hAcceptance->SetMarkerStyle(20);
  hAcceptance->SetMarkerSize(0.5);
  hAcceptance->SetName("Acceptance");  
  hAcceptance->SetStatisticOption(TEfficiency::kFNormal);
  TEfficiency* hHLTEfficiency = new TEfficiency((*hHLTGenDielectronInvMass),(*hGenPassIDdielectronInvMass));
  hHLTEfficiency->SetTitle("HLT Efficiency");
  hHLTEfficiency->SetMarkerStyle(20);
  hHLTEfficiency->SetMarkerSize(0.5);
  hHLTEfficiency->SetName("HLTEfficiency");
  hHLTEfficiency->SetStatisticOption(TEfficiency::kFNormal);
  TEfficiency* hEfficiency = new TEfficiency((*hHLTGenDielectronInvMass),(*hGenDielectronInvMass));
  hEfficiency->SetTitle("Total Efficiency");
  hEfficiency->SetMarkerStyle(20);
  hEfficiency->SetMarkerSize(0.5);
  hEfficiency->SetName("Efficiency");
  hEfficiency->SetStatisticOption(TEfficiency::kFNormal);
  
  //Drawing on canvases
  TCanvas*canvas2 = new TCanvas("cAllEfficiencies","",10,10,1400,700);
  canvas2->Divide(3);
  canvas2->cd(1);
  TVirtualPad*p1 = canvas2->cd(1);
  p1->SetLogx();
  hMatchedEfficiency->Draw();
  gPad->Update();
  auto graph1 = hMatchedEfficiency->GetPaintedGraph();
  graph1->SetMinimum(0);
  graph1->SetMaximum(1);
  graph1->GetXaxis()->SetMoreLogLabels();
  graph1->GetXaxis()->SetNoExponent();
  graph1->Draw("PE");
  gPad->Update();

  canvas2->cd(2);
  TVirtualPad*p2 = canvas2->cd(2);
  p2->SetLogx();
  hIDEfficiency->Draw();
  gPad->Update();
  auto graph2 = hIDEfficiency->GetPaintedGraph();
  graph2->SetMinimum(0);
  graph2->SetMaximum(1);
  graph2->GetXaxis()->SetMoreLogLabels();
  graph2->GetXaxis()->SetNoExponent();
  graph2->Draw("PE");
  gPad->Update();

  canvas2->cd(3);
  TVirtualPad*p3 = canvas2->cd(3);
  p3->SetLogx();
  hHLTEfficiency->Draw();
  gPad->Update();
  auto graph3 = hHLTEfficiency->GetPaintedGraph();
  graph3->SetMinimum(0);
  graph3->SetMaximum(1);
  graph3->GetXaxis()->SetMoreLogLabels();
  graph3->GetXaxis()->SetNoExponent();
  graph3->Draw("PE");
  gPad->Update();
  
  TCanvas*canvas3 = new TCanvas("cEfficiency","",10,10,900,700);
  canvas3->cd();
  TVirtualPad*p4 = canvas3->cd();
  p4->SetLogx();
  hEfficiency->Draw();
  gPad->Update();
  auto graph4 = hEfficiency->GetPaintedGraph();
  graph4->SetMinimum(0);
  graph4->SetMaximum(1);
  graph4->GetXaxis()->SetMoreLogLabels();
  graph4->GetXaxis()->SetNoExponent();
  graph4->Draw("PE");
  gPad->Update();

  TCanvas*canvas4 = new TCanvas("cAcceptance","",10,10,900,700);
  canvas4->cd();
  TVirtualPad*p5 = canvas4->cd();
  p5->SetLogx();
  hAcceptance->Draw();
  gPad->Update();
  auto graph5 = hAcceptance->GetPaintedGraph();
  graph5->SetMinimum(0);
  graph5->SetMaximum(1);
  graph5->GetXaxis()->SetMoreLogLabels();
  graph5->GetXaxis()->SetNoExponent();
  graph5->Draw("PE");
  gPad->Update();
  
  TCanvas*canvas5 = new TCanvas("cPtvsMassProfile","",10,10,900,700);
  canvas5->SetLogx();
  canvas5->SetLogy();
  canvas5->SetGrid();
  TProfile*hpTvsMassProf = hpTvsMass->ProfileX();
  hpTvsMassProf->GetYaxis()->SetMoreLogLabels();
  hpTvsMassProf->GetYaxis()->SetNoExponent();
  hpTvsMassProf->GetYaxis()->SetTitleOffset(1.25);
  hpTvsMassProf->SetTitle("#LTp_{T}#GT vs.Invariant Mass");
  hpTvsMassProf->GetYaxis()->SetTitle("#LTp_{T}#GT  [GeV]");
  hpTvsMassProf->Draw();

  TCanvas*canvas6 = new TCanvas("cMigMatrixGENFSvsGENisHard","",10,10,900,700);
  canvas6->SetLogy();
  canvas6->SetLogx();
  gStyle->SetPalette(1);
  migMatrixGENisHardvsGENFS->Draw("colz");

  TCanvas*canvas7 = new TCanvas("cmigMatrixGENFSvsReco","",10,10,900,700);
  canvas7->SetLogy();
  canvas7->SetLogx();
  gStyle->SetPalette(1);
  migMatrixGENFSvsReco->Draw("colz");

  TCanvas*canvas8 = new TCanvas("cmigMatrixGENisHardvsReco","",10,10,900,700);
  canvas8->SetLogy();
  canvas8->SetLogx();
  gStyle->SetPalette(1);
  migMatrixGENisHardvsReco->Draw("colz");

  TCanvas*canvas9 = new TCanvas("cnEle","",10,10,900,700);
  hnEle->Draw();
    
  canvas1->SaveAs("./plots/hardProcessInvMass.png");
  canvas2->SaveAs("./plots/allEfficiencies.png");
  canvas3->SaveAs("./plots/efficiency.png");
  canvas4->SaveAs("./plots/acceptance.png");
  canvas5->SaveAs("./plots/pTvsMassProf.png");
  canvas6->SaveAs("./plots/migMatrixGENFSvsGENisHard.png");
  canvas7->SaveAs("./plots/migMatrixGENFSvsGReco.png");
  canvas8->SaveAs("./plots/migMatrixGENisHardvsReco.png");
  canvas9->SaveAs("./plots/dielectrons.png");

  rootFile->cd();
  hIDEfficiency->Write();
  hGenMatchedDielectronInvMass->Write();
  hGenDielectronInvMass->Write();
  hGenPassIDdielectronInvMass->Write();
  hGenAllDielectronInvMass->Write();
  hMatchedEfficiency->Write();
  hAcceptance->Write();
  hEfficiency->Write();
  hHLTEfficiency->Write();
  hpTvsMassProf->Write();
  hpTvsMass->Write();
  migMatrixGENisHardvsGENFS->Write();
  migMatrixGENFSvsReco->Write();
  migMatrixGENisHardvsReco->Write();
  canvas1->Write();
  canvas2->Write();
  canvas3->Write();
  canvas4->Write();
  canvas5->Write();
  canvas6->Write();
  canvas7->Write();
  canvas8->Write();
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
