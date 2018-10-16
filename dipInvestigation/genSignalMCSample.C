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

void counter(Long64_t i, Long64_t N);
double calcInvMass(double pt1,double eta1,double phi1,double m1,double pt2,double eta2,double phi2,double m2);
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2);
bool passPromptGenElectron(int ID, int fromfinalstate);
bool passHardProcess(int ID, int hardProces);
bool findGenToRecoMatch(int genIndex,int &recoIndex);

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
int HLT_trigType[MPSIZE],HLT_trigFired[MPSIZE];
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
enum InvMassHist 
{
  KINEMATIC_CUTS,
  RECO_MATCHED,
  ID_CUTS,
  ALL_ELE,
  HLT_CUTS,
  RECO_ELE
};
std::vector<std::string> HLT_trigName;
std::vector<std::string> *pHLT_trigName = &HLT_trigName;

const double massbins[] = {15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,
                            57.5,60,62,64,66,68,70,72,74,76,78.5,81,83.5,86,88.5,91,93.5,96,
                            98.5,101,103.5,106,108.25,108.5,108.75,109,109.25,109.5,109.75,110,                            112.5,115,117.5,120,123,126,129.5,133,
                            137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,210,220,
                            231.5,243,258,273,296.5,320,350,380,410,440,475,510,555,600,650,
                            700,765,830,915,1000,1250,1500,2250,3000};
const double pi=TMath::Pi();
const int numChains = 11;
const int nBins = 2*46;
//Cross sections obtained from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SNUCMSYooDYntuple
const float xSec = 1921.8*3;
const float etaHigh = 2.4;
const float etaGapHigh = 1.566; 
const float etaGapLow = 1.4442;
const float ptHigh = 28;
const float ptLow = 17;
const float eMass = 0.000511;
const float dRMinCut = 0.3;
const int dataLuminosity = 35867;
const TString treeName = "recoTree/DYTree";

const TString pileupRatioName = "/home/hep/wrtabb/git/DY-Analysis/plots/pileup.root";
const TString leg2SFName = "/home/hep/wrtabb/git/DY-Analysis/data/Leg2_SF.root";
const TString medIDSFName = "/home/hep/wrtabb/git/DY-Analysis/data/MediumID_SF.root";
const TString recoSFName = "/home/hep/wrtabb/git/DY-Analysis/data/Reco_SF.root";


const int nSubSamples10to50 = 3;
const int nSubSamples100to200 = 2;
const int ptBinHigh = 499;
const int ptBinLow = 26;
const int nMatrixHistos = 3;
const int nInvMassHistos = 6;
const TString histMatrixNames[nMatrixHistos] = {"hGenHardvsGenFS","hGenFSvsReco",
  "hGenHardvsReco"};
const TString matrixYaxisTitles[nMatrixHistos] = {"Gen-Level Final State Invariant Mass [GeV]",
  "Reco Invariant Mass [GeV]","Reco Invariant Mass [GeV]"};
const TString matrixXaxisTitles[nMatrixHistos] = {"Gen-Level Hard Process Invariant Mass [GeV]"  ,"Gen-Level Final State Invariant Mass [GeV]","Gen-Level Hard Process Invariant Mass [GeV]"};
const TString histInvMassNames[nInvMassHistos] = {"hGenInvMass","hGenMatchedInvMass",
  "hGenPassIDInvMass","hGenAllInvMass","hHLTGenInvMass","hRecoInvMass"};
const TString histInvMassTitles[nInvMassHistos] = {"Only Kinematic Cuts","Reco-Gen Matched",
  "Medium ID Cuts","Final State: No cuts","HLT Cut","Reconstructed"};

void genSignalMCSample()
{
  TTimeStamp ts_start;
  cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
  TStopwatch totaltime;
  totaltime.Start();
  //gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  //Defining branches
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

  //Loading ntuples
  cout << "Loading ntuples" << endl;

    TString dirName = "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE"; 
 TString baseDir = 
   "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/ikrav/DrellYan_13TeV_2016/v2p3/DYJetsToLL_allMasses_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"; 
  TChain*chains;
  TString files;  
  Long64_t nentries = 0;
  chains = new TChain(treeName);
  int nFiles = 0;
  int maxFiles = 100;//maximum number of files to load
  cout << "****************************************************" << endl;
  cout << "Loading ntuples" << endl;
    for(int k=0;k<1000;k++) {
      if(nFiles>=maxFiles) break;
      files = baseDir;
      files += dirName;
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
 
      chains->SetBranchAddress("GENEvt_weight",&GENEvt_weight,&b_GENEvt_weight);
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
      chains->SetBranchAddress("nPileUp", &nPileUp, &b_nPileUp);
      chains->SetBranchAddress("Electron_pT", &Electron_pT, &b_Electron_pT);
      chains->SetBranchAddress("Electron_eta",&Electron_eta, &b_Electron_eta);
      chains->SetBranchAddress("Electron_phi",&Electron_phi, &b_Electron_phi);
      chains->SetBranchAddress
        ("Electron_passMediumID",&Electron_passMediumID,&b_Electron_passMediumID);
      chains->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
      chains->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
      chains->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
      chains->SetBranchAddress("HLT_trigName",&pHLT_trigName);
      
 cout << endl;
 cout << "Total Events Loaded: " << nentries << endl;

  TH1F*histInvMass[nInvMassHistos];
  for(int i=0;i<nInvMassHistos;i++){
    histInvMass[i]=new TH1F(histInvMassNames[i],"",nBins,massbins);
    histInvMass[i]->Sumw2();
    histInvMass[i]->GetXaxis()->SetTitle("invariant mass [GeV]");
    histInvMass[i]->GetXaxis()->SetMoreLogLabels();
    histInvMass[i]->GetXaxis()->SetNoExponent();
    histInvMass[i]->SetTitle(histInvMassTitles[i]);
  }

  TFile*pileupRatioFile  = new TFile(pileupRatioName);
  TH1F*hPileupRatio = (TH1F*)pileupRatioFile->Get("hPileupRatio");
  TFile*fileLeg2SF = new TFile(leg2SFName);
  TH2F*hLeg2SF = (TH2F*)fileLeg2SF->Get("EGamma_SF2D");
  TFile*fileMedIDSF = new TFile(medIDSFName);
  TH2F*hMedIDSF = (TH2F*)fileMedIDSF->Get("EGamma_SF2D");
  TFile*fileRecoSF = new TFile(recoSFName);
  TH2F*hRecoSF = (TH2F*)fileRecoSF->Get("EGamma_SF2D");

 
 //Event Loop
 cout << endl;
 cout << "Starting Event Loop" << endl;
 double invMassFSR, rapidity, dileptonPt, dileptonEta, xSecWeightAlone, xSecWeightGen, 
   weightNoPileup, genWeight, varGenWeight, totalWeight, lumiEffective, nEffective, localEntry,
   sumGenWeight, sumRawGenWeight, pileupWeight, sfReco1, sfReco2, sfID1, sfID2, sfHLT;
 double invMassHardProcess,sfWeight;
 Long64_t count = 0;
 double nEvents = 250000;
 double lumi = dataLuminosity;
 //double lumi = nEvents/xSec[MC50to100];//luminosity of 50to100
 TString HLTname = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
 TString trigName;
 int trigNameSize;
 long int nTooManyDielectrons = 0;
 long int nTooManyDielectronsFS = 0; 
 double eEta1, eEta2, ePt1, ePt2;

      nentries = chains->GetEntries();
      xSecWeightGen=lumi*(xSec/1.0);      
      xSecWeightAlone=lumi*(xSec/nentries);
      sumGenWeight = 0;
      for(Long64_t i=0;i<nentries;i++){
	localEntry = chains->LoadTree(i);
	b_GENEvt_weight->GetEntry(localEntry);
	genWeight = GENEvt_weight/fabs(GENEvt_weight);	//normalized genweight
	sumGenWeight += genWeight;
	varGenWeight += GENEvt_weight*GENEvt_weight; //variance of genweights
	sumRawGenWeight += GENEvt_weight; 	
      }
      nEffective = (sumRawGenWeight*sumRawGenWeight)/varGenWeight;
      lumiEffective = nEffective/xSec;
      pileupWeight = hPileupRatio->GetBinContent(hPileupRatio->FindBin(nPileUp));
      genWeight = GENEvt_weight/fabs(GENEvt_weight);
      genWeight = genWeight/sumGenWeight;

      for(Long64_t i=0;i<nentries;i++)
	{      
	  chains->GetEntry(i);
	  counter(i,nentries);
	  

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
	      if(GENLepton_fromHardProcessFinalState[kLep]==1 && 
                 GENLepton_fromHardProcessFinalState[lLep]==1)
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
		  if(!passDileptonKinematics(Electron_pT[iEle],Electron_pT[jEle],
                    Electron_eta[iEle],Electron_eta[jEle])) continue; 	
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

	  if(idxRecoEle1>=0&&idxRecoEle2>=0)
	    invMassReco=calcInvMass(Electron_pT[idxRecoEle1],Electron_eta[idxRecoEle1],
              Electron_phi[idxRecoEle1],eMass,Electron_pT[idxRecoEle2],
              Electron_eta[idxRecoEle2],Electron_phi[idxRecoEle2],eMass);	  
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
	    GENLepton_phi[idxGenEle1],eMass,GENLepton_pT[idxGenEle2],GENLepton_eta[idxGenEle2],            GENLepton_phi[idxGenEle2],eMass);		  
	  invMassFSR = calcInvMass(GENLepton_pT[idxGenEleFS1],GENLepton_eta[idxGenEleFS1],
            GENLepton_phi[idxGenEleFS1],eMass,GENLepton_pT[idxGenEleFS2],
            GENLepton_eta[idxGenEleFS2],GENLepton_phi[idxGenEleFS2],eMass);		  

          eEta1 = Electron_eta[idxRecoEle1];
          eEta2 = Electron_eta[idxRecoEle2];
          ePt1 = Electron_pT[idxRecoEle1];
          ePt2 = Electron_pT[idxRecoEle2];

          if(ePt1<ptBinLow) ePt1 = ptBinLow;//pull this information from the histograms
          if(ePt2<ptBinLow) ePt2 = ptBinLow;//raise bin
          if(ePt1>ptBinHigh) ePt1 = ptBinHigh;//lower bin
          if(ePt2>ptBinHigh) ePt2 = ptBinHigh;//

	  pileupWeight = hPileupRatio->GetBinContent(hPileupRatio->FindBin(nPileUp));
	  genWeight = GENEvt_weight/fabs(GENEvt_weight);
	  genWeight = genWeight/sumGenWeight;
	  sfReco1=hRecoSF->GetBinContent(hRecoSF->FindBin(eEta1,ePt1));
          sfReco2=hRecoSF->GetBinContent(hRecoSF->FindBin(eEta2,ePt2));
          sfID1=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eEta1,ePt1));
          sfID2=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eEta2,ePt2));
          sfHLT=(hLeg2SF->GetBinContent(hLeg2SF->FindBin(eEta1,ePt1)))*
            (hLeg2SF->GetBinContent(hLeg2SF->FindBin(eEta2,ePt2)));
          sfWeight = sfReco1*sfReco2*sfID1*sfID2*sfHLT;
	  totalWeight = genWeight*xSecWeightGen*pileupWeight*sfWeight;
          totalWeight = sfWeight;

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

	  // First, fill histogram for all dielectrons
	  histInvMass[ALL_ELE]->Fill(invMassFSR,totalWeight);
            
	  // Apply kinematic acceptance criteria
	  if(!passDileptonKinematics(GENLepton_pT[idxGenEleFS1],GENLepton_pT[idxGenEleFS2],
				     GENLepton_eta[idxGenEleFS1], GENLepton_eta[idxGenEleFS2])) 
	    continue;	      
	  // Both electrons are in kinematic acceptance at gen level
	  histInvMass[KINEMATIC_CUTS]->Fill(invMassFSR,totalWeight);

	  // Apply matching to reconstructed electrons requirement
	  int closestTrackLep1, closestTrackLep2;
	  closestTrackLep1 = closestTrackLep2 = -1;
	  bool genToRecoMatchedLep1 = findGenToRecoMatch(idxGenEleFS1,closestTrackLep1);	      
	  bool genToRecoMatchedLep2 = findGenToRecoMatch(idxGenEleFS2,closestTrackLep2);	      
	  if(!(genToRecoMatchedLep1 && genToRecoMatchedLep2)) continue;

	  // Both electrons are reconstructed
	  histInvMass[RECO_MATCHED]->Fill(invMassFSR,totalWeight);
	  
	  // Apply ID criteria:
	  // Dilepton pair at gen level matched to reco and passing ID at reco level
	  if(!Electron_passMediumID[closestTrackLep1]) continue;
	  if(!Electron_passMediumID[closestTrackLep2]) continue;
	  // Both electrons pass ID
	  histInvMass[ID_CUTS]->Fill(invMassFSR,totalWeight);
	  
	  // Apply HLT requirement
	  if(!passHLT) continue;
          

	  histInvMass[HLT_CUTS]->Fill(invMassFSR,totalWeight);
	  histInvMass[RECO_ELE]->Fill(invMassReco,totalWeight);
	}//end event loop   
      
  TCanvas*canvas[nInvMassHistos];
  TString canvasName;
  TString saveName;
  for(int i=0;i<nInvMassHistos;i++){
    canvasName = "canvas";
    canvasName += i;
    saveName = canvasName;
    saveName += "_gen_SFWeight.png";
    canvas[i]=new TCanvas(canvasName,"",10,10,1000,1000);
    canvas[i]->SetLogy();
    canvas[i]->SetGrid();
    histInvMass[i]->GetXaxis()->SetRangeUser(100,125);
    histInvMass[i]->SetFillColor(kRed+2);
    histInvMass[i]->Draw("hist");
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
