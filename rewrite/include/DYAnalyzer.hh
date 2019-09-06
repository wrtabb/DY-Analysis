
#ifndef DYAnalyzer_HH
#define DYAnalyzer_HH
#include <TString.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>
#include "SharedVariableNames.h"
#include "NtuplesV2P3Location.h"


class DYAnalyzer
{
public:
 double xSecWeight;
 double genWeight;
 double invMass;
 double weight;

 Long64_t LoadTrees();
 double CalcInvMass(TLorentzVector v1,TLorentzVector v2);
 Long64_t GetDYEntries(int iChain);
 Long64_t GetDYEntry(int iChain,Long64_t iEntry);
 int GetGenLeptons(LepType lepType,int &idxHardEle1,int &idxHardEle2,
                                   int &idxFSREle1, int &idxFSREle2);
 //-----Cuts-----//
 bool CutOnKinematics(double pt1,double pt2,double eta1,double eta2);
 bool FindGenToRecoMatch(int genIndex,int &recoIndex);

 void Counter(Long64_t i,Long64_t N,TString name);

 //-----Weights-----//
 double GetTotalWeight(int iChain,double genWeight,double xSecWeight,
                   double eta1,double eta2,double pt1,double pt2);
 double GetGenWeight(int iChain);
 double GetXsecWeight(int iChain,bool useGenWeight);

private:
 double sfReco1;
 double sfReco2;
 double sfID1;
 double sfID2;
 double sfHLT;
 double sfWeight;
 double pileupWeight;
 double totalWeight;
 double sumGenWeight;
 double sumRawGenWeight;
 double varGenWeight;

 TChain*chains[numChains];
 TFile*pileupRatioFile;
 TFile*fileLeg2SF; 
 TFile*fileMedIDSF;
 TFile*fileRecoSF;
 void InitBranches();
 void LoadHistograms();
 
 TH1F*hPileupRatio;
 TH2F*hLeg2SF;  
 TH2F*hMedIDSF; 
 TH2F*hRecoSF;  

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
 TBranch*b_GENnPair;
 TBranch*b_GENLepton_eta;
 TBranch*b_GENLepton_phi;
 TBranch*b_GENLepton_pT;
 TBranch*b_GENLepton_ID;
 TBranch*b_GENLepton_isHardProcess;
 TBranch*b_GENLepton_fromHardProcessFinalState;
};

#endif
