
#ifndef DYAnalyzer_HH
#define DYAnalyzer_HH
#include <TString.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TEfficiency.h>
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

 DYAnalyzer(NtupleVersion ntup, LepType lepType, SampleType sampleType);
 double CalcInvMass(double pt1,double eta1,double phi1,double m1,double pt2,
                    double eta2,double phi2,double m2);
 Long64_t GetDYEntries(int iChain);
 Long64_t GetDYEntry(int iChain,Long64_t iEntry);

 int GetGenLeptons(LepType lepType,int &idxHardEle1,int &idxHardEle2,
                   int &idxFSREle1, int &idxFSREle2);
 int GetRecoElectrons(int &leadEle, int &subEle);

 //-----Cuts-----//
 bool AcceptanceCut(double pt1,double pt2,double eta1,double eta2);
 bool GenToRecoMatchCut(int genIndex,int &recoIndex);
 bool MediumIDCut(bool passID1, bool passID2);
 bool HLTCut();

 //-----Functions-----//
 void GetEfficiencies(TH1*hist0,TH1*hist1,TString name);
 void Counter(Long64_t i,Long64_t N,TString name);
 TH1D*DefineMassHist(BinType type,TString histName,int nBins);

 //-----Weights-----//
 double GetTotalWeight(bool isReco,int iChain,double genWeight,double xSecWeight,
                       double eta1,double eta2,double pt1,double pt2);
 double GetGenWeightSum(int iChain);
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

 //Functions
 Long64_t LoadTrees(NtupleVersion ntup,std::vector<TString> dirNames,SampleType sampleType,
                    LepType lepType);
 void InitBranches(const int numChains,bool isMC,bool isReco);
 void LoadHistograms();
 void EventLoop(LepType lepType,TChain*chain,int numchains);
 vector<vector<double>> ReturnAllParametersAllEvents(); 
 vector<double> ReturnAllParameters(LepType lepType,TChain*chain,Long64_t iEvent);
};

#endif
