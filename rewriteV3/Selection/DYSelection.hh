#ifndef Selection_HH
#define Selection_HH
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

class DYSelection
{
public:
 DYSelection(NtupleVersion ntup,LepType lepType,SampleType sampleType);
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

 Long64_t LoadTrees(NtupleVersion ntup,std::vector<TString>dirNames,SampleType sampleType,
                    LepType lepType); 
 void InitBranches(int numChains,bool isMC,bool isReco,LepType lepType);
 void LoadHistograms();
};
#endif
