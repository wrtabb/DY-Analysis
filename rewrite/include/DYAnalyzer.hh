
#include <TString.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <iostream>
#include "SharedVariableNames.h"
#include "NtuplesV2P3Location.h"

#ifndef DYAnalyzer_HH
#define DYAnalyzer_HH

class DYAnalyzer
{
public:
 TChain*chains[numChains];
 Long64_t LoadTrees();
 void InitBranches();
 double CalcInvMass(TLorentzVector v1,TLorentzVector v2);
 Long64_t GetDYEntries(int iChain);
 Long64_t GetDYEntry(int iChain,Long64_t iEntry);
 void Counter(Long64_t i,Long64_t N,TString name);
 int GetGenLeptons(LepType lepType,int &idxHardEle1,int &idxHardEle2,
                                   int &idxFSREle1, int &idxFSREle2);
 void CutOnKinematics(double pt1,double pt2,double eta1,double eta2);
 void FindGenToRecoMatch(int genIndex,int &recoIndex);
};

#endif
