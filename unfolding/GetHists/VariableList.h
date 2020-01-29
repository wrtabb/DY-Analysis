#include "NtuplesV2P6Location.h"

const TString treeName = "recoTree/DYTree";
std::vector<TString> dirNamesLL = {
 DYLL_M10to50,   
 DYLL_M50to100,
 DYLL_M100to200,
 DYLL_M200to400,
 DYLL_M400to500,
 DYLL_M500to700,
 DYLL_M700to800,
 DYLL_M800to1000, 
 DYLL_M1000to1500,
 DYLL_M1500to2000,
 DYLL_M2000to3000
};

std::vector<TString> dirNamesEW = {
 WW_dir,
 WZ_dir,
 ZZ_dir,
 WJets,
 WJets_ext, 
 WJets_ext2v5
};

std::vector<TString> dirNamesTT = {
 ST_tW,
 ST_tbarW,
 ttbar,
 ttbarBackup,
 ttbarBackup_M0to700,
 ttbar_M0to700,
 ttbar_M700to1000,
 ttbar_M1000toInf
};

std::vector<TString> dirNamesData = {
 DoubleEG_RunB,
 DoubleEG_RunC,
 DoubleEG_RunD,
 DoubleEG_RunE,
 DoubleEG_RunF,
 DoubleEG_RunG,
 DoubleEG_RunHver2,
 DoubleEG_RunHver3
};

std::vector<TString> dirNamesSM = {
 SM_2016B,
 SM_2016C,
 SM_2016D,
 SM_2016E,
 SM_2016F,
 SM_2016G,
 SM_2016Hver2,
 SM_2016Hver3
};

std::vector<TString> dirNamesZtoEE = {
 ZToEE_M50to120,
 ZToEE_M120to200,
 ZToEE_M200to400,
 ZToEE_M400to800,
 ZToEE_M800to1400,
 ZToEE_M1400to2300,
 ZToEE_M2300to3500,
 ZToEE_M3500to4500,
 ZToEE_M4500to6000,
 ZToEE_M6000toInf
};

enum SampleType{
 LL,
 EW,
 TT,
 DATA
};

enum LepType{
 ELE,
 MUON,
 TAU
};

enum ChainLL{
 M10to50,   
 M50to100,
 M100to200,
 M200to400,
 M400to500,
 M500to700,
 M700to800,
 M800to1000, 
 M1000to1500,
 M1500to2000,
 M2000to3000
};

enum ChainEW{
 WW,
 WZ,
 ZZ,
 W_PLUS_JETS
};

enum ChainTT{
 TW,
 T_BAR_W,
 TTbar,
 TT0to700,
 TT700to1000,
 TT1000toInf
};

enum ChainData{
 RUN_B,
 RUN_C,
 RUN_D,
 RUN_E,
 RUN_F,
 RUN_G,
 RUN_H
};
