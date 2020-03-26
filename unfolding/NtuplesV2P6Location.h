#ifndef NtuplesV2P6Location_H
#define NtuplesV2P6Location_H

#include <TString.h>
#include <TChain.h>
#include <TTimeStamp.h>
#include <iostream>
#include <TStopwatch.h>

const TString base= "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/wtabb/DrellYan_13TeV_2016/v2p6/";

//-----DYtoLL-----//
const TString DYLL_M10to50 =     base + "DYLL_M10to50/DYLL_M10to50";
const TString DYLL_M50to100 =    base + "DYLL_M50toInf/DYLL_M50toInf_truncated_M50To100";
const TString DYLL_M100to200 =   base + "DYLL_M100to200/DYLL_M100to200";
const TString DYLL_M200to400 =   base + "DYLL_M200to400/DYLL_M200to400";
const TString DYLL_M400to500 =   base + "DYLL_M400to500/DYLL_M400to500";
const TString DYLL_M500to700 =   base + "DYLL_M500to700/DYLL_M500to700";
const TString DYLL_M700to800 =   base + "DYLL_M700to800/DYLL_M700to800";
const TString DYLL_M800to1000 =  base + "DYLL_M800to1000/DYLL_M800to1000";
const TString DYLL_M1000to1500 = base + "DYLL_M1000to1500/DYLL_M1000to1500";
const TString DYLL_M1500to2000 = base + "DYLL_M1500to2000/DYLL_M1500to2000";
const TString DYLL_M2000to3000 = base + "DYLL_M2000to3000/DYLL_M2000to3000";

//-----Bosons-----//
const TString WW_dir =       base + "WW";
const TString WZ_dir =       base + "WZ";
const TString ZZ_dir =       base + "ZZ";
const TString WJets =        base + "WJetsToLNu_amcatnlo";

//-----tops-----//
const TString ST_tW =               base + "ST_tW";
const TString ST_tbarW =            base + "ST_tbarW";
const TString ttbar =               base + "ttbar";
const TString ttbarBackup_M0to700 = base + "ttbarBackup_truncated_M0To700";
const TString ttbar_M0to700 =       base + "ttbar_truncated_M0To700";
const TString ttbar_M700to1000 =    base + "ttbar_M700to1000";
const TString ttbar_M1000toInf =    base + "ttbar_M1000toInf";

//-----DoubleEG-----//
const TString DoubleEG_RunB =     base + "crab_DoubleEG_RunB";
const TString DoubleEG_RunC =     base + "crab_DoubleEG_RunC";
const TString DoubleEG_RunD =     base + "crab_DoubleEG_RunD";
const TString DoubleEG_RunE =     base + "crab_DoubleEG_RunE";
const TString DoubleEG_RunF =     base + "crab_DoubleEG_RunF";
const TString DoubleEG_RunG =     base + "crab_DoubleEG_RunG";
const TString DoubleEG_RunH =     base + "crab_DoubleEG_RunH";

//-----Single Muon-----//
const TString SM_2016B     = base + "SingleMuon_Run2016B";
const TString SM_2016C     = base + "SingleMuon_Run2016C";
const TString SM_2016D     = base + "SingleMuon_Run2016D";
const TString SM_2016E     = base + "SingleMuon_Run2016E";
const TString SM_2016F     = base + "SingleMuon_Run2016F";
const TString SM_2016G     = base + "SingleMuon_Run2016G";
const TString SM_2016H     = base + "SingleMuon_Run2016H";

//-----ZtoEE-----//
const TString ZToEE_M50to120 =    base + "ZToEE_M50to120";
const TString ZToEE_M120to200 =   base + "ZToEE_M120to200";
const TString ZToEE_M200to400 =   base + "ZToEE_M200to400";
const TString ZToEE_M400to800 =   base + "ZToEE_M400to800";
const TString ZToEE_M800to1400 =  base + "ZToEE_M800to1400";
const TString ZToEE_M1400to2300 = base + "ZToEE_M1400to2300";
const TString ZToEE_M2300to3500 = base + "ZToEE_M2300to3500";
const TString ZToEE_M3500to4500 = base + "ZToEE_M3500to4500";
const TString ZToEE_M4500to6000 = base + "ZToEE_M4500to6000";
const TString ZToEE_M6000toInf =  base + "ZToEE_M6000toInf";


#endif
