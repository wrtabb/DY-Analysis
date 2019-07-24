//-----Forard declarations-----//
void getList(TString directory);
void GetFileList(TString directory);
//-----Number of directories to skim-----//
const int nDir = 53;
//-----Location of directories in my storage area-----//
const TString baseDir = "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/wtabb/DrellYan_13TeV_2016/v2p6/";
//-----List of all directories to be skimmed-----//
TString directory[nDir] = {
"DYLL_M1000to1500/",
"DYLL_M100to200/",
"DYLL_M10to50/ext1v1/",
"DYLL_M10to50/v1/",
"DYLL_M10to50/v2/",
"DYLL_M1500to2000/",
"DYLL_M2000to3000/",
"DYLL_M200to400/",
"DYLL_M400to500/",
"DYLL_M500to700/",
"DYLL_M50toInf/base/",
"DYLL_M50toInf/madgraph/",
"DYLL_M50toInf/madgraph_ext/",
"DYLL_M700to800/",
"DYLL_M800to1000/",
"ST_tW/",
"ST_tbarW/",
"SingleMuon_Run2016B/",
"SingleMuon_Run2016C/",
"SingleMuon_Run2016D/",
"SingleMuon_Run2016E/",
"SingleMuon_Run2016F/",
"SingleMuon_Run2016G/",
"SingleMuon_Run2016Hver2/",
"SingleMuon_Run2016Hver3/",
"WJetsToLNu_amcatnlo/",
"WJetsToLNu_amcatnlo_ext/",
"WJetsToLNu_amcatnlo_ext2v5/",
"WW/",
"WZ/",
"ZZ/",
"ZToEE_M120to200/",
"ZToEE_M1400to2300/",
"ZToEE_M200to400/",
"ZToEE_M2300to3500/",
"ZToEE_M3500to4500/",
"ZToEE_M400to800/",
"ZToEE_M4500to6000/",
"ZToEE_M50to120/",
"ZToEE_M6000toInf/",
"ZToEE_M800to1400/",
"crab_DoubleEG_RunB/",
"crab_DoubleEG_RunC/",
"crab_DoubleEG_RunD/",
"crab_DoubleEG_RunE/",
"crab_DoubleEG_RunF/",
"crab_DoubleEG_RunG/",
"crab_DoubleEG_RunHver2/",
"crab_DoubleEG_RunHver3/",
"ttbar/",
"ttbarBackup/",
"ttbar_M1000toInf/",
"ttbar_M700to1000/"
};
void makeDirectoryList()
{
 Long64_t size;
 //-----Loop over all directories-----// 
 for(int i=0;i<nDir;i++){
  directory[i] = baseDir+directory[i];
  cout << "Processing directory: " << directory[i] << endl;
  cout << "Size of skim_file.txt: " << endl; 
  //Get file_list.txt from each directory and save just the name of the file at 
  //the end of each line
  GetFileList(directory[i]);
  //-----System sleep for 2.0 seconds-----//
  //Without this pause between functions, the second function does not work//
  gSystem->Sleep(2000);

  //Open file from last step, add the full directory name to the beginning of each
  //Line and save with the name "skim_file.txt"
  getList(directory[i]);
  //-----Delete temporary file-----//
  //gSystem->Exec("rm "+directory[i]+"temp_file.txt");
  //Double check that the output file, skim_file.txt exists
  //I had some issues with some directories not having it
  gSystem->Exec("stat -c %s "+directory[i]+"skim_file.txt");
 }
}

void getList(TString directory){
 ifstream inFile(directory+"temp_file.txt");
 ofstream outFile(directory+"skim_file.txt");
 TString contents;
 
 while(true){
  inFile >> contents;
  if(inFile.eof()) break;
  outFile << directory+contents << endl;
 }
 inFile.close();
 outFile.close();
}

void GetFileList(TString directory){
 gSystem->Exec("cat "+directory+"file_list.txt | awk -F '/' '{print $NF}' > "+directory+"temp_file.txt");

}
