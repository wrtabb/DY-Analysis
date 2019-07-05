
const TString pileupFileName = "/home/hep/wrtabb/git/DY-Analysis/data/pileup/MyDataPileupHistogram64000.root";
const TString dataFileName = "/home/hep/wrtabb/git/DY-Analysis/data/dataVsMC.root";
const TString treeName = "recoTree/DYTree";
const int nPileupBins = 100;
const float pileupBinLow = 0;
const float pileupBinHigh = nPileupBins;
void counter(Long64_t i, Long64_t N);
//obtained from https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py#L25
//const double pileupMC[nPileupBins] = {1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,
//					0.000140973 , 0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,
//					0.00919534 ,0.0146697 , 0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,
//					0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 , 0.0559937 ,0.0554468 ,
//					0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,
//					0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,
//					0.0142498 ,0.012804 , 0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,
//					0.00829292 ,0.0076195 ,0.0069806 ,0.0062025, 0.00546581 ,0.00484127 ,
//					0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 , 0.00117884 ,
//					0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,
//					9.88128e-05, 6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,
//					1.73032e-05 ,1.435e-05 ,1.36486e-05, 1.35555e-05 ,1.37491e-05 ,
//					0.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 , 1.34177e-05 ,
//					1.32959e-05 ,1.33287e-05};

enum ChainNum {      
  EE10to50,
  EE50to100,
  EE100to200,
  EE200to400,
  EE400to500,
  EE500to700,
  EE700to800,
  EE800to1000,
  EE1000to1500,
  EE1500to2000,
  EE2000to3000
}; 
const int numChains = 11;
int nPileUp;

void pileupRatio()
{
  TTimeStamp ts_start;
  cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
  TStopwatch totaltime;
  totaltime.Start();
  TBranch*b_nPileUp;
  TString dirNames[numChains] = {
    "/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE/crab_DYLL_M10to50",
    "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_truncated_M50To100/EE",
    "/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE/crab_DYLL_M100to200",
    "/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE",
    "/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/EE"
  };
  
  TString baseDirectory =
    "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/ikrav/DrellYan_13TeV_2016/v2p3";
  TString subDirectoryDY = "/DYJetsToLL_allMasses_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8";
  TChain*chains[numChains];
  
  vector <TString> *subFiles[numChains];  
  for(int iChain=0;iChain<numChains;iChain++) {
    
    subFiles[iChain] = new vector<TString>;
    if(iChain==EE10to50) {
      subFiles[iChain]->push_back(dirNames[iChain]+"_ext1v1");
      subFiles[iChain]->push_back(dirNames[iChain]+"_v1");
      subFiles[iChain]->push_back(dirNames[iChain]+"_v2");
    }
    else if(iChain==EE100to200) {
      subFiles[iChain]->push_back(dirNames[iChain]);
      subFiles[iChain]->push_back(dirNames[iChain]+"_ext");
    }
    else subFiles[iChain]->push_back(dirNames[iChain]);      
  } 
  int nFiles;
  const int maxFiles = 10;
  TString files;  
  Long64_t subDirectorySize;
  Long64_t totalentries = 0;
  for(int iChain=0;iChain<numChains;iChain++) {          
    chains[iChain] = new TChain(treeName);
    //if(!(iChain==EE2000to3000))continue;
    nFiles = 0;
    subDirectorySize = subFiles[iChain]->size();
    
    for(int k=0;k<subDirectorySize;k++)
      {	  	      
	TFileCollection filecoll("dum");//Object for creating a list of files in a directory
	files = baseDirectory;
	files+=subDirectoryDY;
	files+=subFiles[iChain]->at(k);
	files+="/*.root";	  
	filecoll.Add(files);
	chains[iChain]->AddFileInfoList(filecoll.GetList());
	cout << files << endl;
	cout << chains[iChain]->GetEntries() << " events loaded" << endl;	
	if(chains[iChain]->GetEntries()==0){
	  cout << endl;
	  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	  cout << "ERROR: Broken files or files not found in: " << endl;
	  cout << files << endl;
	  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	  cout << endl;
	  return;
	}
      }
    
    totalentries+=chains[iChain]->GetEntries(); 
    
    //Setting addresses for branches
    chains[iChain]->SetBranchAddress("nPileUp", &nPileUp, &b_nPileUp);   
    
  }//end iChain loop
  
  cout << "Total Events Loaded: " << totalentries << endl;
  cout << endl;

  float norm = 1.0;
  TFile*filePU = new TFile(pileupFileName);
  TH1D*hPileupData = (TH1D*)filePU->Get("pileup");
  hPileupData->SetName("hPileupData");
  hPileupData->GetXaxis()->SetTitle("Pileup");
  hPileupData->SetLineColor(kBlue);
  hPileupData->SetLineWidth(2);

  TH1D*hPileupMC = new TH1D("hPileupMC","",nPileupBins,pileupBinLow,pileupBinHigh);
  hPileupMC->Sumw2();
  hPileupMC->GetXaxis()->SetTitle("Pileup");
  hPileupMC->SetLineColor(kRed);
  hPileupMC->SetLineWidth(2);

  Long64_t nentries;
  Long64_t count = 0;
  double localEntry;
  
  for(int iChain=0;iChain<numChains;iChain++){
    //Event loop
    nentries = chains[iChain]->GetEntries();   
    for(Long64_t i=0;i<nentries;i++) {      
      counter(count,totalentries);
      count = count+1; 
      localEntry = chains[iChain]->LoadTree(i);
      b_nPileUp->GetEntry(localEntry);
      hPileupMC->Fill(nPileUp);
    }//end event loop 
  }
  
  double scaleData = norm/hPileupData->Integral();
  double scaleMC = norm/hPileupMC->Integral();
  hPileupData->Scale(scaleData);
  hPileupMC->Scale(scaleMC);
  
  TLegend*legend2 = new TLegend(0.65,0.9,0.9,0.75);
  legend2->SetTextSize(0.02);
  legend2->AddEntry(hPileupData,"Data");
  legend2->AddEntry(hPileupMC,"MC");

  TH1D*hPileupRatio = (TH1D*)hPileupData->Clone();
  hPileupRatio->SetName("hPileupRatio");
  hPileupRatio->Divide(hPileupMC);
  hPileupRatio->SetLineColor(kBlack);
  hPileupRatio->SetTitle("Data/MC Ratio");

  TCanvas*canvas3 = new TCanvas("canvas3","",10,10,1400,700);
  canvas3->Divide(2);
  canvas3->cd(1);
  hPileupData->Draw("hist");
  hPileupMC->Draw("hist,same");
  legend2->Draw("same");
  canvas3->cd(2);
  hPileupRatio->Draw("hist");  
  
  canvas3->SaveAs("./plots/pileup.png");
  TFile*pileupSaveFile = new TFile("./plots/pileup.root","RECREATE");
  pileupSaveFile->cd();
  hPileupData->Write();
  hPileupMC->Write();
  hPileupRatio->Write();
  canvas3->Write();
  pileupSaveFile->Write();
  pileupSaveFile->Close();
}

//Counter for tracking program progress
void counter(Long64_t i, Long64_t N)
{
  int P = 100*(i)/(N);  
  TTimeStamp eventTimeStamp;
  if(i%(N/100)==0) {
    cout << "pileupRatio.C " << "[Time: " << eventTimeStamp.AsString("s") << "] " << P << "%" << endl;
  }
  return;
}
