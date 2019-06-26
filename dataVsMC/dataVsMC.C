#include "/home/hep/wrtabb/git/DY-Analysis/headers/headerDataVsMC.h"
#include "/home/hep/wrtabb/git/DY-Analysis/headers/ntupleSkimLocation.h"


void counter(Long64_t i, Long64_t N);
TLorentzVector getDielectronP4(double pt1,double eta1,double phi1,double m1,double pt2,
 double eta2,double phi2,double m2);
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2);

enum ChainNum {
  QCD20to30,
  QCD30to50,
  QCD50to80,
  QCD80to120,
  QCD120to170,
  QCD170to300,
  QCD300toInf,  
  wJets,
  WW,    
  WWTo2L2Nu,
  ZZ,       
  ZZTo4L,
  WZ,
  WZTo3LNu,
  tt0to700,
  tt700to1000,
  tt1000toInf,
  tW,             
  tbarW,         
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
  EE2000to3000,
  TAUTAU10to50,
  TAUTAU50to100,
  TAUTAU100to200,
  TAUTAU200to400,
  TAUTAU400to500,
  TAUTAU500to700,
  TAUTAU700to800,
  TAUTAU800to1000,
  TAUTAU1000to1500,
  TAUTAU1500to2000,
  TAUTAU2000to3000,
  DataRunB,
  DataRunC,
  DataRunD,
  DataRunE,
  DataRunF,
  DataRunG,
  DataRunH
}; 

enum HistBins {
  FAKES,
  EW,
  TOPS,
  EE,
  DATA,
  UNDEF = -1
};
enum HistTypes {
  INV_MASS,
  INV_MASS_LINEAR,
  VERTICES,
  VERTICES_WEIGHTED,
  PT_LEAD,
  PT_SUB,
  PT_DI,
  ETA_LEAD,
  ETA_SUB,
  ETA_DI,
  RAPIDITY,
  INV_MASS0,
  INV_MASS1,
  INV_MASS2,
  INV_MASS3,
  INV_MASS_SCALED
};

void dataVsMC()
{
  TTimeStamp ts_start;
  cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
  TStopwatch totaltime;
  totaltime.Start();
  bool isMC; //is Monte Carlo
  gStyle->SetOptStat(0);
  //Defining branches
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
  
  //Loading ntuples
  cout << "Loading ntuples" << endl;
  //The names of every directory being loaded
  TString dirNames[numChains] = {QCDReco20to30,QCDReco30to50,QCDReco50to80,QCDReco80to120,
   QCDReco120to170,QCDReco170to300,QCDReco300andUp,WJetsReco,WWReco,WWtoLNuReco,ZZReco,
   ZZtoLReco,WZReco,WZtoLNuReco,TT0to700Reco,TT700to1000Reco,TT1000andUpReco,tWReco,tWantiReco,
   EEM10to50Reco,EEM50to100Reco,EEM100to200Reco,EEM200to400Reco,EEM400to500Reco,
   EEM500to700Reco,EEM700to800Reco,EEM800to1000Reco,EEM1000to1500Reco,EEM1500to2000Reco,
   EEM2000to3000Reco,Taus10to50Reco,Taus50to100Reco,Taus100to200Reco,Taus200to400Reco,
   Taus400to500Reco,Taus500to700Reco,Taus700to800Reco,Taus800to1000Reco,Taus1000to1500Reco,
   Taus1500to2000Reco,Taus2000to3000Reco,runB,runC,runD,runE,runF,runG,runH};

  TChain*chains[numChains];

  vector <TString> *subFiles[numChains];  
  for(int iChain=0;iChain<numChains;iChain++) {
    //Not using WWTo2L2Nu, ZZTo4L, or WZTo3LNu
    if(iChain==WWTo2L2Nu||iChain==ZZTo4L||iChain==WZTo3LNu) continue;
    //Skipping QCD samples due to possible problems with them
    if(iChain==QCD20to30||iChain==QCD30to50||iChain==QCD50to80||iChain==QCD80to120||
       iChain==QCD120to170||iChain==QCD170to300||iChain==QCD300toInf) continue;
    
    subFiles[iChain] = new vector<TString>;
    if(iChain==EE10to50||iChain==TAUTAU10to50) {
      subFiles[iChain]->push_back(dirNames[iChain]+"_ext1v1");
      subFiles[iChain]->push_back(dirNames[iChain]+"_v1");
      subFiles[iChain]->push_back(dirNames[iChain]+"_v2");
    }
    else if(iChain==EE100to200||iChain==TAUTAU100to200||iChain==wJets) {
      subFiles[iChain]->push_back(dirNames[iChain]);
      subFiles[iChain]->push_back(dirNames[iChain]+"_ext");
    }
    else if(iChain==DataRunH) {
      subFiles[iChain]->push_back(dirNames[iChain]+"ver2");
      subFiles[iChain]->push_back(dirNames[iChain]+"ver3");
    }
    else if(iChain==tt0to700) {
      subFiles[iChain]->push_back(dirNames[iChain]);
      subFiles[iChain]->push_back(dirNames[iChain]+"Backup");
    }
    else subFiles[iChain]->push_back(dirNames[iChain]);      
  } 
  
  TString files;  
  Long64_t subDirectorySize;
  Long64_t totalentries = 0;
  for(int iChain=0;iChain<numChains;iChain++) {          
    if(iChain==WWTo2L2Nu||iChain==ZZTo4L||iChain==WZTo3LNu) continue;
    if(iChain==QCD20to30||iChain==QCD30to50||iChain==QCD50to80||iChain==QCD80to120||
       iChain==QCD120to170||iChain==QCD170to300||iChain==QCD300toInf) continue;
    chains[iChain] = new TChain(treeName);
    if(iChain==DataRunB||iChain==DataRunC||iChain==DataRunD||iChain==DataRunE||
       iChain==DataRunF||iChain==DataRunG||iChain==DataRunH) isMC = kFALSE;
    else isMC = kTRUE;
    subDirectorySize = subFiles[iChain]->size();
    for(int k=0;k<subDirectorySize;k++) {	  	      
      TFileCollection filecoll("dum");//Object for creating a list of files in a directory
      files=subFiles[iChain]->at(k);
      files+="/skims_0002/*.root";      
      filecoll.Add(files);
      chains[iChain]->AddFileInfoList(filecoll.GetList());
      cout << files << endl;
      cout << chains[iChain]->GetEntries() << " events loaded" << endl;	 
      
      if(chains[iChain]->GetEntries()==0){//error message if no events loaded
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
    chains[iChain]->SetBranchAddress("Nelectrons", &Nelectrons, &b_Nelectrons);
    chains[iChain]->SetBranchAddress("nVertices", &nVertices, &b_nVertices);
    chains[iChain]->SetBranchAddress("nPileUp", &nPileUp, &b_nPileUp);
    chains[iChain]->SetBranchAddress("Electron_pT", &Electron_pT, &b_Electron_pT);
    chains[iChain]->SetBranchAddress("Electron_eta",&Electron_eta, &b_Electron_eta);
    chains[iChain]->SetBranchAddress("Electron_phi",&Electron_phi, &b_Electron_phi);
    chains[iChain]->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,
      &b_Electron_passMediumID);
    chains[iChain]->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
    chains[iChain]->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
    chains[iChain]->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
    chains[iChain]->SetBranchAddress("HLT_trigName",&pHLT_trigName);   
    if(isMC)//No gen weighting for data 
      chains[iChain]->SetBranchAddress("GENEvt_weight",&GENEvt_weight,&b_GENEvt_weight);
    else continue;
  }//end iChain loop
  
  cout << "Total Events Loaded: " << totalentries << endl;
  cout << endl;
  
  //defining histograms
  TH1F*histos[nHistoTypes][nHistos];  
  for(int i=0;i<nHistoTypes;i++){//type of histogram
    for(int j=0;j<nHistos;j++){//histogram within type
      if(i==VERTICES_WEIGHTED&&j==DATA) continue; //data doesn't get weighted
      if(i==INV_MASS||i==INV_MASS0||i==INV_MASS1||i==INV_MASS2||i==INV_MASS3||i==INV_MASS_SCALED) 
        histos[i][j]=new TH1F(histName[j]+histTypeName[i],"",nBinsInvMass,massbins);   
      else histos[i][j]=new TH1F(histName[j]+histTypeName[i],"",nBins[i],binLow[i],binHigh[i]);      
      if(j==DATA) {
	histos[i][j]->SetLineColor(kBlack);
	histos[i][j]->SetMarkerColor(kBlack);
	histos[i][j]->SetMarkerSize(1);
	histos[i][j]->SetMarkerStyle(20);
      }
      else {
	histos[i][j]->SetFillColor(histFillColors[j]);
	histos[i][j]->SetLineColor(histLineColors[j]);
      }
      
      histos[i][j]->GetXaxis()->SetTitle(xAxisLabels[i]);
    }
  }

  TFile*pileupRatioFile  = new TFile(pileupRatioName);
  TH1F*hPileupRatio = (TH1F*)pileupRatioFile->Get("hPileupRatio");
  TFile*fileLeg2SF = new TFile(leg2SFName);
  TH2F*hLeg2SF = (TH2F*)fileLeg2SF->Get("EGamma_SF2D");
  TFile*fileMedIDSF = new TFile(medIDSFName);
  TH2F*hMedIDSF = (TH2F*)fileMedIDSF->Get("EGamma_SF2D");
  TFile*fileRecoSF = new TFile(recoSFName); 
  TH2F*hRecoSF = (TH2F*)fileRecoSF->Get("EGamma_SF2D");
  TFile*pvzFile = new TFile(pvzFileName);
  TH1D*hPVzSF = (TH1D*)pvzFile->Get("PVz_SF");
  
  TH2F*hSFvsInvMassAll = new TH2F("hSFvsInvMassAll","",nBinsInvMass,massbins,150,0,1.5);
  TH2F*hSFvsInvMassHLT = new TH2F("hSFvsInvMassHLT","",nBinsInvMass,massbins,150,0,1.5);
  TH2F*hSFvsInvMassReco = new TH2F("hSFvsInvMassReco","",nBinsInvMass,massbins,150,0,1.5);
  TH2F*hSFvsInvMassID = new TH2F("hSFvsInvMassID","",nBinsInvMass,massbins,150,0,1.5);

  cout << "Starting Event Loop" << endl;
  double varGenWeight, lumiEffective, nEffective, localEntry, sumGenWeight, sumRawGenWeight, 
    totalWeight, sfWeight, weightNoPileup, xSecWeight, genWeight, pileupWeight,pvzWeight, 
    xSecWeightAlone;
  Long64_t nentries;
  Long64_t count = 0;
  int sampleCategory = -1;
  TString compareHLT = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
  TString trigName;
  int trigNameSize;
  double lumi = dataLuminosity;//luminosity for xsec weighting
  ofstream genWeightFile;
  genWeightFile.open("genWeights.txt");
  genWeightFile << "iChain, Nevents, Neffective, LumiEffective" << endl;
  double sfReco1,sfReco2,sfID1,sfID2,sfHLT;//efficiency scale factors
  double eEta1, eEta2, ePt1, ePt2;

  TAxis*xAxis;
  double binx,binWidth;   
  Long64_t nData = 0;
  //Loop over samples
  for(int iChain=0;iChain<numChains;iChain++) {
    if(iChain==WWTo2L2Nu||iChain==ZZTo4L||iChain==WZTo3LNu) continue;
    if(iChain==QCD20to30||iChain==QCD30to50||iChain==QCD50to80||iChain==QCD80to120||
       iChain==QCD120to170||iChain==QCD170to300||iChain==QCD300toInf) continue;
    //determine if chain is MC or Data
    if(iChain==DataRunB||iChain==DataRunC||iChain==DataRunD||iChain==DataRunE||
       iChain==DataRunF||iChain==DataRunG||iChain==DataRunH) isMC = false; 
    else isMC = true;

    cout << endl;
    cout << "Processing chain: " << dirNames[iChain] << endl;
    cout << endl;
    
    //Determine which sample is being looped over
    if(!isMC) sampleCategory = DATA;
    else if (isMC){
      if(iChain==wJets||iChain==QCD20to30||iChain==QCD30to50||iChain==QCD50to80||
         iChain==QCD80to120||iChain==QCD120to170||iChain==QCD170to300||iChain==QCD300toInf) 
	   sampleCategory = FAKES;
      else if(iChain==WW||iChain==ZZ||iChain==WZ||iChain==WWTo2L2Nu||iChain==ZZTo4L||
              iChain==WZTo3LNu||iChain==TAUTAU10to50||iChain==TAUTAU50to100||
              iChain==TAUTAU100to200||iChain==TAUTAU200to400||iChain==TAUTAU400to500||
              iChain==TAUTAU500to700||iChain==TAUTAU700to800||iChain==TAUTAU800to1000||
	      iChain==TAUTAU1000to1500||iChain==TAUTAU1500to2000||iChain==TAUTAU2000to3000)
	        sampleCategory = EW;
      else if(iChain==tt0to700||iChain==tt700to1000||iChain==tt1000toInf||iChain==tW||
              iChain==tbarW)
	        sampleCategory = TOPS;
      else if(iChain==EE10to50||iChain==EE50to100||iChain==EE100to200||iChain==EE200to400||
	      iChain==EE400to500||iChain==EE500to700||iChain==EE700to800||iChain==EE800to1000||
	      iChain==EE1000to1500||iChain==EE1500to2000||iChain==EE2000to3000)
	        sampleCategory = EE;
    }
    if(sampleCategory==UNDEF) {
      cout << "Tried to load a chain which does not exist!" << endl;
      return;
    }
    
    nentries = chains[iChain]->GetEntries();
    cout << "Sample: " << sampleCategory << endl;
    cout << "Number of Entries: " << nentries << endl; 
    //Finding normalized genWeights,sums,variances
    sumGenWeight = 0.0;
    sumRawGenWeight = 0.0;
    varGenWeight = 0.0;
    if(isMC){ //Only MonteCarlo gets gen weights
      for(Long64_t i=0;i<nentries;i++){
	localEntry = chains[iChain]->LoadTree(i);
	b_GENEvt_weight->GetEntry(localEntry);
	genWeight = GENEvt_weight/fabs(GENEvt_weight);	//normalized genweight
	sumGenWeight += genWeight;
	varGenWeight += GENEvt_weight*GENEvt_weight; //variance of genweights
	sumRawGenWeight += GENEvt_weight; 
      }          
      nEffective = (sumRawGenWeight*sumRawGenWeight)/varGenWeight;
      lumiEffective = nEffective/xSec[iChain];
      genWeightFile << iChain << ", " << nentries << ", " << nEffective << ", " << 
        lumiEffective << endl;
    }//end if isMC
    
    //Event loop
    for(Long64_t i=0;i<nentries;i++) {      
      counter(count,totalentries);
      count = count+1; 
      chains[iChain]->GetEntry(i);

      //Weights
      totalWeight = 1.0;
      sfWeight = 1.0;
      weightNoPileup = 1.0;
      xSecWeight = 1.0;
      genWeight = 1.0;
      pileupWeight = 1.0;
      xSecWeightAlone = 1.0; 
      pvzWeight = 1.0;
      //HLT cut
      trigNameSize = pHLT_trigName->size();
      bool passHLT = kFALSE;	  
      for(int iHLT=0;iHLT<trigNameSize;iHLT++) {
	trigName = pHLT_trigName->at(iHLT);	  
	if(trigName.CompareTo(compareHLT)==0) {
	  if(HLT_trigFired[iHLT]==1) {
	    passHLT = kTRUE;	
	  }
	  else {
	    passHLT = kFALSE;
	  }		     
	  break; 
	}
      } 
      if(!passHLT) continue;
     
      bool passNumEle = kFALSE;
      int numDielectrons = 0;
      int subEle = -1;
      int leadEle = -1;
      double invMass = -5000;
      double rapidity = -5000;
      double dielectronPt = -5000;
      double dielectronEta = -5000;
      TLorentzVector dielectronP4;

      //Electron loop
      for(int iEle = 0; iEle < Nelectrons; iEle++) {
	if(!Electron_passMediumID[iEle]) continue;
	for(int jEle = iEle+1; jEle < Nelectrons; jEle++) {
	  if(!Electron_passMediumID[jEle]) continue;
	  if(passDileptonKinematics(Electron_pT[iEle],Electron_pT[jEle],Electron_eta[iEle],
	    Electron_eta[jEle])){
            numDielectrons++;
	    if(Electron_pT[iEle]>Electron_pT[jEle]){
	      leadEle = iEle; subEle = jEle;
	    }
	    else {
	      leadEle = jEle; subEle = iEle;
	    }	  
          }
        }//end jEle loop
      }//end iEle loop
       
     if(numDielectrons==1){
        passNumEle = kTRUE;
     }
     if(!passNumEle) continue;
     if(leadEle<0||subEle<0) continue;
     dielectronP4 = getDielectronP4(Electron_pT[leadEle],Electron_eta[leadEle],
       Electron_phi[leadEle],eMass,Electron_pT[subEle],Electron_eta[subEle],
       Electron_phi[subEle],eMass);
     invMass=dielectronP4.M();
     rapidity=dielectronP4.Rapidity();
     dielectronPt=dielectronP4.Pt();
     dielectronEta=dielectronP4.Eta();  
     if(invMass<-1000||rapidity<-1000||dielectronPt<-1000||dielectronEta<-1000) continue;
     
     eEta1 = Electron_eta[leadEle];
     eEta2 = Electron_eta[subEle];
     ePt1 = Electron_pT[leadEle];
     ePt2 = Electron_pT[subEle];
     
     if(ePt1<ptBinLow) ePt1 = ptBinLow;
     if(ePt2<ptBinLow) ePt2 = ptBinLow;
     if(ePt1>ptBinHigh) ePt1 = ptBinHigh;
     if(ePt2>ptBinHigh) ePt2 = ptBinHigh;

     //Determining weighting factors
     if(isMC){
       xSecWeight=lumi*(xSec[iChain]/1.0);//xSecWeight when used with genWeight 
       xSecWeightAlone = lumi*(xSec[iChain]/nentries);//xSecWeight when used without genWeight
       pileupWeight = hPileupRatio->GetBinContent(hPileupRatio->FindBin(nPileUp));
       genWeight = (GENEvt_weight/fabs(GENEvt_weight))/sumGenWeight;
       //pvzWeight = hPVzSF->GetBinContent(hPVzSF->FindBin(PVz));//need ntuples with PVz 
       if(sampleCategory==FAKES){
         totalWeight = genWeight*xSecWeight*pileupWeight*pvzWeight;
         weightNoPileup = genWeight*xSecWeight;
       }
       else{
         sfReco1=hRecoSF->GetBinContent(hRecoSF->FindBin(eEta1,ePt1));
         sfReco2=hRecoSF->GetBinContent(hRecoSF->FindBin(eEta2,ePt2));
         sfID1=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eEta1,ePt1));
         sfID2=hMedIDSF->GetBinContent(hMedIDSF->FindBin(eEta2,ePt2));
         sfHLT=(hLeg2SF->GetBinContent(hLeg2SF->FindBin(eEta1,ePt1)))*
           (hLeg2SF->GetBinContent(hLeg2SF->FindBin(eEta2,ePt2)));
         sfWeight = sfReco1*sfReco2*sfID1*sfID2*sfHLT;
         totalWeight = sfWeight*genWeight*xSecWeight*pileupWeight;
         weightNoPileup = sfWeight*genWeight*xSecWeight;

         hSFvsInvMassAll->Fill(invMass,sfWeight);
         hSFvsInvMassHLT->Fill(invMass,sfHLT);
         hSFvsInvMassID->Fill(invMass,sfID1*sfID2);
         hSFvsInvMassReco->Fill(invMass,sfReco1*sfReco2);
       }
     }//end if isMC
    
     xAxis = histos[INV_MASS_SCALED][sampleCategory]->GetXaxis();
     binx = xAxis->FindBin(invMass);
     binWidth = histos[INV_MASS_SCALED][sampleCategory]->GetXaxis()->GetBinWidth(binx);
 
     histos[INV_MASS0][sampleCategory]->Fill(invMass,xSecWeightAlone);
     histos[INV_MASS1][sampleCategory]->Fill(invMass,xSecWeight*genWeight);
     histos[INV_MASS2][sampleCategory]->Fill(invMass,xSecWeight*pileupWeight*genWeight);
     histos[INV_MASS3][sampleCategory]->
       Fill(invMass,xSecWeight*pileupWeight*genWeight*sfWeight);
     histos[INV_MASS][sampleCategory]->Fill(invMass,totalWeight);
     histos[INV_MASS_SCALED][sampleCategory]->Fill(invMass,totalWeight/binWidth);
     if(invMass<invMassLow||invMass>invMassHigh) continue;
     histos[INV_MASS_LINEAR][sampleCategory]->Fill(invMass,totalWeight);
     histos[PT_LEAD][sampleCategory]->Fill(Electron_pT[leadEle],totalWeight);
     histos[PT_SUB][sampleCategory]->Fill(Electron_pT[subEle],totalWeight);
     histos[PT_DI][sampleCategory]->Fill(dielectronPt,totalWeight);	    
     histos[ETA_LEAD][sampleCategory]->Fill(Electron_eta[leadEle],totalWeight);
     histos[ETA_SUB][sampleCategory]->Fill(Electron_eta[subEle],totalWeight);
     histos[ETA_DI][sampleCategory]->Fill(dielectronEta,totalWeight);
     histos[RAPIDITY][sampleCategory]->Fill(rapidity,totalWeight);
     histos[VERTICES][sampleCategory]->Fill(nVertices,weightNoPileup);
     if(!isMC){
      nData++;
      continue; //no Vertex weighting for data
     }
     histos[VERTICES_WEIGHTED][sampleCategory]->Fill(nVertices,totalWeight);
    }//end event loop   
  }//end chain loop 
  genWeightFile.close();
  
  TFile *rootFile = new TFile("/home/hep/wrtabb/git/DY-Analysis/data/dataVsMC.root","RECREATE");
  rootFile->cd();
  for(int i=0;i<nHistoTypes;i++) {
    for(int j=0;j<nHistos;j++){
      if(i==VERTICES_WEIGHTED&&j==DATA) continue;
      histos[i][j]->Write();
    }
  }
  
  TProfile*hSFvsInvMassProfAll = hSFvsInvMassAll->ProfileX();
  hSFvsInvMassProfAll->GetXaxis()->SetMoreLogLabels();
  hSFvsInvMassProfAll->GetXaxis()->SetNoExponent();
  hSFvsInvMassProfAll->SetTitle("#LTSF#GT vs.Invariant Mass");
  hSFvsInvMassProfAll->Write();
  TProfile*hSFvsInvMassProfHLT = hSFvsInvMassHLT->ProfileX();
  hSFvsInvMassProfHLT->GetXaxis()->SetMoreLogLabels();
  hSFvsInvMassProfHLT->GetXaxis()->SetNoExponent();
  hSFvsInvMassProfHLT->SetTitle("#LTSF#GT vs.Invariant Mass");
  hSFvsInvMassProfHLT->Write();
  TProfile*hSFvsInvMassProfReco = hSFvsInvMassReco->ProfileX();
  hSFvsInvMassProfReco->GetXaxis()->SetMoreLogLabels();
  hSFvsInvMassProfReco->GetXaxis()->SetNoExponent();
  hSFvsInvMassProfReco->SetTitle("#LTSF#GT vs.Invariant Mass");
  hSFvsInvMassProfReco->Write();
  TProfile*hSFvsInvMassProfID = hSFvsInvMassID->ProfileX();
  hSFvsInvMassProfID->GetXaxis()->SetMoreLogLabels();
  hSFvsInvMassProfID->GetXaxis()->SetNoExponent();
  hSFvsInvMassProfID->SetTitle("#LTSF#GT vs.Invariant Mass");
  hSFvsInvMassProfID->Write();
  rootFile->Write();
  rootFile->Close();
  
  totaltime.Stop();
  Double_t TotalCPURunTime = totaltime.CpuTime();
  Double_t TotalRunTime = totaltime.RealTime();
  TTimeStamp ts_end;
  cout << endl;
  cout << "**************************************************************************" << endl;
  cout << "Total CPU RunTime: " << TotalCPURunTime/60 << " minutes" << endl;
  cout << "Total Real RunTime: " << TotalRunTime/60 << " minutes" << endl;
  cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;   
  cout << "Number of Events Processed: " << count << endl;
  cout << "Number of data events Processed: " << nData << endl;
  cout << "**************************************************************************" << endl;
  cout << endl;
  
}//end main function

//Counter for tracking program progress
void counter(Long64_t i, Long64_t N)
{
  int P = 100*(i)/(N);  
  TTimeStamp eventTimeStamp;
  if(i%(N/100)==0) {
    cout << "dataVsMC.C " << "[Time: " << eventTimeStamp.AsString("s") << "] " << P 
      << "%" << endl;
  }
  return;
}

//Kinematic cuts
bool passDileptonKinematics(double pt1,double pt2,double eta1,double eta2)
{
  if(abs(eta1)>etaGapLow && abs(eta1)<etaGapHigh) return kFALSE;
  if(abs(eta2)>etaGapLow && abs(eta2)<etaGapHigh) return kFALSE; 
  if(abs(eta1)>etaHigh||abs(eta2)>etaHigh) return kFALSE; 
  if(!((pt1>ptLow && pt2>ptHigh)||(pt1>ptHigh && pt2>ptLow))) return kFALSE;
  return kTRUE;
}

TLorentzVector getDielectronP4(double pt1,double eta1,double phi1,double m1,double pt2,
  double eta2,double phi2,double m2)
{
  TLorentzVector vElectron1;
  TLorentzVector vElectron2;
  vElectron1.SetPtEtaPhiM(pt1,eta1,phi1,m1);
  vElectron2.SetPtEtaPhiM(pt2,eta2,phi2,m2);
  return vElectron1+vElectron2;
}
