#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"
#include "/home/hep/wrtabb/git/DY-Analysis/headers/headerFunctions.h"
//forward declaration for counter
//void counter(Long64_t i, Long64_t N);

//strings for naming histgrams and for file locations
const TString mcDist = "/home/hep/wrtabb/git/DY-Analysis/data/efficiencyAndMigration.root";
const TString mcEff = "/home/hep/wrtabb/git/DY-Analysis/plots/efficiencies.root";
const TString mcSF =  "/home/hep/wrtabb/git/DY-Analysis/data/dataVsMC.root";
const TString histSaveName = "toyUnfold.root";
const TString recoName[] = {"hReco","hMCReco","hAltReco"};
const TString trueName[] = {"hTrue","hMCTrue","hAltTrue"};
const TString matrixName[] = {"hMatrixDontUse","hMatrix","hAltMatrix"};
const TString backName[] = {"hBack1","hBack2","hBack3"};

//Number of events to process
const int nEvents = 1e7;
//Number of histograms per array
const int nHists = 3;

//parameters for what to include and how to do unfolding
const bool exactClosure = false;//set exact closure
const bool effInc = false; //include efficiency
const bool backInc = false;//include toy background

void toyModel()
{
  TH1::SetDefaultSumw2();
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0); 
  gROOT->SetBatch(true);

  double massTrue,massMeasured,massBack,smear;
  double massMeasuredMC,massMeasuredData; 
  double backDistInt,totalInt,backWeight,weight;
  double effData,effMC,rand,rho;

  //file to store state of exactClosure, effInc, and BackInc for unfolding
  ofstream parameterFile;
  parameterFile.open("parameters.txt");
  parameterFile << exactClosure << " " << effInc << " " << backInc << endl;
  parameterFile.close();

  //Open files to pull MC distributions for toy models
  TFile*file = new TFile(mcDist);
  TFile*fileSF = new TFile(mcSF);
  TFile*fileEff = new TFile(mcEff);
  //Open files for scale factors and efficiencies
  TProfile*profileSF = (TProfile*)fileSF->Get("hSFvsInvMassAll_pfx");
  TEfficiency*efficiency = (TEfficiency*)fileEff->Get("Efficiency");

  //Mass smearing model
  double resSD = 3;
  double resMean = 0.0;
  double resScale = 1.0;
  TF1*fResolutionModel = new TF1("fPhiResolutionModel","gaus(0)",-20,20);
  fResolutionModel->SetParameters(resScale,resMean,resSD);

  //MC mass distribution for filling toy models
  TH1D*hMassDist = (TH1D*)file->Get("hGenInvMass"); 
   hMassDist->SetName("hMassDist");

  //MC background distribution for filling toy models
  TH1D*hBackDist;
  TH1D*hFakes = (TH1D*)fileSF->Get("hFakesInvMass");
  TH1D*hEW = (TH1D*)fileSF->Get("hEWInvMass");
  TH1D*hTops = (TH1D*)fileSF->Get("hTopsInvMass");
  //Set any negative bin values to zero
  for(int i=0;i<89;i++){
    if(hFakes->GetBinContent(i)<0.0) hFakes->SetBinContent(i,0.0);
    if(hEW->GetBinContent(i)<0.0) hEW->SetBinContent(i,0.0);
    if(hTops->GetBinContent(i)<0.0) hTops->SetBinContent(i,0.0);
  }
  hBackDist = (TH1D*)hFakes->Clone("hBackDist");
   hBackDist->Add(hEW);
   hBackDist->Add(hTops);

  //Determination of background weighting factor to closely match full MC
  if(backInc) backDistInt = hBackDist->Integral();
  else backDistInt = 0;//if background not included, backWeight = 0
  totalInt = backDistInt + hMassDist->Integral();
  weight = 1.0;//weight for signal distributions
  backWeight = weight*backDistInt/totalInt;//ratio of background to total mass 

  TRandom3*random = new TRandom3();

  //Initializing histograms
  TH1D*hReco[nHists];
  TH1D*hTrue[nHists];
  TH1D*hBack[nHists];
  TH2D*hMatrix[nHists];

  //starting loop over histograms
  Long64_t N = 0;
  for(int j=0;j<nHists;j++){
    if(!exactClosure){
      gRandom->SetSeed(j+1);
      random->SetSeed(j+1);
    }
    else{
      gRandom->SetSeed(1);
      random->SetSeed(1);
    }
    //defining histograms
    hReco[j] = new TH1D(recoName[j],"",nLogBins2,massbins2);
    hTrue[j] = new TH1D(trueName[j],"",nLogBins,massbins);
    hMatrix[j] = new TH2D(matrixName[j],"",nLogBins,massbins,nLogBins2,massbins2);
    hBack[j] = new TH1D(backName[j],"",nLogBins2,massbins2);
    
    //loop to fill distributions
    for(Long64_t i=0;i<nEvents;i++){
      //counter for keeping track of program progress
      counter(N,nHists*nEvents);
      N++;

      //parameters for filling histograms
      massTrue = hMassDist->GetRandom();
      smear = fResolutionModel->GetRandom();
      massMeasured = massTrue+smear;
      massBack = hBackDist->GetRandom();

      bool seenInData = true;
      bool seenInMC = true;

      //efficiency factors for filling histograms
      if(effInc){
        rho = profileSF->GetBinContent(profileSF->FindBin(massTrue));
        effMC = efficiency->GetEfficiency(hMassDist->FindBin(massMeasured));
        effData = effMC*rho;
        rand = random->Rndm(); 
        if(rand>effMC) seenInMC = false;
        if(rand>effData) seenInData = false; 
      }
      else rho = 1.0;
      massMeasuredData = massMeasuredMC = massMeasured;
      if(!seenInData) massMeasuredData = 0;
      if(!seenInMC) massMeasuredMC = 0;
      
      //Filling histograms
      //if(exactClosure){
      //  if(seenInMC) hReco[j]->Fill(massMeasuredMC,rho*weight);
      //  else hReco[j]->Fill(massMeasuredMC,weight);          
      //}
      //else hReco[j]->Fill(massMeasuredData,weight);
      hReco[j]->Fill(massMeasuredData,weight);
      hTrue[j]->Fill(massTrue,weight);
      hBack[j]->Fill(massBack,backWeight);
      if(seenInMC){
        hMatrix[j]->Fill(massTrue,massMeasuredMC,weight*rho);
        hMatrix[j]->Fill(massTrue,0.0,weight*(1-rho));
      }
      else hMatrix[j]->Fill(massTrue,massMeasuredMC,weight);
    }
}//end j loop
  //Adding background to reco distribution
  //note: reco[0] is distribution with background included
  //reco[1] and reco[2] are signal only
  hReco[0]->Add(hBack[0]);

  //determination of canvas save name based on parameters used
  TString saveName;
  saveName = 
    "/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase1/step1MigrationMatrix";
  if(effInc) saveName += "_EffInc";
  else saveName += "_NoEff";
  if(exactClosure) saveName += "_ClosureTest";
  else saveName += "_NoClosure";
  if(backInc) saveName += "_BackInc";
  else saveName += "_NoBack";

  TCanvas*canvas=new TCanvas("canvas","",10,10,1000,1000);
  canvas->SetLogy();
  canvas->SetLogx();
  hMatrix[1]->GetXaxis()->SetMoreLogLabels();
  hMatrix[1]->GetXaxis()->SetNoExponent();
  hMatrix[1]->GetYaxis()->SetMoreLogLabels();
  hMatrix[1]->GetYaxis()->SetNoExponent();
  hMatrix[1]->GetXaxis()->SetTitle("true mass [GeV]");
  hMatrix[1]->GetYaxis()->SetTitle("reco mass [GeV]");
  hMatrix[1]->GetYaxis()->SetTitleOffset(1.5);
  hMatrix[1]->Draw("colz");
  saveName += ".png";
  canvas->SaveAs(saveName);
  TFile*file2 = new TFile(histSaveName,"recreate");
  file2->cd();
  for(int i=0;i<nHists;i++){
    hMatrix[i]->Write();
    hReco[i]->Write();
    hTrue[i]->Write();
    hBack[i]->Write();
  }
  hMassDist->Write();
  file2->Write();
  file2->Close();
  }
/*
void counter(Long64_t i, Long64_t N)
{
  
  Long64_t P = 100*(i)/(N);
  TTimeStamp eventTimeStamp;
  if(i%(N/100)==0) {
    cout << "toyModel.C " << "[Time: " << eventTimeStamp.AsString("s") << "] " << P
    << "%" << endl;
  }
return;
}
*/
