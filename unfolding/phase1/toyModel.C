#include "TRandom.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TTimeStamp.h"

//forward declaration for counter
void counter(Long64_t i, Long64_t N);

//defining bins for histograms
const int nLogBins = 43;
const int nLogBins2 = 2*nLogBins;
const float  massbins[] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81,
 86, 91,96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,200, 220, 243, 273,
 320, 380, 440, 510, 600, 700, 830, 1000, 1500,3000};
const float massbins2[] = {15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,
 57.5,60,62,64,66,68,70,72,74,76,78.5,81,83.5,86,88.5,91,93.5,96,
 98.5,101,103.5,106,108,110,112.5,115,117.5,120,123,126,129.5,133,
 137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,210,220,
 231.5,243,258,273,296.5,320,350,380,410,440,475,510,555,600,650,
 700,765,830,915,1000,1250,1500,2250,3000};

//strings for naming histgrams and for file locations
const TString mcDist = "/home/hep/wrtabb/git/DY-Analysis/plots/plotsDY.root";
const TString mcSF =  "/home/hep/wrtabb/git/DY-Analysis/plots/dataVsMC.root";
const TString histSaveName = "toyUnfold.root";
const TString recoName[] = {"hReco","hMCReco","hAltReco"};
const TString trueName[] = {"hTrue","hMCTrue","hAltTrue"};
const TString matrixName[] = {"hMatrixDontUse","hMatrix","hAltMatrix"};
const TString backName[] = {"hBack1","hBack2","hBack3"};

//Number of events to process
const int nEvents = 1e8;
//Number of histograms per array
const int nHists = 3;

//parameters for what to include and how to do unfolding
const bool exactClosure = true;//set exact closure
const bool effInc = true; //include efficiency
const bool backInc = true;//include toy background

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

  //Open files for scale factors and efficiencies
  TProfile*profileSF = (TProfile*)fileSF->Get("hSFvsInvMassAll_pfx");
  TEfficiency*efficiency = (TEfficiency*)file->Get("Efficiency");

  //Mass smearing model
  TF1*fResolutionModel = new TF1("fPhiResolutionModel","gaus(0)",-20,20);
  fResolutionModel->SetParameters(1,0,3);

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
      if(exactClosure){
        if(seenInMC) hReco[j]->Fill(massMeasuredMC,rho*weight);
        else hReco[j]->Fill(massMeasuredMC,weight);          
      }
      else hReco[j]->Fill(massMeasuredData,weight);
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
    "/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase1Plots/step1MigrationMatrix";
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

void counter(Long64_t i, Long64_t N)
{
  int P = 100*(i)/(N);
  TTimeStamp eventTimeStamp;
  if(i%(N/100)==0) {
    cout << "toyModel.C " << "[Time: " << eventTimeStamp.AsString("s") << "] " << P
    << "%" << endl;
  }
return;
}
