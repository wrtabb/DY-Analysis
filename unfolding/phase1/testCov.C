#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"
#include "TLatex.h"
using namespace std;

enum Reglarization {//Strength of regularization
  NO_REG,           //No regularization
  CONST_REG,        //User defined regularization
  VAR_REG           //TUnfoldDensity determines best choice of regularization strength
};

const int nBins = 43;
const int binLow = 15;
const int binHigh = 3000;
//File Names
const TString fileName= "toyUnfold.root";       
const int nSamples = 1000;//number of samples to create
void testCov()
{
  ///////////////////////////////////////
  //  Choose Regularization Type       //
  //  NO_REG for tau=0                 //  
  //  VAR_REG for TUnfold Picks Tau    //
  //  CONST_REG to choose your own tau //
  ///////////////////////////////////////
  
  //int regType = NO_REG;
  int regType = VAR_REG;
  //int regType = CONST_REG;
  
  TH1::SetDefaultSumw2();
  //gROOT->SetBatch(kTRUE);
  //Load the files
  TFile*file= new TFile(fileName);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  //Determine parameters used to create the distributions
  ifstream parameterFile("parameters.txt");
  bool exactClosure,effInc,backInc;
  parameterFile >> exactClosure >> effInc >> backInc;
  double vecSum[nLogBins],vecAvg[nLogBins],vecM[nSamples][nLogBins],vecSum2[nLogBins],
  vecStd[nLogBins],vector[nSamples][nLogBins2];

  //Define hisograms
  TH1D*hReco[nSamples];
  //hReco[0] is the original, unchanged reco mass from toyModel.C
  hReco[0] = (TH1D*)file->Get("hReco");//reconstructed mass
   hReco[0]->SetMarkerStyle(20);
   hReco[0]->SetMarkerColor(kBlack);
   hReco[0]->SetLineColor(kBlack);
  TH1D*hGen = (TH1D*)file->Get("hTrue");//true mass
  TH1D*hBack;
  if(backInc){
    if(exactClosure)
     hBack = (TH1D*)file->Get("hBack1");
    else
     hBack = (TH1D*)file->Get("hBack2");
  }
  TH2D*hMatrix = (TH2D*)file->Get("hMatrix");//migration matrix
  hGen->SetFillColor(kRed+2);
  hGen->SetLineColor(kRed+2);
 
  //Initialize array of histograms and then smear each bin using a gaussian
  //each hReco[i] has its bin counts varied by a gaussian with the width of the bin error
  for(int i=1;i<nSamples;i++){
   TString recoName = "hReco";
   recoName += i;
   hReco[i] = (TH1D*)file->Get("hReco");
   hReco[i]->SetName(recoName);
   hReco[i]->SetMarkerStyle(20);
   hReco[i]->SetMarkerColor(kBlack);
   hReco[i]->SetLineColor(kBlack);
   for(int k=1;k<nLogBins2+1;k++){
    double error = hReco[0]->GetBinError(k);
    TF1*fGaus = new TF1("fGaus","gaus(0)",-5*error,5*error);
    fGaus->SetParameters(1,0,error);
    double smear = fGaus->GetRandom();
    double bin = hReco[0]->GetBinContent(k)+smear;
    hReco[i]->SetBinContent(k,bin);
   }
  }
 
 //Loop to do unfolding on each sample
 for(int j=0;j<nSamples;j++){
  ////////////////////////////
  //  Regularization Modes  //
  ////////////////////////////
  //TUnfold::ERegMode regMode = TUnfold::kRegModeNone; //breaks 
  //TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
  //TUnfold::ERegMode regMode = TUnfold::kRegModeDerivative;
  TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
  //TUnfold::ERegMode regMode = TUnfold::kRegModeMixed; //breaks

  ///////////////////////////
  //  Types of Constraint  //
  ///////////////////////////
  //TUnfold::EConstraint constraintMode = TUnfold::kEConstraintNone;
  TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;
  
  /////////////////////
  //  Density Modes  //
  /////////////////////
  //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeNone;
  TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
  //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeUser;
  //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidthAndUser;
 
  /////////////////////////////////////
  //  Horizontal vs Vertical Output  //
  /////////////////////////////////////
  //TUnfold::EHistMap outputMap = TUnfold::kHistMapOutputVert;
  TUnfold::EHistMap outputMap = TUnfold::kHistMapOutputHoriz;

  //////////////////////////////////////
  //  Constructor for TUnfoldDensity  //
  //////////////////////////////////////
  TUnfoldDensity unfold(hMatrix,outputMap,regMode,constraintMode,densityFlags);
  unfold.SetInput(hReco[j]);//the measured distribution

  //////////////////////////////
  //  Background Subtraction  //
  //////////////////////////////
  double backScale = 1.0;
  double backScaleError = 0.0;//scale error for background
  if(backInc)
   unfold.SubtractBackground(hBack,"background",backScale,backScaleError);

  ////////////////////////////
  //  Add Systematic Error  //
  ////////////////////////////
  //For this part, you need a migration matrix calculated using a different MC generator
  //TUnfoldSys::ESysErrMode sysErrorType = TUnfoldSys::kSysErrModeMatrix;
  //TUnfoldSys::ESysErrMode sysErrorType = TUnfoldSys::kSysErrModeShift;
  //TUnfoldSys::ESysErrMode sysErrorType = TUnfoldSys::kSysErrModeRelative;
  //unfold.AddSysError(hAltMatrix,"systematic error",outputMap,sysErrorType);

  ////////////////////////////
  //  Begin Regularization  //
  ////////////////////////////
  Int_t iBest;
  TSpline *logTauX,*logTauY;
  TGraph *lCurve;
  TGraph*bestLcurve;
  TGraph*bestLogTauLogChi2;
  if(regType == VAR_REG){
    Int_t nScan=30;//This number chosen only because it was given in the tutorial
    Double_t tauMin = 0.0;//If tauMin=tauMax, TUnfold automatically chooses a range
    Double_t tauMax = 0.0;//Not certain how they choose the range
    iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
    cout<< "tau=" << unfold.GetTau() << endl;
    cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()<<" / "<<unfold.GetNdf()<<"\n";
    Double_t t[1],x[1],y[1];
    logTauX->GetKnot(iBest,t[0],x[0]);
    logTauY->GetKnot(iBest,t[0],y[0]);
    bestLcurve=new TGraph(1,x,y);
    bestLogTauLogChi2=new TGraph(1,t,x);
  }
  else if(regType == NO_REG){
    double tau = 0;
    unfold.DoUnfold(tau,hReco[j]);
  }
  else{//user defined
    double tau = 1e-2;//larger tau introduces more MC bias
    unfold.DoUnfold(tau,hReco[j]);
  }

  //The Unfolded Distribution
  TString unfName = "unfolded";
  unfName += j;
  TH1*hUnfolded = unfold.GetOutput(unfName);
   hUnfolded->SetMarkerStyle(25);
   hUnfolded->SetMarkerColor(kBlue+2);
   hUnfolded->SetMarkerSize(1);
  TString unfEName = "unfoldedE";
  unfEName += j;
  TH1F*hUnfoldedE = new TH1F(unfEName,"",nBins,massbins);
   hUnfoldedE->SetMarkerStyle(25);
   hUnfoldedE->SetMarkerColor(kBlue+2);
   hUnfoldedE->SetMarkerSize(1);
  TString eMatName = "errorMatrix";
  eMatName += j;
  TH2*histEmatTotal=unfold.GetEmatrixTotal(eMatName);
  for(int i=0;i<nBins;i++){
    double c = hUnfolded->GetBinContent(i+1);
    hUnfoldedE->SetBinContent(i+1,c);
    hUnfoldedE->SetBinError(i+1,TMath::Sqrt(histEmatTotal->GetBinContent(i+1,i+1)));
  }
  TString recoRebinName = "hRecoRebin";
  recoRebinName += j;
  TH1F*hRecoRebin=(TH1F*)hReco[j]->Clone(recoRebinName);
   hRecoRebin->Rebin(2);
  TString ratioName = "ratio";
  ratioName += j;
  TH1F*ratio = (TH1F*)hUnfoldedE->Clone(ratioName);
   ratio->Divide(hGen);

  double x[nBins],res[nBins];
  double chi = hUnfoldedE->Chi2Test(hGen,"CHI2/NDF",res);//chi2/ndf to print on plot
  double pValues = hUnfoldedE->Chi2Test(hGen,"P",res);//outputs chi2,prob,ndf,igood
  TLatex*chiLabel = new TLatex(500.0,150000,Form("#chi^{2}/ndf = %lg", chi));	

  //Drawing unfolded distributions. 
  //Comment this section out or set batch to true if looping over many samples
/*  
  const float padmargins = 0.03;
  TString canName = "canvas";
  canName += j;
  TCanvas*canvas1 = new TCanvas(canName,"",10,10,1200,1000);
  TPad*pad1 = new TPad("","",0,0.3,1.0,1.0);
  pad1->SetBottomMargin(padmargins);
  pad1->SetGrid();
  pad1->SetLogy();
  pad1->SetLogx();
  pad1->SetTicks(1,1);
  pad1->Draw();
  pad1->cd();
  TLine*line = new TLine(binLow,1,binHigh,1);
  line->SetLineColor(kRed);
  TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
  legend->SetTextSize(0.02);
  legend->AddEntry(hGen,"True Distribution");
  legend->AddEntry(hReco[j],"Measured Distribution");
  legend->AddEntry(hUnfolded,"Unfolded Distribution");
  hGen->SetLabelSize(0);
  hGen->SetTitleSize(0);
  hGen->Draw("hist");
  hRecoRebin->Draw("PE,same");
  hUnfoldedE->Draw("PE,same");
  legend->Draw("same");
  chiLabel->Draw("same");

  canvas1->cd();
  TPad*pad2 = new TPad("","",0,0.05,1,0.3);
  pad2->SetLogx();
  pad2->SetTopMargin(padmargins);
  pad2->SetBottomMargin(0.2);
  pad2->SetGrid();
  pad2->SetTicks(1,1);
  pad2->Draw();
  pad2->cd();  
  ratio->GetYaxis()->SetLabelSize(0.06);
  ratio->GetYaxis()->SetTitleSize(0.08);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->SetTitle("Unfolded/Gen");
  ratio->GetXaxis()->SetTitle("mass [GeV]");
  ratio->GetXaxis()->SetLabelSize(0.1);
  ratio->GetXaxis()->SetTitleSize(0.1);
  ratio->GetXaxis()->SetNoExponent();
  ratio->GetXaxis()->SetMoreLogLabels();
  ratio->GetYaxis()->SetRangeUser(0.8,1.2);
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerColor(kBlack);
  ratio->Draw("PE");
  line->Draw("same");
  TString saveName = "/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase1Plots/";
  saveName += "testCovUnfolded";
  saveName += j;
  saveName += ".png";
  canvas1->SaveAs(saveName);
*/
  //place each unfolded distribution into an array
  for(int k=1;k<nLogBins+1;k++){
    vector[j-1][k-1] = hUnfolded->GetBinContent(k);
  }

 }//end unfolding loop

 //Calculate covariance matrix
 //matrixB is made up of each vector with the average vector subtracted from it
 TMatrixD matrixB(nLogBins,nSamples);
 TMatrixD matrixBT(nSamples,nLogBins);//transpose of matrixB
 TMatrixD covM(nLogBins,nLogBins);//covariance matrix
 TMatrixD corrM(nLogBins,nLogBins);

 for(int i=0;i<nLogBins;i++){
  vecSum[i] = 0;
  vecSum2[i]=0;
 }
 for(int iEle=0;iEle<nLogBins;iEle++){
  for(int jVec=0;jVec<nSamples;jVec++){
   vecSum[iEle] += vector[jVec][iEle];//sum of element iEle summed over vectors
   vecSum2[iEle] += vector[jVec][iEle]*vector[jVec][iEle];//sum of squares of elements
  }
  vecAvg[iEle] = vecSum[iEle]/nSamples;
  vecStd[iEle] = sqrt(vecSum2[iEle]/nSamples-vecSum[iEle]*vecSum[iEle]/(nSamples*nSamples*1.0));
 }

 for(int iEle=0;iEle<nLogBins;iEle++){
  for(int jVec=0;jVec<nSamples;jVec++){
   vecM[jVec][iEle] = vector[jVec][iEle]-vecAvg[iEle];//vectors minus average vector
   matrixB(iEle,jVec) = vecM[jVec][iEle];
  }
 }

 matrixBT.Transpose(matrixB);//the transpose of matrixB calculated above
 covM = matrixB*matrixBT;//covariance matrix without proper weight
 TH2D*hCovM = new TH2D("hCovM","",nLogBins,massbins,nLogBins,massbins);
 TH2D*hCorrM = new TH2D("hCorrM","",nLogBins,massbins,nLogBins,massbins);
 double valCov,valCorr;

 for(int jEle=0;jEle<nLogBins;jEle++){
  for(int iEle=0;iEle<nLogBins;iEle++){
   covM(iEle,jEle) = covM(iEle,jEle)/nSamples;//covariance matrix
   corrM(iEle,jEle) = covM(iEle,jEle)/(vecStd[iEle]*vecStd[jEle]);//correlation matrix
   valCov = covM(iEle,jEle);//value in each bin to place in the histogram
   valCorr = corrM(iEle,jEle);//value in each bin to place in the histogram
   hCovM->SetBinContent(iEle+1,jEle+1,valCov);
   hCorrM->SetBinContent(iEle+1,jEle+1,valCorr);
   
  }
 }

 //Draw covariance and correlation matrices
 TCanvas*can = new TCanvas("can","",1000,1000);
 can->SetLogy();
 can->SetLogx();
 can->SetLogz();
 hCovM->SetTitle("Covariance Matrix");
 hCovM->Draw("colz");
 TString canName = "/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase1Plots/testCovM_";
 canName += nSamples;
 canName += "Samples.png";
 can->SaveAs(canName);

 TCanvas*can1 = new TCanvas("can1","",1000,1000);
 can1->SetLogy();
 can1->SetLogx();
 hCorrM->SetTitle("Correlation Matrix");
 hCorrM->Draw("colz");
 TString can1Name = "/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase1Plots/testCorrM_";
 can1Name += nSamples;
 can1Name += "Samples.png";
 can1->SaveAs(can1Name);

}




















