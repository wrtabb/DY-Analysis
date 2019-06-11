#include "TUnfold.h"
#include "TUnfoldSys.h"
#include "TUnfoldDensity.h"
#include "TUnfoldBinning.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
using namespace std;

enum Reglarization {//Strength of regularization
  NO_REG,           //No regularization
  CONST_REG,        //User defined regularization
  VAR_REG           //TUnfoldDensity determines best choice of regularization strength
};
const int dataLumi = 35867;//data luminosity (pb)
const int nBins = 43;
const int binLow = 15;
const int binHigh = 3000;
const double massbins[] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,
 110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,
 1500,3000};
const double massbins2[] = {15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,
 57.5,60,62,64,66,68,70,72,74,76,78.5,81,83.5,86,88.5,91,93.5,96,98.5,101,103.5,106,108,110,
 112.5,115,117.5,120,123,126,129.5,133,137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,
 210,220,231.5,243,258,273,296.5,320,350,380,410,440,475,510,555,600,650,700,765,830,915,1000,
 1250,1500,2250,3000};
const int nLogBins2 = 86;
//File Names
const TString fileName = "outputDataUnfold.root";       

void doUnfold()
{
  ///////////////////////////////////////
  //  Choose Regularization Type       //
  //  NO_REG for tau=0                 //  
  //  VAR_REG for TUnfold Picks Tau    //
  //  CONST_REG to choose your own tau //
  ///////////////////////////////////////
  
  //gROOT->SetBatch(kTRUE);
  int regType = NO_REG;
  //int regType = VAR_REG;
  //int regType = CONST_REG;
  
  TH1::SetDefaultSumw2();
  //Load the files
  TFile*file= new TFile(fileName);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  //Define hisograms
  TH1D*hReco = (TH1D*)file->Get("hData");//reconstructed mass
  TH1D*hTrue = (TH1D*)file->Get("hTrue");//true mass
  TH1D*hBack = (TH1D*)file->Get("hBack");//background mass
  TH2D*hMatrix = (TH2D*)file->Get("hMatrix");//migration matrix
  TH2D*hCovM = (TH2D*)file->Get("hCovM");
  TH2D*hCovMinv = (TH2D*)file->Get("hCovMinv");
  hReco->SetMarkerStyle(20);
  hReco->SetMarkerColor(kBlack);
  hReco->SetLineColor(kBlack);
  hTrue->SetFillColor(kRed+2);
  hTrue->SetLineColor(kRed+2);
  hTrue->SetTitle("");

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
  //unfold.SetInput(hReco,0.0,0.0,hCovM,hCovMinv);//the measured distribution
  unfold.SetInput(hReco);//the measured distribution

  //////////////////////////////
  //  Background Subtraction  //
  //////////////////////////////
  double backScale = 1.0;
  double backScaleError = 0.0;//scale error for background
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
    unfold.DoUnfold(tau,hReco);
  }
  else{//user defined
    double tau = 1e-8;//larger tau introduces more MC bias
    unfold.DoUnfold(tau,hReco);
  }

  //The Unfolded Distribution
  TH1*hUnfolded = unfold.GetOutput("Unfolded");
   hUnfolded->SetMarkerStyle(25);
   hUnfolded->SetMarkerColor(kBlue+2);
   hUnfolded->SetMarkerSize(1);
  TH2*histEmatStat=unfold.GetEmatrixInput("unfolding stat error matrix");
  TH2*histEmatTotal=unfold.GetEmatrixTotal("unfolding total error matrix");
  TH1F*hUnfoldedE = new TH1F("Unfolded with errors",";(gen)",nBins,massbins);
   hUnfoldedE->SetMarkerStyle(25);
   hUnfoldedE->SetMarkerColor(kBlue+2);
   hUnfoldedE->SetMarkerSize(1);
  for(int i=0;i<nBins;i++){
    double c = hUnfolded->GetBinContent(i+1);
    hUnfoldedE->SetBinContent(i+1,c);
    hUnfoldedE->SetBinError(i+1,TMath::Sqrt(histEmatTotal->GetBinContent(i+1,i+1)));
  }
  TH1F*hRecoRebin=(TH1F*)hReco->Clone("hRecoRebin");
   hRecoRebin->Rebin(2);
  TH1F*ratio = (TH1F*)hUnfoldedE->Clone("ratio");
   ratio->Divide(hTrue);

  double x[nBins],res[nBins];
  double chi = hUnfoldedE->Chi2Test(hTrue,"CHI2/NDF",res);//chi2/ndf to print on plot
  double pValues = hUnfoldedE->Chi2Test(hTrue,"P",res);//outputs chi2,prob,ndf,igood
  TLatex*chiLabel = new TLatex(500.0,150000,Form("#chi^{2}/ndf = %lg", chi));	

  const float padmargins = 0.03;
  TCanvas*canvas1 = new TCanvas("canvas1","",10,10,1200,1000);
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
  legend->AddEntry(hTrue,"True Distribution");
  legend->AddEntry(hReco,"Measured Distribution");
  legend->AddEntry(hUnfolded,"Unfolded Distribution");
  hTrue->SetLabelSize(0);
  hTrue->SetTitleSize(0);
  hTrue->Draw("hist");
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
  ratio->GetYaxis()->SetTitle("Unfolded/Truth");
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

  canvas1->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase3/dataUnfolded_only2Ele.png");
 
  double width,nUnfold;
  TH1D*hCross = new TH1D("hCross","",nBins,massbins);
   hCross->SetMarkerStyle(20);
   hCross->SetMarkerColor(kBlack);
   hCross->SetTitle("Cross Section");
   hCross->GetXaxis()->SetMoreLogLabels();
   hCross->GetXaxis()->SetNoExponent();
   hCross->GetXaxis()->SetTitle("mass [GeV]");
   hCross->GetYaxis()->SetTitle("d#sigma/dm [pb/GeV]");
  for(int k=1;k<nBins+1;k++){
   width = hUnfoldedE->GetXaxis()->GetBinWidth(k);
   nUnfold = hUnfoldedE->GetBinContent(k)/(width*dataLumi);
   hCross->SetBinContent(k,nUnfold);
  }

  TCanvas*canvas2 = new TCanvas("canvas2","",10,10,1400,1000);
  canvas2->SetLogy();
  canvas2->SetLogx();
  canvas2->SetLogz();
  canvas2->SetGrid();
  //histEmatStat->Draw("colz");
  hCross->Draw("PE");
  canvas2->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase3/xsec_only2Ele.png");
  TFile*fileXsec = new TFile("xsec.root","recreate");
  hCross->Write();
  canvas1->Write();
  hUnfoldedE->Write();
  histEmatTotal->Write();
  hTrue->Write();
  hReco->Write();
  fileXsec->Write();
  fileXsec->Close();
}
