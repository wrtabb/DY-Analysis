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

const int nBins = 20;
const double massbins[44] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,
                             110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,
                             380,440,510,600,700,830,1000,1500,3000};
//File Names
//const TString fileName= "sampleHistos.root";         //Toy model
const TString fileName2= "closureFinerDataBinning.root"; //MC InvMass projections
const TString fileName= "dataUnfoldTest.root"; //Data

void testTUnfold()
{
  TH1::SetDefaultSumw2();
  //Load the files
  TFile*file= new TFile(fileName);
  TFile*file2= new TFile(fileName2);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  //Define hisograms
  //TH1F*hReco = (TH1F*)file->Get("hMCMeasured");//Toy Model
  //TH1F*hGen = (TH1F*)file->Get("hMCTrue");
  //TH2F*hMatrix = (TH2F*)file->Get("hMigrationMatrix");
  //TH2F*hMatrixSys = (TH2F*)file->Get("hMigrationMatrix2");
  //TH1F*hReco = (TH1F*)file->Get("migMatrixGENFSvsReco_py");//MC InvMass
  //TH1F*hGen = (TH1F*)file->Get("migMatrixGENFSvsReco_px");
  //TH2F*hMatrix = (TH2F*)file->Get("migMatrixGENFSvsReco");
  //TH2F*hMatrixSys = (TH2F*)file->Get("migMatrixGENFSvsReco");
  TH1F*hBackground = (TH1F*)file->Get("hMCBackground");//Data
  TH1F*hReco = (TH1F*)file->Get("hData");
  TH1F*hGen = (TH1F*)file2->Get("hTrue");
  TH2F*hMatrix = (TH2F*)file->Get("hMatrix");
  hReco->SetMarkerStyle(20);
  hReco->SetMarkerColor(kBlack);
  hReco->SetLineColor(kBlack);
  //hGen->SetFillColor(kRed+2);
  //hGen->SetLineColor(kRed+2);
  ////////////////////////////
  //  Regularization Modes  //
  ////////////////////////////
  //TUnfold::ERegMode regMode = TUnfold::kRegModeNone; //!!!!!Breaks if option selected
  //TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
  //TUnfold::ERegMode regMode = TUnfold::kRegModeDerivative;
  TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
  //TUnfold::ERegMode regMode = TUnfold::kRegModeMixed;  //!!!!!Breaks if option selected

  ///////////////////////////
  //  Types of Constraint  //
  ///////////////////////////
  //TUnfold::EConstraint constraintMode = TUnfold::kEConstraintNone;
  TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea; //!!!!!Breaks if option selected
  
  /////////////////////
  //  Density Modes  //
  /////////////////////
  //!!!!!All of these options seem to work
  //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeNone;
  TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
  //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeUser;
  //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidthAndUser;
 
  /////////////////////////////////////
  //  Horizontal vs Vertical Output  //
  /////////////////////////////////////
  //!!!!!Either of these works depending on the matrix being used
  TUnfold::EHistMap outputMap = TUnfold::kHistMapOutputVert;
  //TUnfold::EHistMap outputMap = TUnfold::kHistMapOutputHoriz;

  //////////////////////////////////////
  //  Constructor for TUnfoldDensity  //
  //////////////////////////////////////
  TUnfoldDensity unfold(hMatrix,outputMap,regMode,constraintMode,densityFlags);
  unfold.SetInput(hReco);//the measured distribution

  //////////////////////////////
  //  Background Subtraction  //
  //////////////////////////////
  unfold.SubtractBackground(hBackground,"background",1.0,0);

  ////////////////////////////
  //  Add Systematic Error  //
  ////////////////////////////
  //For this part, you need a migration matrix calculated using a different MC generator
  //unfold.AddSysError(hMatrixSys,"signalshape_SYS",
  //                   outputMap,TUnfoldSys::kSysErrModeMatrix);

  ////////////////////////////
  //  Begin Regularization  //
  ////////////////////////////
  Int_t nScan=30;//This number chosen only because it was given in the tutorial
  Double_t tauMin = 0.0;//If tauMin=tauMax, TUnfold automatically chooses a range
  Double_t tauMax = 0.0;//Not certain how they choose the range
  Int_t iBest;
  TSpline *logTauX,*logTauY;
  TGraph *lCurve;
  iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
  cout<< "tau=" << unfold.GetTau() << endl;
  cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()<<" / "<<unfold.GetNdf()<<"\n";
  Double_t t[1],x[1],y[1];
  logTauX->GetKnot(iBest,t[0],x[0]);
  logTauY->GetKnot(iBest,t[0],y[0]);
  TGraph *bestLcurve=new TGraph(1,x,y);
  TGraph *bestLogTauLogChi2=new TGraph(1,t,x);

  //The Unfolded Distribution
  TH1*hUnfolded = unfold.GetOutput("Unfolded");
  hUnfolded->SetMarkerStyle(25);
  hUnfolded->SetMarkerColor(kBlue+2);
  hUnfolded->SetMarkerSize(2);
  TH2 *histEmatStat=unfold.GetEmatrixInput("unfolding stat error matrix");
  TH2 *histEmatTotal=unfold.GetEmatrixTotal("unfolding total error matrix");
  //TH1 *histGlobalCorr=unfold.GetRhoItotal("histGlobalCorr",0,0,0,kFALSE);
  //TH2 *histCorrCoeff=unfold.GetRhoIJtotal("histCorrCoeff",0,0,0,kFALSE);
  TCanvas*canvas1 = new TCanvas("canvas1","",10,10,1200,1200);
  canvas1->SetLogy();
  canvas1->SetLogx();
  TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
  legend->SetTextSize(0.02);
  legend->AddEntry(hGen,"Gen Distribution");
  legend->AddEntry(hReco,"Measured Distribution");
  legend->AddEntry(hUnfolded,"Unfolded Distribution");
  hGen->GetXaxis()->SetTitle("Invariant mass [GeV]");
  hGen->Draw("hist");
  hReco->Draw("PE,same");
  hUnfolded->Draw("PE,same");
  legend->Draw("same");

  //TCanvas*canvas2 = new TCanvas("canvas2","",10,10,1200,1200);
  //canvas2->SetLogy();
  //canvas2->SetLogx();
  //hMatrix->GetXaxis()->SetTitle("Measured");
  //hMatrix->GetYaxis()->SetTitle("True");
  //hMatrix->Draw("colz");

  TCanvas*canvas3 = new TCanvas("canvas3","",10,10,1200,1200);
  canvas3->Divide(2);
  canvas3->cd(1);
  lCurve->Draw("AL");
  bestLcurve->SetMarkerColor(kRed);
  bestLcurve->SetMarkerSize(2);
  bestLcurve->Draw("*");
  canvas3->cd(2);
  logTauX->Draw();
  bestLogTauLogChi2->SetMarkerColor(kRed);
  bestLogTauLogChi2->SetMarkerSize(2);
  bestLogTauLogChi2->Draw("*");
  //canvas1->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/testUnfoldData.png");
  //canvas3->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/testUnfoldDataCurves.png");
  
}
