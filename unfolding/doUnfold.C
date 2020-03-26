#include "/home/hep/wrtabb/DY-Analysis/headers/header1.h"
#include "/home/hep/wrtabb/DY-Analysis/headers/drawOptions.h"

const TString fileName = "data/inputMatrix.root";
const TString fileDataName = "data/inputData.root";
const TString fileBackgroundName = "data/inputBackground.root";


//-----Globals-----//
enum RegType {//Strength of regularization
  NO_REG,           //No regularization
  CONST_REG,        //User defined regularization
  VAR_REG           //TUnfoldDensity determines best choice of regularization strength
};
const int binLow = 15;
const int binHigh = 3000;
//-----Forward declarations of functions-----//
TH1F*unfold(RegType regType,bool closure,TH1D*hReco,TH1D*hBack,TH2D*hMatrix);

void doUnfold()
{
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);
 gROOT->SetBatch(kTRUE);
 TH1::SetDefaultSumw2(); 
 //Load the files
 TFile*file= new TFile(fileName);
 TFile*fileData = new TFile(fileDataName);
 TFile*fileBack = new TFile(fileBackgroundName);
 TH1D*hData = (TH1D*)fileData->Get("hReco");
 TH1D*hBack = (TH1D*)fileBack->Get("hBack");
 TH2D*hMatrix = (TH2D*)file->Get("hMatrix");
 
 TH1F*hUnfolded;
 hUnfolded = unfold(VAR_REG,false,hData,hBack,hMatrix);
}

TH1F*unfold(RegType regType,bool closure,TH1D*hReco,TH1D*hBack,TH2D*hMatrix)
{
  if(closure) hReco = hMatrix->ProjectionY();
  hReco->SetMarkerStyle(20);
  hReco->SetMarkerColor(kBlack);
  hReco->SetLineColor(kBlack);
  TH1D*hTrue = hMatrix->ProjectionX();
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
  unfold.SetInput(hReco);//the measured distribution

  //////////////////////////////
  //  Background Subtraction  //
  //////////////////////////////
  double backScale = 1.0;
  double backScaleError = 0.0;//scale error for background
  if(!closure)unfold.SubtractBackground(hBack,"background",backScale,backScaleError);

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
  TH1F*hUnfoldedE = new TH1F("Unfolded with errors",";(gen)",nLogBins,massbins);
   hUnfoldedE->SetMarkerStyle(25);
   hUnfoldedE->SetMarkerColor(kBlue+2);
   hUnfoldedE->SetMarkerSize(1);
  for(int i=0;i<nLogBins;i++){
    double c = hUnfolded->GetBinContent(i+1);
    hUnfoldedE->SetBinContent(i+1,c);
    hUnfoldedE->SetBinError(i+1,TMath::Sqrt(histEmatTotal->GetBinContent(i+1,i+1)));
  }
  TH1F*hRecoRebin=(TH1F*)hReco->Clone("hRecoRebin");
   hRecoRebin->Rebin(2);
  TH1F*ratio = (TH1F*)hUnfoldedE->Clone("ratio");
   ratio->Divide(hTrue);

  double x[nLogBins],res[nLogBins];
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
  //ratio->GetYaxis()->SetRang
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerColor(kBlack);
  ratio->Draw("PE");
  line->Draw("same");

  if(!closure) canvas1->SaveAs("plots/dataUnfolded.png");
  else canvas1->SaveAs("plots/closureTest.png");

  double width,nUnfold;
  TH1D*hCross = new TH1D("hCross","",nLogBins,massbins);
   hCross->SetMarkerStyle(20);
   hCross->SetMarkerColor(kBlack);
   hCross->SetTitle("Cross Section");
   hCross->GetXaxis()->SetMoreLogLabels();
   hCross->GetXaxis()->SetNoExponent();
   hCross->GetXaxis()->SetTitle("mass [GeV]");
   hCross->GetYaxis()->SetTitle("d#sigma/dm [pb/GeV]");
  for(int k=1;k<nLogBins+1;k++){
   width = hUnfoldedE->GetXaxis()->GetBinWidth(k);
   nUnfold = hUnfoldedE->GetBinContent(k)/(width*dataLuminosity);
   hCross->SetBinContent(k,nUnfold);
  }

  TCanvas*canvas2 = new TCanvas("canvas2","",10,10,1400,1000);
  canvas2->SetLogy();
  canvas2->SetLogx();
  canvas2->SetLogz();
  canvas2->SetGrid();
  hCross->Draw("PE");
  if(!closure)canvas2->SaveAs("plots/xsec.png");
  TFile*fileXsec = new TFile("data/unfolded.root","recreate");
  hCross->Write();
  canvas1->Write();
  hUnfoldedE->SetName("hUnfolded");
  hUnfoldedE->Write();
  histEmatTotal->SetName("hCovMatrix");
  histEmatTotal->Write();
  hTrue->Write();
  hReco->Write();
  fileXsec->Write();
  fileXsec->Close();

  return hUnfoldedE;
}
