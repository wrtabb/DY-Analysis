//#include "/home/hep/wrtabb/DY-Analysis/headers/header1.h"
//#include "/home/hep/wrtabb/DY-Analysis/headers/drawOptions.h"
#include "VariableList.h"

const TString fileName = "unfoldingHists.root";

enum RegType {//Strength of regularization
  NO_REG,           //No regularization
  CONST_REG,        //User defined regularization
  VAR_REG_LCURVE,   //TUnfoldDensity determines best choice of regularization strength
  VAR_REG_SCANSURE, //TUnfoldDensity determines best choice of regularization strength
  VAR_REG_SCANTAU   //TUnfoldDensity determines best choice of regularization strength
};
enum Bins{
  DEFAULT,
  ONE_EXTRA,
  TWO_EXTRA
};
const int binLow = 0;
const int binHigh = 50;
//-----Forward declarations of functions-----//
TH1F*unfold(RegType regType,TH1D*hReco,TH1D*hClosure,TH1D*hTrue,TH2D*hMatrix,TString recoName,
            bool closure);

void diffBinningUnfold()
{
 TH1::SetDefaultSumw2();
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);
 //gROOT->SetBatch(true);

 TFile*file = new TFile(fileName);

 std::vector<TH1D*> hReco;
 std::vector<TH1D*> hClosure;
 std::vector<TH1D*> hTrue;
 std::vector<TH2D*> hMatrix;
 int nDistributions = 7;
 for(int i=0;i<nDistributions;i++){
  TString hRecoTitle = "hReco";
  hRecoTitle += i;
  TString hClosureTitle = "hRecoClosure";
  hClosureTitle += i;
  TString hTrueTitle = "hTrue";
  hTrueTitle += i;
  TString hMatrixTitle = "hMatrix";
  hMatrixTitle += i;

  hReco.push_back( (TH1D*)file->Get(hRecoTitle) );
  hClosure.push_back( (TH1D*)file->Get(hClosureTitle) );
  hTrue.push_back( (TH1D*)file->Get(hTrueTitle) );
  hMatrix.push_back( (TH2D*)file->Get(hMatrixTitle) );

  //unfold(NO_REG,hReco.at(i),hClosure.at(i),hTrue.at(i),hMatrix.at(i),hRecoTitle,false);
  //unfold(NO_REG,hReco.at(i),hClosure.at(i),hTrue.at(i),hMatrix.at(i),hRecoTitle,true);
  unfold(VAR_REG_LCURVE,hReco.at(i),hClosure.at(i),hTrue.at(i),hMatrix.at(i),hRecoTitle,false);
  //unfold(VAR_REG_LCURVE,hReco.at(i),hClosure.at(i),hTrue.at(i),hMatrix.at(i),hRecoTitle,true);
 }

}


TH1F*unfold(RegType regType,TH1D*hReco,TH1D*hClosure,TH1D*hTrue,TH2D*hMatrix,TString recoName,             bool closure)
{
 //Do the unfolding
  TH1F*hBlank;
  if(closure){
   hReco = (TH1D*)hClosure->Clone();
   recoName += "Closure";
  }
  if(regType==NO_REG) recoName += "NoReg";
  hReco->SetMarkerStyle(20);
  hReco->SetMarkerColor(kBlack);
  hReco->SetLineColor(kBlack);

  hTrue->SetFillColor(kRed+2);
  hTrue->SetLineColor(kRed+2);
  hTrue->SetTitle("");
  int nBins = nBinsTrue;
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
  if(regType == VAR_REG_LCURVE){
    Int_t nScan=5;//This number chosen only because it was given in the tutorial
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
  else if(regType == VAR_REG_SCANSURE){
   
  }
  else if(regType == VAR_REG_SCANTAU){

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
  TH1F*hUnfoldedE = (TH1F*)hUnfolded->Clone("Unfolded with errors");

   hUnfoldedE->SetMarkerStyle(25);
   hUnfoldedE->SetMarkerColor(kBlue+2);
   hUnfoldedE->SetMarkerSize(1);
  for(int i=0;i<=nBins;i++){
    double c = hUnfolded->GetBinContent(i+1);
    hUnfoldedE->SetBinContent(i+1,c);
    hUnfoldedE->SetBinError(i+1,TMath::Sqrt(histEmatTotal->GetBinContent(i+1,i+1)));
  }
  TH1D*hRecoRebin = (TH1D*)hReco->Rebin(nBinsTrue,"hRecoRebin",binningTrue);
  TH1F*ratio = (TH1F*)hUnfoldedE->Clone("ratio");
   ratio->Divide(hTrue);
  double xChiLabel = 35;
  double yChiLabel = 5e5;
  double x[nBins],res[nBins];
  double chi = hUnfoldedE->Chi2Test(hTrue,"CHI2/NDF",res);//chi2/ndf to print on plot
  double pValues = hUnfoldedE->Chi2Test(hTrue,"P",res);//outputs chi2,prob,ndf,igood
  TLatex*chiLabel = new TLatex(xChiLabel,yChiLabel,Form("#chi^{2}/ndf = %lg", chi));	


  const float padmargins = 0.03;
  TCanvas*canvas1 = new TCanvas("canvas1","",10,10,1200,1000);
  TPad*pad1 = new TPad("","",0,0.3,1.0,1.0);
  pad1->SetBottomMargin(padmargins);
  pad1->SetGrid();
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
  pad2->SetTopMargin(padmargins);
  pad2->SetBottomMargin(0.2);
  pad2->SetGrid();
  pad2->SetTicks(1,1);
  pad2->Draw();
  pad2->cd();  
  ratio->SetMinimum(0.7);
  ratio->SetMaximum(1.3);
  ratio->GetYaxis()->SetLabelSize(0.06);
  ratio->GetYaxis()->SetTitleSize(0.08);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->SetTitle("Unfolded/Truth");
  ratio->GetXaxis()->SetLabelSize(0.1);
  ratio->GetXaxis()->SetTitleSize(0.1);
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerColor(kBlack);
  ratio->Draw("PE");
  line->Draw("same");

  TString saveName = "./plots/unfolded";
  saveName += recoName;
  saveName += ".png";
  canvas1->SaveAs(saveName);
  delete canvas1; 
/*
  TFile*fileXsec = new TFile("unfoldedSaved.root","recreate");
  canvas1->Write();
  hUnfoldedE->SetName("hUnfolded");
  hUnfoldedE->Write();
  histEmatTotal->SetName("hCovMatrix");
  histEmatTotal->Write();
  hTrue->Write();
  hReco->Write();
  fileXsec->Write();
  fileXsec->Close();
*/
  delete hReco;
  delete hTrue;
  delete hMatrix;
  return hUnfoldedE;
}
 
