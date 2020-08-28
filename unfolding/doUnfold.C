//#include "/home/hep/wrtabb/DY-Analysis/headers/header1.h"
//#include "/home/hep/wrtabb/DY-Analysis/headers/drawOptions.h"
#include "VariableList.h"

const TString fileName = "data/migrationMatrix.root";
const TString fileDataName = "data/inputData.root";
const std::vector<TString> backFile = {
 "data/backgroundFAKES.root",
 "data/backgroundEW.root",
 "data/backgroundTT.root"
};
std::vector<TString> variableNames = {
 "Rapidity",
 "Mass"
};

enum RegType {//Strength of regularization
  NO_REG,           //No regularization
  CONST_REG,        //User defined regularization
  VAR_REG_LCURVE,   //TUnfoldDensity determines best choice of regularization strength
  VAR_REG_SCANSURE, //TUnfoldDensity determines best choice of regularization strength
  VAR_REG_SCANTAU   //TUnfoldDensity determines best choice of regularization strength
};
enum VarType{
  RAPIDITY,
  MASS
};
enum Bins{
  DEFAULT,
  ONE_EXTRA,
  TWO_EXTRA
};
const int binLow = 15;
const int binHigh = 3000;
//-----Forward declarations of functions-----//
TH1F*unfold(VarType var,RegType regType,bool closure,TH1D*hReco,TH1D*hBack,TH2D*hMatrix,
            TString saveName
           );
TH1D*GetBackgrounds(VarType var);
TH1D*Rebin1D(TH1D*hist,TString histName,TH1D*hBinning);

void doUnfold()
{
 TH1::SetDefaultSumw2();
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);
 gROOT->SetBatch(true);

 //Load the files
 TFile*file = new TFile(fileName);
 TFile*fileData = new TFile(fileDataName);

 //Load input histograms
 TH1D*hDataM = (TH1D*)fileData->Get("hRecoMass");
 TH1D*hBackM = GetBackgrounds(MASS);
 TH2D*hMatrixM = (TH2D*)file->Get("hMatrixMass");
 
 TH1D*hDataY = (TH1D*)fileData->Get("hRecoRapidity");
 TH1D*hBackY = GetBackgrounds(RAPIDITY);
 TH2D*hMatrixY = (TH2D*)file->Get("hMatrixRapidity");

 //unfold(MASS,NO_REG,false,hDataM,hBackM,hMatrixM,"dataNoRegUnfoldedMass");
 //unfold(MASS,VAR_REG_LCURVE,false,hDataM,hBackM,hMatrixM,"dataLCurveUnfoldedMass");
 unfold(RAPIDITY,NO_REG,false,hDataY,hBackY,hMatrixY,"dataNoRegUnfoldedRapidity");
 //unfold(RAPIDITY,VAR_REG_LCURVE,false,hDataY,hBackY,hMatrixY,"dataLCurveUnfoldedRapidity");
}

TH1F*unfold(VarType var,RegType regType,bool closure,TH1D*hReco,TH1D*hBack,TH2D*hMatrix,
            TString saveName
           )
{
  TH1F*hBlank;
  if(closure) hReco = hMatrix->ProjectionY();
  hReco->SetMarkerStyle(20);
  hReco->SetMarkerColor(kBlack);
  hReco->SetLineColor(kBlack);
  TH1D*hTrue = hMatrix->ProjectionX();
  hTrue->SetFillColor(kRed+2);
  hTrue->SetLineColor(kRed+2);
  hTrue->SetTitle("");
  int nBins;
  nBins = hTrue->GetXaxis()->GetNbins();
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
  if(!closure) unfold.SubtractBackground(hBack,"background",backScale,backScaleError);

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
  TH1F*hUnfoldedE = (TH1F*)hTrue->Clone("Unfolded with errors");
   hUnfoldedE->SetMarkerStyle(25);
   hUnfoldedE->SetMarkerColor(kBlue+2);
   hUnfoldedE->SetMarkerSize(1);

  for(int i=0;i<nBins;i++){
    double c = hUnfolded->GetBinContent(i+1);
    hUnfoldedE->SetBinContent(i+1,c);
    hUnfoldedE->SetBinError(i+1,TMath::Sqrt(histEmatTotal->GetBinContent(i+1,i+1)));
  }

  TH1D*hRecoRebin = Rebin1D(hReco,"hRecoRebin",hTrue);
  hRecoRebin->SetMarkerStyle(20);
  hRecoRebin->SetMarkerColor(kBlack);
  hRecoRebin->SetLineColor(kBlack);
  TH1F*ratio = (TH1F*)hUnfoldedE->Clone("ratio");
   ratio->Divide(hTrue);
  double xChiLabel;
  if(var==MASS) xChiLabel = 500;
  else if(var==RAPIDITY) xChiLabel = 0.5;
  double x[nBins],res[nBins];
  double chi = hUnfoldedE->Chi2Test(hTrue,"CHI2/NDF",res);//chi2/ndf to print on plot
  double pValues = hUnfoldedE->Chi2Test(hTrue,"P",res);//outputs chi2,prob,ndf,igood
  TLatex*chiLabel = new TLatex(xChiLabel,1e6,Form("#chi^{2}/ndf = %lg", chi));	


  const float padmargins = 0.03;
  const float yAxisMinimum = 0.1;
  const float yAxisMaximum= 1e9;
  TCanvas*canvas1 = new TCanvas("canvas1","",10,10,1200,1000);
  TPad*pad1 = new TPad("","",0,0.3,1.0,1.0);
  pad1->SetBottomMargin(padmargins);
  pad1->SetGrid();
  pad1->SetLogy();
  if(var==MASS) pad1->SetLogx();
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
  hTrue->SetMinimum(yAxisMinimum);
  hTrue->SetMaximum(yAxisMaximum);
  hTrue->Draw("hist");
  hRecoRebin->Draw("PE,same");
  hUnfoldedE->Draw("PE,same");
  legend->Draw("same");
  chiLabel->Draw("same");

  canvas1->cd();
  TPad*pad2 = new TPad("","",0,0.05,1,0.3);
  if(var==MASS) pad2->SetLogx();
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
  if(var==MASS) ratio->GetXaxis()->SetTitle("mass [GeV]");
  else if(var==RAPIDITY) ratio->GetXaxis()->SetTitle("rapidity");
  ratio->GetXaxis()->SetLabelSize(0.1);
  ratio->GetXaxis()->SetTitleSize(0.1);
  ratio->GetXaxis()->SetNoExponent();
  ratio->GetXaxis()->SetMoreLogLabels();
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerColor(kBlack);
  ratio->Draw("PE");
  line->Draw("same");

  TString savePlot = "plots/";
  savePlot += saveName;
  savePlot += ".png";
  canvas1->SaveAs(savePlot);
  
  double width,nUnfold;
  TH1D*hCross = (TH1D*)hTrue->Clone("hCross");
   hCross->SetMarkerStyle(20);
   hCross->SetMarkerColor(kBlack);
   hCross->SetTitle("Cross Section");
   hCross->GetXaxis()->SetMoreLogLabels();
   hCross->GetXaxis()->SetNoExponent();
   hCross->GetXaxis()->SetTitle("mass [GeV]");
   hCross->GetYaxis()->SetTitle("d#sigma/dm [pb/GeV]");
  
  
  for(int k=1;k<nBins+1;k++){
   width = hUnfoldedE->GetXaxis()->GetBinWidth(k);
   nUnfold = hUnfoldedE->GetBinContent(k)/(width*dataLuminosity);
   hCross->SetBinContent(k,nUnfold);
  }

  TString fileSaveName = "data/unfolded";
  if(var==MASS) fileSaveName += "Mass";
  else if(var==RAPIDITY) fileSaveName += "Rapidity";
  else {
   cout << "Variable must be MASS or RAPIDITY!" << endl;
   return hBlank;
  }
  if(closure) fileSaveName += "Closure";
  fileSaveName += ".root";

  TFile*fileXsec = new TFile(fileSaveName,"recreate");
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

  delete hReco;
  delete hTrue;
  delete hMatrix;
  return hUnfoldedE;
}
 
TH1D*GetBackgrounds(VarType var)
{
 int nFiles = backFile.size();
 TFile*file[nFiles];
 TH1D*hBack[nFiles];
 TH1D*hBackSum;
 TString histName = "hReco";
 histName += variableNames.at(var);
 for(int i=0;i<nFiles;i++){
  file[i] = new TFile(backFile.at(i));
  hBack[i] = (TH1D*)file[i]->Get(histName);
  if(i==0) hBackSum = (TH1D*)hBack[0]->Clone("hBackSum");
  else hBackSum->Add(hBack[i]);
 }
 return hBackSum;
}

TH1D*Rebin1D(TH1D*hist,TString histName,TH1D*hBinning)
{
 int nBinsOld = hist->GetNbinsX();
 int nBinsNew = hBinning->GetNbinsX();
 TH1D*hBlank;
 if(nBinsNew > nBinsOld){
  cout << "*********************************************************" << endl;
  cout << "ERROR: new binning must have fewer bins than old binning!" << endl;
  cout << "*********************************************************" << endl;
  return hBlank;
 }
 double newbinning[nBinsNew];
 for(int i=0;i<=nBinsNew;i++){
  if(i==0) newbinning[i] = hBinning->GetBinLowEdge(i+1);
  else newbinning[i] = newbinning[i-1]+hBinning->GetBinWidth(i);
 }
 TH1D*hRebin = new TH1D(histName,"",nBinsNew,newbinning);
 double y,x;
 double nEntries;
 double histErrors[nBinsOld+2];
 for(int j=1;j<=nBinsOld;j++){
  histErrors[j] = hist->GetBinError(j);
  x = hist->GetXaxis()->GetBinCenter(j);
  nEntries = hist->GetBinContent(j);
  hRebin->Fill(x,nEntries);
  hRebin->GetBin(x);
 } //end x bin loop
 
 //This part combines the errors from the bins of the old histogram that go into each
 //bin of the new histogram
 double histBinWidth;
 double newBinWidth;
 double nBinsOldInNew[nBinsNew];
 double newBinUpperEdge,oldBinUpperEdge;
 int k = 1;
 //Loop over bins in new histogram
 for(int i=1;i<=nBinsNew;i++){
  //find the upper edge of the current bin
  newBinUpperEdge = hRebin->GetBinLowEdge(i)+hRebin->GetBinWidth(i);
  double newError2 = 0;
  //There must be more bins in the old histogram than the new
  //Loop over bins in the old histogram stopping when the upper edge of the new histogram
  //reaches the upper edge of the old histogram
  //For each bin from the old histogram that goes into the new one,
  //the error is found and the squares of the errors of each bin are added
  for(int j=k;j<=nBinsOld;j++){
   oldBinUpperEdge = hist->GetBinLowEdge(j)+hist->GetBinWidth(j);
   if(oldBinUpperEdge > newBinUpperEdge) continue;
   k = j;
   newError2 += hist->GetBinError(j)*hist->GetBinError(j);
  } 	
  hRebin->SetBinError(i,sqrt(newError2));
 }
 return hRebin;
}

