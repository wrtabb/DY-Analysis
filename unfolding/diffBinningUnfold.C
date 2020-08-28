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
std::vector<TString> histNames = {
 "double bins all",
 "double bins 115-3000",
 "double bins 15-115",
 "double bins 15-60,120-3000",
 "double bins 60-120",
 "split highest bin",
 "split lowest bin"
};
std::vector<TString> plotSave = {
 "unfoldedDoubleBinsAll.png",
 "unfoldedDoubleBins_115_3000.png",
 "unfoldedDoubleBins_15_115.png",
 "unfoldedDoubleBins_15_60_120_3000.png",
 "unfoldedDoubleBins_60_120.png",
 "unfoldedSplitHighestBin.png",
 "unfoldedSplitLowestBin.png"
};

//-----Forward declarations of functions-----//
TH1F*unfold(VarType var,RegType regType,bool closure,TH1D*hReco,TH1D*hBack,TH2D*hMatrix,
            TString saveName);
TH1D*GetBackgrounds(VarType var);
TH2D*Rebin2D(TH2D*hist,TString histName,std::vector<double> binning);
TH2D*RebinTH2(TH2D*hist,TString histName,TH2D*hBinning);
TH1D*Rebin1D(TH1D*hist,TString histName,std::vector<double> binning);
TH1D*Rebin1D(TH1D*hist,TString histName,TH1D*hBinning);

void diffBinningUnfold()
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

 std::vector<TH1D*> hReco;
 hReco.push_back(Rebin1D(hDataM,"1D "+histNames.at(0),massbinsReco));
 hReco.push_back(Rebin1D(hDataM,"1D "+histNames.at(1),massbinsReco0));
 hReco.push_back(Rebin1D(hDataM,"1D "+histNames.at(2),massbinsReco1));
 hReco.push_back(Rebin1D(hDataM,"1D "+histNames.at(3),massbinsReco2));
 hReco.push_back(Rebin1D(hDataM,"1D "+histNames.at(4),massbinsReco3));
 hReco.push_back(Rebin1D(hDataM,"1D "+histNames.at(5),massbinsReco4));
 hReco.push_back(Rebin1D(hDataM,"1D "+histNames.at(6),massbinsReco5));

 std::vector<TH1D*> hBack;
 hBack.push_back(Rebin1D(hBackM,"Back "+histNames.at(0),massbinsReco));
 hBack.push_back(Rebin1D(hBackM,"Back "+histNames.at(1),massbinsReco0));
 hBack.push_back(Rebin1D(hBackM,"Back "+histNames.at(2),massbinsReco1));
 hBack.push_back(Rebin1D(hBackM,"Back "+histNames.at(3),massbinsReco2));
 hBack.push_back(Rebin1D(hBackM,"Back "+histNames.at(4),massbinsReco3));
 hBack.push_back(Rebin1D(hBackM,"Back "+histNames.at(5),massbinsReco4));
 hBack.push_back(Rebin1D(hBackM,"Back "+histNames.at(6),massbinsReco5));
  
 std::vector<TH2D*> hMatrix;
 hMatrix.push_back(Rebin2D(hMatrixM,"2D "+histNames.at(0),massbinsReco));
 hMatrix.push_back(Rebin2D(hMatrixM,"2D "+histNames.at(1),massbinsReco0));
 hMatrix.push_back(Rebin2D(hMatrixM,"2D "+histNames.at(2),massbinsReco1));
 hMatrix.push_back(Rebin2D(hMatrixM,"2D "+histNames.at(3),massbinsReco2));
 hMatrix.push_back(Rebin2D(hMatrixM,"2D "+histNames.at(4),massbinsReco3));
 hMatrix.push_back(Rebin2D(hMatrixM,"2D "+histNames.at(5),massbinsReco4));
 hMatrix.push_back(Rebin2D(hMatrixM,"2D "+histNames.at(6),massbinsReco5));

 //unfold(MASS,VAR_REG_LCURVE,false,hDataM,hBackM,hMatrixM,"testRegular.png");
 //unfold(MASS,VAR_REG_LCURVE,false,hReco.at(0),hBack.at(0),hMatrix.at(0),"testRebinned.png");
 cout << "********************************************" << endl;
 for(int i=0;i<hReco.size();i++)
 {
  cout << "hDataM entries: " << hDataM->Integral() << endl;
  cout << "hReco.at("<< i << ") entries: " << hReco.at(i)->Integral() << endl;
  cout << "********************************************" << endl;
  unfold(MASS,NO_REG,false,hReco.at(i),hBack.at(i),hMatrix.at(i),plotSave.at(i));
 }
  //unfold(MASS,VAR_REG_LCURVE,false,hDataM,hBackM,hMatrixM,plotSave.at(0));
 //unfold(MASS,VAR_REG_LCURVE,true,hDataM,hBackM,hMatrixM);
// unfold(MASS,VAR_REG_LCURVE,false,hDataM,hBackM,hMatrixM);
 //unfold(RAPIDITY,VAR_REG_LCURVE,true,hDataY,hBackY,hMatrixY);
 //unfold(RAPIDITY,VAR_REG_LCURVE,false,hDataY,hBackY,hMatrixY);
}

TH1F*unfold(VarType var,RegType regType,bool closure,TH1D*hReco,TH1D*hBack,TH2D*hMatrix,
            TString saveName)
{
  TH1F*hBlank;
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
   //this is an alternative way to find tau for regularization
   //Not yet implemented 
  }
  else if(regType == VAR_REG_SCANTAU){
   //this is an alternative way to find tau for regularization
   //Not yet implemented 
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

  //TH1D*hRecoRebin = (TH1D*)hReco->Rebin(nBins,"hRecoRebin",massbinsTrue);
  TH1D*hRecoRebin = Rebin1D(hReco,"hRecoRebin",hTrue);
  hRecoRebin->SetMarkerStyle(20);
  hRecoRebin->SetMarkerColor(kBlack);
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

  TString saveNameFinal = "./plots/binningTest/";
  saveNameFinal += saveName;
  canvas1->SaveAs(saveNameFinal);

  delete hReco;
  delete hTrue;
  delete hMatrix;
  delete hUnfolded;
  delete canvas1;
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

TH2D*Rebin2D(TH2D*hist,TString histName,std::vector<double> binning)
{
 //This function takes a 2D histogram as input and gives a new 2D histogram
 //As output with new Y-axis binning defined by input std::vector<double> binning
 
 int nBinsX = hist->GetNbinsX();//number of true bins in input histogram
 int nBinsY = hist->GetNbinsY();//number of reco bins in the input histogram
 int nBins = nBinsX*nBinsY;//total number of bins in input histogram 
 int nBinsReco = binning.size()-1;//number of reco bins for the output rebinned histogram
 double newbinning[nBinsReco];//binning placed in array for defining TH2

 for(int i=0;i<=nBinsReco;i++){
  newbinning[i] = binning.at(i);
 }

 //Initialize the output histogram
 TH2D*hRebin = new TH2D(histName,"",nLogBinsMassTrue,massbinsTrue,nBinsReco,newbinning);

 //Loop over each bin in the histogram
 //For each bin, define the position x and y from the bin center
 //Then do histogram fill for each of the entries in that bin
 double y,x;
 double nEntries;
 int binIn,binOut;
 for(int i=1;i<=nBinsX;i++){
  for(int j=1;j<=nBinsY;j++){
   x = hist->GetXaxis()->GetBinCenter(i);
   y = hist->GetYaxis()->GetBinCenter(j);
   binIn = hist->GetBin(i,j);
   nEntries = hist->GetBinContent(binIn);
   hRebin->Fill(x,y,nEntries);
  }//end y bin loop
 }//end x bin loop
 return hRebin;
}


TH1D*Rebin1D(TH1D*hist,TString histName,std::vector<double> binning)
{
 int nBinsHist = hist->GetNbinsX();
 int nBinsReco = binning.size()-1;
 double newbinning[nBinsReco];
 for(int i=0;i<=nBinsReco;i++){
  newbinning[i] = binning.at(i);
 }
 TH1D*hRebin = new TH1D(histName,"",nBinsReco,newbinning);

 double y,x;
 double nEntries;
 for(int j=1;j<=nBinsHist;j++){
  x = hist->GetXaxis()->GetBinCenter(j);  
  nEntries = hist->GetBinContent(j);
  hRebin->Fill(x,nEntries);
 } //end x bin loop
 return hRebin;
}

TH1D*Rebin1D(TH1D*hist,TString histName,TH1D*hBinning)
{
 int nBinsHist = hist->GetNbinsX();
 int nBinsReco = hBinning->GetNbinsX();
 double newbinning[nBinsReco];
 for(int i=0;i<=nBinsReco;i++){
  if(i==0) newbinning[i] = hBinning->GetBinLowEdge(i+1);
  else newbinning[i] = newbinning[i-1]+hBinning->GetBinWidth(i);
 }
 TH1D*hRebin = new TH1D(histName,"",nBinsReco,newbinning);

 double y,x;
 double nEntries;
 for(int j=1;j<=nBinsHist;j++){
  x = hist->GetXaxis()->GetBinCenter(j);  
  nEntries = hist->GetBinContent(j);
  hRebin->Fill(x,nEntries);
 } //end x bin loop
 return hRebin;
}
