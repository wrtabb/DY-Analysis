#include "VariableList.h"

enum Variable{
 INV_MASS,
 RAPIDITY,
 PT
};
void DoPlots(Variable var,TH2D*hMatrix,TString histName,TString plotTitle);
void GetConditionNumber(TH2D*hResponse);
void GetPercentOnDiagonal(TH2D*hist);
TH2D*MakeResponse(TH2D*hMatrix,TString histName);

const TString fileNameMatrix = "data/migrationMatrix.root";

void makeMigrationPlots()
{
 TH1::SetDefaultSumw2();
 gROOT->SetBatch(true);
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);

 TFile*fileMatrix = new TFile(fileNameMatrix);
 TH2D*hMatrixMass = (TH2D*)fileMatrix->Get("hMatrixMass");
 hMatrixMass->RebinY(2);
 TH2D*hMatrixRapidity = (TH2D*)fileMatrix->Get("hMatrixRapidity");
 hMatrixRapidity->RebinY(2);

 TH2D*hResponseMass = MakeResponse(hMatrixMass,"hResponseMass");
 TH2D*hResponseRapidity = MakeResponse(hMatrixRapidity,"hResponseRapidity");

 DoPlots(INV_MASS,hResponseMass,"responseMatrixMass","Mass Response Matrix");
 DoPlots(RAPIDITY,hResponseRapidity,"responseMatrixRapidity","Rapidity Response Matrix");
 DoPlots(INV_MASS,hMatrixMass,"migrationMatrixMass","Mass Migration Matrix");
 DoPlots(RAPIDITY,hMatrixRapidity,"migrationMatrixRapidity","Rapidity Migration Matrix");

 GetConditionNumber(hResponseMass);
 GetConditionNumber(hResponseRapidity);

 GetPercentOnDiagonal(hResponseMass);
 GetPercentOnDiagonal(hResponseRapidity);
}

void GetPercentOnDiagonal(TH2D*hist)
{
 int nBinsX = hist->GetNbinsX();
 int nBinsY = hist->GetNbinsY();
 cout << endl;
 cout << "Percent on diagonal in each column for " << hist->GetName() << endl;
 for(int i=1;i<nBinsX+1;i++){
  cout << i << ": " << 100*hist->GetBinContent(i,i) << "%" << endl; 
 }
}

void GetConditionNumber(TH2D*hResponse)
{
 TString histName = hResponse->GetName();
 //Turn hist into matrix
 int nBinsX = hResponse->GetNbinsX();
 int nBinsY = hResponse->GetNbinsY();
 TMatrixD matrix(nBinsY,nBinsX);
 for(int i=0;i<nBinsX;i++){
  for(int j=0;j<nBinsY;j++){
   matrix(j,i) = hResponse->GetBinContent(i+1,j+1);
  }
 }
 TDecompSVD decomp(matrix);
 double condition = decomp.Condition();
 cout << "The condition number for " << histName << ": " << condition << endl;

 double determinant;
 TMatrixD mInverse = matrix.Invert(&determinant);
 cout << "The determinant of " << histName << " is " << determinant << endl;
}

TH2D*MakeResponse(TH2D*hMatrix,TString histName)
{
 int nBinsX = hMatrix->GetNbinsX();
 int nBinsY = hMatrix->GetNbinsY();
 for(int i=1;i<nBinsX+1;i++){
  for(int j=1;j<nBinsY+1;j++){
   if(hMatrix->GetBinContent(i,j)<0) hMatrix->SetBinContent(i,j,0.0);
  }
 }
 TH2D*hResponse = (TH2D*)hMatrix->Clone(histName);
 //TH2D*hResponse = new TH2D(histName,"",nBinsX,massbinsTrue,nBinsY,massbinsTrue);
 double sumY;
 for(int i=1;i<nBinsX+1;i++){
  sumY = 0;
  for(int j=1;j<nBinsY+1;j++){
   sumY += hMatrix->GetBinContent(i,j);
  }
  for(int j=1;j<nBinsY+1;j++){
   hResponse->SetBinContent(i,j,hMatrix->GetBinContent(i,j)/sumY);
   if(hMatrix->GetBinContent(i,j) < 0) cout << "(" << i << ", " << j << "): " << hMatrix->GetBinContent(i,j) << endl;
  }
 }
 return hResponse;
}

void DoPlots(Variable var,TH2D*hMatrix,TString histName,TString plotTitle)
{
 TCanvas*canvas = new TCanvas("canvas","",0,0,1200,1000);
 canvas->SetGrid();
 if(var==INV_MASS){
  canvas->SetLogx();
  canvas->SetLogy();
  hMatrix->GetYaxis()->SetNoExponent();
  hMatrix->GetYaxis()->SetMoreLogLabels();
  hMatrix->GetXaxis()->SetNoExponent();
  hMatrix->GetXaxis()->SetMoreLogLabels();
  hMatrix->GetXaxis()->SetTitle("true mass [GeV]");
  hMatrix->GetYaxis()->SetTitle("reco mass [GeV]");
 }
 else if(var==RAPIDITY){
  hMatrix->GetXaxis()->SetTitle("true Y");
  hMatrix->GetYaxis()->SetTitle("reco Y");
 }
 canvas->SetLogz();
 canvas->SetRightMargin(0.13);
 canvas->SetLeftMargin(0.13);
 hMatrix->SetTitle(plotTitle);
 hMatrix->Draw("colz");

 TString saveNameMatrix = "plots/";
 saveNameMatrix += histName;
 saveNameMatrix += "NoNegativeBins";
 saveNameMatrix += ".png";

 canvas->SaveAs(saveNameMatrix);
 //delete canvas;
}
