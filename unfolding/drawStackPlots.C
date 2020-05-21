#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TRatioPlot.h"
#include "VariableList.h"


const std::vector<TString> fileName = {
 "data/backgroundFAKES.root",
 "data/backgroundEW.root",
 "data/backgroundTT.root",
 "data/migrationMatrix.root",
 "data/inputData.root"
};

std::vector<TString> histName = {
 "hFakes",
 "hEW",
 "hTops",
 "hLL",
 "hData"
};
enum VarType{
 MASS,
 RAPIDITY
};
const std::vector<TString> variableName = {
 "Mass",
 "Rapidity"
};
const std::vector<TString> xAxisTitles = {
 "dilepton invariant mass [GeV]",
 "dilepton rapidity"
};
const float yAxisLow = 0.1;
const float yAxisHigh = 1e8;
const int nVariables = 2;
const int nSamples = 5;
const Color_t histFillColors[nSamples] = {
  kViolet+5,//fakes 
  kRed+2,   //EW
  kBlue+2,  //Tops
  kOrange-2,//LL 
  kWhite    //Data
};
const Color_t histLineColors[nSamples] = {
  kViolet+5,//fakes 
  kRed+2,   //EW
  kBlue+2,  //Tops
  kOrange-2,//LL 
  kBlack//Data
};

void drawStackPlots()
{
  gStyle->SetOptStat(0);
  //gROOT->SetBatch(true);
  TFile*file[nSamples];
  for(int i=0;i<nSamples;i++){
   file[i] = new TFile(fileName.at(i));
   cout << fileName.at(i) << endl;
   file[i]->ls();
  }

  //blank histogram to hold properties
  TH1D*hBlank[nVariables];
 
  TH1D*hHist[nVariables][nSamples];
  int count = 0;
  for(int i=0;i<nVariables;i++){
   TString loadName = "hReco";
   if(i==MASS) loadName += "Mass";
   else if(i==RAPIDITY) loadName += "Rapidity";
   cout << loadName << endl;
   for(int j=0;j<nSamples;j++){
    if(i==MASS) histName.at(j) += "Mass";
    else if(i==RAPIDITY) histName.at(j) += "Rapidity";
    hHist[i][j] = (TH1D*)file[j]->Get(loadName);
    hHist[i][j]->SetName(histName.at(j));    
    hHist[i][j]->SetFillColor(histFillColors[j]);
    hHist[i][j]->SetLineColor(histLineColors[j]);
    hHist[i][j]->GetXaxis()->SetTitle(xAxisTitles.at(i));
    hHist[i][j]->GetXaxis()->SetNoExponent();
    hHist[i][j]->GetXaxis()->SetMoreLogLabels();
    //Setting negative bin content to zero
    int maxBin = hHist[i][j]->GetMaximumBin();
    double x = hHist[i][j]->GetXaxis()->GetBinCenter(maxBin);
    for(int k=1;k<x+1;k++){
     if(hHist[i][j]->GetBinContent(k) < 0) hHist[i][j]->SetBinContent(k,0);
    }
    hHist[i][j]->Rebin(2);
   }
  }
 

 THStack*hStack[nVariables];
 for(int i=0;i<nVariables;i++){
  TString stackName = "hStack";
  stackName += variableName.at(i);
  hStack[i] = new THStack(stackName,"");
  for(int j=0;j<nSamples-1;j++){
   hStack[i]->Add(hHist[i][j]);  
  }
 }
 TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
 legend->SetTextSize(0.02);
 legend->AddEntry(hHist[0][DATA],"Data");
 legend->AddEntry(hHist[0][LL],"#gamma^{*}/Z #rightarrow e^{-}e^{+}");
 legend->AddEntry(hHist[0][TT],"t#bar{t}+tW+#bar{t}W");
 legend->AddEntry(hHist[0][EW],"EW (Dibosons, #gamma^{*}/Z #rightarrow #tau^{-}#tau^{+})");
 legend->AddEntry(hHist[0][FAKES],"Fakes (W+Jets)");

 TCanvas*canvas[nVariables];
 for(int i=0;i<nVariables;i++){
  TString cName = "canvas";
  cName += variableName.at(i);
  canvas[i]=new TCanvas(cName,"",0,0,1200,1000);
  canvas[i]->SetGrid();
  canvas[i]->SetLogy();
  if(i==MASS){
   canvas[i]->SetLogx();
  }
  hStack[i]->SetMinimum(yAxisLow);
  hStack[i]->SetMaximum(yAxisHigh);
  hStack[i]->Draw("hist");
  hHist[i][DATA]->SetMarkerStyle(20);
  hHist[i][DATA]->Draw("pe,same");
  legend->Draw("same");
 }
}//end invMassDraw

