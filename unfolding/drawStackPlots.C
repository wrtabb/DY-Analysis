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
 RAPIDITY,
 DI_PT
};
const std::vector<TString> variableName = {
 "Mass",
 "Rapidity",
 "DiPT"
};
const std::vector<TString> xAxisTitles = {
 "dilepton invariant mass [GeV]",
 "dilepton rapidity",
 "dilepton p_{T}"
};
const float yAxisLow = 0.1;
const float yAxisHigh = 1e8;
const int nVariables = 3;
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
  gROOT->SetBatch(true);
  TFile*file[nSamples];
  for(int i=0;i<nSamples;i++){
   file[i] = new TFile(fileName.at(i));
   cout << fileName.at(i) << endl;
   file[i]->ls();
  }

  TH1D*hMCSum[nVariables]; 
  TH1D*histos[nVariables][nSamples];
  TH1D*hDataMCRatio[nVariables];
  const float ratioLow = 0.7;
  const float ratioHigh = 1.3;
  for(int i=0;i<nVariables;i++){
   TString loadName = "hReco";
   TString hRatioName = "hDataMCRatio";
   TString hSumName = "hMCSum";
   loadName += variableName.at(i);
   hRatioName += variableName.at(i);
   hSumName += variableName.at(i);

   for(int j=0;j<nSamples;j++){
    histName.at(j) += variableName.at(i);
    histos[i][j] = (TH1D*)file[j]->Get(loadName);
    histos[i][j]->SetName(histName.at(j));    
    histos[i][j]->SetFillColor(histFillColors[j]);
    histos[i][j]->SetLineColor(histLineColors[j]);
    histos[i][j]->GetXaxis()->SetNoExponent();
    histos[i][j]->GetXaxis()->SetMoreLogLabels();

    //Setting negative bin content to zero
    int maxBin = histos[i][j]->GetMaximumBin();
    double x = histos[i][j]->GetXaxis()->GetBinCenter(maxBin);
    for(int k=1;k<x+1;k++){
     if(histos[i][j]->GetBinContent(k) < 0) histos[i][j]->SetBinContent(k,0);
    }

    histos[i][j]->Rebin(2);//all reco distributions have twice the bins for unfolding purposes
    if(j==0) hMCSum[i] = (TH1D*)histos[i][0]->Clone(hSumName);
    else if(j!=DATA && j!=0) hMCSum[i]->Add(histos[i][j]);
    }
    hMCSum[i]->SetMarkerStyle(22);
    histos[i][DATA]->SetMarkerStyle(20);
    hDataMCRatio[i] = (TH1D*)histos[i][DATA]->Clone(hRatioName);
    hDataMCRatio[i]->Divide(hMCSum[i]);
    hDataMCRatio[i]->SetLineColor(kBlack);
    hDataMCRatio[i]->SetMarkerColor(kBlack);
    hDataMCRatio[i]->SetMarkerStyle(20);
    hDataMCRatio[i]->GetYaxis()->SetRangeUser(ratioLow,ratioHigh);
   }
 

 THStack*hStack[nVariables];
 for(int i=0;i<nVariables;i++){
  TString stackName = "hStack";
  stackName += variableName.at(i);
  hStack[i] = new THStack(stackName,"");
  for(int j=0;j<nSamples-1;j++){
   hStack[i]->Add(histos[i][j]);  
  }
 }
 TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
 legend->SetTextSize(0.02);
 legend->AddEntry(histos[0][DATA],"Data");
 legend->AddEntry(histos[0][LL],"#gamma^{*}/Z #rightarrow e^{-}e^{+}");
 legend->AddEntry(histos[0][TT],"t#bar{t}+tW+#bar{t}W");
 legend->AddEntry(histos[0][EW],"EW (Dibosons, #gamma^{*}/Z #rightarrow #tau^{-}#tau^{+})");
 legend->AddEntry(histos[0][FAKES],"Fakes (W+Jets)");

  TCanvas*canvas[nVariables];
  TPad*pad1[nVariables];
  TPad*pad2[nVariables];
  TLine*line[nVariables];
  TVirtualPad*p1[nVariables];
  TVirtualPad*p2[nVariables];
  const float padmargins = 0.03;
  float x1[nVariables]={0,-2.4};
  float x2[nVariables]={3000,2.5};
  for(int i=0;i<nVariables;i++){
    TString canvasName = "canvas_";
    canvasName+=i;
    canvas[i] = new TCanvas(canvasName,"",10,10,1000,1000);
    line[i]=new TLine(x1[i],1,x2[i],1);
    line[i]->SetLineColor(kRed);
    pad1[i] = new TPad("","", 0,0.3,1,1.0);
    pad1[i]->SetBottomMargin(padmargins);
    pad1[i]->SetGrid();
    pad1[i]->SetLogy();
    pad1[i]->SetTicks(1,1);
    pad1[i]->Draw();
    pad1[i]->cd();
    if(i==MASS) pad1[i]->SetLogx();
    hStack[i]->Draw("hist");
    hStack[i]->GetXaxis()->SetTitle(xAxisTitles.at(i));
    hStack[i]->GetXaxis()->SetNoExponent();
    hStack[i]->GetXaxis()->SetMoreLogLabels();
    hStack[i]->SetMinimum(yAxisLow);
    hStack[i]->SetMaximum(yAxisHigh);
    hStack[i]->GetXaxis()->SetLabelSize(0);
    hStack[i]->GetXaxis()->SetTitleSize(0);
    canvas[i]->Update();
    histos[i][DATA]->Draw("pe,same");
    legend->Draw("same");

    canvas[i]->cd();
    pad2[i] = new TPad("","",0,0.05,1,0.3);
    p2[i] = pad2[i];
    if(i==MASS) pad2[i]->SetLogx();
    pad2[i]->SetTopMargin(padmargins);
    pad2[i]->SetBottomMargin(0.2);
    pad2[i]->SetGrid();
    pad2[i]->SetTicks(1,1);
    pad2[i]->Draw();
    pad2[i]->cd();
    hDataMCRatio[i]->GetYaxis()->SetLabelSize(0.06);
    hDataMCRatio[i]->GetYaxis()->SetTitleSize(0.08);
    hDataMCRatio[i]->GetYaxis()->SetTitleOffset(0.3);
    hDataMCRatio[i]->GetYaxis()->SetTitle("Data/MC");
    hDataMCRatio[i]->GetXaxis()->SetLabelSize(0.1);
    hDataMCRatio[i]->GetXaxis()->SetTitleSize(0.1);
    hDataMCRatio[i]->GetXaxis()->SetNoExponent();
    hDataMCRatio[i]->GetXaxis()->SetMoreLogLabels();
    hDataMCRatio[i]->Draw("PE");
    line[i]->Draw("same");
    canvas[i]->Update();
    TString saveName = "./plots/";
    saveName += "stacked_";
    saveName += variableName.at(i);
    saveName += ".png";
    canvas[i]->SaveAs(saveName);
  }

}//end invMassDraw

