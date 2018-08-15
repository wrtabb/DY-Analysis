#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TRatioPlot.h"

//Histogram parameters
const int nHistoTypes = 11; 
const int nHistos = 5;
const float axisLow = 0.1;
const float axisHigh = 100000000;
const float padmargins = 0.03;
const TString hStackName = "hStack";
const TString dataFileName = "./plots/dataVsMC.root";
const TString histName[nHistos] = {
  "hFakes", "hEW", "hTops", "hMC", "hData"
};
const TString histTypeName[nHistoTypes] = {
  "InvMass", "InvMassLinear", "Vert", "VertWeighted", "pTLead", "pTSub", "pTDi", "EtaLead", "EtaSub", "EtaDi",
  "Rapidity"
};
const Color_t histFillColors[nHistos] = {
  kViolet+5, kRed+2, kBlue+2, kOrange-2, kWhite
};
const Color_t histLineColors[nHistos] = {
  kViolet+3, kRed+4, kBlue+3, kOrange+3, kBlack
};
const TString xAxisLabels[nHistoTypes] = {
  "mass [GeV]",
  "mass [GeV]",
  "Number of vertices",
  "Number of vertices",
  "p_{T} [GeV]",
  "p_{T} [GeV]",
  "p_{T} [GeV]",
  "#eta",
  "#eta",
  "#eta",
  "Y"
};
float binLowInvMass = 0;
const float binHighInvMass = 3000;
//InvMass Linear Plot
const float binLowInvMassLinear = 60;
const float binHighInvMassLinear = 120;
////nVertices
const float binLowVert = 0;
const float binHighVert = 50;
////pT 
const float binLowpT = 0;
const float binHighpT = 500;
////eta
const float binLowEta = -2.5;
const float binHighEta = 2.5;
////rapidity
const float binLowY = -2.5;
const float binHighY = 2.5;

const TString plotTitle[nHistoTypes] = {
  "Dielectron invariant mass 15-3000 GeV",
  "Dielectron invariant mass 60-120 GeV",
  "Number of vertices without PU weighting",
  "Number of vertices with PU weighting",
  "Leading electron p_{T}",
  "Sub-leading electron p_{T}",
  "Dielectron p_{T}",
  "Leading electron #eta",
  "Sub-leading electron #eta",
  "Dielectron #eta",
  "Dielectron Rapidity"
};
enum HistBins {
  FAKES,
  EW,
  TOPS,
  EE,
  DATA
};
enum HistTypes {
  INV_MASS,
  INV_MASS_LINEAR,
  VERTICES,
  VERTICES_WEIGHTED,
  PT_LEAD,
  PT_SUB,
  PT_DI,
  ETA_LEAD,
  ETA_SUB,
  ETA_DI,
  RAPIDITY
};

void drawDataVsMC()
{
  gStyle->SetOptStat(0);
  
  TFile*file = new TFile(dataFileName);
  if(!file){
    cout << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "File failed to load!!!" << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << endl;
  }

  ////////////////////////////////////
  //                                //
  //-----Filling histograms and-----//
  //-----Setting draw options-------//
  //                                //
  ////////////////////////////////////

  TH1F*histos[nHistoTypes][nHistos];  
  TH1F*hDataMCRatio[nHistoTypes];
  TH1F*hMCSum[nHistoTypes];
  TString hSumName;
  TString hRatioName;
  for(int i=0;i<nHistoTypes;i++){//type of histogram
    hSumName = "hMCSum";
    hSumName += i;
    hRatioName = "hDataMCRatio";
    hRatioName += i;
    for(int j=0;j<nHistos;j++){//histogram within type (backgrounds, signal, etc)  
      if(i==VERTICES_WEIGHTED&&j==DATA) histos[i][j]= (TH1F*)file->Get(histName[j]+histTypeName[i-1]);
      else histos[i][j]=(TH1F*)file->Get(histName[j]+histTypeName[i]); 

      //Setting negative bin content to ZERO!!!!!!!!!!!!!!!!!!!!!! 
      int maxBin = histos[i][j]->GetMaximumBin();
      double x = histos[i][j]->GetXaxis()->GetBinCenter(maxBin);
      for(int k=0;k<x;k++){     
	if(histos[i][j]->GetBinContent(k) < 0) histos[i][j]->SetBinContent(k,0);
      }
      
      if(j==DATA) {
	histos[i][j]->SetLineColor(kBlack);
	histos[i][j]->SetMarkerColor(kBlack);
	histos[i][j]->SetMarkerSize(1);
	histos[i][j]->SetMarkerStyle(20);
      }
      else {
	histos[i][j]->SetFillColor(histFillColors[j]);
	histos[i][j]->SetLineColor(histLineColors[j]);
      }      
      histos[i][j]->GetXaxis()->SetTitle(xAxisLabels[i]);
      if(j==0) hMCSum[i] = (TH1F*)histos[i][0]->Clone(hSumName);
      else if(j!=DATA&&j!=0)
	hMCSum[i]->Add(histos[i][j]);
    }
    hDataMCRatio[i] = (TH1F*)histos[i][DATA]->Clone(hRatioName);
    hDataMCRatio[i]->Divide(hMCSum[i]);
    hDataMCRatio[i]->SetLineColor(kBlack);
    hDataMCRatio[i]->SetMarkerColor(kBlack);
    hDataMCRatio[i]->SetMarkerStyle(20);
    hDataMCRatio[i]->GetYaxis()->SetRangeUser(0.5,1.5);    
  }

  //////////////////////////////////////////
  //                                      //
  //-----Place histograms into stacks-----//
  //                                      //
  //////////////////////////////////////////

  THStack*hStack[nHistoTypes];
  for(int i=0;i<nHistoTypes;i++){  
    hStack[i] = new THStack(hStackName+histTypeName[i],"");
    for(int j=0;j<nHistos;j++){
      if(j==DATA) continue; //no data in stacks
      hStack[i]->Add(histos[i][j]);
      hStack[i]->SetTitle(plotTitle[i]);
      hStack[i]->SetMinimum(axisLow);
      hStack[i]->SetMaximum(axisHigh);
    }
  }

  TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
  legend->SetTextSize(0.02);
  legend->AddEntry(histos[0][DATA],"Data");
  legend->AddEntry(histos[0][EE],"#gamma^{*}/Z #rightarrow e^{-}e^{+}");
  legend->AddEntry(histos[0][TOPS],"t#bar{t}+tW+#bar{t}W");
  legend->AddEntry(histos[0][EW],"EW (Dibosons, #gamma^{*}/Z #rightarrow #tau^{-}#tau^{+})");
  legend->AddEntry(histos[0][FAKES],"Fakes (W+Jets)");  

  TCanvas*canvas[nHistoTypes];
  TPad*pad1[nHistoTypes];
  TPad*pad2[nHistoTypes];
  TLine*line[nHistoTypes];
  TVirtualPad*p1[nHistoTypes];
  TVirtualPad*p2[nHistoTypes];
  float x1[nHistoTypes]={0,60,0,0,0,0,0,-2.5,-2.5,-2.5,-2.5};
  float x2[nHistoTypes]={3000,120,50,50,500,500,500,2.5,2.5,2.5,2.5};
  for(int i=0;i<nHistoTypes;i++){ 
    TString canvasName = "canvas";
    canvasName+=i;
    canvas[i] = new TCanvas(canvasName,"",10,10,1000,1000);    
    line[i]=new TLine(x1[i],1,x2[i],1);
    line[i]->SetLineColor(kRed);
    pad1[i] = new TPad("","", 0,0.3,1,1.0);
    pad1[i]->SetBottomMargin(padmargins); 
    pad1[i]->SetGrid();
    pad1[i]->SetLogy();
    pad1[i]->Draw(); 
    pad1[i]->cd(); 
    if(i==INV_MASS) pad1[i]->SetLogx();
    hStack[i]->Draw("hist");
    hStack[i]->GetXaxis()->SetTitle(xAxisLabels[i]);
    hStack[i]->GetXaxis()->SetNoExponent();
    hStack[i]->GetXaxis()->SetMoreLogLabels();
    hStack[i]->SetMinimum(axisLow);
    hStack[i]->GetXaxis()->SetLabelSize(0);
    hStack[i]->GetXaxis()->SetTitleSize(0);
    hStack[i]->GetYaxis()->SetTitle("Events");
    canvas[i]->Update();
    if(i==VERTICES_WEIGHTED) histos[i-1][DATA]->Draw("same,PE");
    else histos[i][DATA]->Draw("same,PE");
    legend->Draw("same");

    canvas[i]->cd();
    pad2[i] = new TPad("","",0,0.05,1,0.3);
    p2[i] = pad2[i];
    if(i==INV_MASS) p2[i]->SetLogx();
    pad2[i]->SetTopMargin(padmargins);
    pad2[i]->SetBottomMargin(0.2);
    pad2[i]->SetGrid();
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
    TString saveName = "./plots/dataVsMC/dataVsMC_";
    saveName+=histTypeName[i];
    saveName+=".png";
    canvas[i]->SaveAs(saveName);     
  }

}//end invMassDraw
