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

  TH1F*histos[nHistoTypes][nHistos];  
  for(int i=0;i<nHistoTypes;i++){//type of histogram
    for(int j=0;j<nHistos;j++){//histogram within type
      if(i==VERTICES_WEIGHTED&&j==DATA) histos[i][j]= (TH1F*)file->Get(histName[j]+histTypeName[i-1]);
      else histos[i][j]=(TH1F*)file->Get(histName[j]+histTypeName[i]); 
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
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //Setting negative bin content found to ZERO!!!!!!!!!!!!!!!!!!!!!!
  for(int i=0;i<3000;i++){
      for(int j=0;j<nHistoTypes;j++){
	if(histos[INV_MASS][j]->GetBinContent(i) < 0) histos[INV_MASS][j]->SetBinContent(i,0);
      }
    }
  //////////////////////////////////////////////////////////////////////////////////////////////////
	
  //Place histograms into stacks
  THStack*hStack[nHistoTypes];
  TRatioPlot*hDataMCRatio[nHistoTypes];
  for(int i=0;i<nHistoTypes;i++){  
    hStack[i] = new THStack(hStackName+histTypeName[i],"");
    for(int j=0;j<nHistos;j++){
      if(j==DATA) continue; //no data in stacks
      hStack[i]->Add(histos[i][j]);
      hStack[i]->SetTitle(plotTitle[i]);
      hStack[i]->SetMinimum(axisLow);
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
  for(int i=0;i<nHistoTypes;i++){
    
    TString canvasName = "canvas";
    canvasName+=i;
    canvas[i] = new TCanvas(canvasName,"",10,10,1000,1000);

    if(i==INV_MASS){
      canvas[i]->SetLogx();
      canvas[i]->SetLogy();
    }
    else canvas[i]->SetLogy();
    canvas[i]->cd();
    
    /*
    hStack[i]->Draw("bar");
    hStack[i]->GetXaxis()->SetTitle(xAxisLabels[i]);
    hStack[i]->GetXaxis()->SetNoExponent();
    hStack[i]->GetXaxis()->SetMoreLogLabels();
    hStack[i]->SetMinimum(axisLow);
    canvas[i]->Update();
    if(i==VERTICES_WEIGHTED) histos[i-1][DATA]->Draw("same,PE");
    else histos[i][DATA]->Draw("same,PE");
    legend->Draw("same");
    */
    
    hDataMCRatio[i]=new TRatioPlot(hStack[i],histos[i][DATA]);
    hDataMCRatio[i]->Draw();
    hDataMCRatio[i]->GetUpperPad()->cd();
    if(i==INV_MASS||i==INV_MASS_LINEAR) legend->Draw("same");
    hDataMCRatio[i]->GetLowerRefGraph()->SetMinimum(0.5);
    hDataMCRatio[i]->GetLowerRefGraph()->SetMaximum(1.5);
    if(i==INV_MASS||i==INV_MASS_LINEAR){
      hDataMCRatio[i]->GetLowerRefXaxis()->SetNoExponent();
      hDataMCRatio[i]->GetLowerRefXaxis()->SetMoreLogLabels();
    }
    hDataMCRatio[i]->GetUpperRefXaxis()->SetTitle(xAxisLabels[i]);
    hDataMCRatio[i]->GetLowerRefYaxis()->SetTitle("MC/Data Ratio");
    hDataMCRatio[i]->GetLowerRefYaxis()->SetLabelSize(0.02);
    hDataMCRatio[i]->GetUpperRefYaxis()->SetTitle("Entries");
    hDataMCRatio[i]->GetUpperPad()->SetGrid();
    hDataMCRatio[i]->GetUpperPad()->SetGrid();
    canvas[i]->Update();  
    TString saveName = "./plots/dataVsMC/dataVsMC_";
    saveName+=histTypeName[i];
    saveName+=".png";
    canvas[i]->SaveAs(saveName);
    
  }

}//end invMassDraw
