#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TRatioPlot.h"


const double massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 
			     106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 
			     380, 440, 510, 600, 700, 830, 1000, 1500, 3000};
const int nMassBins = 43;
const int nCanvas = 6;
const int nHistos = 10;
const float axisLow = 1.0;

const TString dataFileName = "./plots/dataVsMC.root";
const TString histName[nHistos] = {"hFakes", "hFakeslinear", "hEW", "hEWlinear", "hTops", "hTopslinear",   
				   "hMCInvMass", "hMCInvMasslinear", "hDataInvMass", "hDataInvMasslinear"};
const Color_t histFillColors[8] = {kViolet+5, kViolet+5, kRed+2, kRed+2, kBlue+2, kBlue+2, kOrange-2, kOrange-2};
const Color_t histLineColors[8] = {kViolet+3, kViolet+3, kRed+4, kRed+4, kBlue+3, kBlue+3, kOrange-3, kOrange-3};

enum HistBins {
  BINS_FAKES,
  BINS_FAKES_LINEAR,
  BINS_EW,
  BINS_EW_LINEAR,
  BINS_TOPS,
  BINS_TOPS_LINEAR,
  BINS_EE,
  BINS_EE_LINEAR,
  BINS_DATA,
  BINS_DATA_LINEAR    
};

void drawDataVsMC()
{
  gStyle->SetOptStat(0);
  TFile*file = new TFile(dataFileName);
  //Getting histograms
  TH1F*histos[nHistos];
  for(int i=0; i<nHistos;i++) {
    histos[i] = (TH1F*)file->Get(histName[i]);
    histos[i]->SetMinimum(axisLow);
    histos[i]->SetTitle("MC vs. Data");
    if(i==BINS_DATA||i==BINS_DATA_LINEAR) {
      histos[i]->SetLineColor(kBlack);
      histos[i]->SetMarkerColor(kBlack);
      histos[i]->SetMarkerSize(1);
      histos[i]->SetMarkerStyle(20);
    }
    else {
      histos[i]->SetFillColor(histFillColors[i]);
      histos[i]->SetLineColor(histLineColors[i]);
    }
  }//end histogram loop
  
//Place histograms into stacks
  THStack*hStack = new THStack("hStack","");
  THStack*hStacklinear = new THStack("hStacklinear","");
  hStack->SetMinimum(axisLow);
  hStacklinear->SetMinimum(axisLow);
  for(int i=0;i<nHistos;i++) {
    if(i==BINS_DATA||i==BINS_DATA_LINEAR) continue;
    if(i%2!=0) {
      hStacklinear->Add(histos[i]);
    }
    if(i%2==0) {
      hStack->Add(histos[i]);
    }
  }  

  TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
  legend->SetTextSize(0.02);
  legend->AddEntry(histos[BINS_DATA],"Data");
  legend->AddEntry(histos[BINS_EE],"#gamma^{*}/Z #rightarrow e^{-}e^{+}");
  legend->AddEntry(histos[BINS_TOPS],"t#bar{t}+tW+#bar{t}W");
  legend->AddEntry(histos[BINS_EW],"EW (Dibosons, #gamma^{*}/Z #rightarrow #tau^{-}#tau^{+})");
  legend->AddEntry(histos[BINS_FAKES],"Fakes (W+Jets)");

  TCanvas*canvas1 = new TCanvas("canvas3","",10,10,1000,1000);
  canvas2->SetLogx();
  canvas2->SetLogy();

  TCanvas*canvas2 = new TCanvas("canvas4","",10,10,1000,1000);
  canvas2->SetLogy();

  auto hDataMCRatio = new TRatioPlot(hStack,histos[BINS_DATA]);
  hDataMCRatio->GetXaxis()->SetTitle("m_{ee} [GeV]");  
  canvas1->cd();
  hDataMCRatio->Draw();
  hDataMCRatio->GetUpperPad()->cd();
  legend->Draw("same");
  hDataMCRatio->GetLowerRefGraph()->SetMinimum(0.7);
  hDataMCRatio->GetLowerRefGraph()->SetMaximum(1.3);
  hDataMCRatio->GetLowerRefXaxis()->SetNoExponent();
  hDataMCRatio->GetLowerRefXaxis()->SetMoreLogLabels();
  hDataMCRatio->GetUpperRefXaxis()->SetTitle("Dielectron Invariant Mass [GeV]");
  hDataMCRatio->GetLowerRefYaxis()->SetTitle("MC/Data Ratio");
  hDataMCRatio->GetUpperRefYaxis()->SetTitle("Entries");
  canvas1->Update();
  

  auto hDataMCRatiolinear = new TRatioPlot(hStacklinear,histos[BINS_DATA_LINEAR]);
  hDataMCRatiolinear->GetXaxis()->SetTitle("m_{ee} [GeV]");  
  canvas2->cd();
  hDataMCRatiolinear->Draw();
  hDataMCRatiolinear->GetUpperPad()->cd();
  legend->Draw("same");
  hDataMCRatiolinear->GetLowerRefGraph()->SetMinimum(0.7);
  hDataMCRatiolinear->GetLowerRefGraph()->SetMaximum(1.3);
  hDataMCRatiolinear->GetUpperRefXaxis()->SetTitle("Dielectron Invariant Mass [GeV]");
  hDataMCRatiolinear->GetLowerRefYaxis()->SetTitle("MC/Data Ratio");
  hDataMCRatiolinear->GetUpperRefYaxis()->SetTitle("Entries");
  canvas2->Update();

  canvas1->SaveAs("./plots/dataVsMClog.png");
  canvas2->SaveAs("./plots/dataVsMClinear.png");
}//end invMassDraw
