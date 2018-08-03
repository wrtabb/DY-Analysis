#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
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
    histos[i]->GetXaxis()->SetTitle("m_{ee} [GeV]");
    histos[i]->GetXaxis()->SetMoreLogLabels();
    histos[i]->GetXaxis()->SetNoExponent();
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
  hStack->SetMinimum(axisLow);
  THStack*hStacklinear = new THStack("hStacklinear","");
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
  
  TCanvas*canvas1 = new TCanvas("canvas1","",10,10,1000,1000);
  canvas1->SetLogx();
  canvas1->SetLogy();
  TLegend*legend = new TLegend(0.6,0.9,0.9,0.7);
  legend->SetTextSize(0.02);
  legend->AddEntry(histos[BINS_DATA],"Data");
  legend->AddEntry(histos[BINS_EE],"#gamma^{*}/Z #rightarrow e^{-}e^{+}");
  legend->AddEntry(histos[BINS_TOPS],"t#bar{t}+tW+#bar{t}W");
  legend->AddEntry(histos[BINS_EW],"EW (Dibosons, #gamma^{*}/Z #rightarrow #tau^{-}#tau^{+})");
  legend->AddEntry(histos[BINS_FAKES],"Fakes (W+Jets)");
  hStack->Draw("bar");
  histos[BINS_DATA]->Draw("PE,same");
  legend->Draw("same");

  TCanvas*canvas2 = new TCanvas("canvas2","",10,10,1000,1000);
  canvas2->SetLogx();
  canvas2->SetLogy();
  hStacklinear->Draw("bar");
  histos[BINS_DATA_LINEAR]->Draw("PE,same");

  TCanvas*canvas3 = new TCanvas("canvas3","",10,10,1000,1000);
  canvas3->SetLogx();
  canvas3->SetLogy();
  histos[BINS_EE]->Draw("bar");

  canvas1->SaveAs("./plots/dataVsMClog.png");
  canvas2->SaveAs("./plots/dataVsMClinear.png");
  canvas3->SaveAs("./plots/DYEE.png");
}//end invMassDraw
