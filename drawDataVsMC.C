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

const TString pileupFileName = "./plots/MyDataPileupHistogram.root";
const TString dataFileName = "./plots/dataVsMC.root";
const TString histName[nHistos] = {"hFakes", "hFakeslinear", "hEW", "hEWlinear", "hTops", "hTopslinear",   
				   "hMCInvMass", "hMCInvMasslinear", "hDataInvMass", "hDataInvMasslinear"};
const Color_t histFillColors[8] = {kViolet+5, kViolet+5, kRed+2, kRed+2, kBlue+2, kBlue+2, kOrange-2, kOrange-2};
const Color_t histLineColors[8] = {kViolet+3, kViolet+3, kRed+4, kRed+4, kBlue+3, kBlue+3, kOrange-3, kOrange-3};

const int nPileupBins = 75;
//obtained from https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py#L25
const double pileupMC[nPileupBins] = {1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,
					0.000140973 , 0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,
					0.00919534 ,0.0146697 , 0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,
					0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 , 0.0559937 ,0.0554468 ,
					0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,
					0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,
					0.0142498 ,0.012804 , 0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,
					0.00829292 ,0.0076195 ,0.0069806 ,0.0062025, 0.00546581 ,0.00484127 ,
					0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 , 0.00117884 ,
					0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,
					9.88128e-05, 6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,
					1.73032e-05 ,1.435e-05 ,1.36486e-05, 1.35555e-05 ,1.37491e-05 ,
					.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 , 1.34177e-05 ,
					1.32959e-05 ,1.33287e-05};

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
  /*
  TCanvas*canvas1 = new TCanvas("canvas1","",10,10,1000,1000);
  canvas1->SetLogx();
  canvas1->SetLogy();

  TCanvas*canvas2 = new TCanvas("canvas2","",10,10,1000,1000);
  canvas2->SetLogy();

  auto hDataMCRatio = new TRatioPlot(hStack,histos[BINS_DATA]);
  hDataMCRatio->GetXaxis()->SetTitle("m_{ee} [GeV]");  
  canvas1->cd();
  hDataMCRatio->Draw();
  hDataMCRatio->GetUpperPad()->cd();
  legend->Draw("same");
  hDataMCRatio->GetLowerRefGraph()->SetMinimum(0.5);
  hDataMCRatio->GetLowerRefGraph()->SetMaximum(1.5);
  hDataMCRatio->GetLowerRefXaxis()->SetNoExponent();
  hDataMCRatio->GetLowerRefXaxis()->SetMoreLogLabels();
  hDataMCRatio->GetUpperRefXaxis()->SetTitle("Dielectron Invariant Mass [GeV]");
  hDataMCRatio->GetLowerRefYaxis()->SetTitle("MC/Data Ratio");
  hDataMCRatio->GetLowerRefYaxis()->SetLabelSize(0.02);
  hDataMCRatio->GetUpperRefYaxis()->SetTitle("Entries");
  canvas1->Update();  

  auto hDataMCRatiolinear = new TRatioPlot(hStacklinear,histos[BINS_DATA_LINEAR]);
  hDataMCRatiolinear->GetXaxis()->SetTitle("m_{ee} [GeV]");  
  canvas2->cd();
  hDataMCRatiolinear->Draw();
  hDataMCRatiolinear->GetUpperPad()->cd();
  legend->Draw("same");
  hDataMCRatiolinear->GetLowerRefGraph()->SetMinimum(0.5);
  hDataMCRatiolinear->GetLowerRefGraph()->SetMaximum(1.5);
  hDataMCRatiolinear->GetUpperRefXaxis()->SetTitle("Dielectron Invariant Mass [GeV]");
  hDataMCRatiolinear->GetLowerRefYaxis()->SetTitle("MC/Data Ratio");
  hDataMCRatiolinear->GetLowerRefYaxis()->SetLabelSize(0.02);
  hDataMCRatiolinear->GetUpperRefYaxis()->SetTitle("Entries");
  canvas2->Update();

  //canvas1->SaveAs("./plots/dataVsMClog.png");
  //canvas2->SaveAs("./plots/dataVsMClinear.png");
  */

  //Pileup
  float norm = 1.0;
  TFile*filePU = new TFile(pileupFileName);
  TH1F*hPileupData = new TH1F("hPileupData","",100,0,100);
  hPileupData->Sumw2();
  hPileupData->GetXaxis()->SetTitle("Pileup");
  hPileupData->SetLineColor(kBlue+2);
  hPileupData->SetLineWidth(2);

  TH1F*hPileupMC = new TH1F("hPileupMC","",100,0,100);
  hPileupMC->Sumw2();
  hPileupMC->GetXaxis()->SetTitle("Pileup");
  hPileupMC->SetLineColor(kOrange+2);
  hPileupMC->SetLineWidth(2);

  for(int i=0; i<nPileupBins;i++){
    hPileupMC->SetBinContent(i,pileupMC[i]);
    hPileupData = (TH1F*)filePU->Get("pileup");
  }
  hPileupData->Scale(norm/hPileupData->Integral());
  hPileupMC->Scale(norm/hPileupMC->Integral());
  
  TH1F*hPileupRatio = (TH1F*)hPileupData->Clone();
  hPileupRatio->Divide(hPileupMC);
  hPileupRatio->SetMarkerStyle(20);
  hPileupRatio->SetMarkerColor(kBlack);
  hPileupRatio->SetMarkerSize(0.5);
  hPileupRatio->SetTitle("Data/MC Ratio");

  TCanvas*canvas3 = new TCanvas("canvas3","",10,10,1400,700);
  canvas3->Divide(2);
  canvas3->cd(1);
  hPileupData->Draw("hist");
  hPileupMC->Draw("hist,same");
  canvas3->cd(2);
  hPileupRatio->Draw("PE");

  TFile*pileupSaveFile = new TFile("./plots/pileup.root","RECREATE");
  pileupSaveFile->cd();
  hPileupData->Write();
  hPileupMC->Write();
  hPileupRatio->Write();
  canvas3->Write();
  pileupSaveFile->Write();
  pileupSaveFile->Close();
  
  
}//end invMassDraw
