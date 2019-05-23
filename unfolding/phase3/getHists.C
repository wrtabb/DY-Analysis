#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"


const TString file1Name= "/home/hep/wrtabb/git/DY-Analysis/data/dataVsMC.root";
const TString file2Name= "/home/hep/wrtabb/git/DY-Analysis/data/mcHists.root";
const TString fileTestName = "/home/hep/wrtabb/git/DY-Analysis/efficiency/testMC.root";
const TString fileCorrName = "inputCorrelations_10000Samples.root";

void getHists()
{
 TH1::SetDefaultSumw2();
 TFile*file1 = new TFile(file1Name);
 TFile*file2 = new TFile(file2Name);
 TFile*fileTest = new TFile(fileTestName);
 TFile*fileCorr= new TFile(fileCorrName);

 //gROOT->SetBatch(true);
 gStyle->SetOptStat(0);
 gStyle->SetPalette(1);
 //TH1D*hMC = (TH1D*)file2->Get("hRecoInvMass");
 TH1D*hMC = (TH1D*)fileTest->Get("hMC");
  hMC->Rebin(2);
  hMC->SetName("hMC");
  hMC->SetFillColor(kOrange-2);
 TH1D*hData = (TH1D*)file1->Get("hDataInvMass");
  hData->SetName("hData");
 //TH1D*hTrue = (TH1D*)file2->Get("hTrue");
 TH1D*hTrue = (TH1D*)fileTest->Get("hTrue");
  hTrue->SetName("hTrue");
 TH1D*hBack = (TH1D*)file1->Get("hFakesInvMass");
 TH1D*hEW = (TH1D*)file1->Get("hEWInvMass");
 TH1D*hTops = (TH1D*)file1->Get("hTopsInvMass");
 TH2D*hCovM = (TH2D*)fileCorr->Get("hCovM");
 TH2D*hCovMinv = (TH2D*)fileCorr->Get("hCovMinv");
 //Set any negative bin values to 0
 for(int i=1;i<nLogBins2+1;i++){
  if(hBack->GetBinContent(i)<0) hBack->SetBinContent(i,0);
  if(hEW->GetBinContent(i)<0) hEW->SetBinContent(i,0);
  if(hTops->GetBinContent(i)<0) hTops->SetBinContent(i,0);
 }

 hBack->Add(hEW);
 hBack->Add(hTops);
 hBack->SetName("hBack");
 //TH2D*hMatrix = (TH2D*)file2->Get("migMatrixGENisHardvsReco"); 
 TH2D*hMatrix = (TH2D*)fileTest->Get("hMatrix"); 
  hMatrix->SetName("hMatrix");
  hMatrix->GetXaxis()->SetNoExponent();
  hMatrix->GetXaxis()->SetMoreLogLabels();
  hMatrix->GetYaxis()->SetNoExponent();
  hMatrix->GetYaxis()->SetMoreLogLabels();
  hMatrix->GetXaxis()->SetTitle("True mass [GeV]");
  hMatrix->GetYaxis()->SetTitle("Reco mass [GeV]");

 TCanvas*canMatrix = new TCanvas("canMatrix","",10,10,1000,1000);
  canMatrix->SetGrid();
  canMatrix->SetLogy();
  canMatrix->SetLogx();
  canMatrix->SetLogz();
  hMatrix->Draw("colz");
  canMatrix->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase3/matrix.png");
 TH1D*hBackRebin = (TH1D*)hBack->Clone("hBackRebin");
  hBackRebin->Rebin(2);
 TH1D*hDataRebin = (TH1D*)hData->Clone("hDataRebin");
  hDataRebin->Rebin(2);
 THStack*hStack = new THStack("hStack","");
  hStack->Add(hBackRebin);
  hStack->Add(hMC);
 TH1D*hStackHist = new TH1D("hStackHist","",nLogBins,massbins);
 double stackContent;
 for(int i=1;i<nLogBins+1;i++){
  stackContent = hBackRebin->GetBinContent(i) + hMC->GetBinContent(i); 
  hStackHist->SetBinContent(i,stackContent); 
 }

 TH1D*hRatio = (TH1D*)hDataRebin->Clone("ratio");
 float ratioRange = 0.3;
 float ratioUpper = 1.0 + ratioRange;
 float ratioLower = 1.0 - ratioRange;
 hRatio->Divide(hStackHist);
 hRatio->SetLineColor(kBlack);
 hRatio->SetMarkerColor(kBlack);
 hRatio->SetMarkerStyle(20);
 hRatio->GetYaxis()->SetLabelSize(0.06);
 hRatio->GetYaxis()->SetTitleSize(0.08);
 hRatio->GetYaxis()->SetTitleOffset(0.3);
 hRatio->GetYaxis()->SetTitle("data/MC");
 hRatio->GetXaxis()->SetLabelSize(0.1);
 hRatio->GetXaxis()->SetTitleSize(0.1);
 hRatio->GetXaxis()->SetNoExponent();
 hRatio->GetXaxis()->SetMoreLogLabels();
 hRatio->GetYaxis()->SetRangeUser(ratioLower,ratioUpper);

 TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
  legend->SetTextSize(0.02);
  legend->AddEntry(hData,"Data");
  legend->AddEntry(hMC,"MC Signal");
  legend->AddEntry(hBack,"MC Background (tops, EW, Fakes)");

 //TStack doesn't allow some histogram options to be definined
 //hBlank is a blank histogram to hold these options
 TH1D*hBlank = new TH1D("hBlank","",nLogBins,massbins);
  hBlank->GetXaxis()->SetLabelSize(0);
  hBlank->GetYaxis()->SetRangeUser(0.1,2.5e8);
  hBlank->SetTitle("Monte Carlo and Data Before Unfolding");
 float marg = 0.03;
 float bot = 0.2;
 TCanvas*canvas = new TCanvas("canvas","",10,10,1000,1000);
 TVirtualPad*p;
 TLine*line = new TLine(0,1,3000,1);
 line->SetLineColor(kRed);
 TPad*pad1 = new TPad("pad1","",0,0.3,1,1.0);
 pad1->SetBottomMargin(marg);
 pad1->SetGrid();
 pad1->SetLogx();
 pad1->SetLogy();
 pad1->SetTicks(1,1);
 pad1->Draw();
 pad1->cd();
 hBlank->Draw();
 hStack->Draw("hist,same");
 hDataRebin->Draw("PE,same");
 legend->Draw("same");
 canvas->Update();

 canvas->cd();
 TPad*pad2 = new TPad("pad2","",0,0.05,1,0.3);
 p = pad2;
 pad2->SetTopMargin(marg);
 pad2->SetBottomMargin(bot);
 pad2->SetGrid();
 pad2->SetTicks(1,1);
 pad2->SetLogx();
 pad2->Draw();
 pad2->cd();
 hRatio->Draw("PE");
 line->Draw("same");
 canvas->Update();

 canvas->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase3/dataVsMC.png");
 TFile*outputFile = new TFile("outputDataUnfold.root","recreate");
 hMC->Write();
 hData->Write();
 hTrue->Write();
 hBack->Write();
 hMatrix->Write();
 outputFile->Write();
 outputFile->Close();
}
