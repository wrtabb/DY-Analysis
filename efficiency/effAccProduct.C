#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"

const TString effLoc = "/home/hep/wrtabb/git/DY-Analysis/data/plotsDY.root";
const TString fileSave = "/home/hep/wrtabb/git/DY-Analysis/data/effAcc.root";
const int rangeBinNum = 11;
const double rangeBinning[] = {10,50,100,200,400,500,700,800,1000,1500,2000,3000};

void effAccProduct()
{
 TFile*effFile = new TFile(effLoc);
 
 TEfficiency*hEff = (TEfficiency*)effFile->Get("Efficiency");
 TEfficiency*hAcc = (TEfficiency*)effFile->Get("Acceptance");
 TH1D*hEffAcc = new TH1D("effAcc","",nLogBins,massbins);
 double eff, acc,avg;
 int j;
 for(int i=1;i<nLogBins+1;i++){
  eff = hEff->GetEfficiency(i);
  acc = hAcc->GetEfficiency(i);
  hEffAcc->SetBinContent(i,eff*acc);
 }

 TLine*lineLow = new TLine(50,0,50,1);
 TLine*lineHigh = new TLine(100,0,100,1);
 TCanvas*canvas = new TCanvas("canvas","",10,10,1400,1000);
  canvas->SetGrid();
  canvas->SetLogx();
 hEffAcc->GetYaxis()->SetRangeUser(0,1);
 hEffAcc->GetXaxis()->SetMoreLogLabels();
 hEffAcc->GetXaxis()->SetNoExponent();
 hEffAcc->Draw();
 lineLow->Draw("same");
 lineHigh->Draw("same");

 double avgBin,sumBin;
 int binLow = 8;
 int binHigh = 19;
 sumBin = 0;
 for(int i=binLow;i<binHigh;i++){
  sumBin = sumBin + hEffAcc->GetBinContent(i);
  cout << hEffAcc->GetBinContent(i) << endl;
 }
 avgBin = sumBin/(binHigh-binLow);
 cout << "Average: " << avgBin << endl;
 TFile*saveFile = new TFile(fileSave,"recreate");
 hEff->Write();
 hAcc->Write();
 hEffAcc->Write();
 saveFile->Write();
 saveFile->Close();
}
