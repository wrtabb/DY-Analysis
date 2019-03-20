#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"

const TString fileName1 = "/home/hep/wrtabb/git/DY-Analysis/data/efficiencyAndMigration.root";
const TString fileName2 = "/home/hep/wrtabb/git/DY-Analysis/data/effMassRangeNoWeight.root";
//const int nMassBins= 11;
//const double massbins[] = {10,50,100,200,400,500,700,800,1000,1500,2000,3000};

void massBinEff()
{
 gStyle->SetOptStat(0);
 TFile*file1 = new TFile(fileName1);
 TFile*file2 = new TFile(fileName2);
 
 TEfficiency*hAccW = (TEfficiency*)file1->Get("Acceptance");
 TEfficiency*hEffW = (TEfficiency*)file1->Get("Efficiency");
 TH1D*hTotEffW = new TH1D("hTotEffW","",nLogBins,massbins);
 
 TEfficiency*hAcc = (TEfficiency*)file2->Get("Acceptance");
 TEfficiency*hEff = (TEfficiency*)file2->Get("Efficiency");
 TH1D*hTotEff = new TH1D("hTotEff","",nLogBins,massbins);

 double eff,acc;
 for(int i=1;i<nLogBins+1;i++){
  eff = hEffW->GetEfficiency(i);
  acc = hAccW->GetEfficiency(i);
  hTotEffW->SetBinContent(i,eff*acc);

  eff = hEff->GetEfficiency(i);
  acc = hAcc->GetEfficiency(i);
  hTotEff->SetBinContent(i,eff*acc);
 }

 hTotEffW->GetXaxis()->SetMoreLogLabels();
 hTotEffW->GetXaxis()->SetNoExponent();
 hTotEffW->GetYaxis()->SetRangeUser(0,1);
 hTotEffW->SetTitle("Efficiency with Weights");
 TCanvas*canvasW = new TCanvas("canvasW","",10,10,1400,1000);
 canvasW->SetGrid();
 canvasW->SetLogx(); 
 hTotEffW->Draw();

 hTotEff->GetXaxis()->SetMoreLogLabels();
 hTotEff->GetXaxis()->SetNoExponent();
 hTotEff->GetYaxis()->SetRangeUser(0,1);
 hTotEff->SetTitle("Efficiency without Weights");
 TCanvas*canvas = new TCanvas("canvas","",10,10,1400,1000);
 canvas->SetGrid();
 canvas->SetLogx(); 
 hTotEff->Draw();
}
