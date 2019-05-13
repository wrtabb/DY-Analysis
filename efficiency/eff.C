#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"
const TString fileName = "/home/hep/wrtabb/git/DY-Analysis/data/efficiencyAndMigration.root";
const int nHists = 5;
const TString effName[nHists] = {
 "Acceptance",
 "Reco-Gen matched efficiency",
 "ID efficiency",
 "HLT efficiency",
 "Total efficiency"
};
void eff()
{
 gStyle->SetOptStat(0);
 TFile*inputFile = new TFile(fileName);
 TH1D*hist[nHists];
 TEfficiency*eff[nHists];

 hist[0] = (TH1D*)inputFile->Get("hGenAllInvMass");
 hist[1] = (TH1D*)inputFile->Get("hGenInvMass");
 hist[2] = (TH1D*)inputFile->Get("hGenMatchedInvMass");
 hist[3] = (TH1D*)inputFile->Get("hGenPassIDInvMass");
 hist[4] = (TH1D*)inputFile->Get("hHLTGenInvMass");

 for(int i=0;i<nHists;i++){
  hist[i]->SetBinContent(0,0);
  if(i>0) {
   eff[i-1] = new TEfficiency((*hist[i]),(*hist[i-1]));
   eff[i-1]->SetName(effName[i]);
  }
 }
 eff[4] = new TEfficiency((*hist[4]),(*hist[1]));
 eff[4]->SetName("Total efficiency");
 eff[4]->SetTitle("Total efficiency");
 
 TCanvas*canvas = new TCanvas("canvas","",10,10,1200,1000);
 canvas->Divide(2,2);
 canvas->SetLogx();
 TH1D*hBlank[nHists];

 for(int i=0;i<nHists;i++){
  hBlank[i] = new TH1D("","",1000,15,3000);
  hBlank[i]->SetMinimum(0);
  hBlank[i]->SetMaximum(1);
  hBlank[i]->SetTitle(effName[i]);
  hBlank[i]->GetXaxis()->SetNoExponent();
  hBlank[i]->GetXaxis()->SetMoreLogLabels();
  hBlank[i]->GetXaxis()->SetTitle("mass [GeV]");
  if(i>nHists-2) continue;
  canvas->cd(i+1);
  gPad->SetLogx();
  gPad->SetGrid();
  hBlank[i]->Draw();
  eff[i]->Draw("same");
 }
 TCanvas*canvas2 = new TCanvas("canvas2","",10,10,1200,1000);
 canvas2->SetGrid();
 canvas2->SetLogx();
 hBlank[4]->Draw();
 eff[4]->Draw("same"); 

 canvas->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/efficiency/efficiencies.png");
 canvas2->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/efficiency/totalEfficiency.png");
}
