#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"

const TString effLoc = "/home/hep/wrtabb/git/DY-Analysis/plots/plotsDY.root";
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

 hEffAcc->Draw();
 TFile*saveFile = new TFile("effAcc.root","recreate");
 hEff->Write();
 hAcc->Write();
 hEffAcc->Write();
 //hEffAccRebin->Write();
 saveFile->Write();
 saveFile->Close();
}
