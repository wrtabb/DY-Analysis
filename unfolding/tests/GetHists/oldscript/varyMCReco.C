
TString fileName = "histograms_for_unfolding.root";
TString outputName = "mc_varied_for_unfolding.root";
void varyMCReco()
{
 TH1::SetDefaultSumw2();
 TFile*file = new TFile(fileName);
 TH1D*hist = (TH1D*)file->Get("hMC");
 hist->SetLineColor(kRed);
 TH1D*newHist = (TH1D*)hist->Clone();
 newHist->SetLineColor(kBlack); 
 newHist->SetMarkerColor(kBlack);
 newHist->SetMarkerStyle(20);

 TRandom3 gen;
 double variation;
 double newBinContent;
 double binContent;
 double mean;
 double sigma;

 for(int i=1;i<=hist->GetNbinsX();i++){
  binContent = hist->GetBinContent(i);
  mean = binContent;
  sigma = binContent*0.01;
  newBinContent = gen.Gaus(mean,sigma);
  newHist->SetBinContent(i,newBinContent);
 }
 TCanvas*canvas = new TCanvas("canvas","",0,0,1200,1000);
 canvas->SetGrid();
 canvas->SetLogx();
 canvas->SetLogy();
 hist->GetXaxis()->SetMoreLogLabels();
 hist->GetXaxis()->SetNoExponent();
 hist->Draw("hist");
 newHist->Draw("pe,same");

 TFile*outFile = new TFile(outputName,"recreate");
 newHist->Write();
 outFile->Write();
 outFile->Close();
}
