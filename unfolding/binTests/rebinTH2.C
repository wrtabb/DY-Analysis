#include "VariableList.h"

const TString fileName = "unfoldingHists.root";
TH2D*RebinTH2(TH2D*hist,TString histName,std::vector<double> binning);
TH2D*RebinTH2(TH2D*hist,TString histName,TH2D*hBinning);

void rebinTH2()
{
 //gROOT->SetBatch(true);
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);
 TFile*file = new TFile(fileName);
 TH2D*hMatrix0 = (TH2D*)file->Get("hMatrix0");
 TH2D*hMatrix1 = (TH2D*)file->Get("hMatrix1");
 TH2D*hRebin = RebinTH2(hMatrix0,"histRebin",binningReco1);

 TH1D*proj0 = hMatrix0->ProjectionY();
 TH1D*proj1 = hMatrix1->ProjectionY();
 TH1D*proj0Rebin = hRebin->ProjectionY();
 int nBins = proj1->GetNbinsX();
 double diff;
 for(int i=1;i<=nBins;i++){
  diff = abs(proj1->GetBinContent(i)-proj0Rebin->GetBinContent(i));
  cout << i << ", " << diff << endl;
 }

 TCanvas*canvas=new TCanvas("canvas","",0,0,1000,1000);
 canvas->SetGrid();
 proj0->SetMarkerStyle(20);
 proj0->SetMarkerColor(kBlue);
 proj0Rebin->SetMarkerStyle(20);
 proj0Rebin->SetMarkerColor(kRed);
 proj1->Draw("hist");
 proj0->Draw("pe,same");
 proj0Rebin->Draw("pe,same");
}

TH2D*RebinTH2(TH2D*hist,TString histName,std::vector<double> binning)
{
 int nBinsHist = hist->GetNbinsY();
 int nBinsReco = binning.size()-1;
 double newbinning[nBinsReco];
 for(int i=0;i<=nBinsReco;i++){
  newbinning[i] = binning.at(i);
 }
 TH2D*hRebin = new TH2D(histName,"",nBinsTrue,binningTrue,nBinsReco,newbinning);

 double y,x;
 for(int i=1;i<=nBinsTrue;i++){ 
  for(int j=1;j<=nBinsHist;j++){
   x = hist->GetXaxis()->GetBinCenter(i);  
   y = hist->GetYaxis()->GetBinCenter(j);  
   int nEntries = hist->GetBinContent(i,j);
   for(int k=0;k<nEntries;k++){
    hRebin->Fill(x,y);
   }//end filling hRebin
  }//end y bin loop
 } //end x bin loop
 return hRebin;
}

TH2D*RebinTH2(TH2D*hist,TString histName,TH2D*hBinning)
{
 int nBinsHist = hist->GetNbinsY();
 int nBinsReco = hBinning->GetNbinsY();
 double newbinning[nBinsReco];
 TH1D*hBinningProjection = hBinning->ProjectionY();
 for(int i=0;i<=nBinsReco;i++){
  if(i==0) newbinning[i] = hBinningProjection->GetBinLowEdge(i+1);
  else newbinning[i] = newbinning[i-1]+hBinningProjection->GetBinWidth(i);
 }
 TH2D*hRebin = new TH2D(histName,"",nBinsTrue,binningTrue,nBinsReco,newbinning);

 double y,x;
 for(int i=1;i<=nBinsTrue;i++){ 
  for(int j=1;j<=nBinsHist;j++){
   x = hist->GetXaxis()->GetBinCenter(i);  
   y = hist->GetYaxis()->GetBinCenter(j);  
   int nEntries = hist->GetBinContent(i,j);
   for(int k=0;k<nEntries;k++){
    hRebin->Fill(x,y);
   }//end filling hRebin
  }//end y bin loop
 } //end x bin loop
 return hRebin;
}
