#include "VariableList.h"

const TString fileName = "unfoldingHists.root";
TH1D*RebinTH1(TH1D*hist,TString histName,std::vector<double> binning);

void rebinTH1()
{
 //gROOT->SetBatch(true);
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);
 //Get files and histogram
 TFile*file = new TFile(fileName);
 TH1D*hReco = (TH1D*)file->Get("hReco0");
 //Define binning (binnings defined in VariableList.h)
 std::vector<double> binningReco = binningTrueVec;
 //place binning into an array for the Root function
 int nBins = binningReco.size()-1;
 double newbinning[nBins];
 for(int i=0;i<=nBins;i++){
  newbinning[i] = binningReco.at(i);
 }

 //Create rebinned version of input histogram using my function
 TH1D*hRebinFunc = RebinTH1(hReco,"histRebin",binningReco);
 //Create rebinned version of input histogram using ROOT function
 TH1D*hRebinRoot=(TH1D*)hReco->Rebin(nBins,"histRebin",newbinning);

 TCanvas*c1 = new TCanvas("c1","",0,0,1000,1000);

 hRebinFunc->SetMarkerStyle(20);
 hRebinRoot->Draw("hist");
 hRebinFunc->Draw("pe,same");
 
}

TH1D*RebinTH1(TH1D*hist,TString histName,std::vector<double> binning)
{
 int nBinsHist = hist->GetNbinsX();
 int nBinsReco = binning.size()-1;
 double newbinning[nBinsReco];
 for(int i=0;i<=nBinsReco;i++){
  newbinning[i] = binning.at(i);
 }
 TH1D*hRebin = new TH1D(histName,"",nBinsReco,newbinning);

 int bin;
 double y,x;
 for(int j=1;j<=nBinsHist;j++){
  x = hist->GetXaxis()->GetBinCenter(j);  
  int nEntries = hist->GetBinContent(j);
  for(int k=0;k<nEntries;k++){
   hRebin->Fill(x);
  }//end filling hRebin
 } //end x bin loop
 return hRebin;
}

//int maxBin = histos[i][j]->GetMaximumBin();
//double x = histos[i][j]->GetXaxis()->GetBinCenter(maxBin);
