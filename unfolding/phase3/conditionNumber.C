#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"
#include "/home/hep/wrtabb/git/DY-Analysis/headers/drawOptions.h"

const TString file1Name= "/home/hep/wrtabb/git/DY-Analysis/data/dataVsMC.root";
const TString file2Name= "outputDataUnfold.root";

void conditionNumber()
{
 //gROOT->SetBatch(true);
 gStyle->SetOptStat(0);
 gStyle->SetPalette(1);
 TFile*file1 = new TFile(file1Name);
 TFile*file2 = new TFile(file2Name);

 TH2D*hMatrix = (TH2D*)file2->Get("hMatrix");
  hMatrix->SetName("hMatrix");
  hMatrix->RebinY(2);
 hist2DPlot(0,hMatrix,"colz",true,true,true); 
 TMatrixD matrix(nLogBins,nLogBins);
 for(int i=0;i<nLogBins;i++){
  for(int j=0;j<nLogBins;j++){
   matrix(i,j) = hMatrix->GetBinContent(i,j);
  }
 }
 TDecompSVD svd(matrix);
 TVectorD sig = svd.GetSig();
 double condN = sig.Max()/sig.Min();
 sig.Print();
 cout << "Condition number = " << condN << endl;
}
