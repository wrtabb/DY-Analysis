#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"


const TString file1Name= "/home/hep/wrtabb/git/DY-Analysis/data/dataVsMC.root";
const TString file2Name= "outputDataUnfold.root";


void conditionNumber()
{
 gStyle->SetOptStat(0);
 gStyle->SetPalette(1);
 TFile*file1 = new TFile(file1Name);
 TFile*file2 = new TFile(file2Name);

 TH2D*hMatrix = (TH2D*)file2->Get("hMatrix");
  hMatrix->SetName("hMatrix");
 TMatrixD matrix(nLogBins2,nLogBins);
 for(int i=1;i<nLogBins+1;i++){
  for(int j=1;j<nLogBins2+1;j++){
   matrix(j-1,i-1) = hMatrix->GetBinContent(i,j);
  }
 }
 //matrix.Draw("colz");
 TDecompSVD svd(matrix);
 TVectorD sig = svd.GetSig();
 double condN = sig.Max()/sig.Min();
 sig.Print();
 cout << "Condition number = " << condN << endl;
}
