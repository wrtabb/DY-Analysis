#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"
#include "/home/hep/wrtabb/git/DY-Analysis/headers/drawOptions.h"

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
  hMatrix->RebinY(2);
 TH2D*hResponse = new TH2D("hResponse","",nLogBins,massbins,nLogBins,massbins);
 TMatrixD matrix(nLogBins,nLogBins);
 TMatrixD response(nLogBins,nLogBins);
 for(int i=0;i<nLogBins;i++){
  double sum = 0;
  for(int j=0;j<nLogBins;j++){
   matrix(i,j) = hMatrix->GetBinContent(i+1,j+1);
   sum += matrix(i,j);
  }
  for(int j=0;j<nLogBins;j++){
   response(i,j) = matrix(i,j)/sum;
   hResponse->SetBinContent(i+1,j+1,response(i,j));
  }
 }
 TMatrixD responseInv = response;
 double det;
 responseInv.Invert(&det);
 TDecompSVD svd(response);
 TVectorD sig = svd.GetSig();
 double condN = sig.Max()/sig.Min();
 hist2DPlot(0,hResponse,"colz",true,true,true); 
 sig.Print();
 cout << "Condition number = " << condN << endl;
 cout << "Determinant = " << det << endl;
}
