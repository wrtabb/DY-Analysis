#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"


const TString file1Name= "/home/hep/wrtabb/git/DY-Analysis/data/dataVsMC.root";
const TString file2Name= "/home/hep/wrtabb/git/DY-Analysis/data/efficiencyAndMigration.root";


void conditionNumber()
{
 TFile*file1 = new TFile(file1Name);
 TFile*file2 = new TFile(file2Name);

 TH2D*hMatrix = (TH2D*)file2->Get("migMatrixGENisHardvsReco");
  hMatrix->SetName("hMatrix");
 TMatrixD matrix(nLogBins2,nLogBins2);

 for(int i=1;i<nLogBins+1;i++){
  for(int j=1;j<nLogBins2+1;j++){
   if(j>nLogBins+1) matrix(j-1,i-1) = 0;
   else matrix(j-1,i-1) = hMatrix->GetBinContent(i,j);
  }
 }
 
/*
 TMatrixD A(3,3);
 A(0,0) = 1;
 A(0,1) = 0;
 A(1,2) = 3;
 A(1,0) = 0;
 A(1,1) = -3;
 A(1,2) = 9;
 A(2,0) = 1;
 A(2,1) = 2;
 A(2,2) = 4;
*/

 TDecompSVD svd(matrix);
 TVectorD sig = svd.GetSig();
 double condN = sqrt(sig(0)/sig(42));
 sig.Print();
 cout << "Condition number = " << condN << endl;
 cout << "Alternate condition number:" << svd.Condition() << endl;
}
