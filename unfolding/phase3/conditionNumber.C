#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"


const TString file1Name= "/home/hep/wrtabb/git/DY-Analysis/data/dataVsMC.root";
const TString file2Name= "/home/hep/wrtabb/git/DY-Analysis/data/efficiencyAndMigration.root";


void conditionNumber()
{
 TFile*file1 = new TFile(file1Name);
 TFile*file2 = new TFile(file2Name);

 TH2D*hMatrix = (TH2D*)file2->Get("migMatrixGENisHardvsReco");
  hMatrix->SetName("hMatrix");
 TMatrixD matrix(nLogBins2,nLogBins);

 for(int i=1;i<nLogBins+1;i++){
  for(int j=1;j<nLogBins2+1;j++){
   matrix(j-1,i-1) = hMatrix->GetBinContent(i,j);
  }
 }
 
/*
 TRandom3*rand = new TRandom3();
 TMatrixD matrix(3,3);
 for(int i=0;i<3;i++){
  for(int j=0;j<3;j++){
   double r = rand->Rndm()+r;
   
   if(j==1&&i==0) matrix(i,j) = 1;
   else matrix(i,j) = r;
  }
 }
*/
 //matrix.Print();
 TDecompSVD svdMatrix(matrix);
 svdMatrix.Print();

 cout << svdMatrix.Condition() << endl; 
 
}
