#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"
#include "/home/hep/wrtabb/git/DY-Analysis/headers/drawOptions.h"
const TString inputFileName = "/home/hep/wrtabb/git/DY-Analysis/unfolding/phase3/outputDataUnfold.root";

void inversionTest()
{
 gStyle->SetOptStat(0);
 gStyle->SetPalette(1);

 //Define file
 TFile*inputFile = new TFile(inputFileName);

 //Initialize histograms
 TH2D*hMatrix = (TH2D*)inputFile->Get("hMatrix");
 TH2D*hResponse = new TH2D("hResponse","",nLogBins,massbins,nLogBins,massbins);
 TH1D*hMC = (TH1D*)inputFile->Get("hMC");
 TH1D*hData = (TH1D*)inputFile->Get("hData");
 TH1D*hTrue = (TH1D*)inputFile->Get("hTrue");
 TH1D*hBack = (TH1D*)inputFile->Get("hBack");
 TH1D*hUnfolded = (TH1D*)hTrue->Clone("hUnfolded");

 //Rebin so that the matrix will be square
 hData->Rebin(2);
 hBack->Rebin(2);
 hMatrix->RebinY(2);

 //Subtract background from data
 hData->Add(hBack,-1);
 
 //Initialize matrix and vectors
 TMatrixD matrix(nLogBins,nLogBins);
 TMatrixD response(nLogBins,nLogBins);
 TMatrixD unfold(nLogBins,nLogBins);
 TVectorD vData(nLogBins);

 //Definte matrix and vectors
 for(int i=0;i<nLogBins;i++){
  for(int j=0;j<nLogBins;j++){
   matrix(i,j) = hMatrix->GetBinContent(i+1,j+1);
   if(i==0) vData(j) = hData->GetBinContent(j+1);
  }
 }

 //Normalize matrix to make response matrix
 //Each column adds up to 1
 for(int i=0;i<nLogBins;i++){
  double sum = 0;
  for(int j=0;j<nLogBins;j++){
   sum += matrix(i,j);
  }
  for(int j=0;j<nLogBins;j++){
   if(sum!=0) response(i,j) = matrix(i,j)/sum;
   else cout << "Could not normalize! i = " << i << ", j = " << j << endl;
  }
 }
 
/*
 //Check that response matrix is properly normalized
 for(int i=0;i<nLogBins;i++){
  double testSum = 0;
  for(int j=0;j<nLogBins;j++){
   testSum += response(i,j);
  }
  //each testSum must equal 1
  cout << "The sum of column i = " << i << " is: " << testSum << endl;
 }
*/

 //invert response matrix to get unfolding matrix
 unfold = response;
 unfold.Invert();

 //multiply tranposed unfolding matrix by input vector to get unfolded vector
 TVectorD vUnfolded = (unfold.T())*vData;
 
 //Place vUnfolded and response in histograms for plotting
 for(int i=0;i<nLogBins;i++){
  hUnfolded->SetBinContent(i+1,vUnfolded(i));
  for(int j=0;j<nLogBins;j++){
   hResponse->SetBinContent(i+1,j+1,response(i,j));
  }
 }

 //draw results
 TCanvas*c2 = new TCanvas("c2","",0,0,1000,1000);
 hResponse->SetTitle("Response matrix");
 hResponse->GetXaxis()->SetTitle("True mass [GeV]");
 hResponse->GetYaxis()->SetTitle("Reco mass [GeV]");
 hist2DPlot(c2,hResponse,"colz",true,true,true);
 c2->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase3/responseMatrix.png");

 TCanvas*canvas = new TCanvas("canvas","",0,0,1200,1000);
 hUnfolded->SetFillColor(kRed+2);
 hUnfolded->SetTitle("Unfolding: Inversion method");
 hTrue->SetMarkerStyle(20);
 histPlot(canvas,hUnfolded,"hist",hTrue,"PE,same",true,true);
 canvas->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase3/unfoldInversionTest.png");

 TCanvas*canMatrix = new TCanvas("canMatrix","",0,0,1000,1000);
 hMatrix->SetTitle("Migration matrix");
 hist2DPlot(canMatrix,hMatrix,"colz",true,true,true);
 canMatrix->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase3/migrationMatrix.png");
 
 TFile saveFile("tempSave.root","recreate");
 hResponse->Write();
 hUnfolded->Write();
 hTrue->Write();
 hMC->Write();
 hMatrix->Write();
 hData->Write();
 hBack->Write();
 saveFile.Close();
}

