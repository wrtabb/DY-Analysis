#include "/home/hep/wrtabb/DY-Analysis/headers/header1.h"
#include "/home/hep/wrtabb/DY-Analysis/headers/drawOptions.h"

const TString fileName = "/home/hep/wrtabb/DY-Analysis/unfolding/phase3/outputDataUnfold.root";
//Unfold MC or Data
//MC is a closure test
const bool isMC = true;

void inversionTest()
{
 gStyle->SetOptStat(0);
 gStyle->SetPalette(1);

 //Define file
 TFile*inputFile = new TFile(fileName);

 //Initialize histograms
 TH2D*hMatrix = (TH2D*)inputFile->Get("hMatrix");
 TH2D*hResponse = new TH2D("hResponse","",nLogBins,massbins,nLogBins,massbins);
 TH1D*hMC = (TH1D*)inputFile->Get("hMC");
 TH1D*hData = (TH1D*)inputFile->Get("hData");
 TH1D*hTrue = (TH1D*)inputFile->Get("hTrue");
 TH1D*hBack = (TH1D*)inputFile->Get("hBack");
 TH1D*hUnfolded = new TH1D("hUnfolded","",nLogBins,massbins); 

 //Rebin so that the matrix will be square
 hMC->Rebin(2);
 hData->Rebin(2);
 hBack->Rebin(2);
 hMatrix->RebinY(2);

 //Subtract background from data
 hData->Add(hBack,-1);
 
 nLogBins = nLogBins + 1;
 //Initialize matrix and vectors
 TMatrixD matrix(nLogBins,nLogBins);
 TMatrixD response(nLogBins,nLogBins);
 TMatrixD unfold(nLogBins,nLogBins);
 TVectorD vData(nLogBins);

 //Define matrix and vectors
 for(int i=0;i<nLogBins;i++){
  for(int j=0;j<nLogBins;j++){
   matrix(i,j) = hMatrix->GetBinContent(i,j);
  }
  if(isMC) vData(i) = hMC->GetBinContent(i);
  else vData(i) = hData->GetBinContent(i);
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

 //invert response matrix to get unfolding matrix
 unfold = response;
 double det;
 unfold.Invert(&det);
 cout << "Determinant of unfolding matrix: " << det << endl;

 //multiply tranposed unfolding matrix by input vector to get unfolded vector
 TVectorD vUnfolded = (unfold.T())*vData;
 
 //Place vUnfolded and response in histograms for plotting
 for(int i=1;i<nLogBins;i++){
  hUnfolded->SetBinContent(i,vUnfolded(i));
  for(int j=1;j<nLogBins;j++){
   hResponse->SetBinContent(i,j,response(i,j));
  }
 }

 //draw results
 TCanvas*c3 = new TCanvas("c3","",0,0,1000,1000);
 c3->SetGrid();
 unfold.Draw("colz");
 c3->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase3/unfoldingMatrix.png");

 TCanvas*c2 = new TCanvas("c2","",0,0,1000,1000);
 hResponse->SetTitle("Response matrix");
 hResponse->GetXaxis()->SetTitle("True mass [GeV]");
 hResponse->GetYaxis()->SetTitle("Reco mass [GeV]");
 hist2DPlot(c2,hResponse,"colz",true,true,true);
 c2->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase3/responseMatrix.png");

 TCanvas*canvas = new TCanvas("canvas","",0,0,1200,1000);
 hUnfolded->SetTitle("Unfolding: Inversion method");
 hTrue->SetMarkerStyle(20);
 ratioPlot(canvas,hTrue,hData,hUnfolded);
 TString unfSaveName = "/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase3/unfoldInversionTest";
 if(isMC) unfSaveName += "_MC.png";
 else unfSaveName += "_Data.png";
 canvas->SaveAs(unfSaveName);

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

