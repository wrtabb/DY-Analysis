#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF2.h"
#include "TMath.h"
#include "Math/PdfFuncMathCore.h"
#include "TFile.h"

Double_t trueDistribution(Double_t *val, Double_t *par){
  
  double norm = par[0];
  double expDecay = par[1];
  double peakNormRel = par[2];
  double peakMean = par[3];
  double peakWidth = par[4];

  double x = val[0];

  double result = norm * exp( -x/expDecay)
    + norm * peakNormRel * ROOT::Math::breitwigner_pdf( x, peakWidth, peakMean );

  return result;
}

Double_t resolutionFunction(Double_t *val, Double_t *par){
  
  double mean = par[0];
  double sigma = par[1];

  double x = val[0];

  bool isNormalized = kTRUE;
  
  double result = TMath::Gaus(x, mean, sigma, isNormalized);
  return result;
}

Double_t integrand2D(Double_t *val, Double_t *par){

  double norm = par[0];
  double expDecay = par[1];
  double peakNormRel = par[2];
  double peakMean = par[3];
  double peakWidth = par[4];
  //
  double x = val[0];
  double y = val[1];
  //
  double parTrue[5] = {norm, expDecay, peakNormRel, peakMean, peakWidth};
  double valTrue[1] = {y};

  //
  const double resMean = 0.0;
  double resWidth = par[5];
  double parRes[2] = {resMean, resWidth};
  double valRes[1] = {x-y};

  double result = trueDistribution(valTrue, parTrue) * resolutionFunction(valRes, parRes);
  return result;
}

Double_t integrand1D(Double_t *val, Double_t *par){

  //
  // Convolution integral: F(x) = integral( trueFunc(y) * resolFunc(x-y) ) dy,
  //    so the integrand is trueFunc(y) * resolFunc(x-y)
  //    where y is the variable and x is a parameter
  //
  double x = par[6];

  //
  double norm = par[0];
  double expDecay = par[1];
  double peakNormRel = par[2];
  double peakMean = par[3];
  double peakWidth = par[4];
  //
  double y = val[0];
  //
  double parTrue[5] = {norm, expDecay, peakNormRel, peakMean, peakWidth};
  double valTrue[1] = {y};

  //
  const double resMean = 0.0;
  double resWidth = par[5];
  double parRes[2] = {resMean, resWidth};
  double valRes[1] = {x-y};

  double result = trueDistribution(valTrue, parTrue) * resolutionFunction(valRes, parRes);
  return result;
}



Double_t observedDistribution(Double_t *val, Double_t *par){

  double norm = par[0];
  double expDecay = par[1];
  double peakNormRel = par[2];
  double peakMean = par[3];
  double peakWidth = par[4];
  double resWidth  = par[5];
  //
  double x = val[0];

  TF1 *integrand1DFunc = new TF1("integrand1DFunc",integrand1D, -10, 150, 7);
  integrand1DFunc->SetParameters(norm, expDecay, peakNormRel, peakMean, peakWidth, resWidth, x);

  // Integrate from -7*sigma to +7 sigma, should be close enough to from -infinity to infinity for a Gaussian
  double result = integrand1DFunc->Integral(x-7*resWidth, x+7*resWidth);
  return result;
}

Double_t migrationMatrix(Double_t *val, Double_t *par){

  double norm = par[0];
  double expDecay = par[1];
  double peakNormRel = par[2];
  double peakMean = par[3];
  double peakWidth = par[4];
  double resWidth  = par[5];
  //
  double x = val[0];
  double y = val[1];
  TF2 *integrand2DFunc = new TF2("integrand2DFunc",integrand2D, -10, 150, 0,100,6);
  integrand2DFunc->SetParameters(norm, expDecay, peakNormRel, peakMean, peakWidth, resWidth);

  // Integrate from -7*sigma to +7 sigma, should be close enough to from -infinity to infinity for a Gaussian
  double result = integrand2DFunc->Integral(x-7*resWidth, x+7*resWidth,y-7*resWidth, y+7*resWidth);
  return result;
}


void integrate(){

  // Draw functions
  TF1 *trueDist = new TF1("trueDist", trueDistribution, 0, 100, 5);
  trueDist->SetParameters(100, 20, 10, 60, 5);

//  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,800);
//  c1->cd();

  trueDist->SetNpx(1000);
//  trueDist->Draw();

  TF1 *observedDist = new TF1("observedDist", observedDistribution, 0, 100, 6);
  observedDist->SetParameters(100, 20, 10, 60, 5, 3);
  observedDist->SetNpx(1000);
  observedDist->SetLineColor(kBlue);
//  observedDist->Draw("same");

  // Draw histograms
  const int nbins = 100;
  const double xmin = 0;
  const double xmax = 100;
  TH1F *histTrue = new TH1F("histTrue", "histTrue", nbins, xmin, xmax);
  TH1F *histObs  = new TH1F("histObs", "histObs", nbins, xmin, xmax);
  for(int i=1; i<=nbins; i++){
    float yieldTrue = trueDist->Integral( histTrue->GetXaxis()->GetBinLowEdge(i),
				      histTrue->GetXaxis()->GetBinUpEdge(i));
    histTrue->SetBinContent(i,yieldTrue);

    float yieldObs = observedDist->Integral( histTrue->GetXaxis()->GetBinLowEdge(i),
					     histTrue->GetXaxis()->GetBinUpEdge(i));
    histObs->SetBinContent(i,yieldObs);
  }

//  TCanvas *c2 = new TCanvas("c2","c2",50,50,800,800);
//  c2->cd();
//  histTrue->Draw();
//  histObs->SetLineColor(kRed);
//  histObs->Draw("same");

  // 2D histogram: migration matrix
  TF2 *migrationMatrixFunc = new TF2("migrationMatrixFunc", migrationMatrix, -10, 150, 0, 100, 6);
  migrationMatrixFunc->SetParameters(100, 20, 10, 60, 5, 3);
  migrationMatrixFunc->SetNpx(1000);	
  migrationMatrixFunc->SetNpy(1000);	
  migrationMatrixFunc->Draw();

  TH2F *migrationHist = new TH2F("migrationHist", "migrationHist", nbins, xmin, xmax, nbins, xmin, xmax);
  
  for(int iBin = 1; iBin <= nbins; iBin++){
    for(int jBin = 1; jBin <= nbins; jBin++){
      double xlow = migrationHist->GetXaxis()->GetBinLowEdge(iBin);
      double xhi  = migrationHist->GetXaxis()->GetBinUpEdge (iBin);
      double ylow = migrationHist->GetYaxis()->GetBinLowEdge(jBin);
      double yhi  = migrationHist->GetYaxis()->GetBinUpEdge (jBin);
      double yield = migrationMatrixFunc->Integral(xlow, xhi, ylow, yhi);
      migrationHist->SetBinContent(iBin, jBin, yield);
    }
  }

  TCanvas *c3 = new TCanvas("c3","c3",100,100,800,800);
  c3->cd();
  TH1F*projX = (TH1F*)migrationHist->ProjectionX();
  TH1F*projY = (TH1F*)migrationHist->ProjectionY();
  projY->SetMarkerStyle(20);
  //projX->Draw("hist");
  //projY->Draw("pe,same");
  migrationHist->Draw("colz");
 
  TFile*saveFile = new TFile("unfoldingHistsIntegrationMethod.root","recreate");
  migrationHist->Write();
  histTrue->Write();
  histObs->Write();
  saveFile->Close();

}
