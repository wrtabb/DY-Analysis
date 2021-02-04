#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF2.h"
#include "TFile.h"
#include "TMath.h"
#include "Math/PdfFuncMathCore.h"

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

void DrawProjections(TH1F*hTrue,TH1F*hReco,TH2F*hMatrix)
{
	TH1F*projX = (TH1F*)hMatrix->ProjectionX();
	TH1F*projY = (TH1F*)hMatrix->ProjectionY();
	TH1F*ratioY = (TH1F*)projY->Clone("ratioY");
	ratioY->Divide(hTrue);
	TH1F*ratioX = (TH1F*)projX->Clone("ratioX");
	ratioX->Divide(hReco);

	hTrue->SetLineColor(kRed);
	hReco->SetLineColor(kBlue);
	projX->SetMarkerColor(kBlue);
	projX->SetMarkerStyle(25);
	projY->SetMarkerColor(kRed);
	projY->SetMarkerStyle(25);

	int lastBin = hTrue->GetNbinsX();
	double binLow = hTrue->GetXaxis()->GetBinLowEdge(1);
	double binHigh = hTrue->GetXaxis()->GetBinUpEdge(lastBin);
	const float padmargins = 0.03;

	TCanvas*canvas1 = new TCanvas("canvas1","",10,10,1200,1000);
	TPad*pad1 = new TPad("","",0,0.3,1.0,1.0);
	pad1->SetBottomMargin(padmargins);
	pad1->SetGrid();
	pad1->SetTicks(1,1);
	pad1->Draw();
	pad1->cd();
	TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
	legend->SetTextSize(0.02);
	legend->AddEntry(hTrue,"True Distribution");
	legend->AddEntry(hReco,"Measured Distribution");
	legend->AddEntry(projX,"x-projection");
	legend->AddEntry(projY,"y-projection");

	hTrue->SetLabelSize(0);
	hTrue->SetTitleSize(0);
	hTrue->SetTitle("projections vs. histograms");
	hTrue->Draw("hist");
	hReco->Draw("hist,same");
	projX->Draw("pe,same");
	projY->Draw("pe,same");
	legend->Draw("same");

	canvas1->cd();
	TPad*pad2 = new TPad("","",0,0.05,1,0.3);
	pad2->SetTopMargin(padmargins);
	pad2->SetBottomMargin(0.2);
	pad2->SetGrid();
	pad2->SetTicks(1,1);
	pad2->Draw();
	pad2->cd();

	ratioX->SetMinimum(0.7);
	ratioX->SetMaximum(1.3);
	ratioX->GetYaxis()->SetLabelSize(0.06);
	ratioX->GetYaxis()->SetTitleSize(0.08);
	ratioX->GetYaxis()->SetTitleOffset(0.3);
	ratioX->GetYaxis()->SetTitle("proj/hist");
	ratioX->GetXaxis()->SetLabelSize(0.1);
	ratioX->GetXaxis()->SetTitleSize(0.1);
	ratioX->SetMarkerStyle(20);
	ratioX->SetMarkerColor(kBlue);

	ratioY->SetMinimum(0.7);
	ratioY->SetMaximum(1.3);
	ratioY->GetYaxis()->SetLabelSize(0.06);
	ratioY->GetYaxis()->SetTitleSize(0.08);
	ratioY->GetYaxis()->SetTitleOffset(0.3);
	ratioY->GetXaxis()->SetLabelSize(0.1);
	ratioY->GetXaxis()->SetTitleSize(0.1);
	ratioY->SetMarkerStyle(20);
	ratioY->SetMarkerColor(kRed);

	ratioX->Draw("pe");
	ratioY->Draw("pe,same");
	TLine*line = new TLine(binLow,1,binHigh,1);
	line->SetLineColor(kBlack);
	line->Draw("same");

	TString saveName = "plots/projectionsVsHistograms_";
	saveName += binLow;
	saveName += "to";
	saveName += binHigh;
	saveName += ".png";
	canvas1->SaveAs(saveName);
}

void integrate()
{
	TH1::SetDefaultSumw2();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	//Define functions from which to fill histograms
	TF1 *trueDist = new TF1("trueDist", trueDistribution, -10, 150, 5);
	trueDist->SetParameters(100, 20, 10, 60, 5);
	TF1 *observedDist = new TF1("observedDist", observedDistribution, -10, 150, 6);
	observedDist->SetParameters(100, 20, 10, 60, 5, 3);

	//Histogram parameters
	const double xmin = 0;
	const double xmax = 100;
	const int nbins = xmax-xmin;

	
	//Define histograms
	TH1F *histTrue = new TH1F("histTrue", "histTrue", nbins, xmin, xmax);
	TH1F *histObs  = new TH1F("histObs", "histObs", nbins, xmin, xmax);
	
	//Fill histograms
	for(int i=1; i<=nbins; i++){
		double binLowEdge = histTrue->GetXaxis()->GetBinLowEdge(i);
		double binHighEdge = histTrue->GetXaxis()->GetBinUpEdge(i);
		float yieldTrue = trueDist->Integral(binLowEdge,binHighEdge);
		float yieldObs = observedDist->Integral(binLowEdge,binHighEdge);

		histTrue->SetBinContent(i,yieldTrue);
		histObs->SetBinContent(i,yieldObs);
	}

	// 2D histogram: migration matrix
	TF2 *integrand2DFunc = new TF2("integrand2DFunc", integrand2D, 0, 100, 0, 100, 6);
	integrand2DFunc->SetParameters(100, 20, 10, 60, 5, 3);

	TH2F *migrationHist = new TH2F("migrationHist", "migrationHist",nbins,xmin,xmax,nbins, 
				       xmin,xmax);

	double xlow,xhi,ylow,yhi,yield;
	for(int iBin = 0; iBin <= nbins+1; iBin++){//loop over columns
		for(int jBin = 0; jBin <= nbins+1; jBin++){//looop over rows
			if(iBin==0) xlow = -10;
			else xlow = migrationHist->GetXaxis()->GetBinLowEdge(iBin);

			if(iBin==nbins+1) xhi = 150;
			else xhi  = migrationHist->GetXaxis()->GetBinUpEdge (iBin);

			if(jBin==0) ylow = -10;
			else ylow = migrationHist->GetYaxis()->GetBinLowEdge(jBin);

			if(jBin==nbins+1) yhi = 150;
			else  yhi  = migrationHist->GetYaxis()->GetBinUpEdge (jBin);

			double yield = integrand2DFunc->Integral(xlow, xhi, ylow, yhi);
			migrationHist->SetBinContent(iBin, jBin, yield);
		}//end loop over rows
	}//end loop over columns

	//Draw all histograms alongside migration matrix projections
	//The projections should match the true and observed histograms
	DrawProjections(histTrue,histObs,migrationHist);
}

