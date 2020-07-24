#include "VariableList.h"

void DrawPlots(TH1D*hTrue,TH1D*hReco,TH2D*hMatrix);

void makeToyModels()
{ 
 TH1::SetDefaultSumw2();
 gStyle->SetOptStat(0);
 gStyle->SetPalette(1);

 TH1D*hReco; 
 TH1D*hTrue;
 TH2D*hMatrix;
 
 const double mean = 20;
 const double sigma = 5;
 const double mean_smeared = 0;
 const double sigma_smeared = 0.8;
 const Long64_t nEntries = 1e8;

 TF1*func = new TF1("func","1/(x+1)+gaus(0)",0,50);
 func->SetParameters(1.0,mean,sigma);
 TRandom3 gen;

 std::vector< std::vector<double> > recoBinning;
 recoBinning.push_back(binningReco0);
 recoBinning.push_back(binningReco1);
 recoBinning.push_back(binningReco2);
 recoBinning.push_back(binningReco3);
 recoBinning.push_back(binningReco4);
 recoBinning.push_back(binningReco5);
 recoBinning.push_back(binningReco6);

 int nDistributions = recoBinning.size();
 TFile*saveFile = new TFile("unfoldingHists.root","recreate");
 for(int k=0;k<nDistributions;k++){
  const int nBinsReco = (recoBinning.at(k)).size()-1;
  double binningReco[nBinsReco+1];
  for(int j=0;j<nBinsReco+1;j++){
   binningReco[j]=(recoBinning.at(k)).at(j);
  }
  TString hRecoTitle = "hReco";
  hRecoTitle += k;
  TString hTrueTitle = "hTrue";
  hTrueTitle += k;
  TString hMatrixTitle = "hMatrix";
  hMatrixTitle += k;

  hReco = new TH1D(hRecoTitle,"",nBinsReco,binningReco);
  hTrue = new TH1D(hTrueTitle,"",nBinsTrue,binningTrue);
  hMatrix = new TH2D(hMatrixTitle,"",nBinsTrue,binningTrue,nBinsReco,binningReco);
  double peak,peak_smeared;
  double slope,slope_smeared;
  for(int i=0;i<nEntries;i++){
   peak = func->GetRandom();;
   peak_smeared = peak+gen.Gaus(mean_smeared,sigma_smeared);
 
   hReco->Fill(peak_smeared);
   hTrue->Fill(peak);
   hMatrix->Fill(peak,peak_smeared);  
  }
  //---Randomize bins of hReco---//
  double binContent,vary;
  for(int i=1;i<=nBinsReco;i++){
   binContent = hReco->GetBinContent(i);
   vary = binContent*0.005;
   hReco->SetBinContent(i,binContent+gen.Gaus(0,vary));
  }
 saveFile->cd();
 hReco->Write();
 hTrue->Write();
 hMatrix->Write(); 
 }//end loop over distributions
 saveFile->Close();
 //DrawPlots(hTrue,hReco,hMatrix);
}

void DrawPlots(TH1D*hTrue,TH1D*hReco,TH2D*hMatrix)
{
 hTrue->SetFillColor(kRed+2);
 hTrue->SetLineColor(kRed+2);
 hReco->SetMarkerStyle(20);
 
 TString cName = "canvas";
 TCanvas*canvas = new TCanvas(cName,"",0,0,1000,1000);
 canvas->SetGrid();

 hTrue->Draw("hist");
 TH1D*hRecoRebin = (TH1D*)hReco->Rebin(nBinsTrue,"hRecoRebin",binningTrue);

 hRecoRebin->Draw("pesame");
 TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
 legend->SetTextSize(0.02);
 legend->AddEntry(hTrue,"true");
 legend->AddEntry(hReco,"reco");
 legend->Draw("same");

 TCanvas*canvas2 = new TCanvas("canvas2","",0,0,1000,1000);
 canvas2->SetGrid();
 canvas2->SetLogz();
 hMatrix->GetXaxis()->SetTitle("true");
 hMatrix->GetYaxis()->SetTitle("reco");
 hMatrix->Draw("colz");

 canvas->SaveAs("./plots/oneDimDistributions.png");
 canvas2->SaveAs("./plots/migrationMatrix.png");
 TFile*file = new TFile("unfoldingHists.root","recreate");
 hTrue->Write();
 hReco->Write();
 hMatrix->Write();
 file->Close();
}

