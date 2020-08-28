#include "VariableList.h"

void DrawPlots(TH1D*hTrue,TH1D*hReco,TH2D*hMatrix);

void makeToyModels()
{ 
 TH1::SetDefaultSumw2();
 gStyle->SetOptStat(0);
 gStyle->SetPalette(1);

 TH1D*hReco; 
 TH1D*hRecoClosure; 
 TH1D*hTrue;
 TH2D*hMatrix;
 
 const double mean = 20;
 const double sigma = 5;
 const double mean_smeared = 0;
 const double sigma_smeared = 1.0;
 const Long64_t nEntries = 1e7;

 TF1*func = new TF1("func","1/(x+1)+gaus(0)",0,50);
 func->SetParameters(1.0,mean,sigma);
 TRandom3 gen1;
 TRandom3 gen2;
 gen1.SetSeed(82);
 gen2.SetSeed(1981);

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
  cout << k+1 << " out of " << nDistributions << endl;
  const int nBinsReco = (recoBinning.at(k)).size()-1;
  double binningReco[nBinsReco+1];
  for(int j=0;j<nBinsReco+1;j++){
   binningReco[j]=(recoBinning.at(k)).at(j);
  }
  TString hRecoTitle = "hReco";
  hRecoTitle += k;
  TString hClosureTitle = "hRecoClosure";
  hClosureTitle += k;
  TString hTrueTitle = "hTrue";
  hTrueTitle += k;
  TString hMatrixTitle = "hMatrix";
  hMatrixTitle += k;

  hReco = new TH1D(hRecoTitle,"",nBinsReco,binningReco);
  hRecoClosure = new TH1D(hClosureTitle,"",nBinsReco,binningReco);

  hTrue = new TH1D(hTrueTitle,"",nBinsTrue,binningTrue);
  hMatrix = new TH2D(hMatrixTitle,"",nBinsTrue,binningTrue,nBinsReco,binningReco);
  double peak,peakReco,peak_smeared,peakReco_smeared;
  double slope,slope_smeared;
  for(int i=0;i<nEntries;i++){
   peak = func->GetRandom();
   peakReco = func->GetRandom();
   peak_smeared = peak+gen1.Gaus(mean_smeared,sigma_smeared);
   peakReco_smeared = peakReco+gen2.Gaus(mean_smeared,sigma_smeared);
 
   hRecoClosure->Fill(peak_smeared);
   hReco->Fill(peakReco_smeared);
   hTrue->Fill(peak);
   hMatrix->Fill(peak,peak_smeared);  

  }
  saveFile->cd();
  hReco->Write();
  hRecoClosure->Write();
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


