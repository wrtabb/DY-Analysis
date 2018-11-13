



void testPlots()
{
  TFile*rootfile=new TFile("dataUnfoldTest.root","recreate");
  TFile*fileData=new TFile("tempData.root");
  TFile*fileMatrix=new TFile("plotsDY.root");

  //TH1F*hBackgroundFakes=(TH1F*)fileData->Get("hFakesInvMass");
  //TH1F*hBackgroundEW=(TH1F*)fileData->Get("hEWInvMass");
  //TH1F*hBackgroundTops=(TH1F*)fileData->Get("hTopsInvMass");
  //int maxBin = hBackgroundFakes->GetMaximumBin();
  //double x = hBackgroundFakes->GetXaxis()->GetBinCenter(maxBin);
  //for(int k=0;k<2*43;k++){
  //  if(hBackgroundFakes->GetBinContent(k+1) < 0) hBackgroundFakes->SetBinContent(k+1,0);
  //  if(hBackgroundEW->GetBinContent(k+1) < 0) hBackgroundEW->SetBinContent(k+1,0);
  //  if(hBackgroundTops->GetBinContent(k+1) < 0) hBackgroundTops->SetBinContent(k+1,0);
 // }

  //TH1F*hBackground=(TH1F*)hBackgroundFakes->Clone();
  //hBackground->Add(hBackgroundEW);
  //hBackground->Add(hBackgroundTops);
  //hBackground->SetName("hMCBackground");  
  TH1F*hBackground=(TH1F*)fileData->Get("hMCBackground"); 
  TH1F*hData=(TH1F*)fileData->Get("hData");

  TH2F*hMatrix=(TH2F*)fileMatrix->Get("migMatrixGENFSvsReco");
  hMatrix->SetName("hMatrix");
  rootfile->cd();
  hBackground->Write();
  hData->Write();
  hMatrix->Write();
  rootfile->Write();
  rootfile->Close(); 
}
