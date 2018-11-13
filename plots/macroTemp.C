


void macroTemp()
{
  TFile*file=new TFile("dataUnfoldTest.root");
  
  TH2F*hMatrix = (TH2F*)file->Get("hMatrix");
  TH1F*hData = (TH1F*)file->Get("hData");
  TH1F*hBackground = (TH1F*)file->Get("hMCBackground");
  TH1D*hProjY = (TH1D*)hMatrix->ProjectionY();
  TH1D*hProjX = (TH1D*)hMatrix->ProjectionX();
  
  hData->SetMarkerStyle(25);
  hProjY->SetMarkerStyle(20);
  hProjX->SetMarkerStyle(20);
  hProjY->SetMarkerColor(kRed);
  hProjX->SetMarkerColor(kBlue);
  const int maxBin = hBackground->GetMaximumBin();  
  const int nBins = hBackground->GetXaxis()->GetBinCenter(maxBin);
  for(int i=0;i<nBins;i++){
    if(hBackground->GetBinContent(i+1)<0) hBackground->SetBinContent(i+1,0);
  }
  hData->Add(hBackground,-1.0);
  TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
  legend->SetTextSize(0.02);
  legend->AddEntry(hData,"Data");
  legend->AddEntry(hProjY,"Y-Projection, Reco");
  legend->AddEntry(hProjX,"X-Projection, Gen");

  TH1F*ratio = (TH1F*)hData->Clone("ratio");
  ratio->Divide(hProjY);
  ratio->GetYaxis()->SetRangeUser(0.7,1.3);
  const float padmargins = 0.03;
  TCanvas*canvas1 = new TCanvas("canvas1","",10,10,1200,1200);
  TPad*pad1 = new TPad("","",0,0.3,1.0,1.0);
  pad1->SetBottomMargin(padmargins);
  pad1->SetGrid();
  pad1->SetLogy();
  pad1->SetLogx();
  pad1->SetTicks(1,1);
  pad1->Draw();
  pad1->cd();
  TLine*line = new TLine(15,1,3000,1);
  line->SetLineColor(kRed);
  hData->SetLabelSize(0);
  hData->SetTitleSize(0);
  hData->Draw("hist");
  hProjY->Draw("PE,same");
  hProjX->Draw("PE,same");
  legend->Draw("same");
  canvas1->Update();

  canvas1->cd();
  TPad*pad2 = new TPad("","",0,0.05,1,0.3);
  pad2->SetLogx();  
  pad2->SetTopMargin(padmargins);
  pad2->SetBottomMargin(0.2);
  pad2->SetGrid();
  pad2->SetTicks(1,1);
  pad2->Draw();
  pad2->cd();
  ratio->GetYaxis()->SetLabelSize(0.06);
  ratio->GetYaxis()->SetTitleSize(0.08);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->SetTitle("Unfolded/Gen");
  ratio->GetXaxis()->SetLabelSize(0.1);
  ratio->GetXaxis()->SetTitleSize(0.1);
  ratio->GetXaxis()->SetNoExponent();
  ratio->GetXaxis()->SetMoreLogLabels();
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerColor(kBlack);
  ratio->Draw("PE");
  line->Draw("same");
  canvas1->Update();

  
  canvas1->SaveAs("dataRecoComparison.png");
}

