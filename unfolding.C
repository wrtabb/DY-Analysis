#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

const double massbins[44] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106, 
  110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380, 440, 510, 600, 700, 830, 
  1000, 1500, 3000};
const int nMassBins = 43;
const int nCanvas = 8;
const TString matrixFileName = "./plots/plotsDY.root";
const TString dataFileName = "./plots/dataVsMC.root";
void unfolding()
{
  //gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  //Getting Data and Backgrounds from root files
  TFile*fData = new TFile(dataFileName); 
  TFile*fMigrationMatrix = new TFile(matrixFileName);
  TH1F*hEW = (TH1F*)fData->Get("hEWInvMass");
  TH1F*hFakes = (TH1F*)fData->Get("hFakesInvMass");
  TH1F*hTops = (TH1F*)fData->Get("hTopsInvMass");
  TH1F*hMC = (TH1F*)fData->Get("hMCInvMass");
  TH1F*hBackground = (TH1F*)hEW->Clone("hBackground");
  hBackground->Add(hFakes);
  hBackground->Add(hTops);
  TH1F*hData = (TH1F*)fData->Get("hDataInvMass");
  hData->SetLineColor(kBlack);
  hData->SetMarkerColor(kBlack);
  hData->SetMarkerStyle(20);
  hData->Add(hBackground,-1);
  hMC->SetLineColor(kRed+2);
  hMC->SetMarkerStyle(21);
  hMC->SetMarkerColor(kRed+2);
  TH1F*hMCpostFSR = (TH1F*)fMigrationMatrix->Get("hHLTGenDielectronInvMass");
  hMCpostFSR->SetLineColor(kBlue+2);
  hMCpostFSR->SetMarkerColor(kBlue+2);
  hMCpostFSR->SetMarkerStyle(21);
  int maxBin = hFakes->GetMaximumBin();
      double x = hFakes->GetXaxis()->GetBinCenter(maxBin);
      for(int k=1;k<x+1;k++){
        if(hFakes->GetBinContent(k) < 0) hFakes->SetBinContent(k,0);
      }
   
  //Getting histograms of migration matrices from files
  TH2F*migMatrixGENisHardvsGENFS = (TH2F*)fMigrationMatrix->Get("migMatrixGENisHardvsGENFS");
  TH2F*migMatrixGENFSvsReco = (TH2F*)fMigrationMatrix->Get("migMatrixGENFSvsReco");
  TH2F*migMatrixGENisHardvsReco = (TH2F*)fMigrationMatrix->Get("migMatrixGENisHardvsReco");
  TH1F*hGenFS = (TH1F*)fMigrationMatrix->Get("hHLTGenDielectronInvMass");  
  hGenFS->SetLineColor(kBlack);
  TH1F*hReco = (TH1F*)fMigrationMatrix->Get("hRecoInvMass");
  hReco->SetTitle("Gen-Level Final State vs. Reconstructed");
  hReco->SetMarkerStyle(20);
  hReco->SetMarkerColor(kRed);
  hReco->SetLineColor(kRed);
  
  //Defininig migration matrix objects
  TMatrixD migrationGENvsGEN(nMassBins,nMassBins);
  TMatrixD migrationFSvsReco(nMassBins,nMassBins);
  TMatrixD migrationHardvsReco(nMassBins,nMassBins); 
 
  for(int i=1; i<=nMassBins; i++)
    {
      for(int j=1;j<=nMassBins;j++)
	{
	  migrationGENvsGEN(i-1,j-1) = migMatrixGENisHardvsGENFS->GetBinContent(i,j);
	  migrationFSvsReco(i-1,j-1) = migMatrixGENFSvsReco->GetBinContent(i,j);
	  migrationHardvsReco(i-1,j-1) = migMatrixGENisHardvsReco->GetBinContent(i,j);
	}
    }

  //Defining response matrix objects (normalized columns from migration matrix)
  TMatrixD responseGENvsGEN(nMassBins,nMassBins);
  TMatrixD responseFSvsReco(nMassBins,nMassBins);
  TMatrixD responseHardvsReco(nMassBins,nMassBins);

  TVectorD yieldsMeasured(nMassBins);
  TVectorD yieldsTrue(nMassBins);
  for(int i=1; i<=nMassBins; i++)
    {
      yieldsMeasured(i-1) = hData->GetBinContent(i);
    }
  
  for(int i=0; i<nMassBins; i++)
    {
      float sumGENvsGEN=0;
      float sumFSvsReco=0;
      float sumHardvsReco=0;
      for(int j=0; j<nMassBins; j++)
	{
	  sumGENvsGEN += migrationGENvsGEN(i,j);
	  sumFSvsReco += migrationFSvsReco(i,j);
	  sumHardvsReco += migrationHardvsReco(i,j);
	}      
      for(int j=0; j<nMassBins; j++)
	{
	  if(sumGENvsGEN!=0) responseGENvsGEN(i,j) = migrationGENvsGEN(i,j)/sumGENvsGEN;
	  else
	    {
	      responseGENvsGEN = 0;
	      cout << "GENvsGEN not properly normalized!" << endl;
	    }
	  if(sumFSvsReco!=0) responseFSvsReco(i,j) = migrationFSvsReco(i,j)/sumFSvsReco;
	  else
	    {
	      responseFSvsReco = 0;
	      cout << "FSvsReco not properly normalized!" << endl;
	    }
	  if(sumHardvsReco!=0) responseHardvsReco(i,j) = migrationHardvsReco(i,j)/sumHardvsReco;
	  else
	    {
	      responseHardvsReco = 0;
	      cout << "HardvsReco not properly normalized!" << endl;
	    }
	}
    }  

  //Defining unfolding matrices (the inverse of response matrices)
  TMatrixD unfoldingGENvsGEN = responseGENvsGEN;
  TMatrixD unfoldingFSvsReco = responseFSvsReco;
  TMatrixD unfoldingHardvsReco = responseHardvsReco;
  Double_t det;
  unfoldingGENvsGEN.Invert(&det);
  unfoldingFSvsReco.Invert(&det);
  unfoldingHardvsReco.Invert(&det);
  const float axisLabelOffset = 1.65;
  const float axisLabelSize = 0.03;
  const int nMatrices = 6;
  //Placing all matrices in histograms for plotting
  TH2F*hresponseGENvsGEN = (TH2F*)migMatrixGENisHardvsGENFS->Clone("hresponseGENvsGEN");
  hresponseGENvsGEN->
    SetTitle("Gen-Level Final State vs. Gen-Level Hard Process Response Matrix");
  hresponseGENvsGEN->GetYaxis()->SetTitleOffset(axisLabelOffset);
  hresponseGENvsGEN->GetXaxis()->SetTitleOffset(axisLabelOffset);
  hresponseGENvsGEN->GetYaxis()->SetTitleSize(axisLabelSize);
  hresponseGENvsGEN->GetXaxis()->SetTitleSize(axisLabelSize);
  TH2F*hresponseFSvsReco = (TH2F*)migMatrixGENFSvsReco->Clone("hresponseFSvsReco");
  hresponseFSvsReco->SetTitle("Reconstructed vs. Gen-Level Final State Response Matrix");
  hresponseFSvsReco->GetYaxis()->SetTitleOffset(axisLabelOffset);
  hresponseFSvsReco->GetXaxis()->SetTitleOffset(axisLabelOffset);
  hresponseFSvsReco->GetYaxis()->SetTitleSize(axisLabelSize);
  hresponseFSvsReco->GetXaxis()->SetTitleSize(axisLabelSize);
  TH2F*hresponseHardvsReco = (TH2F*)migMatrixGENisHardvsReco->Clone("hresponseGENvsGEN");
  hresponseHardvsReco->SetTitle("Reconstructed vs. Gen-Level Hard Process Response Matrix");
  hresponseHardvsReco->GetYaxis()->SetTitleOffset(axisLabelOffset);
  hresponseHardvsReco->GetXaxis()->SetTitleOffset(axisLabelOffset);
  hresponseHardvsReco->GetYaxis()->SetTitleSize(axisLabelSize);
  hresponseHardvsReco->GetXaxis()->SetTitleSize(axisLabelSize);
  TH2F*hunfoldingGENvsGEN = (TH2F*)migMatrixGENisHardvsGENFS->Clone("hunfoldingGENvsGEN");
  hunfoldingGENvsGEN->
    SetTitle("Gen-Level Final State vs. Gen-Level Hard Process Unfolding Matrix");
  hunfoldingGENvsGEN->GetYaxis()->SetTitleOffset(axisLabelOffset);
  hunfoldingGENvsGEN->GetXaxis()->SetTitleOffset(axisLabelOffset);
  hunfoldingGENvsGEN->GetYaxis()->SetTitleSize(axisLabelSize);
  hunfoldingGENvsGEN->GetXaxis()->SetTitleSize(axisLabelSize);
  TH2F*hunfoldingFSvsReco = (TH2F*)migMatrixGENFSvsReco->Clone("hunfoldingFSvsReco");
  hunfoldingFSvsReco->SetTitle("Reconstructed vs. Gen-Level Final State Unfolding Matrix");
  hunfoldingFSvsReco->GetYaxis()->SetTitleOffset(axisLabelOffset);
  hunfoldingFSvsReco->GetXaxis()->SetTitleOffset(axisLabelOffset);
  hunfoldingFSvsReco->GetYaxis()->SetTitleSize(axisLabelSize);
  hunfoldingFSvsReco->GetXaxis()->SetTitleSize(axisLabelSize);
  TH2F*hunfoldingHardvsReco = (TH2F*)migMatrixGENisHardvsReco->Clone("hunfoldingHardvsReco");
  hunfoldingHardvsReco->SetTitle("Reconstructed vs. Gen-Level Hard Process Unfolding Matrix");
  hunfoldingHardvsReco->GetYaxis()->SetTitleOffset(axisLabelOffset);
  hunfoldingHardvsReco->GetXaxis()->SetTitleOffset(axisLabelOffset);
  hunfoldingHardvsReco->GetYaxis()->SetTitleSize(axisLabelSize);
  hunfoldingHardvsReco->GetXaxis()->SetTitleSize(axisLabelSize);

  for(int i=0; i<nMassBins; i++)
    {
      for(int j=0; j<nMassBins; j++)
	{
	  hresponseGENvsGEN->SetBinContent(i+1,j+1,responseGENvsGEN(i,j));
	  hresponseFSvsReco->SetBinContent(i+1,j+1,responseFSvsReco(i,j));
	  hresponseHardvsReco->SetBinContent(i+1,j+1,responseHardvsReco(i,j));
	  hunfoldingGENvsGEN->SetBinContent(i+1,j+1,unfoldingGENvsGEN(i,j));
	  hunfoldingFSvsReco->SetBinContent(i+1,j+1,unfoldingFSvsReco(i,j));
	  hunfoldingHardvsReco->SetBinContent(i+1,j+1,unfoldingHardvsReco(i,j));	  
	}      
    }

  TVectorD yieldsUnfolded = (unfoldingFSvsReco.T())*yieldsMeasured;
  TH1F *hMassUnfolded = new TH1F("hMassUnfolded","",nMassBins,massbins);
  hMassUnfolded->SetTitle("Gen-Level Final State vs. Unfolded");
  hMassUnfolded->GetXaxis()->SetTitle("m_{ee} [GeV]");
  //hMassUnfolded->SetLineColor(kRed);
  hMassUnfolded->SetMarkerColor(kRed+2);
  hMassUnfolded->SetLineColor(kRed+2);
  hMassUnfolded->SetMarkerStyle(20);
  hMassUnfolded->GetXaxis()->SetNoExponent();
  hMassUnfolded->GetXaxis()->SetMoreLogLabels();
  for(int i=0; i<nMassBins; i++)
    {
      hMassUnfolded->SetBinContent(i+1,yieldsUnfolded(i));
    }
  /* //plotting slices of a response matrix
  TCanvas*canvas2[nMassBins];  
  TCanvas*canvas3 = new TCanvas("canvas3","",10,10,900,700);
  TLatex*label = new TLatex();
  label->SetTextFont(42);
  label->SetTextSize(0.02);
  const int firstbin = 1;
  const int lastbin = 44;
  TH2F*hprojMat[nMassBins];
  canvas3->SetLogx();
  canvas3->SetLogy();
  hresponseFSvsReco->Draw("colz");
  TH1D*proj[nMassBins];
  for(int i=0;i<nMassBins;i++)
    {
      TString canvasName2 = "canvas2_";
      TString projName = "proj";
      TString cloneName = "clone";
      TString labelText = "Bin Range: ";
      labelText+=massbins[i];
      labelText+="to";
      labelText+=massbins[i+1];
      
      cloneName+=i;
      projName+=i;
      canvasName2+=i;
      hprojMat[i] = (TH2F*)hresponseFSvsReco->Clone(cloneName); 
      proj[i] = hprojMat[i]->ProjectionY(Form("bin%d",i+1),i+1,i+2);
      proj[i]->GetXaxis()->SetTitleOffset(1);
      canvas2[i] = new TCanvas(canvasName2,"",10,10,900,700);
      canvas2[i]->SetLogx();
      canvas2[i]->SetLogy(); 
      canvas2[i]->cd();
      proj[i]->Draw("hist");
      label->DrawLatex(0.6,0.88,labelText);
      cout << labelText << endl;
    }
  */
  //Plotting and saving histograms
  TCanvas*canvas[nCanvas];
  TString canvasName = "canvas";
  for(int i=0;i<nCanvas;i++)
    {
      canvasName+=i;
      canvas[i] = new TCanvas(canvasName,"",10,10,1200,1200);
      canvas[i]->SetLogx();
      canvas[i]->SetLogy();
    }
  
  canvas[0]->cd();
  hresponseGENvsGEN->Draw("colz");
  canvas[1]->cd();
  hresponseFSvsReco->Draw("colz");
  canvas[2]->cd();
  hresponseHardvsReco->Draw("colz");
  canvas[3]->cd();
  hunfoldingGENvsGEN->Draw("colz");
  canvas[4]->cd();
  hunfoldingFSvsReco->Draw("colz");
  canvas[5]->cd();
  hunfoldingHardvsReco->Draw("colz");
  canvas[6]->cd();
  canvas[6]->SetGrid();
  hData->Draw("p");
  hMC->Draw("p,same");
  TLegend*legend1 = new TLegend(0.55,0.9,0.9,0.75);
  legend1->SetTextSize(0.02);
  legend1->AddEntry(hData,"Data with MC Background Subtracted");
  legend1->AddEntry(hMC,"MC Signal Only");
  legend1->Draw("same");
  canvas[7]->cd();
  canvas[7]->SetGrid();
  hMassUnfolded->Draw("pe");
  hMCpostFSR->Draw("pe,same");
  TLegend*legend2 = new TLegend(0.65,0.9,0.9,0.75);
  legend2->SetTextSize(0.02);
  legend2->AddEntry(hMassUnfolded,"Unfolded data");
  legend2->AddEntry(hMCpostFSR,"MC post-FSR");
  legend2->Draw("same");
    
  
  const TString canvasSaveName[nCanvas] = 
    {
      "./plots/unfolding/responseGENvsGEN.png",
      "./plots/unfolding/responseFSvsReco.png",
      "./plots/unfolding/responseHardvsReco.png",
      "./plots/unfolding/unfoldingGENvsGEN.png",
      "./plots/unfolding/unfoldingFSvsReco.png",
      "./plots/unfolding/unfoldingHardvsReco.png",
      "./plots/unfolding/signalDataVsMC.png",
      "./plots/unfolding/unfoldedDataVsMC.png"
      
    };
   
  TFile *rootFile = new TFile("./plots/unfoldingMatrices.root","RECREATE");
  rootFile->cd();
  hresponseGENvsGEN->Write();
  hresponseFSvsReco->Write();
  hresponseHardvsReco->Write();  
  hunfoldingGENvsGEN->Write();
  hunfoldingFSvsReco->Write();
  hunfoldingHardvsReco->Write();
  hMassUnfolded->Write();
  hGenFS->Write();
  hReco->Write();
  for(int i=0;i<nCanvas;i++)
    {
      canvas[i]->Write();
      canvas[i]->SaveAs(canvasSaveName[i]);
    }
  rootFile->Write();
  rootFile->Close();
 
}//end main
