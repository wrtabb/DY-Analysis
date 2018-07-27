#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

const double massbins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 
			     106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 
			     380, 440, 510, 600, 700, 830, 1000, 1500, 3000};
const int nMassBins = 43;
const int nHistos = 3;
const int nCanvas = 6;
const TString dataFileName = "./plots/plotsDY.root";

void unfolding()
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  TFile*fMigrationMatrix = new TFile(dataFileName);

  TH2F*migMatrixGENisHardvsGENFS = (TH2F*)fMigrationMatrix->Get("migMatrixGENisHardvsGENFS");
  TH2F*migMatrixGENFSvsReco = (TH2F*)fMigrationMatrix->Get("migMatrixGENFSvsReco");
  TH2F*migMatrixGENisHardvsReco = (TH2F*)fMigrationMatrix->Get("migMatrixGENisHardvsReco");
  
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

  TMatrixD responseGENvsGEN(nMassBins,nMassBins);
  TMatrixD responseFSvsReco(nMassBins,nMassBins);
  TMatrixD responseHardvsReco(nMassBins,nMassBins);
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
      
      // Fill the response matrices
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


  // "Unfolding" matrix: inverted response matrix
  TMatrixD unfoldingGENvsGEN = responseGENvsGEN;
  TMatrixD unfoldingFSvsReco = responseFSvsReco;
  TMatrixD unfoldingHardvsReco = responseHardvsReco;

  Double_t det;
  unfoldingGENvsGEN.Invert(&det);
  unfoldingFSvsReco.Invert(&det);
  unfoldingHardvsReco.Invert(&det);
  // Repack the response and unfolding matrices into 2D histograms for plotting
  TH2F*hresponseGENvsGEN = (TH2F*)migMatrixGENisHardvsGENFS->Clone("hresponseGENvsGEN");
  hresponseGENvsGEN->SetTitle("Gen vs Gen response matrix");
  TH2F*hresponseFSvsReco = (TH2F*)migMatrixGENFSvsReco->Clone("hresponseFSvsReco");
  hresponseFSvsReco->SetTitle("FS vs Reco response matrix");
  TH2F*hresponseHardvsReco = (TH2F*)migMatrixGENisHardvsReco->Clone("hresponseGENvsGEN");
  hresponseHardvsReco->SetTitle("Hard vs Reco response matrix");

  TH2F*hunfoldingGENvsGEN = (TH2F*)migMatrixGENisHardvsGENFS->Clone("hunfoldingGENvsGEN");
  hunfoldingGENvsGEN->SetTitle("Gen vs Gen unfolding matrix");
  TH2F*hunfoldingFSvsReco = (TH2F*)migMatrixGENFSvsReco->Clone("hunfoldingFSvsReco");
  hunfoldingFSvsReco->SetTitle("FS vs Reco unfolding matrix");
  TH2F*hunfoldingHardvsReco = (TH2F*)migMatrixGENisHardvsReco->Clone("hunfoldingHardvsReco");
  hunfoldingHardvsReco->SetTitle("Hard vs Reco unfolding matrix");

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
  TCanvas*canvas[nCanvas];
  TString canvasName = "canvas";

  for(int i=0;i<nCanvas;i++)
    {
      canvasName+=i;
      canvas[i] = new TCanvas(canvasName,"",10,10,1000,700);
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

TString canvasSaveName[nCanvas] = 
  {
    "./plots/responseGENvsGEN.png",
    "./plots/responseFSvsReco.png",
    "./plots/responseHardvsReco.png",
    "./plots/unfoldingGENvsGEN.png",
    "./plots/unfoldingFSvsReco.png",
    "./plots/unfoldingHardvsReco.png"
  };

  TFile *rootFile = new TFile("./plots/unfoldingMatrices.root","RECREATE");
  rootFile->cd();
  hresponseGENvsGEN->Write();
  hresponseFSvsReco->Write();
  hresponseHardvsReco->Write();  
  hunfoldingGENvsGEN->Write();
  hunfoldingFSvsReco->Write();
  hunfoldingHardvsReco->Write();
  for(int i=0;i<nCanvas;i++)
    {
      canvas[i]->Write();
      canvas[i]->SaveAs(canvasSaveName[i]);
    }
  rootFile->Write();
  rootFile->Close();
}//end main
