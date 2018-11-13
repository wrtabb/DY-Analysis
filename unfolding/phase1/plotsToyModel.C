
void counter(Long64_t i, Long64_t N);
const int nLogBins = 43;
const int nLogBins2 = 2*nLogBins;
const float  massbins[] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 
  86, 91,96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,200, 220, 243, 273, 
  320, 380, 440, 510, 600, 700, 830, 1000, 1500,3000};
const float massbins2[] = {15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,
                            57.5,60,62,64,66,68,70,72,74,76,78.5,81,83.5,86,88.5,91,93.5,96,
                            98.5,101,103.5,106,108,110,112.5,115,117.5,120,123,126,129.5,133,
                            137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,210,220,
                            231.5,243,258,273,296.5,320,350,380,410,440,475,510,555,600,650,
                            700,765,830,915,1000,1250,1500,2250,3000};
const int nBins = 40;
const int nBins2 = 80;
const float binLow = 20;
const float binHigh = 200;

const double pi = TMath::Pi();
const TString toyModelName = "toyData.root";
const TString toyTreeName = "toyData";
const TString histSaveName = "toyUnfold.root";

void plotsToyModel()
{
  //Histogram style
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  TH1::SetDefaultSumw2();
  //Define variables and branches from tree
  float massTrue,massMeasured;
  TBranch*b_massTrue;
  TBranch*b_massMeasured;

  //Define file and tree for loading data
  TFile*fToyData = new TFile(toyModelName);
  TTree*tree=(TTree*)fToyData->Get(toyTreeName);
  tree->SetBranchAddress("massTrue",&massTrue,&b_massTrue);
  tree->SetBranchAddress("massMeasured",&massMeasured,&b_massMeasured);

  //Define histograms
  TH1D*hMassTrue = new TH1D("hMassTrue","",nBins,binLow,binHigh);
  hMassTrue->SetLineColor(kBlue+2);
  hMassTrue->SetLineWidth(2);
  TH1D*hMassMeasured = new TH1D("hMassMeasured","",nBins2,binLow,binHigh);
  hMassMeasured->SetMarkerColor(kRed+2);
  hMassMeasured->SetMarkerStyle(21);
  hMassMeasured->SetMarkerSize(1);
  hMassMeasured->SetLineColor(kRed+2);
  hMassMeasured->SetLineWidth(2);
  TH2D*hMatrix = new TH2D("hMatrix","",nBins2,binLow,binHigh,nBins,binLow,binHigh);

  //Begin event loop
  Long64_t nEntries = tree->GetEntries();
  for(Long64_t i=0;i<nEntries;i++){
    counter(i,nEntries);
    tree->GetEntry(i);
    hMassTrue->Fill(massTrue);
    hMassMeasured->Fill(massMeasured);
    hMatrix->Fill(massTrue,massMeasured);
  }
  TH1D*hMassMeasuredRebin = (TH1D*)hMassMeasured->Clone();
  hMassMeasuredRebin->Rebin(2);
  TCanvas*canvas = new TCanvas("canvas","",10,10,1000,1000);
  //canvas->SetLogy();
  //canvas->SetLogx();
  canvas->SetGrid();
  TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
  legend->SetTextSize(0.02);
  legend->AddEntry(hMassTrue,"True Mass");
  legend->AddEntry(hMassMeasuredRebin,"Measured Mass");
  //hMassTrue->GetXaxis()->SetMoreLogLabels();
  //hMassTrue->GetXaxis()->SetNoExponent();
  hMassTrue->GetXaxis()->SetTitle("mass [GeV]");
  //hMassMeasured->GetXaxis()->SetMoreLogLabels();
  //hMassMeasured->GetXaxis()->SetNoExponent();
  //hMassMeasuredRebin->GetXaxis()->SetMoreLogLabels();
  //hMassMeasuredRebin->GetXaxis()->SetNoExponent();
  hMassTrue->Draw("hist");
  hMassMeasuredRebin->Draw("PE,same");
  legend->Draw("same");

  TCanvas*cMatrix = new TCanvas("cMatrix","",10,10,1000,1000);
  //cMatrix->SetLogy();
  //cMatrix->SetLogx();
  hMatrix->GetYaxis()->SetTitleOffset(1.5);
  //hMatrix->GetXaxis()->SetMoreLogLabels();
  //hMatrix->GetYaxis()->SetMoreLogLabels();
  //hMatrix->GetXaxis()->SetNoExponent();
  //hMatrix->GetYaxis()->SetNoExponent();
  hMatrix->GetXaxis()->SetTitle("True Mass [GeV]");
  hMatrix->GetYaxis()->SetTitle("Measured Mass [GeV]");
  hMatrix->Draw("colz");
  
  canvas->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase1Plots/step1TruthandReco_RecoInMassRange.png");
  cMatrix->SaveAs("/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase1Plots/step1MigrationMatrix_RecoInMassRange.png");
  cout << "Number of entries processed: " << nEntries << endl;
  TFile*file = new TFile(histSaveName,"recreate");
  file->cd();
  hMatrix->Write();
  hMassMeasured->Write();
  hMassTrue->Write();
  file->Write();
  file->Close();
}

//Counter for tracking program progress
void counter(Long64_t i, Long64_t N)
{
  int P = 100*(i)/(N);
  TTimeStamp eventTimeStamp;
  if(i%(N/100)==0) {
    cout << "plotsToyModel.C " << "[Time: " << eventTimeStamp.AsString("s") << "] " << P
    << "%" << endl;
  }
return;
}
   
