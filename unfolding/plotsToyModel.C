
void counter(Long64_t i, Long64_t N);
const int nBins = 43;
const double massbins[] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86,
                           91,96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171,
                           185,200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000,
                           1500,3000};
const double pi = TMath::Pi();
const TString toyModelName = "toyData.root";
const TString toyTreeName = "toyData";


void plotsToyModel()
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  float massTrue,massMeasured;
  TBranch*b_massTrue;
  TBranch*b_massMeasured;

  TFile*fToyData = new TFile(toyModelName);
  TTree*tree=(TTree*)fToyData->Get(toyTreeName);
  tree->SetBranchAddress("massTrue",&massTrue,&b_massTrue);
  tree->SetBranchAddress("massMeasured",&massMeasured,&b_massMeasured);

  TH1F*hMassTrue = new TH1F("hMassTrue","",nBins,massbins);
  hMassTrue->SetLineColor(kBlue+2);
  hMassTrue->SetLineWidth(2);
  TH1F*hMassMeasured = new TH1F("hMassMeasured","",nBins,massbins);
  hMassMeasured->SetLineColor(kRed+2);
  hMassMeasured->SetLineWidth(2);
  
  TH2F*hMatrix = new TH2F("hMatrix","",nBins,massbins,nBins,massbins);
  Long64_t nEntries = tree->GetEntries();
  for(Long64_t i=0;i<nEntries;i++){
    counter(i,nEntries);
    tree->GetEntry(i);
    hMassTrue->Fill(massTrue);
    hMassMeasured->Fill(massMeasured);
    hMatrix->Fill(massTrue,massMeasured);
  }

  TCanvas*canvas = new TCanvas("canvas","",10,10,1000,1000);
  canvas->SetLogy();
  canvas->SetLogx();
  canvas->SetGrid();
  TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
  legend->SetTextSize(0.02);
  legend->AddEntry(hMassTrue,"True Mass");
  legend->AddEntry(hMassMeasured,"Measured Mass");
  hMassTrue->GetXaxis()->SetMoreLogLabels();
  hMassTrue->GetXaxis()->SetNoExponent();
  hMassMeasured->GetXaxis()->SetMoreLogLabels();
  hMassMeasured->GetXaxis()->SetNoExponent();
  hMassTrue->Draw("hist");
  hMassMeasured->Draw("hist,same");
  legend->Draw("same");

  TCanvas*cMatrix = new TCanvas("cMatrix","",10,10,1000,1000);
  cMatrix->SetLogy();
  cMatrix->SetLogx();
  hMatrix->GetYaxis()->SetTitleOffset(1.5);
  hMatrix->GetXaxis()->SetMoreLogLabels();
  hMatrix->GetYaxis()->SetMoreLogLabels();
  hMatrix->GetXaxis()->SetNoExponent();
  hMatrix->GetYaxis()->SetNoExponent();
  hMatrix->GetXaxis()->SetTitle("True Mass [GeV]");
  hMatrix->GetYaxis()->SetTitle("Measured Mass [GeV]");
  hMatrix->Draw("colz");
  
  canvas->SaveAs("./plots/step1TruthandReco.png");
  cMatrix->SaveAs("./plots/step1MigrationMatrix.png");
  cout << "Number of entries processed: " << nEntries << endl;
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
   
