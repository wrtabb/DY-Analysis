
std::vector<double> chi2Values = {
 6.87245,
 22.6519,
 16.7891,
 26.9312,
 12.0057,
 39.2215,
 51.7758
};
std::vector<TString> binLabels = {
 "double bins all",
 "double bins 25-50",
 "double bins 0-15",
 "double bins 0-15, 25-50",
 "double bins 15-25",
 "split highest bin",
 "split lowest bin"
};
double GetUpperBound(std::vector<double> vec);

void plotChi2()
{
 gStyle->SetOptStat(0);
 int nBins = chi2Values.size();
 TH1D*hChi2 = new TH1D("Chi2","",nBins,0,nBins); 
 hChi2->SetMarkerStyle(20);
 hChi2->SetMarkerColor(kBlue+2);
 hChi2->SetLineColor(kBlue+2);
 hChi2->GetYaxis()->SetTitle("#chi^{2}/ndf");
 double upperBound = GetUpperBound(chi2Values);
 hChi2->GetYaxis()->SetRangeUser(0,upperBound);
 for(int i=1;i<=nBins;i++){
  hChi2->SetBinContent(i,chi2Values.at(i-1));
  hChi2->GetXaxis()->SetBinLabel(i,binLabels.at(i-1));
  hChi2->SetBinError(i,sqrt(chi2Values.at(i-1)));
 }

 TCanvas*canvas = new TCanvas("canvas","",0,0,1000,1000);
 canvas->SetGrid();
 hChi2->Draw("pe"); 
 canvas->SaveAs("./plots/chi2FromBinningToy.png");
}

double GetUpperBound(std::vector<double> vec)
{
 double upperValue = 0;
 for(int i=0;i<vec.size();i++){
  if(upperValue < vec.at(i)) upperValue = vec.at(i);
 }
 double finalValue = upperValue + sqrt(upperValue)*1.1;
 return finalValue;
}
