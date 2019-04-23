
void counter(Long64_t i, Long64_t N);


void covTest()
{
 gStyle->SetOptStat(0);
 Long64_t nPoints = 10000;
 Long64_t nSampleStart = 100;
 double mean1 = 100;
 double mean2 = 500; 
 int stepSize = nSampleStart;
 TF1*distX = new TF1("distX","TMath::Poisson(x,[0])",0,1000);
 distX->SetParameter(0,mean1); 
 distX->SetLineColor(kBlue);
 TF1*distY = new TF1("distY","TMath::Poisson(x,[0])",0,1000);
 distY->SetParameter(0,mean2); 
 distY->SetLineColor(kRed);
 TH1D*hGraph = new TH1D("hGraph","",nPoints,nSampleStart,stepSize*nPoints);
 int nSamples = nSampleStart;
 double x,y,xy,xAvg,yAvg,xyAvg,xStd,yStd,xyStd,covXY;

 distX->Draw();
 distY->Draw("same");
 Long64_t count = 0;
for(Long64_t j=0;j<nPoints;j++){
 double xSum = 0.0;
 double ySum = 0.0;
 double xySum = 0.0;
 double xSum2 = 0.0;
 double ySum2 = 0.0;
 double xySum2 = 0.0;
 
 for(Long64_t i=0;i<nSamples;i++){
  x = distX->GetRandom();
  y = distY->GetRandom();
  xy = x*y;

  xSum += x;
  ySum += y;
  xySum += xy;

  xSum2 += x*x;
  ySum2 += y*y;
  xySum2 += xy*xy;
 }

 xAvg = xSum/nSamples;
 yAvg = ySum/nSamples;
 xyAvg = xySum/nSamples;
 xStd = sqrt(xSum2/nSamples-xSum*xSum/(nSamples*nSamples*1.0));
 yStd = sqrt(ySum2/nSamples-ySum*ySum/(nSamples*nSamples*1.0));
 xyStd = sqrt(xySum2/nSamples-xySum*xySum/(nSamples*nSamples*1.0));
 //covXY = 0.5*(xyStd*xyStd-xStd*xStd-yStd*yStd);
 covXY = (xyAvg-xAvg*yAvg)*(xyAvg-xAvg*yAvg);
 hGraph->SetBinContent(j,covXY);
/*
 cout << "****************************************************************" << endl;
 cout << "Number of samples: " << nSamples << endl;
 cout << "Averages (x, y): " << "(" << xAvg << ", " << yAvg << ")" << endl;
 cout << "Std Devs (x, y): " << "(" << xStd << ", " << yStd << ")" << endl;
 cout << "Covariance: " << covXY << endl;
 cout << "nSamples: " << nSamples << endl;
*/
 nSamples = nSamples + stepSize;
}
TCanvas*canvas = new TCanvas("canvas","",10,10,1200,1000);
canvas->SetGrid();
canvas->SetLogy();
hGraph->GetXaxis()->SetTitle("Number of samples");
hGraph->GetYaxis()->SetTitle("Covariance^{2}");
hGraph->SetMarkerStyle(20);
hGraph->Draw("P");
TString saveName = "/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/covariance/";
saveName += "testCovariance";
saveName += nSamples;
saveName += "Samples";
saveName += "_Mean";
saveName += mean1;
saveName += "_";
saveName += mean2;
saveName += ".png";
canvas->SaveAs(saveName);
}

void counter(Long64_t i, Long64_t N)
{
 int P = 100*(i)/(N);  
 TTimeStamp eventTimeStamp;
 if(i%(N/100)==0)
  cout << "covTest.C " << "[Time: " << eventTimeStamp.AsString("s") << "] " << P << "%" << endl;
 return;
}

