
void counter(Long64_t i, Long64_t N);

void covTest()
{
 TH1::SetDefaultSumw2();
 gStyle->SetOptStat(0);
 Long64_t nPoints = 1000;
 Long64_t nSampleStart = 1000;
 double mean1 = 200;
 double mean2 = 500; 
 int stepSize = nSampleStart;
 const bool corrTest = true;//do correlation matrix instead of covariance
 
 TF1*fitFunc = new TF1("fitFunc","[0]*x^[1]",0,nPoints*nSampleStart);
  fitFunc->SetParameter(0,10);
  fitFunc->SetParameter(1,0.5);

 int nSamples = nSampleStart;
 double x,y,xy,xAvg,yAvg,xyAvg,xAvg2,yAvg2,xyAvg2,xStd,yStd,xyStd,covXY,corrXY;
 double xSum,ySum,xySum,xSum2,ySum2,xySum2,binError;
 TH1D*hGraph = new TH1D("hGraph","",nPoints,0,nPoints*nSampleStart);
 TRandom3 gen;

 for(Long64_t j=1;j<nPoints+1;j++){
  xSum = 0.0;
  ySum = 0.0;
  xySum = 0.0;
  xSum2 = 0.0;
  ySum2 = 0.0;
  xySum2 = 0.0;
  counter(j,nPoints);

  for(Long64_t i=0;i<nSamples;i++){
   x = gen.Poisson(mean1);
   y = gen.Poisson(mean2);
   xy = x*y;

   xSum += x;
   ySum += y;
   xySum += xy;

   xSum2 += x*x;
   ySum2 += y*y;
   xySum2 += xy*xy;
  }

  xAvg = xSum/(1.0*nSamples);
  yAvg = ySum/(1.0*nSamples);
  xyAvg = xySum/(1.0*nSamples);

  xAvg2 = xSum2/(1.0*nSamples);
  yAvg2 = ySum2/(1.0*nSamples);
  xyAvg2 = xySum2/(1.0*nSamples);

  xStd = sqrt(xAvg2-xAvg*xAvg);
  yStd = sqrt(yAvg2-yAvg*yAvg);
  xyStd = sqrt(xyAvg2-xyAvg*xyAvg);

  //must set error bars in order to do fitting
  binError = 1/sqrt(nSamples);
  covXY = abs(xyAvg-xAvg*yAvg);
  corrXY = covXY/((xStd*yStd));
  hGraph->SetBinError(j,binError);
  if(corrTest) hGraph->SetBinContent(j,corrXY);
  else hGraph->SetBinContent(j,covXY);
  nSamples = nSamples + stepSize;
 }

 //Draw the graph
 TCanvas*canvas = new TCanvas("canvas","",10,10,1200,1000);
 canvas->SetGrid();
 canvas->SetLogy();
 hGraph->GetXaxis()->SetTitle("Number of samples");
 if(corrTest) hGraph->GetYaxis()->SetTitle("correlation");
 else hGraph->GetYaxis()->SetTitle("covariance");
 hGraph->SetMarkerStyle(20);
 hGraph->Fit("fitFunc");
 hGraph->Draw("P");

/*
 TString saveName = "/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/covariance/";
 if(corrTest) saveName += "testCorrelations";
 else saveName += "testCovariance";
 saveName += "_Mean";
 saveName += mean1;
 saveName += "_";
 saveName += mean2;
 saveName += "_Range";
 saveName += nSampleStart*nPoints;
 saveName += ".png";
 canvas->SaveAs(saveName);
*/
}

void counter(Long64_t i, Long64_t N)
{
 int P = 100*(i)/(N);  
 TTimeStamp eventTimeStamp;
 if(i%(N/100)==0)
  cout << "covTest.C " << "[Time: " << eventTimeStamp.AsString("s") << "] " << P << "%" << endl;
 return;
}

