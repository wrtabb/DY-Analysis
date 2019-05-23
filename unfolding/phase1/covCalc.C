#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"

const Long64_t nSamples = 10000;//number of input vectors to create
const TString fileName = "toyUnfold.root";//location of toy model distributions

void covCalc()
{
 TTimeStamp ts_start;
 cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
 TStopwatch totaltime;
 totaltime.Start();

 //initialize variables, arrays, TFile, and set histogram options
 double error,binContent,smear;
 vector <double>*inVec[nSamples];
 TFile*file = new TFile(fileName);
 TH1::SetDefaultSumw2();
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);
 gROOT->SetBatch(true);
 //Get histograms from files
 TH1D*hTrue = (TH1D*)file->Get("hTrue");
  hTrue->SetFillColor(kRed+2);
  hTrue->SetLineColor(kRed+2);
 TH1D*hBack = (TH1D*)file->Get("hBack1");
 TH2D*hMatrix = (TH2D*)file->Get("hMatrix");

 //Get reco distribution
 TH1D*hReco = (TH1D*)file->Get("hReco");
  hReco->SetMarkerStyle(20);
  hReco->SetMarkerColor(kBlack);
  hReco->SetLineColor(kBlack);

 //Array of reco distributions smeared with a gaussian
 //These will make up the input vectors for unfolding
 TRandom3 gen;
 //Begin loop over samples doing unfold on each one and placing it in an array
 for(Long64_t i=1;i<nSamples+1;i++){
  cout << i << "/" << nSamples << " (" << 100*(i)/(nSamples) << "%)" << endl;
  gRandom->SetSeed(i);
  ////////////////////////////////////////////////////////////////////////
  // Part 1: smearing each histogram to create various input histograms //
  ////////////////////////////////////////////////////////////////////////
  TString recoName = "hRecoSmeared";
  recoName += i;
  TH1D*hRecoSmeared = (TH1D*)hReco->Clone("hRecoSmeared");
  //vary each bin with a gaussian of width equal to the bin error
  for(int j=1;j<nLogBins2+1;j++){
   error = hReco->GetBinError(j);
   smear = gen.Gaus(0,error);
   binContent = hReco->GetBinContent(j)+smear;
   hRecoSmeared->SetBinContent(j,binContent);
   //Now hReco is varied
  }
  /////////////////////////////////////////////////////
  // Part 2: Doing the unfolding with TUnfoldDensity //
  /////////////////////////////////////////////////////
/*
  //Set unfolding options
  TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
  TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;
  TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
  TUnfold::EHistMap outputMap = TUnfold::kHistMapOutputHoriz;

  //Do the unfolding
  TUnfoldDensity unfold(hMatrix,outputMap,regMode,constraintMode,densityFlags);
  unfold.SetInput(hRecoSmeared);
  Int_t iBest;
  TSpline *logTauX,*logTauY;
  TGraph *lCurve;
  TGraph*bestLcurve;
  TGraph*bestLogTauLogChi2;
  Int_t nScan=30;
  Double_t tauMin = 0.0;
  Double_t tauMax = 0.0;
  iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
  Double_t t[1],x[1],y[1];
  logTauX->GetKnot(iBest,t[0],x[0]);
  logTauY->GetKnot(iBest,t[0],y[0]);
  bestLcurve=new TGraph(1,x,y);
  bestLogTauLogChi2=new TGraph(1,t,x);
  TString unfName = "unfolded";
  unfName += i;
  TH1*hUnfolded = unfold.GetOutput(unfName);
*/
  //Now that unfolding has been done, place results into an array
  //hRecoSmeared->Rebin(2);
  inVec[i-1] = new vector<double>;
  for(int j=1;j<nLogBins2+1;j++){
   //inVec[i-1]->push_back(hUnfolded->GetBinContent(j));
   inVec[i-1]->push_back(hRecoSmeared->GetBinContent(j));
  }
 }//end loop over samples

 ///////////////////////////////////////////
 // Part 3: Calculating covariance matrix //
 ///////////////////////////////////////////

 TH2D*hCovM = new TH2D("hCovM","",nLogBins2,massbins2,nLogBins2,massbins2);
 TH2D*hCovMinv = new TH2D("hCovMinv","",nLogBins2,massbins2,nLogBins2,massbins2);
 TH2D*hCorrM = new TH2D("hCorrM","",nLogBins2,massbins2,nLogBins2,massbins2);
 double x,y,xy,xAvg,yAvg,xyAvg,xAvg2,yAvg2,xyAvg2,xStd,yStd,xyStd,covXY,corrXY;
 double xSum,ySum,xySum,xSum2,ySum2,xySum2;
 TMatrixD fCovM(nLogBins2,nLogBins2);

  for(int i=1;i<nLogBins2+1;i++){
   for(int j=1;j<nLogBins2+1;j++){
    xSum = 0;
    ySum = 0;
    xySum = 0;
    xSum2 = 0;
    ySum2 = 0;
    xySum2 = 0; 
    for(int iSample=0;iSample<nSamples;iSample++){
    x = inVec[iSample]->at(i-1);
    y = inVec[iSample]->at(j-1);
    xy = x*y;
    xSum += x;
    ySum += y;
    xySum += xy;

    xSum2 += x*x;
    ySum2 += y*y;
    xySum2 += xy*xy;
   }//end loop over samples
  
   xAvg = xSum/(1.0*nSamples);
   yAvg = ySum/(1.0*nSamples);
   xyAvg = xySum/(1.0*nSamples);

   xAvg2 = xSum2/(1.0*nSamples);
   yAvg2 = ySum2/(1.0*nSamples);
   xyAvg2 = xySum2/(1.0*nSamples);

   xStd = sqrt(xAvg2-xAvg*xAvg);
   yStd = sqrt(yAvg2-yAvg*yAvg);
   xyStd = sqrt(xyAvg2-xyAvg*xyAvg);

   covXY = (xyAvg-xAvg*yAvg);
   corrXY = covXY/(xStd*yStd);
   hCovM->SetBinContent(i,j,covXY);
   hCorrM->SetBinContent(i,j,corrXY);
   fCovM(i-1,j-1) = covXY;
  }//end j loop
 }//end i loop
 TMatrixD fCovMinv = fCovM;
 fCovMinv.Invert();
 for(int i=0;i<nLogBins2;i++){
  for(int j=0;j<nLogBins2;j++){
   hCovMinv->SetBinContent(i+1,j+1,fCovMinv(i,j));
  }
 } 
 TCanvas*canvas = new TCanvas("canvas","",10,10,1200,1200);
 canvas->SetLogx();
 canvas->SetLogy();
 canvas->SetLogz();
 hCovM->GetYaxis()->SetNoExponent();
 hCovM->GetYaxis()->SetMoreLogLabels();
 hCovM->GetXaxis()->SetNoExponent();
 hCovM->GetXaxis()->SetMoreLogLabels();
 hCovM->SetTitle("Unfolding covariance matrix");
 hCovM->GetXaxis()->SetTitle("mass [GeV]");
 hCovM->GetYaxis()->SetTitle("mass [GeV]");
 hCovM->Draw("colz");
 TString saveName1 = "/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase1/";
 saveName1 += "covariance_";
 saveName1 += nSamples;
 saveName1 += "Samples.png";
 canvas->SaveAs(saveName1);
 TCanvas*canvas2 = new TCanvas("canvas2","",10,10,1200,1200);
 canvas2->SetLogx();
 canvas2->SetLogy();
 canvas2->SetLogz();
 hCorrM->GetYaxis()->SetNoExponent();
 hCorrM->GetYaxis()->SetMoreLogLabels();
 hCorrM->GetXaxis()->SetNoExponent();
 hCorrM->GetXaxis()->SetMoreLogLabels();
 hCorrM->SetTitle("Unfolding correlations matrix");
 hCorrM->GetXaxis()->SetTitle("mass [GeV]");
 hCorrM->GetYaxis()->SetTitle("mass [GeV]");
 hCorrM->Draw("colz");
 TString saveName2 = "/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase1/";
 saveName2 += "inputCorrelations_";
 saveName2 += nSamples;
 saveName2 += "Samples.png";
 canvas2->SaveAs(saveName2);

 TString saveFileName = "./data/inputCorrelations_";//location to save correlations
 saveFileName += nSamples;
 saveFileName += "Samples.root";
 TFile*saveFile = new TFile(saveFileName,"recreate");
 saveFile->cd();
 hCovM->Write();
 hCovMinv->Write();
 hCorrM->Write();
 canvas->Write();
 canvas2->Write();
 saveFile->Write();
 saveFile->Close();

 totaltime.Stop();
 Double_t TotalCPURunTime = totaltime.CpuTime();
 Double_t TotalRunTime = totaltime.RealTime();
 TTimeStamp ts_end;
 cout << endl;
 cout << "**************************************************************************" << endl;
 cout << "Total CPU RunTime: " << TotalCPURunTime/60 << " minutes" << endl;
 cout << "Total Real RunTime: " << TotalRunTime/60 << " minutes" << endl;
 cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
 cout << "Number of Samples Processed: " << nSamples << endl;
 cout << "**************************************************************************" << endl;
 cout << endl;

}//end main()



























