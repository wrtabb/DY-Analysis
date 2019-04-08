#include "/home/hep/wrtabb/git/DY-Analysis/headers/header1.h"

const Long64_t nSamples = 1000;//number of input vectors to create
const TString fileName = "toyUnfold.root";//location of toy model distributions

void calcCorr()
{
 //initialize variables, arrays, TFile, and set histogram options
 double error,binContent,smear;
 double vector[nSamples][nLogBins];
 TFile*file = new TFile(fileName);
 TH1::SetDefaultSumw2();
 TCanvas*c[nSamples];
 TCanvas*c1[nSamples];
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);

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
 TH1D*hRecoSmeared[nSamples];
 //Begin loop over samples doing unfold on each one and placing it in an array
 for(Long64_t i=1;i<nSamples+1;i++){
  ////////////////////////////////////////////////////////////////////////
  // Part 1: smearing each histogram to create various input histograms //
  ////////////////////////////////////////////////////////////////////////
  //gRandom->SetSeed(time(0));
  TString recoName = "hRecoSmeared";
  recoName += i;
  hRecoSmeared[i-1] = (TH1D*)hReco->Clone("hRecoSmeared");
  //vary each bin with a gaussian of width equal to the bin error
  for(int j=1;j<nLogBins2+1;j++){
   error = hReco->GetBinError(j);
   //not certain what range to place on this TF1
   TF1*fGaus = new TF1("fgaus","gaus(0)",-100*error,100*error);
   fGaus->SetParameters(1,0,error);
   smear = fGaus->GetRandom();
   binContent = hReco->GetBinContent(j)+smear;
   hRecoSmeared[i-1]->SetBinContent(j,binContent);
   //Now hReco is varied
  }

  /////////////////////////////////////////////////////
  // Part 2: Doing the unfolding with TUnfoldDensity //
  /////////////////////////////////////////////////////

  //Set unfolding options
  TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
  TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;
  TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
  TUnfold::EHistMap outputMap = TUnfold::kHistMapOutputHoriz;

  //Do the unfolding
  TUnfoldDensity unfold(hMatrix,outputMap,regMode,constraintMode,densityFlags);
  unfold.SetInput(hRecoSmeared[i-1]);
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

  //Drawing for each sample
  //Comment out this part if you run on lots of events
/*
  TString cName = "canvas";
  cName += i;
  TString hName = "hist";
  hName += i;
  c[i] = new TCanvas(cName,"",10,10,1000,1000);
  c[i]->SetLogx();
  c[i]->SetLogy();
  c[i]->SetGrid();
  c[i]->cd();
  TH1D*hRecoRebin = (TH1D*)hRecoSmeared[i-1]->Clone();
  hRecoRebin->Rebin(2);
  TH1D*hUnfDraw = (TH1D*)hUnfolded->Clone(); 
   hUnfDraw->SetMarkerStyle(25);
   hUnfDraw->SetMarkerColor(kBlue+2);
   hUnfDraw->SetMarkerSize(1);

  hTrue->Draw("hist"); 
  hUnfDraw->Draw("PE,same");
  hRecoRebin->Draw("same");
*/
/*
  TString ematName = "eMatrix";
  ematName += i;
  TH2*histEmatTotal=unfold.GetEmatrixTotal(ematName);
  TString cName1 = "canvas1";
  cName1 += i;
  c1[i] = new TCanvas(cName1,"",10,10,1000,1000);
  c1[i]->SetLogx();
  c1[i]->SetLogy();
  c1[i]->SetLogz();
  c1[i]->SetGrid();
  c1[i]->cd();
  histEmatTotal->Draw("colz");
  */
  //Now that unfolding has been done, place results into an array
  for(int j=1;j<nLogBins+1;j++){
    vector[i-1][j-1] = hUnfolded->GetBinContent(j);
  }
 }//end loop over samples

 ///////////////////////////////////////////
 // Part 3: Calculating covariance matrix //
 ///////////////////////////////////////////
 
 //Initializing matrices
 //matrixB is made up of each vector with the average vector subtracted from it
 TMatrixD matrixB(nLogBins,nSamples);
 TMatrixD matrixBT(nSamples,nLogBins);//transpose of matrixB
 TMatrixD covM(nLogBins,nLogBins);//covariance matrix
 TMatrixD corrM(nLogBins,nLogBins);//correlation matrix
 
 //arrays needed for covariance calculation
 double vecSum[nLogBins],vecAvg[nLogBins],vecM[nSamples][nLogBins],vecSum2[nLogBins],
  vecStd[nLogBins];
 
 //set vector sums such that all components are zero
 for(int i=0;i<nLogBins;i++){
  vecSum[i] = 0;
  vecSum2[i]=0;
 }

 //loop over all vectors and components to calculate the average of each component
 for(int iEle=0;iEle<nLogBins;iEle++){
  for(int jVec=0;jVec<nSamples;jVec++){
   vecSum[iEle] += vector[jVec][iEle];//sum of element iEle summed over vectors
   vecSum2[iEle] += vector[jVec][iEle]*vector[jVec][iEle];//sum of squares of elements
  }
  vecAvg[iEle] = vecSum[iEle]/nSamples;//vector of averages
  //vector of standard deviations
  vecStd[iEle] = sqrt(vecSum2[iEle]/nSamples-vecSum[iEle]*vecSum[iEle]/(nSamples*nSamples*1.0));
 }
 
 //calculation of matrixB
 for(int iEle=0;iEle<nLogBins;iEle++){
  for(int jVec=0;jVec<nSamples;jVec++){
   vecM[jVec][iEle] = vector[jVec][iEle]-vecAvg[iEle];//vectors minus average vector
   matrixB(iEle,jVec) = vecM[jVec][iEle];
  }
 }

 matrixBT.Transpose(matrixB);//the transpose of matrixB calculated above
 covM = matrixB*matrixBT;//covariance matrix without proper weight
 //Define a histogram to hold covariance matrix so it can be plotted
 TH2D*hCovM = new TH2D("hCovM","",nLogBins,massbins,nLogBins,massbins);
 double val;
 //calculate covariance and correlation matrices
 for(int jEle=0;jEle<nLogBins;jEle++){
  for(int iEle=0;iEle<nLogBins;iEle++){
   covM(iEle,jEle) = covM(iEle,jEle)/nSamples;//covariance matrix
   corrM(iEle,jEle) = covM(iEle,jEle)/(vecStd[iEle]*vecStd[jEle]);//correlation matrix
   val = covM(iEle,jEle);//value in each bin to place in the histogram
   hCovM->SetBinContent(iEle+1,jEle+1,val);
  }
 }

 //Drawing the covariance matrix
 TCanvas*canvas = new TCanvas("canvas","",1000,1000);
 canvas->SetLogy();
 canvas->SetLogx();
 canvas->SetLogz();
 hCovM->Draw("colz");
 TString saveName = "/home/hep/wrtabb/git/DY-Analysis/plots/unfolding/phase1Plots/covariance_";
 saveName += nSamples;
 saveName += "Samples.png";
 canvas->SaveAs(saveName);

 //Print out of vectors for sanity check
 /*
 TString vecName;
 TString avgMName;
 for(int j=1;j<nSamples+1;j++){
  vecName = "Vector";
  vecName += j;
  vecName += ": (";
  cout << vecName;
  for(int i=0;i<nLogBins;i++){
   cout << vector[j-1][i];
   if(i!=nLogBins-1) cout << ", ";
  } 
  cout << ")" << endl;
 }
 cout << "Avg Vector: (";
 for(int i=0;i<nLogBins;i++){
  cout << vecAvg[i];
  if(i!=nLogBins-1) cout << ", ";
 }
 cout <<")" << endl;

 for(int j=1;j<nSamples+1;j++){
  avgMName= "Vector-Avg";
  avgMName+= j;
  avgMName+= ": (";
  cout << avgMName;
  for(int i=0;i<nLogBins;i++){
   cout << vecM[j-1][i];
   if(i!=nLogBins-1) cout << ", ";
  } 
  cout << ")" << endl;
 }
 */
 //Output bin errors to a .txt file
 ofstream errorFile;
 errorFile.open("myCovErrors.txt");
 double binError;
 for(int i=1;i<nLogBins+1;i++){
  binError = TMath::Sqrt(hCovM->GetBinContent(i,i));
  errorFile << "Bin: " << i << ", Error: " << binError << endl;
 } 
 errorFile.close();
}//end main()



























