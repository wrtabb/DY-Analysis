enum RegType {      //Strength of regularization
  NO_REG,           //No regularization
  CONST_REG,        //User defined regularization
  VAR_REG_LCURVE,   //TUnfoldDensity determines best choice of regularization strength
  VAR_REG_SCANSURE, //TUnfoldDensity determines best choice of regularization strength
  VAR_REG_SCANTAU   //TUnfoldDensity determines best choice of regularization strength
};
enum UnfoldType {
  TUNFOLD,
  INVERSION
};
//true binning
double binningTrue[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50};
//reco binning, each true bin split in half
//I'm in the process of testing different potential binning schemes
//But they all exhibit the same oscillatory behavior
double binningReco[] = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34,34.5,35,35.5,36,36.5,37,37.5,38,38.5,39,39.5,40,40.5,41,41.5,42,42.5,43,43.5,44,44.5,45,45.5,46,46.5,47,47.5,48,48.5,49,49.5,50};
const int nBinsTrue = size(binningTrue)-1;
const int nBinsReco = size(binningReco)-1;

//Forward declarations of functions
void makeToyModels();
void plotUnfolded(TH1D*hReco,TH1D*hTrue,TH1F*hUnfoldedE,UnfoldType unfoldType,bool closure);
void unfoldTUnfold(RegType regType,TH1D*hReco,TH1D*hClosure,TH1D*hTrue,TH2D*hMatrix,bool closure);
void unfoldInversion(TH1D*hReco,TH1D*hClosure,TH1D*hTrue,TH2D*hResponse,bool closure);
double GetConditionNumber(TH2D*hResponse);
TMatrixD makeMatrixFromHist(TH2D*hist);
TVectorD makeVectorFromHist(TH1D*hist);
TH2D*makeResponseMatrix(TH2D*hist);
TH1F*makeHistFromVector(TVectorD vec,TH1D*hist);

double mean = 50;
int nBinsR = 100;
int nBinsT = nBinsR/2;

//Main Function
void test()
{
 TH1::SetDefaultSumw2();
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);
 gROOT->SetBatch(true);

 const TString fileName = "unfoldingHists.root";
 TFile*file = new TFile(fileName);

 //Function to make the toy models
 //Saves to a root file
 makeToyModels();

 //Get toy models saved to root file to perform unfolding on
 TH1D*hReco =     (TH1D*)file->Get("hReco");
 TH1D*hClosure =  (TH1D*)file->Get("hClosure");
 TH1D*hTrue =     (TH1D*)file->Get("hTrue");
 TH2D*hMatrix =   (TH2D*)file->Get("hMatrix");
 TH2D*hResponse = (TH2D*)file->Get("hResponse");

 //Perform the unfolding in two ways:
 // 1. TUnfold
 unfoldTUnfold(NO_REG,hReco,hClosure,hTrue,hMatrix,false);
 // 2. Matrix inversion
 unfoldInversion(hReco,hClosure,hTrue,hResponse,false);
}

void unfoldTUnfold(RegType regType,TH1D*hReco,TH1D*hClosure,TH1D*hTrue,TH2D*hMatrix,bool closure)
{
 cout << endl;
 cout << "*****************************************" << endl;
 cout << "Beginning unfolding process using TUnfold" << endl;
 cout << "*****************************************" << endl;
 cout << endl;

 if(closure){
  hReco = (TH1D*)hClosure->Clone();
 }

 ////////////////////////////
 //  Regularization Modes  //
 ////////////////////////////
 //TUnfold::ERegMode regMode = TUnfold::kRegModeNone; //breaks 
 //TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
 //TUnfold::ERegMode regMode = TUnfold::kRegModeDerivative;
 TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
 //TUnfold::ERegMode regMode = TUnfold::kRegModeMixed; //breaks

 ///////////////////////////
 //  Types of Constraint  //
 ///////////////////////////
 //TUnfold::EConstraint constraintMode = TUnfold::kEConstraintNone;
 TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;

 /////////////////////
 //  Density Modes  //
 /////////////////////
 //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeNone;
 TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
 //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeUser;
 //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidthAndUser;

 /////////////////////////////////////
 //  Horizontal vs Vertical Output  //
 /////////////////////////////////////
 //TUnfold::EHistMap outputMap = TUnfold::kHistMapOutputVert;
 TUnfold::EHistMap outputMap = TUnfold::kHistMapOutputHoriz;

 //////////////////////////////////////
 //  Constructor for TUnfoldDensity  //
 //////////////////////////////////////
 TUnfoldDensity unfold(hMatrix,outputMap,regMode,constraintMode,densityFlags);
 unfold.SetInput(hReco);//the measured distribution

 ////////////////////////////
 //  Begin Regularization  //
 ////////////////////////////
 Int_t iBest;
 TSpline *logTauX,*logTauY;
 TGraph *lCurve;
 TGraph*bestLcurve;
 TGraph*bestLogTauLogChi2;
 if(regType == VAR_REG_LCURVE){
   Int_t nScan=30;//This number chosen only because it was given in the tutorial
   Double_t tauMin = 0.0;//If tauMin=tauMax, TUnfold automatically chooses a range
   Double_t tauMax = 0.0;//Not certain how they choose the range
   iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
   cout<< "tau=" << unfold.GetTau() << endl;
   cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()<<" / "<<unfold.GetNdf()<<"\n";
   Double_t t[1],x[1],y[1];
   logTauX->GetKnot(iBest,t[0],x[0]);
   logTauY->GetKnot(iBest,t[0],y[0]);
   bestLcurve=new TGraph(1,x,y);
   bestLogTauLogChi2=new TGraph(1,t,x);
 }
 else if(regType == VAR_REG_SCANSURE){
  cout << "Option VAR_REG_SCANSURE not yet implemented" << endl;
  return;
 }
 else if(regType == VAR_REG_SCANTAU){
  cout << "Option VAR_REG_SCANTAU not yet implemented" << endl;
  return;
 }
 else if(regType == NO_REG){
   double tau = 0;
   unfold.DoUnfold(tau,hReco);
 }
 else{//user defined
   double tau = 1e-8;//larger tau introduces more MC bias
   unfold.DoUnfold(tau,hReco);
 }

 //Create unfolded histogram
 TH1*hUnfolded = unfold.GetOutput("Unfolded");
 
 //Create error matrices
 TH2*histEmatStat=unfold.GetEmatrixInput("unfolding stat error matrix");
 TH2*histEmatTotal=unfold.GetEmatrixTotal("unfolding total error matrix");

 //Create unfolding histogram with errors
 TH1F*hUnfoldedE = (TH1F*)hUnfolded->Clone("Unfolded with errors");

 //loop over unfolded histogram bins and assign errors to each one
 for(int i=0;i<=nBinsT;i++){
  double c = hUnfolded->GetBinContent(i+1);
  hUnfoldedE->SetBinContent(i+1,c);
  hUnfoldedE->SetBinError(i+1,TMath::Sqrt(histEmatTotal->GetBinContent(i+1,i+1)));
 }
 plotUnfolded(hReco,hTrue,hUnfoldedE,TUNFOLD,closure);
}

void plotUnfolded(TH1D*hReco,TH1D*hTrue,TH1F*hUnfoldedE,UnfoldType unfoldType,bool closure)
{
 //Rebin the reco distribution to plot against the true distribution for easier comparison
 TH1D*hRecoRebin = (TH1D*)hReco->Clone();
 if(unfoldType==TUNFOLD) hRecoRebin->Rebin(2);

 hUnfoldedE->SetMarkerStyle(25);
 hUnfoldedE->SetMarkerColor(kBlue+2);
 hUnfoldedE->SetMarkerSize(1);
 hUnfoldedE->SetFillColor(kWhite);
 hRecoRebin->SetMarkerStyle(20);
 hRecoRebin->SetMarkerColor(kBlack);
 hRecoRebin->SetLineColor(kBlack);
 hRecoRebin->SetFillColor(kWhite);
 hTrue->SetFillColor(kRed+2);
 hTrue->SetLineColor(kRed+2);
 hTrue->SetTitle("");
 hTrue->SetLabelSize(0);
 hTrue->SetTitleSize(0);
 
 //Ratio of unfolded histogram over true histogram
 TH1F*ratio = (TH1F*)hUnfoldedE->Clone("ratio");
 ratio->Divide(hTrue);

 //postion of chi2 label to place on plots
 double xChiLabel = 35;
 double yChiLabel = 5e5;
 double x[nBinsT],res[nBinsT];
 //create chi2 label
 double chi = hUnfoldedE->Chi2Test(hTrue,"CHI2/NDF",res);//chi2/ndf to print on plot
 TLatex*chiLabel = new TLatex(xChiLabel,yChiLabel,Form("#chi^{2}/ndf = %lg", chi));

 //Make canvas and pads and set all options
 const float padmargins = 0.03;
 TCanvas*canvas1 = new TCanvas("canvas1","",10,10,1200,1000);
 TPad*pad1 = new TPad("","",0,0.3,1.0,1.0);
 pad1->SetBottomMargin(padmargins);
 pad1->SetGrid();
 pad1->SetTicks(1,1);
 pad1->Draw();
 pad1->cd();
 double binLow = 0.0;
 double binHigh = 50.0;
 TLine*line = new TLine(binLow,1,binHigh,1);
 line->SetLineColor(kRed);
 TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
 legend->SetTextSize(0.02);
 legend->AddEntry(hTrue,"True Distribution");
 legend->AddEntry(hRecoRebin,"Measured Distribution");
 legend->AddEntry(hUnfoldedE,"Unfolded Distribution");
 hTrue->Draw("hist");
 hRecoRebin->Draw("PE,same");
 hUnfoldedE->Draw("PE,same");
 legend->Draw("same");
 //chiLabel->Draw("same");

 canvas1->cd();
 TPad*pad2 = new TPad("","",0,0.05,1,0.3);
 pad2->SetTopMargin(padmargins);
 pad2->SetBottomMargin(0.2);
 pad2->SetGrid();
 pad2->SetTicks(1,1);
 pad2->Draw();
 pad2->cd();
 ratio->SetMinimum(0.7);
 ratio->SetMaximum(1.3);
 ratio->GetYaxis()->SetLabelSize(0.06);
 ratio->GetYaxis()->SetTitleSize(0.08);
 ratio->GetYaxis()->SetTitleOffset(0.3);
 ratio->GetYaxis()->SetTitle("Unfolded/Truth");
 ratio->GetXaxis()->SetLabelSize(0.1);
 ratio->GetXaxis()->SetTitleSize(0.1);
 ratio->SetMarkerStyle(20);
 ratio->SetMarkerColor(kBlack);
 ratio->Draw("PE");
 line->Draw("same");
 
 TString saveName = "testPlots/unfolded";
 if(unfoldType==TUNFOLD) saveName += "TUnfold";
 else if(unfoldType==INVERSION) saveName += "Inversion";
 if(closure) saveName += "ClosureTest";
 saveName += "_NoRegularization_Mean";
 saveName += mean;
 saveName += "_RecoBins";
 saveName += nBinsR;
 saveName += ".png";
 canvas1->SaveAs(saveName);

}

void makeToyModels()
{//Make toy models to perform unfolding on

 cout << endl;
 cout << "*************************************" << endl;
 cout << "Creating toy models to test unfolding" << endl;
 cout << "*************************************" << endl;
 cout << endl;

 TH1::SetDefaultSumw2();
 gStyle->SetOptStat(0);
 gStyle->SetPalette(1);

 int nBinsT = nBinsR/2;
 
 TH1D*hReco =        new TH1D("hReco","",nBinsR,0,50);//measured distribution
 TH1D*hRecoClosure = new TH1D("hClosure","",nBinsR,0,50);//measured distribution using same seed as migration matrix
 TH1D*hTrue =        new TH1D("hTrue","",nBinsT,0,50);//true distribution
 TH2D*hMatrix =      new TH2D("hMatrix","",nBinsT,0,50,nBinsR,0,50);//2D matrix of true versus reco distributions
 //-----Parameters for toy model-----//
 const double sigma = 5;//Sigma of gaussian part of distribution
 const double mean_smeared = 0;//mean of smearing used to make reco
 const double sigma_smeared = 1.0;//sigma of smearing used to make reco
 const Long64_t nEntries = 1e7;//number of entries

 //Functional model to create distributions
 TF1*func = new TF1("func","1/(x+1)+gaus(0)",0,50);
 func->SetParameters(1.0,mean,sigma);

 TRandom3 gen1;
 TRandom3 gen2;
 gen1.SetSeed(82);
 gen2.SetSeed(1981);

 //Save distributions for later retrieval
 TFile*saveFile = new TFile("unfoldingHists.root","recreate");

 //Fill all distributions
 double peak,peakReco,peak_smeared,peakReco_smeared;
 for(int i=0;i<nEntries;i++){
  peak = func->GetRandom();
  peakReco = func->GetRandom();
  peak_smeared = peak+gen1.Gaus(mean_smeared,sigma_smeared);
  peakReco_smeared = peakReco+gen2.Gaus(mean_smeared,sigma_smeared);

  hRecoClosure->Fill(peak_smeared);
  hReco->Fill(peakReco_smeared);
  hTrue->Fill(peak);
  hMatrix->Fill(peak,peak_smeared);
 }//end loop over entries
 TH2D*hMatrixRebin = (TH2D*)hMatrix->Clone("hMatrixRebin");
 TH2D*hResponse = makeResponseMatrix(hMatrixRebin);
 double conditionNumber = GetConditionNumber(hResponse);
 double xPosition = 15;
 double yPosition = 5;
 TLatex*conditionLabel = new TLatex(xPosition,yPosition,
                                    Form("condition number = %lg", conditionNumber));
 //Draw the response matrix
 TCanvas*canvas2 = new TCanvas("canvas2","",0,0,1000,1000);
 canvas2->SetGrid();
 hResponse->Draw("colz");
 conditionLabel->Draw("same");
 TString responseSaveName = "testPlots/responseMatrix_Mean";
 responseSaveName += mean;
 responseSaveName += "_RecoBins";
 responseSaveName += nBinsR;
 responseSaveName += ".png";
 canvas2->SaveAs(responseSaveName);
 
 //Draw the migration matrix
 TCanvas*canvas = new TCanvas("canvas","",0,0,1000,1000);
 canvas->SetGrid();
 hMatrix->Draw("colz");
 TString plotSaveName = "testPlots/migrationMatrix_Mean";
 plotSaveName += mean;
 plotSaveName += "_RecoBins";
 plotSaveName += nBinsR;
 plotSaveName += ".png";
 canvas->SaveAs(plotSaveName);

 saveFile->cd();
 hReco->Write();
 hRecoClosure->Write();
 hTrue->Write();
 hMatrix->Write();
 hResponse->Write();
 saveFile->Close();
}

double GetConditionNumber(TH2D*hResponse)
{
 TString histName = hResponse->GetName();
 int nBinsX = hResponse->GetNbinsX();
 int nBinsY = hResponse->GetNbinsY();
 TMatrixD matrix = makeMatrixFromHist(hResponse);

 TDecompSVD decomp(matrix);
 double condition = decomp.Condition();
 cout << "The condition number for " << histName << ": " << condition << endl;

 double determinant;
 TMatrixD mInverse = matrix.Invert(&determinant);
 cout << "The determinant of " << histName << " is " << determinant << endl;
 return condition;
}

TH2D*makeResponseMatrix(TH2D*hist)
{
 hist->RebinY(2);
 TH2D*hResponse = (TH2D*)hist->Clone("hResponse");
 int nBinsX = hist->GetNbinsX();
 int nBinsY = hist->GetNbinsY();
 for(int i=1;i<=nBinsX;i++){
  double nEntriesX = 0;
  for(int j=1;j<=nBinsY;j++){
   nEntriesX += hist->GetBinContent(i,j);
  }
  double sum = 0;
  for(int j=1;j<=nBinsY;j++){
   double scaledContent = hist->GetBinContent(i,j)/nEntriesX;
   hResponse->SetBinContent(i,j,scaledContent);
   sum += scaledContent;
  }
 }
 return hResponse;
}
TMatrixD makeMatrixFromHist(TH2D*hist)
{
 int nBinsX = hist->GetNbinsX();
 int nBinsY = hist->GetNbinsY();
 TMatrixD matrix(nBinsY,nBinsX);
 if(nBinsX != nBinsT || nBinsY != nBinsT){
  cout << "makeMatrixFromHist: Matrix bins don't match histogram bins" << endl;
  return matrix;
 }
 for(int i=1;i<=nBinsX;i++){
  for(int j=1;j<=nBinsY;j++){
   matrix(j-1,i-1) = hist->GetBinContent(i,j);
  }
 } 
 return matrix;
}

TVectorD makeVectorFromHist(TH1D*hist)
{
 int nBins = hist->GetNbinsX();
 TVectorD vec(nBins);
 for(int i=1;i<=nBins;i++){
  vec(i-1) = hist->GetBinContent(i);
 } 
 return vec;
}

TH1F*makeHistFromVector(TVectorD vec,TH1D*hist)
{
 TH1F*hReturn = (TH1F*)hist->Clone("hUnfolded");
 int nBins = vec.GetNrows();
 for(int i=0;i<nBins;i++){
  hReturn->SetBinContent(i+1,vec(i));
 }
 return hReturn;
}

void unfoldInversion(TH1D*hReco,TH1D*hClosure,TH1D*hTrue,TH2D*hResponse,bool closure)
{
 cout << endl;
 cout << "**************************************************" << endl;
 cout << "Beginning unfolding process using Matrix inversion" << endl;
 cout << "**************************************************" << endl;
 cout << endl;

 TString unfoldType = "Inversion";
 if(closure){
  hReco = (TH1D*)hClosure->Clone("hReco");
  unfoldType += "Closure";
 }
 //Make response matrix histogram from input histogram matrix
 //Meaning we normalize each column (true bins) 
 //Here I am assuming the reco has twice as many bins as tru
 //This is because this is what I'm using for this specific task
 //This would have to be treated differently if this wasn't the case
 hReco->Rebin(2);

 //Turn histograms into matrices and vectors
 TMatrixD responseM = makeMatrixFromHist(hResponse); 
 TVectorD trueV = makeVectorFromHist(hTrue);
 TVectorD recoV = makeVectorFromHist(hReco);

 //Invert
 TMatrixD invertedM = responseM.Invert();
 TVectorD unfoldedV = invertedM*recoV;

 //Get covariance (assuming Vy = identity)
 TMatrixD invertedMT = invertedM.T();
 TMatrixD Vx = invertedM*invertedMT;

 TH1F*hUnfolded = makeHistFromVector(unfoldedV,hTrue);
 TH1F*hUnfoldedE = (TH1F*)hUnfolded->Clone("unfolded with errors");
 cout << "Number of bins = " << nBinsT << endl;
 
 //Get error bars from diagonal of covariance matrix for unfolded histogram
 for(int i=0;i<nBinsT;i++){
  hUnfoldedE->SetBinError(i+1,TMath::Sqrt(Vx(i,i)));
 }

 
 plotUnfolded(hReco,hTrue,hUnfoldedE,INVERSION,closure);
 
 TCanvas*c = new TCanvas("c","",0,0,1000,1000);
 c->SetGrid();
 hUnfoldedE->Draw("hist");
 hUnfolded->Draw("pe,same");
 c->SaveAs("temp.png");
}

