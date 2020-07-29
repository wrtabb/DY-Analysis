
#include "VariableList.h"

//forward declaration
TH2D*Rebin2D(TH2D*hist,TString histName,std::vector<double> binning);
TCanvas*MakeCanvas(TString cName,bool log);
void OutputEntries(TH2D*h1,TH2D*h2=0,TH2D*h3=0,TH2D*h4=0,TH2D*h5=0);
void CheckBinByBin(TH2D*h1,TH2D*h2);

TString fileName = "data/migrationMatrix.root";

double binning0[] = {0,1,2,3,4,5,10,15,20,25,30,35,40,45,50};
double binning1[] = {0,10,20,30,40,50};
std::vector<double> binningVec0 = {0,1,2,3,4,5,10,15,20,25,30,35,40,45,50};
std::vector<double> binningVec1 = {0,10,20,30,40,50};
int nBins0 = binningVec0.size()-1;
int nBins1 = binningVec1.size()-1;

void rebin2D()
{
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);

 //-----Testing with toy models-----//
 //Initialize histograms with different binnings on the y-axis
 TH2D*h1 = new TH2D("h1","",nBins0,binning0,nBins1,binning1);
 h1->SetTitle("h1");
 TH2D*h2 = new TH2D("h2","",nBins0,binning0,nBins0,binning0);
 h2->SetTitle("h2");

 //Initialize random number generator to fill histograms
 TRandom3 genX,genY;
 genX.SetSeed(1827);
 genY.SetSeed(8282);

 //Set parameters for gaussian shapes for toy histograms 
 double meanX = 25;
 double meanY = 25;
 double sigmaX = 5;
 double sigmaY = 5;
 double x,y;
 int nEntries = 1e6;

 //Fill histograms with same data
 for(int i=0;i<nEntries;i++){
  x = genX.Gaus(meanX,sigmaX);
  y = genY.Gaus(meanY,sigmaY);
  h1->Fill(x,y);
  h2->Fill(x,y);
 }

 //Create rebinned version of h2 that matches h1
 TH2D*h2Rebin = Rebin2D(h2,"h2Rebin",binningVec1);
 h2Rebin->SetTitle("h2Rebin");

 //-----Drell-Yan Distributions-----//
 //File for opening Drell-Yan distributions
 TFile*file = new TFile(fileName);

 //Original migration matrix
 TH2D*hMatrixM = (TH2D*)file->Get("hMatrixMass");
 hMatrixM->SetTitle("Original histogram");

 //Rebinned migration matrix, in this case using the same binning as the original
 //Should return the exact same matrix
 TH2D*hMatrixRebin = Rebin2D(hMatrixM,"hMatrixRebin",massbinsReco);
 hMatrixRebin->SetTitle("Rebinned histogram");

 //Drawing the matrices on separate histograms for visual comparison
 /*
 TCanvas*c1 = MakeCanvas("c1",false);
 h1->Draw("colz");
 TCanvas*c2 = MakeCanvas("c2",false);
 h2Rebin->Draw("colz");
 TCanvas*c3 = MakeCanvas("c3",false);
 h2->Draw("colz");
 */
 //Output of the number of entries in each histogram as a numerical comparison
 //OutputEntries(h1,h2Rebin);
 CheckBinByBin(hMatrixM,hMatrixRebin);
 CheckBinByBin(h1,h2Rebin);
}

void CheckBinByBin(TH2D*h1,TH2D*h2)
{
 //Function to compare the bin contents of each bin in two 2D histograms
 int nBinsX = h1->GetNbinsX();
 int nBinsY = h1->GetNbinsY();
 if( (nBinsX != h2->GetNbinsX()) || (nBinsY != h2->GetNbinsY())){
  //If the histograms have different bins they cannot be compared
  cout << "---------------------------------" << endl;
  cout << "ERROR: CheckBinByBin():" << endl;
  cout << "Histogram bins do not match" << endl;
  cout << "Cannot compare" << endl;
  cout << "---------------------------------" << endl;
  return;
 }

 //loop over all bins x,y and check if bin contents are the same
 int notMatchingCount = 0;
 for(int i=1;i<=nBinsX;i++){
  for(int j=1;j<=nBinsY;j++){
   if(h1->GetBinContent(i,j) - h2->GetBinContent(i,j) != 0){
    //If the bin contents of any bin do not match
    //Output to screen which bin it is
    cout << "Bin " << i << ", " << j << " does not match!" << endl;
    notMatchingCount++;
   }//end if
  }//end loop over Y
 }//end loop over X
 cout << "---------------------------------------------------" << endl;
 cout << "Total number of bins: " << nBinsX*nBinsY << endl;
 cout << "Number of bins not matching: " << notMatchingCount << endl;
 cout << "Percent not matching: " << 100.0*notMatchingCount/(nBinsX*nBinsY) << "%" << endl;
 cout << "---------------------------------------------------" << endl;
}

void OutputEntries(TH2D*h1,TH2D*h2=0,TH2D*h3=0,TH2D*h4=0,TH2D*h5=0)
{
 //Outputs number of entries for up to five histograms
 //Useful check to make sure no entries go missing during rebinning
 cout << "---------------------------" << endl;
 cout << "Entries:" << endl;
 cout << "h1: " << h1->Integral() << endl;
 if(h2) cout << "h2: " << h2->Integral() << endl;
 if(h3) cout << "h3: " << h3->Integral() << endl;
 if(h4) cout << "h4: " << h4->Integral() << endl;
 if(h5) cout << "h5: " << h5->Integral() << endl;
 cout << "---------------------------" << endl;
}

TCanvas*MakeCanvas(TString cName,bool log)
{
 //Makes canvases to draw on
 TCanvas*canvas=new TCanvas(cName,"",0,0,1000,1000);
 canvas->SetGrid();
 if(log){
  canvas->SetLogx();
  canvas->SetLogy();
 }
 canvas->SetLogz();
 return canvas;
}

TH2D*Rebin2D(TH2D*hist,TString histName,std::vector<double> binning)
{
 //This function takes a 2D histogram as input and gives a new 2D histogram
 //As output with new Y-axis binning defined by input std::vector<double> binning
 
 int nBinsHistY = hist->GetNbinsY();//number of reco bins in the input histogram
 int nBinsHistX = hist->GetNbinsX();//numner of true bins in the input histogram
 double binningX[nBinsHistX];//true binning
 int nBinsReco = binning.size()-1;//number of reco bins for the output rebinned histogram

 TH1D*histX = hist->ProjectionX();
 for(int i=0;i<=nBinsHistX+1;i++){
  if(i==0) binningX[i] = histX->GetBinLowEdge(i+1);
  else binningX[i] = binningX[i-1]+histX->GetBinWidth(i);
 }
 
 double newbinning[nBinsReco];//binning placed in array for defining TH2
 for(int i=0;i<=nBinsReco;i++){
  newbinning[i] = binning.at(i);
 }

 //Initialize the output histogram
 TH2D*hRebin = new TH2D(histName,"",nBinsHistX,binningX,nBinsReco,newbinning);

 //Loop over each bin in the histogram
 //For each bin, define the position x and y from the bin center
 //Of the input histogram
 //Then do histogram fill for each of the entries in that bin 
 //For the output histogram
 for(int i=1;i<=nBinsHistX;i++){
  for(int j=1;j<=nBinsHistY;j++){
   double y,x;
   x = hist->GetXaxis()->GetBinCenter(i);
   y = hist->GetYaxis()->GetBinCenter(j);
   int nEntries = hist->GetBinContent(i,j);
   for(int k=0;k<nEntries;k++){
    hRebin->Fill(x,y);
   }//end filling hRebin
  }//end y bin loop
 }//end x bin loop
 return hRebin;
}

