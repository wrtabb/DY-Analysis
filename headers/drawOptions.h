#include "TCanvas.h"
#include "TH1.h"
enum plotType{
 POINT,
 HIST
};
void hist2DPlot(TCanvas*canvas,TH2*hist,Option_t*draw,bool logX,bool logY,bool logZ);
void histPlot(TCanvas*canvas,TH1*hist1,Option_t*draw1,TH1*hist2,Option_t*draw2,bool logX,bool logY);
void ratioPlot(TCanvas*canvas,TH1*hist1,TH1*hist2,TH1*hist3);

void ratioPlot(TCanvas*canvas,TH1*hist1,TH1*hist2,TH1*hist3)
{
 const float padmargins = 0.03;
 const int nBins = 43;
 const int binLow = 15;
 const int binHigh = 3000;
 double x[nBins],res[nBins];
 double pValues = hist3->Chi2Test(hist1,"P",res);
 double chi = hist3->Chi2Test(hist1,"CHI2/NDF",res);
 TLatex*chiLabel = new TLatex(500.0,150000,Form("#chi^{2}/ndf = %lg", chi));

 TPad*pad1 = new TPad("","",0,0.3,1.0,1.0);
 pad1->SetBottomMargin(padmargins);
 pad1->SetGrid();
 pad1->SetLogy();
 pad1->SetLogx();
 pad1->SetTicks(1,1);
 pad1->Draw();
 pad1->cd();
 TLine*line = new TLine(binLow,1,binHigh,1);
 line->SetLineColor(kRed);
 hist1->SetLabelSize(0);
 hist1->SetTitleSize(0);
 hist1->Draw("hist");
 hist2->Draw("PE,same");
 hist3->Draw("PE,same");
 chiLabel->Draw("same");

 canvas->cd();
 TPad*pad2 = new TPad("","",0,0.05,1,0.3);
 pad2->SetLogx();
 pad2->SetTopMargin(padmargins);
 pad2->SetBottomMargin(0.2);
 pad2->SetGrid();
 pad2->SetTicks(1,1);
 pad2->Draw();
 pad2->cd();
 TH1D*ratio = (TH1D*)hist3->Clone("ratio");
 ratio->Divide(hist1);

 ratio->GetYaxis()->SetLabelSize(0.06);
 ratio->GetYaxis()->SetTitleSize(0.08);
 ratio->GetYaxis()->SetTitleOffset(0.3);
 ratio->GetYaxis()->SetTitle("True/Unfolded");
 ratio->GetXaxis()->SetTitle("mass [GeV]");
 ratio->GetXaxis()->SetLabelSize(0.1);
 ratio->GetXaxis()->SetTitleSize(0.1);
 ratio->GetXaxis()->SetNoExponent();
 ratio->GetXaxis()->SetMoreLogLabels();
 ratio->GetYaxis()->SetRangeUser(0.8,1.2);
 ratio->SetMarkerStyle(20);
 ratio->SetMarkerColor(kBlack);
 ratio->Draw("PE");
 line->Draw("same");
}

void histPlot(TCanvas*canvas,TH1*hist1,Option_t*draw1,TH1*hist2,Option_t*draw2,bool logX,bool logY)
{
 canvas->SetGrid();
 if(logX){
  canvas->SetLogx();
  hist1->GetXaxis()->SetNoExponent();
  hist1->GetXaxis()->SetMoreLogLabels();
 }
 if(logY){
  canvas->SetLogy();
 }  
 hist1->Draw(draw1);
 if(hist2!=NULL) hist2->Draw(draw2);
}

void hist2DPlot(TCanvas*canvas,TH2*hist,Option_t*draw,bool logX,bool logY,bool logZ)
{
 canvas->SetLeftMargin(0.12);
 canvas->SetRightMargin(0.16);
 canvas->SetGrid();
 hist->GetYaxis()->SetTitleOffset(1.8);
 if(logX){
  canvas->SetLogx();
  hist->GetXaxis()->SetNoExponent();
  hist->GetXaxis()->SetMoreLogLabels();
 }
 if(logY){
  canvas->SetLogy();
  hist->GetYaxis()->SetNoExponent();
  hist->GetYaxis()->SetMoreLogLabels();
 }
 if(logZ){
  canvas->SetLogz();
 }
 hist->Draw(draw);
}
