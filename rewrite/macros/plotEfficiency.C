
TGraphAsymmErrors* setOptions(TEfficiency*eff,TPad*pad);

const TString fileName = "/home/hep/wrtabb/DY-Analysis/rewrite/data/efficiencies.root";
const TString saveLocation = "/home/hep/wrtabb/DY-Analysis/rewrite/plots/efficiency/";
void plotEfficiency()
{
 TFile*file = new TFile(fileName);

 TEfficiency*acceptance = (TEfficiency*)file->Get("Acceptance");
 TEfficiency*recoGenEfficiency = (TEfficiency*)file->Get("RecoGenEfficiency");
 TEfficiency*mediumIDEfficiency = (TEfficiency*)file->Get("MediumIDEfficiency");

 TF1*func = new TF1("func","x",0,1);
 TGaxis*midAxis = new TGaxis(3000,0,3000,1,"func",510,"+");
 midAxis->SetLabelSize(0);
 TCanvas*canvas = new TCanvas("canvas","",0,0,1200,1000);
 
 TPad*pad1 = new TPad("","",0,0.5,0.5,1,kWhite,0,0);
 pad1->SetFrameBorderMode(0);
 pad1->SetBorderMode(0);
 pad1->SetBorderSize(0);
 pad1->SetRightMargin(0);
 pad1->Draw();
 pad1->cd();
 auto gAcc = setOptions(acceptance,pad1);
 gAcc->GetXaxis()->SetLabelSize(0);
 gAcc->Draw("PE");
 midAxis->Draw("same");

 canvas->cd();
 TPad*pad2 = new TPad("","",0.5,0.5,1,1,kWhite,0,0);
 pad2->SetFrameBorderMode(0);
 pad2->SetBorderMode(0);
 pad2->SetBorderSize(0);
 pad2->SetLeftMargin(0);
 pad2->Draw();
 pad2->cd();
 auto gRecoGen = setOptions(recoGenEfficiency,pad2);
 gRecoGen->GetYaxis()->SetLabelSize(0);
 gRecoGen->GetXaxis()->SetLabelSize(0);
 gRecoGen->Draw("PE");
 midAxis->SetLabelSize(0.035);
 midAxis->SetLabelFont(42);
 midAxis->SetLabelOffset(0.05);
 midAxis->ChangeLabel(1,-1,-1,-1,-1,-1);
 midAxis->Draw("same");


 canvas->cd();
 TPad*pad3 = new TPad("","",0,0,0.5,0.5,kWhite,0,0);
 pad3->SetFrameBorderMode(0);
 pad3->SetBorderMode(0);
 pad3->SetBorderSize(0);
 pad3->SetRightMargin(0);
 pad3->Draw();
 pad3->cd();
 auto gMedID = setOptions(mediumIDEfficiency,pad3);
 gMedID->GetXaxis()->SetTitle("mass [GeV]");
 gMedID->Draw("PE");
 //midAxis->Draw("same");

 canvas->SaveAs(saveLocation+"efficiencies.png");
}

TGraphAsymmErrors* setOptions(TEfficiency*eff,TPad*pad){
 float markSize = 0.6;
 int markStyle = 20;
 pad->SetGrid();
 pad->SetLogx();
 eff->SetMarkerStyle(markStyle);
 eff->SetMarkerSize(markSize);
 eff->Draw();
 gPad->Update();
 auto graph = eff->GetPaintedGraph();
 graph->SetMinimum(0);
 graph->SetMaximum(1); 
 graph->GetXaxis()->SetRangeUser(15,3000);
 graph->GetYaxis()->SetRangeUser(0,1);
 graph->GetXaxis()->SetNoExponent();
 graph->GetXaxis()->SetMoreLogLabels();
 return graph;
}
