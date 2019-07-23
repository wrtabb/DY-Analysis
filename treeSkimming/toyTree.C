
void makeToy();
void readToy();
void skimToy();
void setHistOptions(TH1D*hist,int color);

const TString fileName = "toyTree.root";
const TString writeFileName = "toyTreeSkim.root";
const TString treeName = "toyTree";
const int nEvents = 1e4;

void toyTree()
{
 //makeToy(); 
 //readToy();
 skimToy();
}

void makeToy()
{
 TFile*file = new TFile(fileName,"recreate");
 TTree*tree = new TTree(treeName,"Toy Tree");

 float branch1,branch2,branch3;
 tree->Branch("branch1",&branch1,"branch1/F");
 tree->Branch("branch2",&branch2,"branch2/F");
 tree->Branch("branch3",&branch3,"branch3/F");

 TRandom3 gen;
 for(int i=0;i<nEvents;i++){
  branch1 = gen.Gaus(5,2);
  branch2 = gen.Gaus(10,3);
  branch3 = gen.Gaus(9,4);
  tree->Fill();
 }

 tree->Write();
 file->Close();
}

void readToy()
{
 TFile*openFile = new TFile(fileName);
 TTree*openTree = (TTree*)openFile->Get(treeName);
 TBranch*b_branch1;
 TBranch*b_branch2;
 TBranch*b_branch3;

 float branch1,branch2,branch3;
 openTree->SetBranchAddress("branch1",&branch1,&b_branch1);
 openTree->SetBranchAddress("branch2",&branch2,&b_branch2);
 openTree->SetBranchAddress("branch3",&branch3,&b_branch3);

 TH1D*hBranch1 = new TH1D("hBranch1","",100,0,20);
 TH1D*hBranch2 = new TH1D("hBranch2","",100,0,20);
 TH1D*hBranch3 = new TH1D("hBranch3","",100,0,20);
 setHistOptions(hBranch1,2);
 setHistOptions(hBranch2,4);
 setHistOptions(hBranch3,8);

 int nentries = openTree->GetEntries();
 for(int i=0;i<nentries;i++){
  openTree->GetEntry(i);
  hBranch1->Fill(branch1);
  hBranch2->Fill(branch2);
  hBranch3->Fill(branch3); 
 }
 TCanvas*canvas = new TCanvas("canvas","",0,0,1000,1000);
 canvas->SetGrid();
 hBranch1->Draw("pe");
 hBranch2->Draw("pe,same");
 hBranch3->Draw("pe,same");
}

void skimToy()
{
 TFile*openFile = new TFile(fileName);
 TTree*openTree = (TTree*)openFile->Get(treeName);
 TBranch*b_branch1;
 TBranch*b_branch2;
 TBranch*b_branch3;

 float branch1,branch2,branch3;
 openTree->SetBranchAddress("branch3",&branch3,&b_branch3);
 openTree->SetBranchAddress("branch1",&branch1,&b_branch1);
 openTree->SetBranchAddress("branch2",&branch2,&b_branch2);

 TFile*writeFile = new TFile(writeFileName,"recreate");
 TTree*writeTree = openTree->CloneTree(0); 


 int nentries = openTree->GetEntries();
 for(int i=0;i<nentries;i++){
  openTree->GetEntry(i);
  if(branch2>15||branch1>15) continue;
  writeTree->Fill();
 }
 TBranch*b = writeTree->GetBranch("branch3"); 
 writeTree->GetListOfBranches()->Remove(b);
 TLeaf*l = writeTree->GetLeaf("branch3");
 writeTree->GetListOfLeaves()->Remove(l);
 writeTree->Write();
 writeFile->Close();
}

void setHistOptions(TH1D*hist,int color)
{
 hist->SetMarkerStyle(20);
 hist->SetLineColor(color);
 hist->SetMarkerColor(color);
}
