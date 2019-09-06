#include "../include/DYAnalyzer.hh"


void efficiency(LepType lepType)
{
 DYAnalyzer*drellYan = new DYAnalyzer();
 Long64_t totalentries = drellYan->LoadTrees();
 TH1D*hist = new TH1D("hMass","",nLogBins,massbins);
 hist->SetFillColor(kYellow+2);

 //-----Loop over samples-----//
 Long64_t count = 0;
 bool useGenWeights = true; //this determines how xSecWeights are calculated
 for(int iChain=0;iChain<numChains;iChain++){
  Long64_t nentries = drellYan->GetDYEntries(iChain);
  double genWeight = 1.0;
  if(useGenWeights) genWeight = drellYan->GetGenWeight(iChain);  
  double xSecWeight = drellYan->GetXsecWeight(iChain,useGenWeights);

  //-----Loop over events-----//
  for(Long64_t iEntry=0;iEntry<nentries;iEntry++){
   drellYan->GetDYEntry(iChain,iEntry);
   drellYan->Counter(count,totalentries,"Looping over events: ");
   count++;

   int idxHardEle1 = -1;
   int idxHardEle2 = -1;
   int idxFSREle1  = -1;
   int idxFSREle2  = -1;
   int nDileptons = 
    drellYan->GetGenLeptons(lepType,idxHardEle1,idxHardEle2,idxFSREle1,idxFSREle2);
   double invMassHard = 0;
   double totalWeight = 1.0;

   //-----If dileptons found, put in vectors and calculate invariant mass-----//
   if(nDileptons==1){
    TLorentzVector vHardEle1;
    TLorentzVector vHardEle2;
    vHardEle1.SetPtEtaPhiM(GENLepton_pT[idxHardEle1],GENLepton_eta[idxHardEle1],
                           GENLepton_phi[idxHardEle1],eMass);
    vHardEle2.SetPtEtaPhiM(GENLepton_pT[idxHardEle2],GENLepton_eta[idxHardEle2],
                           GENLepton_phi[idxHardEle2],eMass);
    invMassHard = drellYan->CalcInvMass(vHardEle1,vHardEle2);
    totalWeight = drellYan->AddWeights(iChain,genWeight,xSecWeight,GENLepton_eta[idxHardEle1],
                                       GENLepton_eta[idxHardEle2],GENLepton_pT[idxHardEle1],
                                       GENLepton_pT[idxHardEle2]);
    hist->Fill(invMassHard,totalWeight);
   }
  }//end loop over events
 }//end loop over chains
 TCanvas*canvas = new TCanvas("canvas","",0,0,1200,1000);
 canvas->SetLogx();
 canvas->SetLogy();
 canvas->SetGrid();
 hist->Draw("hist");
}//end efficiency()


