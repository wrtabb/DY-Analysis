#include "../include/DYAnalyzer.hh"

void efficiency()
{
 DYAnalyzer*drellYan = new DYAnalyzer();
 Long64_t totalentries = drellYan->LoadTrees();
 drellYan->InitBranches();
 TH1D*hHard=new TH1D("hHard","",nLogBins,massbins);
 
 //-----Set lepton type-----//
 LepType lepType = ELE;

 Long64_t count = 0;
 //-----Loop over samples-----//
 for(int iChain=0;iChain<numChains;iChain++){
  Long64_t nentries = drellYan->GetDYEntries(iChain);
  
  //-----Loop over events-----//
  for(Long64_t iEntry=0;iEntry<nentries;iEntry++){
   drellYan->GetDYEntry(iChain,iEntry);
   drellYan->Counter(count,totalentries,"Looping over events: ");
   count++;

   int idxHardEle1 = -1;
   int idxHardEle2 = -1;
   int idxFSREle1  = -1;
   int idxFSREle2  = -1;
   int nDielectrons = 
    drellYan->GetGenLeptons(lepType,idxHardEle1,idxHardEle2,idxFSREle1,idxFSREle2);
   double invMassHard = 0;

   //-----If dileptons found, put in vectors and calculate invariant mass-----//
   if(nDielectrons==1){
    TLorentzVector vHardEle1;
    TLorentzVector vHardEle2;
    vHardEle1.SetPtEtaPhiM(GENLepton_pT[idxHardEle1],GENLepton_eta[idxHardEle1],
                       GENLepton_phi[idxHardEle1],eMass);
    vHardEle2.SetPtEtaPhiM(GENLepton_pT[idxHardEle2],GENLepton_eta[idxHardEle2],
                       GENLepton_phi[idxHardEle2],eMass);
    invMassHard = drellYan->CalcInvMass(vHardEle1,vHardEle2);
   }
   hHard->Fill(invMassHard);
  }//end loop over events
 }//end loop over chains
 hHard->Draw("hist");
}//end efficiency()


