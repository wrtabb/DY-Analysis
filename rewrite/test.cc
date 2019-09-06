#include "include/DYAnalyzer.hh"

void test()
{
 DYAnalyzer*dy = new DYAnalyzer();
 Long64_t totalentries = dy->LoadTrees();
 Long64_t nentries;
 TH1D*hist=new TH1D("hist","",nLogBins,massbins);
 double xSecWeight = 1.0;
 double genWeight = 1.0;
 double weight = 1.0;
 double mass;
 Long64_t count = 0;
 for(int iChain=0;iChain<numChains;iChain++){
  nentries = dy->GetDYEntries(iChain);
  
  //-----Chain-level weights-----//
  xSecWeight = dy->GetXsecWeight(iChain,true);
  genWeight = dy->GetGenWeight(iChain);

  //-----Indices for leptons-----//
  int iHard1 = -1;
  int iHard2 = -1;
  int iFSR1 = -1;
  int iFSR2 = -1;
  int nDileptons = -1;

  //-----Event loop-----//
  for(Long64_t i=0;i<nentries;i++){
   dy->GetDYEntry(iChain,i);
   dy->Counter(count,totalentries,"Looping over events: ");
   count++;
   //if(count%100==0) cout << count << endl;
   nDileptons = dy->GetGenLeptons(ELE,iHard1,iHard2,iFSR1,iFSR2);   
   if(nDileptons!=1) continue;

   TLorentzVector vHardEle1;
   TLorentzVector vHardEle2;
   vHardEle1.SetPtEtaPhiM(GENLepton_pT[iHard1],GENLepton_eta[iHard1],
                          GENLepton_phi[iHard1],eMass);
   vHardEle2.SetPtEtaPhiM(GENLepton_pT[iHard2],GENLepton_eta[iHard2],
                          GENLepton_phi[iHard2],eMass);
   mass = dy->CalcInvMass(vHardEle1,vHardEle2);
  
   //-----Get event-level weights-----//
   weight = dy->GetTotalWeight(iChain,genWeight,xSecWeight,GENLepton_eta[iHard1],
                           GENLepton_eta[iHard2],GENLepton_pT[iHard1],GENLepton_pT[iHard2]); 
   hist->Fill(mass,weight);
  }//end event loop
 }//end chain loop
hist->Draw("hist");
 
}
