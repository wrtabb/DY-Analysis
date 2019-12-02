#include "include/DYAnalyzer.hh"

void test()
{
 DYAnalyzer*dy = new DYAnalyzer();
 Long64_t totalentries = dy->LoadTrees();
 Long64_t nentries;
 TH1D*hMassHardProcess=new TH1D("hist","",nLogBins,massbins);
 TH1D*hMass0 = new TH1D("hMassNoCuts","",nLogBins,massbins);
 TH1D*hMass1 = new TH1D("hMassInAcceptance","",nLogBins,massbins);
 TH1D*hMass2 = new TH1D();

 double xSecWeight = 1.0;
 double genWeight = 1.0;
 double weight = 1.0;
 double mass;
  bool useGenWeights = true;
 Long64_t count = 0;
 ofstream saveWeight("data/genWeights.txt");
 for(int iChain=0;iChain<numChains;iChain++){
  nentries = dy->GetDYEntries(iChain);
  //nentries = 1000; 
  //-----Chain-level weights-----//
  //genWeight = dy->GetGenWeight(iChain);
  
  //cross section weights are different if gen weights are being used
  if(genWeight==1.0) useGenWeights = false;
  xSecWeight = dy->GetXsecWeight(iChain,useGenWeights);
  saveWeight << genWeight << endl;

  Long64_t event;
  //-----Event loop-----//
  for(Long64_t i=0;i<nentries;i++){
   event = dy->GetDYEntry(iChain,i);
   //-----Indices for leptons-----//
   int iHard1 = -1;
   int iHard2 = -1;
   int iFSR1 = -1;
   int iFSR2 = -1;
   int nDileptons = -1;
   dy->Counter(count,totalentries,"Looping over events: ");
   count++;
   nDileptons = dy->GetGenLeptons(ELE,iHard1,iHard2,iFSR1,iFSR2);   
   //if(nDileptons!=1) continue;

   bool passAcceptance = dy->AcceptanceCut(GENLepton_pT[iHard1],GENLepton_pT[iHard2],
                                           GENLepton_eta[iHard2],GENLepton_eta[iHard2]);
   TLorentzVector vHardEle1;
   TLorentzVector vHardEle2;
   vHardEle1.SetPtEtaPhiM(GENLepton_pT[iHard1],GENLepton_eta[iHard1],
                          GENLepton_phi[iHard1],eMass);
   vHardEle2.SetPtEtaPhiM(GENLepton_pT[iHard2],GENLepton_eta[iHard2],
                          GENLepton_phi[iHard2],eMass);
   mass = dy->CalcInvMass(vHardEle1,vHardEle2);
   if(mass < 0){
    cout << "Mass < 0: (idx1,idx2) = " << iHard1 << ", " << iHard2 << 
     ", event = " << event << endl;
   continue;
   }

   //-----Get event-level weights-----//
   weight = dy->GetTotalWeight(iChain,genWeight,xSecWeight,GENLepton_eta[iHard1],
                           GENLepton_eta[iHard2],GENLepton_pT[iHard1],GENLepton_pT[iHard2]); 
   hMass0->Fill(mass,weight);
   if(!passAcceptance) mass = 0;
   hMass1->Fill(mass,weight);
   hMassHardProcess->Fill(mass,weight);
  }//end event loop
 }//end chain loop

 saveWeight.close();
 TCanvas*canvas = new TCanvas("camvas","",0,0,1200,1000);
 canvas->SetGrid();
 canvas->SetLogx();
 canvas->SetLogy();
 hMassHardProcess->GetXaxis()->SetNoExponent();
 hMassHardProcess->GetXaxis()->SetMoreLogLabels();
 hMassHardProcess->SetFillColor(kYellow-2);
 hMassHardProcess->Draw("hist");
 dy->GetEfficiencies(hMass0,hMass1);
 canvas->SaveAs("data/invMassHardProcess.png");
}
