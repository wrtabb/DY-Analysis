#include "include/DYAnalyzer.hh"

void test()
{
 DYAnalyzer*dy = new DYAnalyzer();
 Long64_t totalentries = dy->LoadTrees();
 Long64_t nentries;
 TH1D*hMassHardProcess=new TH1D("hMassHardProcess","",nLogBins,massbins);
 TH1D*hMassReco=new TH1D("hMassReco","",nLogBins,massbins);
 TH1D*hist0=new TH1D("hist0","",nLogBins,massbins);
 TH1D*hist1=new TH1D("hist1","",nLogBins,massbins);
 TH1D*hist2=new TH1D("hist2","",nLogBins,massbins);

 double xSecWeight = 1.0;
 double genWeight = 1.0;
 double genWeightSum;
 double weight = 1.0;
 double weightReco = 1.0;
 double massHard;
 double massReco;
 bool useGenWeights = true;
 Long64_t count = 0;
 for(int iChain=0;iChain<numChains;iChain++){
  nentries = dy->GetDYEntries(iChain);

  //genWeightSum is the same for all events in a chain, but the gen weight must be calculated
  //Separately for each event, so the sum is calculated here, but the actual gen weight
  //is calculated within the event loop
  genWeightSum = dy->GetGenWeightSum(iChain);
  
  //xSecWeight is different if gen weights are being used
  xSecWeight = dy->GetXsecWeight(iChain,useGenWeights);

  Long64_t event;
  //-----Event loop-----//
  for(Long64_t i=0;i<nentries;i++){
   event = dy->GetDYEntry(iChain,i);
   genWeight = (GENEvt_weight/fabs(GENEvt_weight))/genWeightSum;
   //-----Indices for leptons-----//
   int iHard1 = -1;
   int iHard2 = -1;
   int iFSR1 = -1;
   int iFSR2 = -1;
   int nGenDileptons,nRecoElectrons;
   int leadEle = -1;
   int subEle = -1;
   int closestTrack1 = -1;
   int closestTrack2 = -1;

   dy->Counter(count,totalentries,"Looping over events: ");
   count++;

   //Function to find two leptons from hard process
   nGenDileptons = dy->GetGenLeptons(ELE,iHard1,iHard2,iFSR1,iFSR2);   

   //Function to find two reconstructed electrons
   //I need to redo this to make it more general so that one function can handle
   //electrons or muons
   nRecoElectrons = dy->GetRecoElectrons(leadEle,subEle);

   //Determines if the selected gen leptons pass acceptance cuts
   bool passAcceptance = dy->AcceptanceCut(GENLepton_pT[iHard1],GENLepton_pT[iHard2],
                                           GENLepton_eta[iHard2],GENLepton_eta[iHard2]);
   bool passRecoMatch1 = dy->GenToRecoMatchCut(iFSR1,closestTrack1);
   bool passRecoMatch2 = dy->GenToRecoMatchCut(iFSR2,closestTrack2);

   //Determine if both leptons pass reco gen match cut
   //So far this is only electrons
   bool passRecoMatch = passRecoMatch1 && passRecoMatch2;

   //Calculate invariant mass from hard process
   massHard = dy->CalcInvMass(GENLepton_pT[iHard1],GENLepton_eta[iHard1],
                              GENLepton_phi[iHard1],eMass,GENLepton_pT[iHard2],
                              GENLepton_eta[iHard2],GENLepton_phi[iHard2],eMass);
   //Calculate invariant mass for reconstructed electrons
   massReco = dy->CalcInvMass(Electron_pT[leadEle],Electron_eta[leadEle],
                              Electron_phi[leadEle],eMass,Electron_pT[subEle],
                              Electron_eta[subEle],Electron_phi[subEle],eMass);

   //-----Get weights-----//
   //The first argument in the function is whether the weight is for reco or not
   //This distinction is important because reco events get scale factor weighting
   //but gen level events do not
   weight = dy->GetTotalWeight(false,iChain,genWeight,xSecWeight,GENLepton_eta[iHard1],
                               GENLepton_eta[iHard2],GENLepton_pT[iHard1],
                               GENLepton_pT[iHard2]); 
   weightReco = dy->GetTotalWeight(true,iChain,genWeight,xSecWeight,GENLepton_eta[iHard1],
                                   GENLepton_eta[iHard2],GENLepton_pT[iHard1],
                                   GENLepton_pT[iHard2]); 

   //Fill histograms
   hMassHardProcess->Fill(massHard,weight);
   hMassReco->Fill(massReco,weightReco);

   //This structure is temporary with these generic looking names just to test how the cuts 
   //are working
   //hist0: no cuts
   //hist1: after acceptance cut
   //hist2: after reco gen matching cut
   hist0->Fill(massHard,weight);
   if(!passAcceptance) massHard = 0;
   hist1->Fill(massHard,weight);
   if(!passRecoMatch) massHard = 0;
   hist2->Fill(massHard,weight);

  }//end event loop
 }//end chain loop

 TCanvas*canvas = new TCanvas("camvas","",0,0,1200,1000);
 canvas->SetGrid();
 canvas->SetLogx();
 canvas->SetLogy();
 hMassHardProcess->GetXaxis()->SetNoExponent();
 hMassHardProcess->GetXaxis()->SetMoreLogLabels();
 hMassHardProcess->SetFillColor(kYellow-2);
 hMassHardProcess->Draw("hist");
 canvas->SaveAs("data/invMassHardProcess.png");

 TFile*histSave = new TFile("data/histograms.root","recreate");
 hist0->Write();
 hist1->Write();
 hist2->Write();
 hMassReco->Write();
 hMassHardProcess->Write();
 histSave->Close();

 dy->GetEfficiencies(hist0,hist1,"Acceptance");
 dy->GetEfficiencies(hist1,hist2,"RecoGenEfficiency");
}
