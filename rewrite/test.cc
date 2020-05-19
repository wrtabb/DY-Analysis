#include "include/DYAnalyzer.hh"

void test()
{
 //initialize DYAnalyzer class
 DYAnalyzer*dy = new DYAnalyzer(V2P6,ELE,EE);
 //Function for loading trees used for analysis
 //The trees being loaded are specified in src/DYAnalysis.cc
 //May change structure, see notes in src/DYAnalysis.cc
 //Function also returns total number of events for all trees
 Long64_t totalentries = dy->LoadTrees(V2P6,dirNames,EE,ELE);
 //nentries will be number of entries in a given sample
 Long64_t nentries;
 //Define histograms
 TH1D*hHardProcess[numChains];
 TString histbasename = "hHardProcess";
 TString histname;
 for(int jChain=0;jChain<numChains;jChain++){
  histname = histbasename;
  histname+=jChain;
  hHardProcess[jChain] = dy->DefineMassHist(LINEAR,histname);
  hHardProcess[jChain]->SetFillColor(jChain+1);
  hHardProcess[jChain]->GetXaxis()->
   SetTitle("Gen-Level Dielectron Mass (isHardProcess) [GeV]"); 
  hHardProcess[jChain]->GetYaxis()->SetRangeUser(0.000001,1000000000);
  hHardProcess[jChain]->GetXaxis()->SetNoExponent();
  hHardProcess[jChain]->GetXaxis()->SetMoreLogLabels();
 }//end loop to define hard process histograms

 TH1D*hMassReco = dy->DefineMassHist(LOG,"hMassReco");
 TH1D*hMassFSR = dy->DefineMassHist(LOG,"hMassFSR");

 TH1D*hist0 = dy->DefineMassHist(LOG,"hist0");
 TH1D*hist1 = dy->DefineMassHist(LOG,"hist1");
 TH1D*hist2 = dy->DefineMassHist(LOG,"hist2");
 TH1D*hist3 = dy->DefineMassHist(LOG,"hist3");
 TH1D*hist4 = dy->DefineMassHist(LOG,"hist4");

 //Weighting factors are all initialized to 1
 //If a weighting factor is to be used, it is redefined later
 double xSecWeight = 1.0;
 double genWeight = 1.0;
 double weight = 1.0;
 double weightReco = 1.0;

 double genWeightSum;
 

 //Specify lepton type
 LepType lepType = ELE;
 double lepMass = -5000;
 if(lepType == ELE) lepMass = eMass;
 else if(lepType == MUON) lepMass = muMass;
 else if(lepType == TAU) lepMass = tauMass;
 else {
  cout << "ERROR: Lepton type not correctly specified!" << endl;
 }

 //Specify if gen weights are to be used
 //cross section weighting is different if gen weights are used
 bool useGenWeights = true;
 Long64_t count = 0;

 TCanvas*canvas1 = new TCanvas("cInvMassHardProcess","",10,10,1000,1000);
 canvas1->SetLogx();
 canvas1->SetLogy();
 canvas1->cd();

 //Begin looping over samples
 //Each separate sample gets placed in an array of chains
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
   //get event i from chain chains[iChain]
   event = dy->GetDYEntry(iChain,i);
   double massHard = -1;
   double massFSR = -1;
   double massReco = -1;
   
   //gen weights are calculated
   genWeight = (GENEvt_weight/fabs(GENEvt_weight))/genWeightSum;

   //-----Indices for leptons-----//
   //hard process leptons
   int iHard1 = -1;
   int iHard2 = -1;
   //Final state radiation leptons
   int iFSR1 = -1;
   int iFSR2 = -1;
   //leading and subleading leptons
   int leadEle = -1;
   int subEle = -1;
   //These are used in the gen to reco match cut
   int closestTrack1 = -1;
   int closestTrack2 = -1;
   //counters for numbers of leptons
   int nGenDileptons,nRecoElectrons;
   
   //Simple counting function to keep track of program progress as it is running
   dy->Counter(count,totalentries,"Looping over events: ");
   count++;

   //Function to find two leptons from hard process
   nGenDileptons = dy->GetGenLeptons(lepType,iHard1,iHard2,iFSR1,iFSR2);   
   if(nGenDileptons!=1) continue;

   //Function to find two reconstructed electrons
   //I need to redo this to make it more general so that one function can handle
   //electrons or muons
   nRecoElectrons = dy->GetRecoElectrons(leadEle,subEle);

   //Determines if the selected gen leptons pass acceptance cuts
   bool passAcceptance = dy->AcceptanceCut(GENLepton_pT[iHard1],GENLepton_pT[iHard2],
                                           GENLepton_eta[iHard1],GENLepton_eta[iHard2]);

   //Takes the final state leptons and tries to find the closest reconstructed lepton
   //If there is not a reco lepton close enough to the fsr lepton, passRecoMatch = false
   bool passRecoMatch1 = dy->GenToRecoMatchCut(iFSR1,closestTrack1);
   bool passRecoMatch2 = dy->GenToRecoMatchCut(iFSR2,closestTrack2);

   //Do reco electrons which match up with gen electrons pass medium ID cuts?
   bool passMediumID = dy->MediumIDCut(Electron_passMediumID[closestTrack1],
                                       Electron_passMediumID[closestTrack2]);

   //Determine if both leptons pass reco gen match cut
   //It only passes the cut if both pass
   bool passRecoMatch = passRecoMatch1 && passRecoMatch2;

   //Determine if event passes HLT cut
   bool passHLT = dy->HLTCut();

   //Calculate invariant mass from hard process
   if(iHard1>=0 && iHard2>=0){
    massHard = dy->CalcInvMass(GENLepton_pT[iHard1],GENLepton_eta[iHard1],
                               GENLepton_phi[iHard1],lepMass,GENLepton_pT[iHard2],
                               GENLepton_eta[iHard2],GENLepton_phi[iHard2],lepMass);
   }
   //Calculate invariant mass from FSR
   if(iFSR1>=0 && iFSR2>=0){
    massFSR = dy->CalcInvMass(GENLepton_pT[iFSR1],GENLepton_eta[iFSR1],
                               GENLepton_phi[iFSR1],lepMass,GENLepton_pT[iFSR2],
                               GENLepton_eta[iFSR2],GENLepton_phi[iFSR2],lepMass);
   }
   //Calculate invariant mass for reconstructed electrons
   if(leadEle>=0 && subEle>=0){
    massReco = dy->CalcInvMass(Electron_pT[leadEle],Electron_eta[leadEle],
                               Electron_phi[leadEle],lepMass,Electron_pT[subEle],
                               Electron_eta[subEle],Electron_phi[subEle],lepMass);
   }
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
   if(massReco>=0) hMassReco->Fill(massReco,weightReco);
   if(massFSR>=0)  hMassFSR->Fill(massFSR,weight);
   if(massHard>=0) hHardProcess[iChain]->Fill(massHard,weight);

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
   if(!passMediumID) massHard = 0;
   hist3->Fill(massHard,weight);
   if(!passHLT) massHard = 0;
   hist4->Fill(massHard,weight);
  }//end event loop
  if(iChain==0) hHardProcess[iChain]->Draw("Bar");      
  else hHardProcess[iChain]->Draw("Barsame");
 }//end chain loop

 canvas1->SaveAs("plots/hardProcessBySample.png");
 TFile*histSave = new TFile("data/histograms.root","recreate");
 hist0->Write();
 hist1->Write();
 hist2->Write();
 hist3->Write();
 canvas1->Write();
 hMassReco->Write();
 hMassFSR->Write();
 histSave->Close();

 dy->GetEfficiencies(hist0,hist1,"Acceptance");
 dy->GetEfficiencies(hist1,hist2,"RecoGenEfficiency");
 dy->GetEfficiencies(hist2,hist3,"MediumIDEfficiency");
 dy->GetEfficiencies(hist3,hist4,"HLTEfficiency");
}
