


void testUnfold1()
{
  TFile *signalFile=new TFile("testUnfold_signal.root","recreate");
  TTree *signalTree=new TTree("signal","event");
  signalTree->Branch("etarec",&etaRec,"etarec/F");
  signalTree->Branch("ptrec",&ptRec,"ptrec/F");
  signalTree->Branch("discr",&discr,"discr/F");
  signalTree->Branch("istriggered",&istriggered,"istriggered/I");
  signalTree->Branch("etagen",&etaGen,"etagen/F");
  signalTree->Branch("ptgen",&ptGen,"ptgen/F");
  cout<<"fill signal tree\n";


}
