#include "DYSelection.hh"

void runDYSelection()
{
 NtupleVersion ntup = V2P3;
 LepType lep = ELE;
 SampleType sample = LL;
 DYSelection*dy = new DYSelection(ntup,lep,sample);
}
