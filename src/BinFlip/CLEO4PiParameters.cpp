#include <Mint/CLEO4PiParameters.h>
#include <FourPiUtils.h>

using namespace std;

CLEO4PiBinning::CLEO4PiBinning(const HyperHistogram& hist):
  m_hist(hist)
{}

string CLEO4PiBinning::type() const {
  return "CLEO4PiBinning";
}

HadronicParameters::EventBinInfo CLEO4PiBinning::binInfo(IDalitzEvent& evt) const {
  HadronicParameters::EventBinInfo binNo;
  binNo.evt = &evt;
  
}
