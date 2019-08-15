#include <Mint/BinFlipChi2.h>

using namespace std ;
using MINT::FitParameter ;

BinFlipParSet::BinFlipParSet(double _zcp_Re, double zcp_Re_step, double _zcp_Im, double zcp_Im_step,
			     double _deltaz_Re, double deltaz_Re_step, double _deltaz_Im, double deltaz_Im_step) :
  zcp_Re("zcp_Re", FitParameter::FIT, _zcp_Re, zcp_Re_step, 0., 0., this),
  zcp_Im("zcp_Im", FitParameter::FIT, _zcp_Im, zcp_Im_step, 0., 0., this),
  deltaz_Re("deltaz_Re", FitParameter::FIT, _deltaz_Re, deltaz_Re_step, 0., 0., this),
  deltaz_Im("deltaz_Im", FitParameter::FIT, _deltaz_Im, deltaz_Im_step, 0., 0., this)
{
}

BinFlipParSet::BinFlipParSet(const string& fname) :
  zcp_Re("zcp_Re", fname.c_str(), this, FitParameter::FIT),
  zcp_Im("zcp_Im", fname.c_str(), this, FitParameter::FIT),
  deltaz_Re("deltaz_Re", fname.c_str(), this, FitParameter::FIT),
  deltaz_Im("deltaz_Im", fname.c_str(), this, FitParameter::FIT)
{}

void BinFlipParSet::fixZeroCPV() {
  deltaz_Re.setCurrentFitVal(0.) ;
  deltaz_Re.fixAndHide() ;
  deltaz_Im.setCurrentFitVal(0.) ;
  deltaz_Im.fixAndHide() ;
}

complex<double> BinFlipParSet::zcp() const {
  return complex<double>(zcp_Re.getCurrentFitVal(), zcp_Im.getCurrentFitVal()) ;
}

complex<double> BinFlipParSet::deltaz() const {
  return complex<double>(deltaz_Re.getCurrentFitVal(), deltaz_Im.getCurrentFitVal()) ;
}

pair<complex<double>, complex<double> > BinFlipParSet::fromXY(double x, double y, double magqoverp, double phi) {
  complex<double> z(-y, -x) ;
  complex<double> qoverp = polar(magqoverp, phi) ;
  complex<double> poverq = pow(qoverp, -1) ;
  complex<double> zcp = z * (qoverp + poverq) / 2. ;
  complex<double> deltaz = z * (qoverp - poverq) / 2. ;
  return make_pair(zcp, deltaz) ;
}

BinFlipChi2Base::BinFlipChi2Base(BinFlipParSet* fitPars) :
  Minimisable(fitPars),
  m_fitPars(fitPars),
  m_zcp(fitPars->zcp()),
  m_dz(fitPars->deltaz())
{}

void BinFlipChi2Base::parametersChanged() {
  m_zcp = m_fitPars->zcp() ;
  m_dz = m_fitPars->deltaz() ;
}

TH2D BinFlipChi2Base::scan2D(unsigned ip1, unsigned ip2, float nSigmaRange, unsigned nBins, bool zeroCentre,
			     float scale) {
  FitParameter* p1 = (FitParameter*)m_fitPars->getParPtr(ip1) ;
  FitParameter* p2 = (FitParameter*)m_fitPars->getParPtr(ip2) ;
  const double v1 = p1->getCurrentFitVal() ;
  const double v2 = p2->getCurrentFitVal() ;
  const double s1 = p1->err() ;
  const double s2 = p2->err() ;
  double v1min = v1 - s1 * nSigmaRange ;
  double v1max = v1 + s1 * nSigmaRange ;
  double v2min = v2 - s2 * nSigmaRange ;
  double v2max = v2 + s2 * nSigmaRange ;
  if(zeroCentre){
    v1min -= v1 ;
    v1max -= v1 ;
    v2min -= v2 ;
    v2max -= v2 ;
  }
  TH2D hscan(("deltachi2_" + p2->name() + "_vs_" + p1->name()).c_str(), "",
	     nBins, v1min*scale, v1max*scale, nBins, v2min*scale, v2max*scale) ;
  const double chi2 = getVal() ;
  for(unsigned i1 = 1 ; i1 < nBins + 1 ; ++i1){
    double iv1 = hscan.GetXaxis()->GetBinCenter(i1) ;
    for(unsigned i2 = 1 ; i2 < nBins + 1 ; ++i2){
      double iv2 = hscan.GetYaxis()->GetBinCenter(i2) ;
      p1->setCurrentFitVal(iv1) ;
      p2->setCurrentFitVal(iv2) ;
      parametersChanged() ;
      double ichi2 = getVal() ;
      double dsigma = sqrt(ichi2 - chi2) ;
      if(dsigma <= nSigmaRange)
	hscan.SetBinContent(i1, i2, dsigma) ;
      else
	hscan.SetBinContent(i1, i2, 0.) ;
    }
  }

  vector<string> titles(1, "Re(z_{CP})") ;
  titles.push_back("Im(z_{CP})") ;
  titles.push_back("Re(#Delta z)") ;
  titles.push_back("Im(#Delta z)") ;
  string xtitle = titles.at(ip1) ;
  string ytitle = titles.at(ip2) ;
  if(zeroCentre){
    xtitle = "#delta " + xtitle ;
    ytitle = "#delta " + ytitle ;
  }
  hscan.SetXTitle(xtitle.c_str()) ;
  hscan.SetYTitle(ytitle.c_str()) ;
  hscan.SetZTitle("N. #sigma") ;
  hscan.SetContour(nSigmaRange) ;
  hscan.SetStats(false) ;

  p1->setCurrentFitVal(v1) ;
  p2->setCurrentFitVal(v2) ;
  parametersChanged() ;
  return hscan ;
}

BinFlipParSet* BinFlipChi2Base::getBinFlipParSet() {
  return m_fitPars ;
}

BinFlipChi2::BinFlipChi2(BinFlipParSet* fitPars, double lifetime, const HadronicParameters& hadronicPars,
			 const TimeBinning& timeBins) :
  BinFlipChi2Base(fitPars),
  m_hadronicPars(hadronicPars),
  m_timeBinning(timeBins),
  m_lifetime(lifetime)
{
  init() ;
}

BinFlipChi2::BinFlipChi2(BinFlipParSet* fitPars, double lifetime,
			 const string& hadronicParsName, const string& timeBinsName,
			 const string& fname) :
  BinFlipChi2Base(fitPars),
  m_hadronicPars(hadronicParsName, fname),
  m_timeBinning(timeBinsName, fname),
  m_lifetime(lifetime)
{
  init() ;
}

void BinFlipChi2::init() {
  checkPhaseBinning() ;
}

void BinFlipChi2::checkPhaseBinning() const {
  if(m_hadronicPars.binning() != *m_timeBinning.phaseBinning())
    throw invalid_argument("The phase binning of the hadronic parameters and binned data don't match!") ;
}

double BinFlipChi2::getVal() {
  double chi2 = 0. ;
  for(unsigned iTimeBin = 0 ; iTimeBin < m_timeBinning.nBinsTime() ; ++iTimeBin){
    for(unsigned iPhaseBin = 1 ; iPhaseBin < m_timeBinning.nBinsPhase() + 1; ++iPhaseBin){
      // Could use the t & t2 for each phase space bin here.
      double Rplus = m_hadronicPars.bin(iPhaseBin).R(m_timeBinning.integratedBin(iTimeBin).t(),
						     m_timeBinning.integratedBin(iTimeBin).t2(),
						     m_lifetime, m_zcp, m_dz) ;
      double Rminus = m_hadronicPars.bin(iPhaseBin).Rbar(m_timeBinning.integratedBin(iTimeBin).t(),
							 m_timeBinning.integratedBin(iTimeBin).t2(),
							 m_lifetime, m_zcp, m_dz) ;
      chi2 += m_timeBinning.chiSquared(iTimeBin, iPhaseBin, Rplus, Rminus) ;
    }
  }
  return chi2 ;
}

BinFlipChi2Simul::BinFlipChi2Simul(BinFlipChi2* flip1, BinFlipChi2* flip2) :
  BinFlipChi2Base((BinFlipParSet*)flip1->getParSet()),
  m_flippers(1, flip1)
{
  m_flippers.push_back(flip2) ;
}

BinFlipChi2Simul::BinFlipChi2Simul(const BinFlipChi2Simul::BinFlippers& flippers) :
  BinFlipChi2Base((BinFlipParSet*)flippers.at(0)->getParSet()),
  m_flippers(flippers)
{}

double BinFlipChi2Simul::getVal() {
  double chi2 = 0. ;
  for(auto* flipper : m_flippers)
    chi2 += flipper->getVal() ;
  return chi2 ;
}

void BinFlipChi2Simul::parametersChanged() {
  BinFlipChi2Base::parametersChanged() ;
  for(auto* flipper : m_flippers)
    flipper->parametersChanged() ;
}
