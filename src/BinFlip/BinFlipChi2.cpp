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

BinFlipChi2::BinFlipChi2(BinFlipParSet* fitPars, const HadronicParameters& hadronicPars,
			 const TimeBinning& timeBins) :
  Minimisable(fitPars),
  m_fitPars(fitPars),
  m_hadronicPars(hadronicPars),
  m_timeBinning(timeBins),
  m_zcp(fitPars->zcp()),
  m_dz(fitPars->deltaz())
{
  init() ;
}

BinFlipChi2::BinFlipChi2(BinFlipParSet* fitPars, const string& hadronicParsName, const string& timeBinsName,
			 const string& fname) :
  Minimisable(fitPars),
  m_fitPars(fitPars),
  m_hadronicPars(hadronicParsName, fname),
  m_timeBinning(timeBinsName, fname),
  m_zcp(fitPars->zcp()),
  m_dz(fitPars->deltaz())
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

void BinFlipChi2::parametersChanged() {
  m_zcp = m_fitPars->zcp() ;
  m_dz = m_fitPars->deltaz() ;
}

double BinFlipChi2::getVal() {
  double chi2 = 0. ;
  for(unsigned iTimeBin = 0 ; iTimeBin < m_timeBinning.nBinsTime() ; ++iTimeBin){
    for(unsigned iPhaseBin = 1 ; iPhaseBin < m_timeBinning.nBinsPhase() ; ++iPhaseBin){
      // Could use the t & t2 for each phase space bin here.
      double Rplus = m_hadronicPars.bin(iPhaseBin).R(m_timeBinning.integratedBin(iTimeBin).t(),
						     m_timeBinning.integratedBin(iTimeBin).t2(),
						     m_zcp, m_dz) ;
      double Rminus = m_hadronicPars.bin(iPhaseBin).Rbar(m_timeBinning.integratedBin(iTimeBin).t(),
							 m_timeBinning.integratedBin(iTimeBin).t2(),
							 m_zcp, m_dz) ;
      chi2 += m_timeBinning.chiSquared(iTimeBin, iPhaseBin, Rplus, Rminus) ;
    }
  }
  return chi2 ;
}
