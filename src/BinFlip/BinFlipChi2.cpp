#include <Mint/BinFlipChi2.h>
#include <Mint/NamedParameterBase.h>

using namespace std ;
using MINT::FitParameter ;
using MINT::NamedParameterBase ;

vector<double> _blindingPars(unsigned long seed, double range, unsigned long offset){
  if(seed == 0)
    return vector<double>();
  vector<double> pars(3, seed + offset);
  pars[1] = -range;
  pars[2] = range;
  return pars;
}

BinFlipParSet::BinFlipParSet(double _zcp_Re, double zcp_Re_step, double _zcp_Im, double zcp_Im_step,
			     double _deltaz_Re, double deltaz_Re_step, double _deltaz_Im, double deltaz_Im_step,
			     unsigned long blindingSeed, double zBlindRange, double deltazBlindRange) :
  zcp_Re("zcp_Re", FitParameter::FIT, _zcp_Re, zcp_Re_step, 0., 0., this, 
	 NamedParameterBase::VERBOSE, nullptr, _blindingPars(blindingSeed, zBlindRange, 0)),
  zcp_Im("zcp_Im", FitParameter::FIT, _zcp_Im, zcp_Im_step, 0., 0., this, 
	 NamedParameterBase::VERBOSE, nullptr, _blindingPars(blindingSeed, zBlindRange, 1)),
  deltaz_Re("deltaz_Re", FitParameter::FIT, _deltaz_Re, deltaz_Re_step, 0., 0., this, 
	    NamedParameterBase::VERBOSE, nullptr, _blindingPars(blindingSeed, deltazBlindRange, 2)),
  deltaz_Im("deltaz_Im", FitParameter::FIT, _deltaz_Im, deltaz_Im_step, 0., 0., this, 
	    NamedParameterBase::VERBOSE, nullptr, _blindingPars(blindingSeed, deltazBlindRange, 3))
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

void BinFlipParSet::floatOnlyDeltaY() {
  zcp_Re.fixAndHide();
  zcp_Im.fixAndHide();
  deltaz_Im.fixAndHide();
}

complex<double> BinFlipParSet::zcp() const {
  return complex<double>(zcp_Re.blindedMean(), zcp_Im.blindedMean()) ;
}

complex<double> BinFlipParSet::deltaz() const {
  return complex<double>(deltaz_Re.blindedMean(), deltaz_Im.blindedMean()) ;
}

complex<double> BinFlipParSet::unblindZcp() const {
  return complex<double>(zcp_Re.getCurrentFitVal(), zcp_Im.getCurrentFitVal()) ;
}

complex<double> BinFlipParSet::unblindDeltaz() const {
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
  m_zcp(fitPars->unblindZcp()),
  m_dz(fitPars->unblindDeltaz())
{}

void BinFlipChi2Base::parametersChanged() {
  m_zcp = m_fitPars->unblindZcp() ;
  m_dz = m_fitPars->unblindDeltaz() ;
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

BinFlipChi2::BinFlipChi2(BinFlipParSet* fitPars, const HadronicParameters& hadronicPars,
			 const TimeBinning& timeBins) :
  BinFlipChi2Base(fitPars),
  m_hadronicPars(hadronicPars),
  m_timeBinning(timeBins)
{
  init() ;
}

BinFlipChi2::BinFlipChi2(BinFlipParSet* fitPars,
			 const string& hadronicParsName, const string& timeBinsName,
			 const string& fname) :
  BinFlipChi2Base(fitPars),
  m_hadronicPars(hadronicParsName, fname),
  m_timeBinning(timeBinsName, fname)
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

double BinFlipChi2::binChiSquared(unsigned iTimeBin, unsigned iPhaseBin) const {
  // Needs the unmixed time.
  double Rplus = m_hadronicPars.bin(iPhaseBin).R(m_timeBinning.meanUnmixedTime(iTimeBin),
						 m_timeBinning.meanUnmixedTime2(iTimeBin),
						 m_timeBinning.getLifetime(), m_zcp, m_dz) ;
  double Rminus = m_hadronicPars.bin(iPhaseBin).Rbar(m_timeBinning.meanUnmixedTime(iTimeBin),
						     m_timeBinning.meanUnmixedTime2(iTimeBin),
						     m_timeBinning.getLifetime(), m_zcp, m_dz) ;
  return m_timeBinning.chiSquared(iTimeBin, iPhaseBin, Rplus, Rminus) ;
}

double BinFlipChi2::getVal() {
  double chi2 = 0. ;
  for(unsigned iTimeBin = 0 ; iTimeBin < m_timeBinning.nBinsTime() ; ++iTimeBin){
    for(unsigned iPhaseBin = 1 ; iPhaseBin < m_timeBinning.nBinsPhase() + 1; ++iPhaseBin){
      chi2 += binChiSquared(iTimeBin, iPhaseBin);
    }
  }
  return chi2 ;
}

deque<deque<TH1F> > BinFlipChi2::plotsVsTime(const string& name) const {
  deque<deque<TH1F> > histos = m_timeBinning.plotsVsTime(name) ;
  unsigned iPhaseBin = 1 ;
  // Loop over histos for each phase bin.
  for(auto& binhistos : histos) {
    // Histos are Nbar+, Nbar-, Rbar, N+, N-, R, R-Rbar, ACP.
    TH1F& hrminus = binhistos[2] ;
    TH1F& hrplus = binhistos[5] ;
    TH1F& hrdiff = binhistos[6] ;
    TH1F& hacp = binhistos[7] ;
    TH1F rminusfit = TH1F(hrminus) ;
    rminusfit.SetName((string(hrminus.GetName()) + "_fit").c_str()) ;
    TH1F rplusfit = TH1F(hrplus) ;
    rplusfit.SetName((string(hrplus.GetName()) + "_fit").c_str()) ;
    TH1F rdifffit = TH1F(hrdiff) ;
    rdifffit.SetName((string(hrdiff.GetName()) + "_fit").c_str()) ;
    TH1F acpfit = TH1F(hacp);
    acpfit.SetName((string(hacp.GetName()) + "_fit").c_str());
    for(unsigned iTimeBin = 0 ; iTimeBin < m_timeBinning.nBinsTime() ; ++iTimeBin){
      double t = m_timeBinning.meanUnmixedTime(iTimeBin) ;
      double t2 = m_timeBinning.meanUnmixedTime2(iTimeBin) ;
      double rplus = m_hadronicPars.bin(iPhaseBin).R(t, t2, m_timeBinning.getLifetime(), m_zcp, m_dz) ;
      double rminus = m_hadronicPars.bin(iPhaseBin).Rbar(t, t2, m_timeBinning.getLifetime(), m_zcp, m_dz) ;
      double rdiff = rplus-rminus ;
      double acp = m_hadronicPars.bin(iPhaseBin).asymmetry(t, t2, m_timeBinning.getLifetime(), m_zcp, m_dz);
      rminusfit.SetBinContent(iTimeBin+1, rminus) ;
      rminusfit.SetBinError(iTimeBin+1, 0.) ;
      rplusfit.SetBinContent(iTimeBin+1, rplus) ;
      rplusfit.SetBinError(iTimeBin+1, 0.) ;
      rdifffit.SetBinContent(iTimeBin+1, rdiff) ;
      rdifffit.SetBinError(iTimeBin+1, 0.) ;
      acpfit.SetBinContent(iTimeBin+1, acp);
      acpfit.SetBinError(iTimeBin+1, 0.);
    }
    binhistos.push_back(rminusfit) ;
    binhistos.push_back(rplusfit) ;
    binhistos.push_back(rdifffit) ;
    binhistos.push_back(acpfit);
    ++iPhaseBin ;
  }
  return histos ;
}

void BinFlipChi2::savePlotsVsTime(const string& name, TFile& outputfile) const {
  outputfile.cd() ;
  deque<deque<TH1F> > histos = plotsVsTime(name) ;
  for(auto& binhistos : histos)
    for(auto& histo : binhistos)
      histo.Write() ;
}

AsymmetryChi2::AsymmetryChi2(BinFlipParSet* fitPars, const HadronicParameters& hadronicPars,
			     const TimeBinning& timeBins) :
  BinFlipChi2(fitPars, hadronicPars, timeBins)
{}

AsymmetryChi2::AsymmetryChi2(BinFlipParSet* fitPars,
			     const string& hadronicParsName, const string& timeBinsName,
			     const string& fname) :
  BinFlipChi2(fitPars, hadronicParsName, timeBinsName, fname)
{}

double AsymmetryChi2::binChiSquared(unsigned iTimeBin, unsigned iPhaseBin) const {
  double expectedA = m_hadronicPars.bin(iPhaseBin).asymmetry(m_timeBinning.meanUnmixedTime(iTimeBin),
							     m_timeBinning.meanUnmixedTime2(iTimeBin),
							     m_timeBinning.getLifetime(), m_zcp, m_dz);
  return m_timeBinning.asymmetryChiSquared(iTimeBin, iPhaseBin, expectedA);
}

BinFlipChi2Simul::BinFlipChi2Simul(BinFlipChi2Base* flip1, BinFlipChi2Base* flip2) :
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
