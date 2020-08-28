#include <Mint/TimeBinning.h>
#include <Mint/NamedParameter.h>
#include <fstream>
#include <TFile.h>
#include <stdexcept>
#include <RooHistError.h>
#include <TSpline.h>
#include <sys/stat.h>

using namespace std ;
using MINT::NamedParameter ;

bool fileExists( const string fileName ) {

  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(fileName.c_str(),&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    //cout << "ConfigurationFileList: WARNING: File \"" << fileName << "\" cannot be accessed!" << endl ;
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }

  return(blnReturn);
}

TimeBinning::Bin::Bin() :
  m_t(0.),
  m_t2(0.),
  m_sumw(0.),
  m_sumw2(0.),
  m_sumwplus(0.),
  m_sumw2plus(0.),
  m_sumwminus(0.),
  m_sumw2minus(0.)
{}

TimeBinning::Bin::Bin(const string& _name, unsigned number, const string& fname) :
  m_t(0.),
  m_t2(0.),
  m_sumw(0.),
  m_sumw2(0.),
  m_sumwplus(0.),
  m_sumw2plus(0.),
  m_sumwminus(0.),
  m_sumw2minus(0.)
{
  string name = getName(_name, number) ;
  NamedParameter<double> sumw(name + "sumw", fname.c_str()) ;
  m_sumw = sumw ;
  NamedParameter<double> sumw2(name + "sumw2", fname.c_str()) ;
  m_sumw2 = sumw2 ;

  NamedParameter<double> _t(name + "t", fname.c_str()) ;
  m_t = _t * m_sumw ;
  NamedParameter<double> _t2(name + "t2", fname.c_str()) ;
  m_t2 = _t2 * m_sumw2 ;
  
  NamedParameter<double> sumwplus(name + "sumwplus", fname.c_str()) ;
  m_sumwplus = sumwplus ;
  NamedParameter<double> sumw2plus(name + "sumw2plus", fname.c_str()) ;
  m_sumw2plus = sumw2plus ;

  NamedParameter<double> sumwminus(name + "sumwminus", fname.c_str()) ;
  m_sumwminus = sumwminus ;
  NamedParameter<double> sumw2minus(name + "sumw2minus", fname.c_str()) ;
  m_sumw2minus = sumw2minus ;
}

void TimeBinning::Bin::Print(const std::string& _name, unsigned number, ostream& os) const {
  os.precision(16) ;
  string name = getName(_name, number) ;
  os << name << "t\t" << t() << endl ;
  os << name << "t2\t" << t2() << endl ;
  os << name << "sumw\t" << m_sumw << endl ;
  os << name << "sumw2\t" << m_sumw2 << endl ;
  os << name << "sumwplus\t" << m_sumwplus << endl ;
  os << name << "sumw2plus\t" << m_sumw2plus << endl ;
  os << name << "sumwminus\t" << m_sumwminus << endl ;
  os << name << "sumw2minus\t" << m_sumw2minus << endl ;
}

string TimeBinning::Bin::getName(const std::string& name, unsigned number) {
  ostringstream sstr ;
  sstr << name << "_bin_" << number << "_" ;
  return sstr.str() ;
}

void TimeBinning::Bin::add(double t, bool isplus, double weight) {
  double weight2 = weight*weight ;
  m_t += t ;
  m_t2 += t*t ;
  m_sumw += weight ;
  m_sumw2 += weight2 ;
  if(isplus){
    m_sumwplus += weight ;
    m_sumw2plus += weight2 ;
  }
  else {
    m_sumwminus += weight ;
    m_sumw2minus += weight2 ;
  }
}

double TimeBinning::Bin::t() const {
  return m_sumw > 0 ? m_t/m_sumw : 0. ;
}

double TimeBinning::Bin::t2() const {
  return m_sumw > 0 ? m_t2/m_sumw : 0. ;
}

double TimeBinning::Bin::nPlus() const {
  return m_sumwplus ;
}

double TimeBinning::Bin::nPlusErr() const {
  return sqrt(m_sumw2plus) ;
}

double TimeBinning::Bin::nMinus() const {
  return m_sumwminus ;
}

double TimeBinning::Bin::nMinusErr() const {
  return sqrt(m_sumw2minus) ;
}

/// Whether to use Poisson errors (weights are all 1).
bool TimeBinning::Bin::usePoissonErrs() const {
  return m_sumw == m_sumw2 ;
}

/// Get the low and high Poisson errors squared on a count.
void TimeBinning::Bin::poissonErrors2(const double n, double& err2low, double& err2high) {
  RooHistError::instance().getPoissonInterval(n, err2low, err2high);
  err2low = n - err2low;
  err2high -= n;
  err2low *= err2low;
  err2high *= err2high;
}

/// Get the low and high errors squared on n. plus and n. minus.
void TimeBinning::Bin::errors2(double& err2pluslow, double& err2plushigh,
			       double& err2minuslow, double& err2minushigh) const {
  if(!usePoissonErrs()){
    err2pluslow = err2plushigh = m_sumw2plus;
    err2minuslow = err2minushigh = m_sumw2minus;
    return;
  }
  poissonErrors2(m_sumwplus, err2pluslow, err2plushigh);
  poissonErrors2(m_sumwminus, err2minuslow, err2minushigh);
}

double TimeBinning::Bin::chiSquared(double R) const {
  if(nMinus() == 0. or nPlus() == 0.)
    return 0. ;
  double numerator = nMinus() - nPlus() * R ;
  numerator *= numerator ;
  double err2pluslow(0.), err2plushigh(0.), err2minuslow(0.), err2minushigh(0.);
  errors2(err2pluslow, err2plushigh, err2minuslow, err2minushigh);

  // There's probably a more correct way to do this, but this should do.
  // It only matters in the case of Poisson errors, which should only occur
  // in pure signal samples (toys or MC).
  double denominator(0.);
  if(nMinus() > nPlus() * R)
    denominator = err2minuslow + err2plushigh * R*R ;
  else
    denominator = err2minushigh + err2pluslow * R*R ;
    
  if(denominator == 0.)
    return 0. ;
  return numerator/denominator ;
}

TimeBinning::Bin TimeBinning::Bin::operator+(TimeBinning::Bin other) const {
  return other += *this;
}

TimeBinning::Bin& TimeBinning::Bin::operator+=(const TimeBinning::Bin& other) {
  m_t += other.m_t;
  m_t2 += other.m_t2;
  m_sumw += other.m_sumw;
  m_sumw2 += other.m_sumw2;
  m_sumwplus += other.m_sumwplus;
  m_sumw2plus += other.m_sumw2plus;
  m_sumwminus += other.m_sumwminus;
  m_sumw2minus += other.m_sumw2minus;
  return *this;
}

TimeBinning::TimeBinning(const vector<double>& timeBins, HadronicParameters::BinningPtr phaseBinning,
			 double lifetime, const TH1* hefficiency) :
  m_timeBins(timeBins),
  m_meant(nBinsTime(), 0.),
  m_meant2(nBinsTime(), 0.),
  m_lifetime(lifetime),
  m_bins(nBinsTime(), Bins(phaseBinning->nBins)),
  m_binsBar(nBinsTime(), Bins(phaseBinning->nBins)),
  m_binsInt(nBinsTime()),
  m_phaseBinning(phaseBinning),
  m_efficiencySpline(hefficiency != nullptr ? new TSpline3(hefficiency) : 0)
{
  setLifetime(lifetime) ;
}

TimeBinning::TimeBinning(const string& name, const string& fname) :
  m_timeBins(),
  m_meant(),
  m_meant2(),
  m_lifetime(),
  m_bins(),
  m_binsBar(),
  m_binsInt(),
  m_phaseBinning(0),
  m_efficiencySpline(nullptr)
{
  NamedParameter<double> timeBins(name + "_timeBins", fname.c_str()) ;
  m_timeBins = timeBins.getVector() ;
  m_meant = vector<double>(nBinsTime(), 0.) ;
  m_meant2 = vector<double>(nBinsTime(), 0.) ;
  NamedParameter<double> lifetime(name + "_lifetime", fname.c_str()) ;
  setLifetime(lifetime) ;
  m_phaseBinning = HadronicParameters::getPhaseBinning(name, fname) ;
  string nameint = name + "_int" ;
  string nameplus = name + "_plus" ;
  string nameminus = name + "_minus" ;
  for(unsigned iTimeBin = 0 ; iTimeBin < nBinsTime() ; ++iTimeBin){
    m_binsInt.push_back(Bin(nameint, iTimeBin, fname)) ;
    m_bins.push_back(Bins()) ;
    m_binsBar.push_back(Bins()) ;
    string nameplusbin = Bin::getName(nameplus + "_time", iTimeBin) + "phase" ;
    string nameminusbin = Bin::getName(nameminus + "_time", iTimeBin) + "phase" ;
    for(unsigned iPhaseBin = 1 ; iPhaseBin < m_phaseBinning->nBins + 1 ; ++iPhaseBin){
      m_bins.back().push_back(Bin(nameplusbin, iPhaseBin, fname)) ;
      m_binsBar.back().push_back(Bin(nameminusbin, iPhaseBin, fname)) ;
    }
  }
  string rootfile = fname + ".root";
  if(fileExists(rootfile)){
    TFile feff(rootfile.c_str());
    string splinename = name + "_efficiencySpline";
    TSpline3* efficiencySpline = (TSpline3*)feff.Get(splinename.c_str());
    if(efficiencySpline)
      m_efficiencySpline = MINT::counted_ptr<TSpline3>(efficiencySpline);
  }
}
      
void TimeBinning::Print(const string& name, ostream& os) const {
  os.precision(16) ;
  os << name << "_timeBins" ;
  for(const auto& bin : m_timeBins)
    os << "\t" << bin ;
  os << endl ;
  os << name << "_lifetime\t" << m_lifetime << endl ;
  string nameint = name + "_int" ;
  string nameplus = name + "_plus" ;
  string nameminus = name + "_minus" ;
  for(unsigned iTimeBin = 0 ; iTimeBin < nBinsTime() ; ++iTimeBin){
    m_binsInt[iTimeBin].Print(nameint, iTimeBin, os) ;
    string nameplusbin = Bin::getName(nameplus + "_time", iTimeBin) + "phase" ;
    string nameminusbin = Bin::getName(nameminus + "_time", iTimeBin) + "phase" ;
    for(unsigned iPhaseBin = 1 ; iPhaseBin < m_phaseBinning->nBins + 1; ++iPhaseBin){
      bin(iTimeBin, iPhaseBin).Print(nameplusbin, iPhaseBin, os) ;
      binBar(iTimeBin, iPhaseBin).Print(nameminusbin, iPhaseBin, os) ;
    }
  }
  m_phaseBinning->Print(name, os) ;
}

void TimeBinning::write(const string& name, const string& fname) const {
  ofstream outfile ;
  outfile.open(fname) ;
  Print(name, outfile) ;
  outfile.close() ;
  if(m_efficiencySpline.get() == nullptr)
    return;
  string rootfile = fname + ".root";
  TFile feff(rootfile.c_str(), "recreate");
  string savename = name + "_efficiencySpline";
  m_efficiencySpline->Write(savename.c_str());
  feff.Close();
}

HadronicParameters::BinningPtr TimeBinning::phaseBinning() const {
  return m_phaseBinning ;
}

int TimeBinning::timeBin(double t) const {
  for(unsigned i = 0 ; i < nBinsTime() ; ++i)
    if(m_timeBins[i] <= t && t < m_timeBins[i+1])
      return i ;
  return -1 ;
}

void TimeBinning::add(IDalitzEvent& evt, int tag, double t, double weight) {
  int timeBinNo = timeBin(t) ;
  if(timeBinNo < 0)
    return ;
  int phaseBin = m_phaseBinning->binNumber(evt) ;
  bool isPlus = phaseBin > 0 ;
  if(tag > 0){
    _bin(timeBinNo, abs(phaseBin)).add(t, isPlus, weight) ;
  }
  else{
    isPlus = !isPlus ;
    _binBar(timeBinNo, abs(phaseBin)).add(t, isPlus, weight) ;
  }
  _integratedBin(timeBinNo).add(t, isPlus, weight) ;
}

double TimeBinning::chiSquared(unsigned iTimeBin, unsigned iPhaseBin, double Rplus, double Rminus) const {
  return bin(iTimeBin, iPhaseBin).chiSquared(Rplus)
    + binBar(iTimeBin, iPhaseBin).chiSquared(Rminus) ;
}

unsigned TimeBinning::nBinsTime() const {
  return m_timeBins.size()-1 ;
}

unsigned TimeBinning::nBinsPhase() const {
  return m_phaseBinning->nBins ;
}

const TimeBinning::Bin& TimeBinning::integratedBin(unsigned i) const {
  return m_binsInt.at(i) ;
}

const TimeBinning::Bin& TimeBinning::bin(unsigned iTimeBin, unsigned iPhaseBin) const {
  return m_bins.at(iTimeBin).at(iPhaseBin-1) ;
}

const TimeBinning::Bin& TimeBinning::binBar(unsigned iTimeBin, unsigned iPhaseBin) const {
  return m_binsBar.at(iTimeBin).at(iPhaseBin-1) ;
}

TimeBinning::Bin& TimeBinning::_integratedBin(unsigned i) {
  return m_binsInt.at(i) ;
}

TimeBinning::Bin& TimeBinning::_bin(unsigned iTimeBin, unsigned iPhaseBin) {
  return m_bins.at(iTimeBin).at(iPhaseBin-1) ;
}

TimeBinning::Bin& TimeBinning::_binBar(unsigned iTimeBin, unsigned iPhaseBin) {
  return m_binsBar.at(iTimeBin).at(iPhaseBin-1) ;
}

deque<TH1F> TimeBinning::plotVsTime(const string& _name, unsigned phaseBin, int tag) const {
  ostringstream stag ;
  stag << tag ;
  string name = Bin::getName(_name + "_" + stag.str() + "_phase", phaseBin) ;
  
  TH1F hplus((name + "nPlus").c_str(), "", nBinsTime(), &m_timeBins[0]) ;
  TH1F hminus((name + "nMinus").c_str(), "", nBinsTime(), &m_timeBins[0]) ;
  TH1F hr((name + "ratio").c_str(), "", nBinsTime(), &m_timeBins[0]) ;
  for(unsigned i = 1 ; i < nBinsTime()+1 ; ++i) {
    const Bin& pbin = tag > 0 ? bin(i-1, phaseBin) : binBar(i-1, phaseBin) ;
    hplus.SetBinContent(i, pbin.nPlus()) ;
    hplus.SetBinError(i, pbin.nPlusErr()) ;
    hminus.SetBinContent(i, pbin.nMinus()) ;
    hminus.SetBinError(i, pbin.nMinusErr()) ;
    double r = pbin.nPlus() > 0. ? pbin.nMinus()/pbin.nPlus() : 0. ;
    double rerr = r > 0 && pbin.nMinus() > 0 ? r * sqrt(pow(pbin.nPlusErr()/pbin.nPlus(), 2) + pow(pbin.nMinusErr()/pbin.nMinus(), 2)) : 0. ;
    hr.SetBinContent(i, r) ;
    hr.SetBinError(i, rerr) ;
  }
  deque<TH1F> histos ;
  histos.push_back(hplus) ;
  histos.push_back(hminus) ;
  histos.push_back(hr) ;
  return histos ;
}

deque<deque<TH1F> > TimeBinning::plotsVsTime(const string& name) const {
  deque<deque<TH1F> > histos ;
  for(unsigned iphase = 1 ; iphase < m_phaseBinning->nBins + 1 ; ++iphase){
    histos.push_back(deque<TH1F>()) ;
    deque<TH1F>& binhistos = histos.back() ;
    for(int tag = -1 ; tag < 2 ; tag += 2){
      deque<TH1F> taghistos = plotVsTime(name, iphase, tag) ;
      binhistos.insert(binhistos.end(), taghistos.begin(), taghistos.end()) ;
    }
    TH1F rdiff(binhistos.back()) ;
    string rdiffname(Bin::getName(name + "_phase", iphase) + "ratiodiff") ;
    rdiff.SetName(rdiffname.c_str()) ;
    rdiff.Add(&binhistos[2], -1.) ;
    binhistos.push_back(rdiff) ;
  }
  return histos ;
}

void TimeBinning::savePlotsVsTime(const string& name, TFile& outputfile) const {
  outputfile.cd() ;
  deque<deque<TH1F> > histos = plotsVsTime(name) ;
  for(auto& binhistos : histos)
    for(auto& histo : binhistos)
      histo.Write() ;
}

double TimeBinning::unmixedTimeMoment(unsigned ibin, double lifetime, int moment) const {
  const double& tmin = m_timeBins.at(ibin) ;
  const double& tmax = m_timeBins.at(ibin+1) ;
  if(m_efficiencySpline.get() == 0){
    double norm = exp(-tmin/lifetime) - exp(-tmax/lifetime) ;
    double timeMoment = (pow(tmin, moment) * exp(-tmin/lifetime) - pow(tmax, moment) * exp(-tmax/lifetime))/norm ;
    if(moment == 0)
      return timeMoment ;
    timeMoment += moment * lifetime * unmixedTimeMoment(ibin, lifetime, moment-1) ;
    return timeMoment ;
  }
  double timeMoment = 0.;
  double norm = 0.;
  unsigned npoints = 1000;
  double deltat = (tmax-tmin)/npoints;
  double t = deltat/2;
  for(unsigned i = 0; i < npoints; ++i){
    double val = exp(-t/lifetime) * m_efficiencySpline->Eval(t);
    timeMoment += val * pow(t, moment);
    norm += val;
    t += deltat;
  }
  timeMoment /= norm;
  return timeMoment;
}

void TimeBinning::setLifetime(double lifetime) {
  m_lifetime = lifetime ;
  for(unsigned ibin = 0 ; ibin < nBinsTime() ; ++ibin){
    m_meant.at(ibin) = unmixedTimeMoment(ibin, lifetime, 1) ;
    m_meant2.at(ibin) = unmixedTimeMoment(ibin, lifetime, 2) ;
  }
}

double TimeBinning::getLifetime() const {
  return m_lifetime ;
}

double TimeBinning::meanUnmixedTime(unsigned ibin) const {
  return m_meant.at(ibin) ;
}

double TimeBinning::meanUnmixedTime2(unsigned ibin) const {
  return m_meant2.at(ibin) ;
}

/// Check that two TimeBinnings use the same binning scheme.
bool TimeBinning::isConsistent(const TimeBinning& other) const {
  return m_timeBins == other.m_timeBins && m_meant == other.m_meant && m_meant2 == other.m_meant2 \
    && m_lifetime == other.m_lifetime && *m_phaseBinning == *other.m_phaseBinning;
}

TimeBinning TimeBinning::operator+(TimeBinning other) const {
  return other += *this;
}

/*
I feel like this should work, but it doesn't build, giving
error: need ‘typename’ before ‘std::deque<T>::const_iterator’ because ‘std::deque<T>’ is a dependent scope

template <typename T>
std::deque<T>& operator+=(std::deque<T>& self, const std::deque<T>& other) {
  if(self.size() != other.size())
    throw std::invalid_argument("Attempt to add deques of different length!");
  std::deque<T>::const_iterator iother = other.begin();
  for(std::deque<T>::iterator iself = self.begin(); iself != self.end(); ++iself)
    *iself += *iother++;
  return self;
}
*/

TimeBinning& TimeBinning::operator+=(const TimeBinning& other) {
  if(!isConsistent(other))
    throw std::invalid_argument("Attempt to combine incompatible TimeBinning instances!");
  //m_bins += other.m_bins;
  //m_binsBar += other.m_binsBar;
  //m_binsInt += other.m_binsInt;
  for(unsigned iTimeBin = 0; iTimeBin < nBinsTime(); ++iTimeBin){
    m_binsInt.at(iTimeBin) += other.m_binsInt.at(iTimeBin);
    for(unsigned iPhaseBin = 0; iPhaseBin < nBinsPhase(); ++iPhaseBin){
      m_bins.at(iTimeBin).at(iPhaseBin) += other.m_bins.at(iTimeBin).at(iPhaseBin);
      m_binsBar.at(iTimeBin).at(iPhaseBin) += other.m_binsBar.at(iTimeBin).at(iPhaseBin);
    }
  }
  return *this;
}
