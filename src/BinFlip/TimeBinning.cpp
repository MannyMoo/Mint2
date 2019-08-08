#include <Mint/TimeBinning.h>
#include <Mint/NamedParameter.h>
#include <fstream>
#include <TFile.h>

using namespace std ;
using MINT::NamedParameter ;

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

bool TimeBinning::Bin::usePoissonErrs() const {
  return m_sumw == m_sumw2 ;
}

double TimeBinning::Bin::chiSquared(double R) const {
  if(nMinus() == 0. or nPlus() == 0.)
    return 0. ;
  double numerator = nMinus() - nPlus() * R ;
  numerator *= numerator ;
  double denominator = m_sumw2minus + m_sumw2plus * R*R ;
  if(denominator == 0.)
    return 0. ;
  return numerator/denominator ;
}

TimeBinning::TimeBinning(const vector<double>& timeBins, HadronicParameters::BinningPtr phaseBinning) :
  m_timeBins(timeBins),
  m_bins(timeBins.size(), Bins(phaseBinning->nBins)),
  m_binsBar(timeBins.size(), Bins(phaseBinning->nBins)),
  m_binsInt(timeBins.size()),
  m_phaseBinning(phaseBinning)
{}

TimeBinning::TimeBinning(const string& name, const string& fname) :
  m_timeBins(),
  m_bins(),
  m_binsBar(),
  m_binsInt(),
  m_phaseBinning(0)
{
  NamedParameter<double> timeBins(name + "_timeBins", fname.c_str()) ;
  m_timeBins = timeBins.getVector() ;
  m_phaseBinning = HadronicParameters::getPhaseBinning(name, fname) ;
  string nameint = name + "_int" ;
  string nameplus = name + "_plus" ;
  string nameminus = name + "_minus" ;
  for(unsigned iTimeBin = 0 ; iTimeBin < m_timeBins.size() ; ++iTimeBin){
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
}
      
void TimeBinning::Print(const string& name, ostream& os) const {
  os.precision(16) ;
  os << name << "_timeBins" ;
  for(const auto& bin : m_timeBins)
    os << "\t" << bin ;
  os << endl ;
  string nameint = name + "_int" ;
  string nameplus = name + "_plus" ;
  string nameminus = name + "_minus" ;
  for(unsigned iTimeBin = 0 ; iTimeBin < m_timeBins.size() ; ++iTimeBin){
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
}

HadronicParameters::BinningPtr TimeBinning::phaseBinning() const {
  return m_phaseBinning ;
}

int TimeBinning::timeBin(double t) const {
  for(unsigned i = 0 ; i < m_timeBins.size()-1 ; ++i)
    if(t < m_timeBins[i+1])
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
  return m_timeBins.size() ;
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
  
  TH1F hplus((name + "nPlus").c_str(), "", m_timeBins.size()-1, &m_timeBins[0]) ;
  TH1F hminus((name + "nMinus").c_str(), "", m_timeBins.size()-1, &m_timeBins[0]) ;
  TH1F hr((name + "ratio").c_str(), "", m_timeBins.size()-1, &m_timeBins[0]) ;
  for(unsigned i = 0 ; i < nBinsTime() ; ++i) {
    const Bin& pbin = tag > 0 ? bin(i, phaseBin) : binBar(i, phaseBin) ;
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

void TimeBinning::savePlotsVsTime(const string& name, TFile& outputfile) const {
  outputfile.cd() ;
  for(int tag = -1 ; tag < 2 ; tag += 2){
    for(unsigned iphase = 1 ; iphase < m_phaseBinning->nBins + 1 ; ++iphase){
      deque<TH1F> histos = plotVsTime(name, iphase, tag) ;
      for(const auto& h : histos)
	h.Write() ;
    }
  }
}