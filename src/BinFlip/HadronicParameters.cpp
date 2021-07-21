#include <Mint/HadronicParameters.h>
#include <TMath.h>
#include <Mint/NamedParameter.h>
#include <fstream>
#include <Mint/MinuitParameterSet.h>
#include <Mint/GaussianConstraintChi2.h>
#include <stdexcept>

using namespace std ;
using MINT::NamedParameter ;
using MINT::FitParameter;
using MINT::NamedParameterBase;
using MINT::MinuitParameterSet;
using MINT::IMinimisable;

int HadronicParameters::EventBinInfo::binSign() const {
  return binNumber > 0 ? 1 : -1 ;
}

bool HadronicParameters::EventBinInfo::isPlus() const {
  return binNumber > 0 ;
}

HadronicParameters::PhaseBinningBase::PhaseBinningBase(unsigned nBins) :
  nBins(nBins)
{}

int HadronicParameters::PhaseBinningBase::binNumber(IDalitzEvent& evt) const {
  return binInfo(evt).binNumber ;
}

bool HadronicParameters::PhaseBinningBase::operator==(const HadronicParameters::PhaseBinningBase& other) const {
  return other.nBins == nBins && other.type() == type() ;
}

bool HadronicParameters::PhaseBinningBase::operator!=(const HadronicParameters::PhaseBinningBase& other) const {
  return !(operator==(other)) ;
}

HadronicParameters::ModelPhaseBinning::ModelPhaseBinning(HadronicParameters::ModelPtr model,
							 HadronicParameters::ModelPtr cpmodel,
							 unsigned nBins) :
  PhaseBinningBase(nBins),
  m_model(model),
  m_cpmodel(cpmodel)
{}

string HadronicParameters::ModelPhaseBinning::type() const {
  return "model" ;
}

bool HadronicParameters::ModelPhaseBinning::isFavoured(const double F, const double Fbar, IDalitzEvent&) const {
  return F > Fbar ;
}

HadronicParameters::EventBinInfo HadronicParameters::ModelPhaseBinning::binInfo(IDalitzEvent& evt) const {
  EventBinInfo binNo ;
  binNo.evt = &evt ;
  binNo.amp = m_model->ComplexVal(evt) ;
  binNo.ampBar = m_cpmodel->ComplexVal(evt) ;
  binNo.F = norm(binNo.amp) ;
  binNo.Fbar = norm(binNo.ampBar) ;
  double phasediff = arg(binNo.amp/binNo.ampBar) ;
  // Invert the phase difference in the suppressed region, effectively using the phase difference
  // from the favoured CP-conjugate point to determine the bin. This assumes no direct CPV, ie, that
  // arg(Fminus/Fbarplus) = -arg(Fplus/Fbarminus). Could use the CP-conjugate event to determine
  // the phase difference if we need to allow for direct CPV.
  if(!isFavoured(binNo.F, binNo.Fbar, evt))
    phasediff *= -1 ;
  if(phasediff < 0.)
    phasediff += 2. * TMath::Pi() ;
  binNo.binNumber = (int(phasediff * nBins/2./TMath::Pi() + 0.5) % nBins + 1) ;
  // Possible issue if F == Fbar ?
  if(binNo.Fbar > binNo.F)
    binNo.binNumber *= -1 ;
  return binNo ;
}

void HadronicParameters::ModelPhaseBinning::Print(const std::string& _name, std::ostream& os) const {
  string name = _name + "_binning_" ;
  os << name << "type\t" << type() << endl ;
  os << name << "nBins\t" << nBins << endl ;
  os << name << "eventPattern" ;
  for(const auto& ipat : m_model->getAmpPtr(0)->getTreePattern().getVectorOfInts())
    os << "\t" << ipat ;
  os << endl ;
  os << name << "cpEventPattern" ;
  for(const auto& ipat : m_cpmodel->getAmpPtr(0)->getTreePattern().getVectorOfInts())
    os << "\t" << ipat ;
  os << endl ;
}

HadronicParameters::BinningPtr HadronicParameters::ModelPhaseBinning::fromConfig(const string& _name,
										 const string& fname) {
  string name = _name + "_binning_" ;
  if(string(NamedParameter<string>(name + "type", fname.c_str())) != "model")
    return BinningPtr(0) ;
  NamedParameter<int> nBins(name + "nBins", fname.c_str()) ;
  NamedParameter<int> pat(name + "eventPattern", fname.c_str()) ;
  NamedParameter<int> cpPat(name + "cpEventPattern", fname.c_str()) ;
  ModelPtr model(new FitAmpSum(DalitzEventPattern(pat), fname.c_str())) ;
  ModelPtr cpModel(new FitAmpSum(DalitzEventPattern(cpPat), fname.c_str())) ;
  return BinningPtr(new ModelPhaseBinning(model, cpModel, nBins)) ;
}

/*bool HadronicParameters::ModelPhaseBinning::operator==(const PhaseBinningBase& other) const {
  if(!PhaseBiningBase::operator==(other))
    return false ;
  const ModelPhaseBinning* ptr = (ModelPhaseBinning*)(&other) ;
  }
*/

HadronicParameters::ModelBinning3Body::ModelBinning3Body(HadronicParameters::ModelPtr model,
							 HadronicParameters::ModelPtr cpmodel,
							 unsigned nBins) :
  ModelPhaseBinning(model, cpmodel, nBins)
{}

string HadronicParameters::ModelBinning3Body::type() const {
  return "model3body" ;
}

bool HadronicParameters::ModelBinning3Body::isFavoured(const double, const double, IDalitzEvent& evt) const {
  const DalitzEventPattern& pat = evt.eventPattern() ;
  int i1(1), i2(2), i3(3) ;
  // Assuming we have a decay of the form P->h+h-h0 (in some order).
  for(unsigned i = 1 ; i < pat.size() ; ++i){
    if(pat[i].charge() > 0)
      i1 = i ;
    else if(pat[i].charge() < 0)
      i2 = i ;
    else 
      i3 = i ;
  }
  double s13 = evt.s(i1, i3) ;
  double s23 = evt.s(i2, i3) ;
  return s13 > s23 ;
}

HadronicParameters::BinningPtr HadronicParameters::ModelBinning3Body::fromConfig(const string& _name,
										 const string& fname) {
  string name = _name + "_binning_" ;
  if(string(NamedParameter<string>(name + "type", fname.c_str())) != "model3body")
    return BinningPtr(0) ;
  NamedParameter<int> nBins(name + "nBins", fname.c_str()) ;
  NamedParameter<int> pat(name + "eventPattern", fname.c_str()) ;
  NamedParameter<int> cpPat(name + "cpEventPattern", fname.c_str()) ;
  ModelPtr model(new FitAmpSum(DalitzEventPattern(pat), fname.c_str())) ;
  ModelPtr cpModel(new FitAmpSum(DalitzEventPattern(cpPat), fname.c_str())) ;
  return BinningPtr(new ModelBinning3Body(model, cpModel, nBins)) ;
}
  

HadronicParameters::Bin::Bin(const string& name, unsigned ibin,
			     double _Fplus, double _Fminus, const complex<double>& X,
			     double _Fbarplus, double _Fbarminus, const complex<double>& Xbar,
			     double norm, double normBar, double sumw, double sumw2) :
  m_name(getName(name, ibin)),
  m_Fplus(new FitParameter(m_name + "_Fplus", FitParameter::HIDE, _Fplus, 0.)),
  m_Fminus(new FitParameter(m_name + "_Fminus", FitParameter::HIDE, _Fminus, 0.)),
  m_C(new FitParameter(m_name + "_C", FitParameter::HIDE, X.real(), 0.)),
  m_S(new FitParameter(m_name + "_S", FitParameter::HIDE, X.imag(), 0.)),
  m_Fbarplus(new FitParameter(m_name + "_Fbarplus", FitParameter::HIDE, _Fbarplus, 0.)),
  m_Fbarminus(new FitParameter(m_name + "_Fbarminus", FitParameter::HIDE, _Fbarminus, 0.)),
  m_Cbar(new FitParameter(m_name + "_Cbar", FitParameter::HIDE, Xbar.real(), 0.)),
  m_Sbar(new FitParameter(m_name + "_Sbar", FitParameter::HIDE, Xbar.real(), 0.)),
  m_sumw(sumw),
  m_sumw2(sumw2),
  m_norm(norm),
  m_normBar(normBar)
{}

HadronicParameters::Bin::Bin(const string& name, unsigned ibin) :
  m_name(getName(name, ibin)),
  m_Fplus(new FitParameter(m_name + "_Fplus", FitParameter::HIDE, 0., 0.)),
  m_Fminus(new FitParameter(m_name + "_Fminus", FitParameter::HIDE, 0., 0.)),
  m_C(new FitParameter(m_name + "_C", FitParameter::HIDE, 0., 0.)),
  m_S(new FitParameter(m_name + "_S", FitParameter::HIDE, 0., 0.)),
  m_Fbarplus(new FitParameter(m_name + "_Fbarplus", FitParameter::HIDE, 0., 0.)),
  m_Fbarminus(new FitParameter(m_name + "_Fbarminus", FitParameter::HIDE, 0., 0.)),
  m_Cbar(new FitParameter(m_name + "_Cbar", FitParameter::HIDE, 0., 0.)),
  m_Sbar(new FitParameter(m_name + "_Sbar", FitParameter::HIDE, 0., 0.)),
  m_sumw(0.),
  m_sumw2(0.),
  m_norm(1.),
  m_normBar(1.)
{}

void HadronicParameters::Bin::initFitPar(FitParPtr& fitPar, FitParPtr alt,
					 const string& fname,
					 const string& parName, const string& altParName) {
  NamedParameter<double> par(m_name + parName, fname.c_str(),
			     NamedParameterBase::QUIET);
  if(!par.gotInitialised() && altParName.size() > 0)
    par = NamedParameter<double>(m_name + altParName, fname.c_str(),
				 NamedParameterBase::QUIET);
  if(!par.gotInitialised()){
    fitPar = alt;
    return;
  }
  int fow = FitParameter::HIDE;
  double mean(0.), step(0.), mi(0.), ma(0.);
  if(par.size() == 1)
    mean = par.getVal(0);
  else{
    fow = int(par.getVal(0));
    mean = par.getVal(1);
    step = par.getVal(2);
  }
  if(par.size() == 5){
    mi = par.getVal(3);
    ma = par.getVal(4);
  }
  fitPar = FitParPtr(new FitParameter(m_name + parName, fow, mean, step, mi, ma));
}
  
HadronicParameters::Bin::Bin(const string& _name, unsigned number, const string& fname) :
  m_name(getName(_name, number))
{
  initFitPar(m_Fplus, nullptr, fname, "_Fplus");
  initFitPar(m_Fminus, nullptr, fname, "_Fminus");
  initFitPar(m_C, nullptr, fname, "_C", "_Xplus_Re");
  initFitPar(m_S, nullptr, fname, "_S", "_Xplus_Im");

  // In case of no CPV, Abar == A
  initFitPar(m_Fbarplus, m_Fplus, fname, "_Fbarplus");
  initFitPar(m_Fbarminus, m_Fminus, fname, "_Fbarminus");
  initFitPar(m_Cbar, m_C, fname, "_Cbar", "_Xbarplus_Re");
  initFitPar(m_Sbar, m_S, fname, "_Sbar", "_Xbarplus_Im");

  NamedParameter<double> sumw(m_name + "_sumw", fname.c_str()) ;
  NamedParameter<double> sumw2(m_name + "_sumw2", fname.c_str()) ;
  NamedParameter<double> norm(m_name + "_norm", fname.c_str()) ;
  NamedParameter<double> normBar(m_name + "_normBar", fname.c_str(),
				 NamedParameterBase::QUIET) ;

  m_norm = norm ;
  if(normBar.gotInitialised())
    m_normBar = normBar;
  else
    m_normBar = norm;
  m_sumw = sumw ;
  m_sumw2 = sumw2 ;
}

void HadronicParameters::Bin::setNoCPV() {
  m_Fbarplus = m_Fplus;
  m_Fbarminus = m_Fminus;
  m_Cbar = m_C;
  m_Sbar = m_S;
}

bool HadronicParameters::Bin::allowsCPV() const {
  return m_Fbarplus.get() != m_Fplus.get();
}

vector<FitParameter*> HadronicParameters::Bin::getFitPars(){
  vector<FitParameter*> pars;
  pars.push_back(m_Fplus.get());
  pars.push_back(m_Fminus.get());
  pars.push_back(m_C.get());
  pars.push_back(m_S.get());
  if(!allowsCPV())
    return pars;
  pars.push_back(m_Fbarplus.get());
  pars.push_back(m_Fbarminus.get());
  pars.push_back(m_Cbar.get());
  pars.push_back(m_Sbar.get());
  return pars;
}

vector<const FitParameter*> HadronicParameters::Bin::getFitPars() const {
  vector<const FitParameter*> pars;
  pars.push_back(m_Fplus.get());
  pars.push_back(m_Fminus.get());
  pars.push_back(m_C.get());
  pars.push_back(m_S.get());
  if(!allowsCPV())
    return pars;
  pars.push_back(m_Fbarplus.get());
  pars.push_back(m_Fbarminus.get());
  pars.push_back(m_Cbar.get());
  pars.push_back(m_Sbar.get());
  return pars;
}

void HadronicParameters::Bin::setParSet(MinuitParameterSet* pSet){
  vector<FitParameter*> pars = getFitPars();
  for(auto& par : pars)
    par->addToParSet(pSet);
}

void HadronicParameters::Bin::add(const HadronicParameters::EventBinInfo& evtPlus,
				  const HadronicParameters::EventBinInfo& evtMinus,
				  double weight) {
  *m_Fplus = *m_Fplus + evtPlus.F * weight ;
  *m_Fminus = *m_Fminus + evtMinus.F * weight ;
  complex<double> dX = conj(evtPlus.amp) * evtPlus.ampBar * weight;
  *m_C = *m_C + dX.real();
  *m_S = *m_S + dX.imag();
  if(allowsCPV()){
    *m_Fbarplus = *m_Fbarplus + evtMinus.Fbar * weight ;
    *m_Fbarminus = *m_Fbarminus + evtPlus.Fbar * weight ;
    complex<double> dXbar = conj(evtMinus.ampBar) * evtMinus.amp * weight ;
    *m_Cbar = *m_C + dXbar.real();
    *m_Sbar = *m_S + dXbar.imag();
  }
  m_sumw += weight ;
  m_sumw2 += weight*weight ;
}

void HadronicParameters::Bin::finaliseSum() {
  double denom = sqrt((*m_Fplus)*(*m_Fbarminus));
  double denomBar = sqrt((*m_Fbarplus)*(*m_Fminus));
  *m_C = *m_C/denom;
  *m_S = *m_S/denom;
  *m_Fplus = *m_Fplus/m_sumw;
  *m_Fminus = *m_Fminus/m_sumw;
  if(allowsCPV()){
    *m_Cbar = *m_Cbar/denomBar;
    *m_Sbar = *m_Sbar/denomBar;
    *m_Fbarplus = *m_Fbarplus/m_sumw;
    *m_Fbarminus = *m_Fbarminus/m_sumw;
  }
}

void HadronicParameters::Bin::normalise(double _norm, double _normBar) {
  setNorm(_norm);
  setNormBar(_normBar);
}

double HadronicParameters::Bin::Fplus() const {
  return *m_Fplus;
}

double HadronicParameters::Bin::Fminus() const {
  return *m_Fminus;
}

complex<double> HadronicParameters::Bin::Xplus() const {
  return complex<double>(*m_C, *m_S);
}

complex<double> HadronicParameters::Bin::Xminus() const {
  return conj(Xplus()) ;
}

double HadronicParameters::Bin::Fbarplus() const {
  return *m_Fbarplus;
}

double HadronicParameters::Bin::Fbarminus() const {
  return *m_Fbarminus;
}

complex<double> HadronicParameters::Bin::Xbarplus() const {
  return complex<double>(*m_Cbar, *m_Sbar);
}

complex<double> HadronicParameters::Bin::Xbarminus() const {
  return conj(Xbarplus()) ;
}

double HadronicParameters::Bin::getNorm() const {
  return m_norm ;
}

void HadronicParameters::Bin::setNorm(double norm) {
  *m_Fplus = *m_Fplus * m_norm;
  *m_Fminus = *m_Fminus * m_norm;
  m_norm = norm ;
  *m_Fplus = *m_Fplus/norm;
  *m_Fminus = *m_Fminus/norm;
}

double HadronicParameters::Bin::getNormBar() const {
  return m_normBar ;
}

void HadronicParameters::Bin::setNormBar(double norm) {
  if(!allowsCPV())
    return;
  *m_Fbarplus = *m_Fbarplus * m_normBar;
  *m_Fbarminus = *m_Fbarminus * m_normBar;
  m_normBar = norm ;
  *m_Fbarplus = *m_Fbarplus/norm;
  *m_Fbarminus = *m_Fbarminus/norm;
}

void HadronicParameters::Bin::Print(const string& _name, unsigned number, ostream& os) const {
  string name = getName(_name, number) ;
  os.precision(16) ;
  os << *m_Fplus << endl;
  os << *m_Fminus << endl;
  os << *m_C << endl;
  os << *m_S << endl;
  if(allowsCPV()){
    os << *m_Fbarplus << endl;
    os << *m_Fbarminus << endl;
    os << *m_Cbar << endl;
    os << *m_Sbar << endl;
  }
  os << name << "_sumw\t" << m_sumw << endl ;
  os << name << "_sumw2\t" << m_sumw2 << endl ;
  os << name << "_norm\t" << m_norm << endl ;
  if(allowsCPV())
    os << name << "_normBar\t" << m_normBar << endl ;
}

string HadronicParameters::Bin::getName(const string& name, unsigned number) {
  ostringstream sstr ;
  sstr << name << "_bin_" << number ;
  return sstr.str() ;
}

complex<double> HadronicParameters::Bin::sumz(const complex<double>& zcp, const complex<double>& dz) const {
  return zcp + dz;
}

complex<double> HadronicParameters::Bin::sumzBar(const complex<double>& zcp, const complex<double>& dz) const {
  return zcp - dz;
}

double HadronicParameters::Bin::R(double t, double t2, double lifetime,
				  const complex<double>& zcp, const complex<double>& dz) const {
  return _R(t, t2, lifetime, zcp, dz, Fplus(), Fminus(), Xplus(), sumz(zcp, dz)) ;
}

double HadronicParameters::Bin::Rbar(double t, double t2, double lifetime,
				     const complex<double>& zcp, const complex<double>& dz) const {
  return _R(t, t2, lifetime, zcp, dz, Fbarplus(), Fbarminus(), Xbarplus(), sumzBar(zcp, dz)) ;
}

pair<double, double> HadronicParameters::Bin::_N(double t, double t2, double lifetime,
						 const complex<double>& zcp, const complex<double>& dz,
						 double _Fplus, double _Fminus,
						 const complex<double>& X, const complex<double>& sumz) const {
  t /= lifetime ;
  t2 /= lifetime * lifetime ;
  double r = _Fminus/_Fplus ;
  double term1 = (1 + 0.25 * t2 * (zcp*zcp - dz*dz).real()) ;
  double term2 = 0.25 * t2 * norm(sumz) ;
  double term3 = sqrt(r) * t * (conj(X) * sumz).real() ;
  double numerator = r * term1 + term2 + term3 ;
  double term4 = sqrt(r) * t * (X * sumz).real() ;
  double denominator = term1 + r * term2 + term4 ;
  return pair(numerator, denominator) ;
}

double HadronicParameters::Bin::_R(double t, double t2, double lifetime,
				   const complex<double>& zcp, const complex<double>& dz,
				   double _Fplus, double _Fminus,
				   const complex<double>& X, const complex<double>& _sumz) const {
  auto n = _N(t, t2, lifetime, zcp, dz, _Fplus, _Fminus, X, _sumz);
  return n.first/n.second;
}

pair<double, double> HadronicParameters::Bin::N(double t, double t2, double lifetime,
						const complex<double>& zcp, const complex<double>& dz) const {
  return _N(t, t2, lifetime, zcp, dz, Fplus(), Fminus(), Xplus(), sumz(zcp, dz)) ;
}

pair<double, double> HadronicParameters::Bin::Nbar(double t, double t2, double lifetime,
						   const complex<double>& zcp, const complex<double>& dz) const {
  return _N(t, t2, lifetime, zcp, dz, Fbarplus(), Fbarminus(), Xbarplus(), sumzBar(zcp, dz)) ;
}

double HadronicParameters::Bin::asymmetry(double t, double t2, double lifetime,
					  const complex<double>& zcp, const complex<double>& dz) const {
  auto _n = N(t, t2, lifetime, zcp, dz);
  auto _nbar = Nbar(t, t2, lifetime, zcp, dz);
  double n = _n.first + _n.second;
  double nbar = _nbar.first + _nbar.second;
  return (n - nbar)/(n + nbar);
}

HadronicParameters::HadronicParameters(const string& name,
				       const HadronicParameters::Bins& bins, 
				       BinningPtr phaseBinning) :
  m_name(name),
  m_bins(bins),
  m_phaseBinning(phaseBinning)
{
  if(m_bins.size() != m_phaseBinning->nBins)
    throw length_error("HadronicParameters: Number of bins given doesn't match \
number of bins in the binning scheme!");
}

HadronicParameters::HadronicParameters(const string& name, BinningPtr phaseBinning) :
  m_name(name),
  m_bins(),
  m_phaseBinning(phaseBinning)
{
  for(unsigned i = 1 ; i < m_phaseBinning->nBins + 1 ; ++i)
    m_bins.push_back(Bin(name, i)) ;
}

HadronicParameters::HadronicParameters(const string& name, const string& fname) :
  m_name(name),
  m_bins(),
  m_phaseBinning(HadronicParameters::getPhaseBinning(name, fname))
{
  for(unsigned i = 1 ; i < m_phaseBinning->nBins + 1 ; ++i)
    m_bins.push_back(Bin(name, i, fname)) ;
}

const HadronicParameters::Bin& HadronicParameters::bin(IDalitzEvent& evt) const {
  int binNo = abs(m_phaseBinning->binNumber(evt)) ;
  return m_bins.at(binNo-1) ;
}

const HadronicParameters::Bin& HadronicParameters::bin(unsigned i) const {
  return m_bins.at(i-1) ;
}

void HadronicParameters::add(IDalitzEvent& evt, double weight) {
  EventBinInfo binNo = m_phaseBinning->binInfo(evt) ;
  DalitzEvent cpEvt(evt) ;
  cpEvt.CP_conjugateYourself() ;
  EventBinInfo cpBinNo = m_phaseBinning->binInfo(cpEvt) ;
  if(binNo.isPlus())
    m_bins.at(binNo.binNumber-1).add(binNo, cpBinNo, weight) ;
  else 
    m_bins.at(cpBinNo.binNumber-1).add(cpBinNo, binNo, weight) ;
}

void HadronicParameters::add(const DalitzEventPattern& pat, TRandom3& rndm, unsigned nevt) {
  for(unsigned i = 0 ; i < nevt ; ++i)
    add(pat, rndm) ;
}

void HadronicParameters::add(const DalitzEventPattern& pat, TRandom3& rndm) {
  DalitzEvent evt(pat, (TRandom*)&rndm) ;
  add(evt) ;
}

pair<double, double> HadronicParameters::integral() const {
  double norm = 0. ;
  double normBar = 0. ;
  for(const auto& bin : m_bins){
    norm += bin.Fplus() + bin.Fminus() ;
    normBar += bin.Fbarplus() + bin.Fbarminus() ;
  }
  return pair<double, double>(norm, normBar) ;
}

void HadronicParameters::finaliseSum() {
  for(auto& bin : m_bins)
    bin.finaliseSum();
  normalise();
}

pair<double, double> HadronicParameters::normalise(double _norm, double _normBar) {
  for(auto& bin : m_bins){
    bin.setNorm(1.) ;
    bin.setNormBar(1.) ;
  }
  pair<double, double> norm = integral() ;
  norm.first /= _norm ;
  norm.second /= _normBar ;
  for(auto& bin : m_bins){
    bin.setNorm(norm.first) ;
    bin.setNormBar(norm.second) ;
  }
  return norm ;
}

void HadronicParameters::Print(ostream& os) const {
  m_phaseBinning->Print(m_name, os) ;
  unsigned i = 1 ;
  for(const auto& bin : m_bins){
    bin.Print(m_name, i, os) ;
    i += 1 ;
  }
}

void HadronicParameters::write(const string& fname) const {
  ofstream fout ;
  fout.open(fname) ;
  Print(fout) ;
  fout.close() ;
}
  
HadronicParameters::BinningPtr
HadronicParameters::getPhaseBinning(const string& name, const string& fname) {
  NamedParameter<string> type(name + "_binning_type", fname.c_str()) ;
  string strtype(type) ;
  if(strtype == string("model"))
    return BinningPtr(ModelPhaseBinning::fromConfig(name, fname)) ;
  if(strtype == string("model3body"))
    return BinningPtr(ModelBinning3Body::fromConfig(name, fname)) ;
  throw invalid_argument("Unknown binning type: " + strtype) ;
}

const HadronicParameters::PhaseBinningBase& HadronicParameters::binning() const {
  return *m_phaseBinning ;
}

HadronicParameters::BinningPtr HadronicParameters::binningPtr() const {
  return m_phaseBinning ;
}

void HadronicParameters::setParSet(MinuitParameterSet* pSet){
  for(auto& bin : m_bins)
    bin.setParSet(pSet);
}

bool HadronicParameters::allowsCPV() const {
  return bin(1).allowsCPV();
}

IMinimisable* HadronicParameters::getConstraints(MinuitParameterSet* pSet) {
  setParSet(pSet);

  auto pars = getConstrainedPars();
  if(pars.size() > 0)
    return new GaussianConstraintChi2(pSet, pars, m_covMatrix);
  return nullptr;
}

vector<const FitParameter*> HadronicParameters::getFloatingPars() const {
  // Get floating parameters from the bins
  vector<const FitParameter*> pars;
  for(auto& bin : m_bins){
    auto binpars = bin.getFitPars();
    for(auto& par : binpars)
      if(par->iFixInit() == FitParameter::FIT)
	pars.push_back(par);
  }
  return pars;
}

vector<FitParameter*> HadronicParameters::getFloatingPars() {
  // Get floating parameters from the bins
  vector<FitParameter*> pars;
  for(auto& bin : m_bins){
    auto binpars = bin.getFitPars();
    for(auto& par : binpars)
      if(par->iFixInit() == FitParameter::FIT)
	pars.push_back(par);
  }
  return pars;
}

vector<const FitParameter*> HadronicParameters::getConstrainedPars() const {
  return getFloatingPars();
}

vector<FitParameter*> HadronicParameters::getConstrainedPars() {
  return getFloatingPars();
}

void HadronicParameters::setCovMatrix(const CovMatrix& cov) {
  if(cov.GetNcols() != getConstrainedPars().size())
    throw length_error("HadronicParameters::setCovMatrix: Inconsistent size\
 of covariance matrix and number of floating parameters!");
  m_covMatrix = cov;
}
		       
