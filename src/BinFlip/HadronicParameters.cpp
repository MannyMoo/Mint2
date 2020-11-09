#include <Mint/HadronicParameters.h>
#include <TMath.h>
#include <Mint/NamedParameter.h>
#include <fstream>

using namespace std ;
using MINT::NamedParameter ;

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
  

HadronicParameters::Bin::Bin(double _Fplus, double _Fminus, const complex<double>& X,
			     double _Fbarplus, double _Fbarminus, const complex<double>& Xbar,
			     double norm, double normBar, double sumw, double sumw2) :
  m_Fplus(_Fplus * sumw * norm),
  m_Fminus(_Fminus * sumw * norm),
  m_X(X),
  m_Fbarplus(_Fbarplus * sumw * normBar),
  m_Fbarminus(_Fbarminus * sumw * normBar),
  m_Xbar(Xbar),
  m_sumw(sumw),
  m_sumw2(sumw2),
  m_norm(norm),
  m_normBar(normBar)
{
  m_X *= sqrt(Fplus() * Fminus()) * sumw * norm ;
  m_Xbar *= sqrt(Fbarplus() * Fbarminus()) * sumw * normBar ;
}

HadronicParameters::Bin::Bin() :
  m_Fplus(0.),
  m_Fminus(0.),
  m_X(0., 0.),
  m_Fbarplus(0.),
  m_Fbarminus(0.),
  m_Xbar(0., 0.),
  m_sumw(0.),
  m_sumw2(0.),
  m_norm(1.),
  m_normBar(1.)
{}

HadronicParameters::Bin::Bin(const string& _name, unsigned number, const string& fname) :
  m_Fplus(0.),
  m_Fminus(0.),
  m_X(0., 0.),
  m_Fbarplus(0.),
  m_Fbarminus(0.),
  m_Xbar(0., 0.),
  m_sumw(0.),
  m_sumw2(0.),
  m_norm(0.),
  m_normBar(1.)
{
  string name = getName(_name, number) ;
  NamedParameter<double> _Fplus(name + "_Fplus", fname.c_str()) ;
  NamedParameter<double> _Fminus(name + "_Fminus", fname.c_str()) ;
  NamedParameter<double> Xplus_Re(name + "_Xplus_Re", fname.c_str()) ;
  NamedParameter<double> Xplus_Im(name + "_Xplus_Im", fname.c_str()) ;

  NamedParameter<double> _Fbarplus(name + "_Fbarplus", fname.c_str()) ;
  NamedParameter<double> _Fbarminus(name + "_Fbarminus", fname.c_str()) ;
  NamedParameter<double> Xbarplus_Re(name + "_Xbarplus_Re", fname.c_str()) ;
  NamedParameter<double> Xbarplus_Im(name + "_Xbarplus_Im", fname.c_str()) ;

  NamedParameter<double> sumw(name + "_sumw", fname.c_str()) ;
  NamedParameter<double> sumw2(name + "_sumw2", fname.c_str()) ;
  NamedParameter<double> norm(name + "_norm", fname.c_str()) ;
  NamedParameter<double> normBar(name + "_normBar", fname.c_str()) ;

  m_norm = norm ;
  m_normBar = normBar ;
  m_sumw = sumw ;
  m_sumw2 = sumw2 ;

  m_Fplus = _Fplus * m_sumw * m_norm ;
  m_Fminus = _Fminus * m_sumw * m_norm ;
  m_X = complex<double>(Xplus_Re, Xplus_Im) * sqrt(Fplus() * Fminus()) * m_sumw * m_norm ;

  m_Fbarplus = _Fbarplus * m_sumw * m_normBar ;
  m_Fbarminus = _Fbarminus * m_sumw * m_normBar ;
  m_Xbar = complex<double>(Xbarplus_Re, Xbarplus_Im) * sqrt(Fbarplus() * Fbarminus()) * m_sumw * m_normBar ;
}


void HadronicParameters::Bin::add(const HadronicParameters::EventBinInfo& evtPlus,
				  const HadronicParameters::EventBinInfo& evtMinus,
				  double weight) {
  m_Fplus += evtPlus.F * weight ;
  m_Fminus += evtMinus.F * weight ;
  m_Fbarplus += evtMinus.Fbar * weight ;
  m_Fbarminus += evtPlus.Fbar * weight ;
  m_X += conj(evtPlus.amp) * evtPlus.ampBar * weight ;
  m_Xbar += conj(evtMinus.ampBar) * evtMinus.amp * weight ;
  m_sumw += weight ;
  m_sumw2 += weight*weight ;
}

double HadronicParameters::Bin::Fplus() const {
  return m_Fplus / m_sumw / m_norm ;
}

double HadronicParameters::Bin::Fminus() const {
  return m_Fminus / m_sumw / m_norm ;
}

complex<double> HadronicParameters::Bin::Xplus() const {
  return m_X / m_sumw / sqrt(Fplus() * Fbarminus()) / m_norm ;
}

complex<double> HadronicParameters::Bin::Xminus() const {
  return conj(Xplus()) ;
}

double HadronicParameters::Bin::Fbarplus() const {
  return m_Fbarplus / m_sumw / m_normBar ;
}

double HadronicParameters::Bin::Fbarminus() const {
  return m_Fbarminus / m_sumw / m_normBar ;
}

complex<double> HadronicParameters::Bin::Xbarplus() const {
  return m_Xbar / m_sumw / sqrt(Fbarplus() * Fminus()) / m_normBar ;
}

complex<double> HadronicParameters::Bin::Xbarminus() const {
  return conj(Xbarplus()) ;
}

double HadronicParameters::Bin::getNorm() const {
  return m_norm ;
}

void HadronicParameters::Bin::setNorm(double norm) {
  m_norm = norm ;
}

double HadronicParameters::Bin::getNormBar() const {
  return m_normBar ;
}

void HadronicParameters::Bin::setNormBar(double norm) {
  m_normBar = norm ;
}

void HadronicParameters::Bin::Print(const string& _name, unsigned number, ostream& os) const {
  string name = getName(_name, number) ;
  complex<double> Xp = Xplus() ;
  complex<double> Xbarp = Xbarplus() ;
  os.precision(16) ;
  os << name << "_Fplus\t" << Fplus() << endl ;
  os << name << "_Fminus\t" << Fminus() << endl ;
  os << name << "_Xplus_Re\t" << Xp.real() << endl ;
  os << name << "_Xplus_Im\t" << Xp.imag() << endl ;
  os << name << "_Fbarplus\t" << Fbarplus() << endl ;
  os << name << "_Fbarminus\t" << Fbarminus() << endl ;
  os << name << "_Xbarplus_Re\t" << Xbarp.real() << endl ;
  os << name << "_Xbarplus_Im\t" << Xbarp.imag() << endl ;
  os << name << "_sumw\t" << m_sumw << endl ;
  os << name << "_sumw2\t" << m_sumw2 << endl ;
  os << name << "_norm\t" << m_norm << endl ;
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

HadronicParameters::HadronicParameters(const HadronicParameters::Bins& bins, 
				       BinningPtr phaseBinning) :
  m_bins(bins),
  m_phaseBinning(phaseBinning)
{}

HadronicParameters::HadronicParameters(BinningPtr phaseBinning) :
  m_bins(phaseBinning->nBins, HadronicParameters::Bin()),
  m_phaseBinning(phaseBinning)
{}

HadronicParameters::HadronicParameters(const string& name, const string& fname) :
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

void HadronicParameters::Print(const string& name, ostream& os) const {
  m_phaseBinning->Print(name, os) ;
  unsigned i = 1 ;
  for(const auto& bin : m_bins){
    bin.Print(name, i, os) ;
    i += 1 ;
  }
}

void HadronicParameters::write(const string& name, const string& fname) const {
  ofstream fout ;
  fout.open(fname) ;
  Print(name, fout) ;
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
