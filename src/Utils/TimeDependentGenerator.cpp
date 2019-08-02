#include <Mint/TimeDependentGenerator.h>
#include <Mint/SplineGenerator.h>
#include <TRandom3.h>
#include <iostream>
#include <sys/stat.h>
#include <Mint/DalitzPdfSaveInteg.h>
#include <stdexcept>

using namespace std;
using namespace MINT;

//
map<string, unsigned> TimeDependentGenerator::GenTimeEvent::infoNames = \
  TimeDependentGenerator::GenTimeEvent::makeInfoNames() ;

map<string, unsigned> TimeDependentGenerator::GenTimeEvent::makeInfoNames() {
  map<string, unsigned> infoNames ;
  infoNames["tag"] = TimeDependentGenerator::GenTimeEvent::ITAG ;
  infoNames["decaytime"] = TimeDependentGenerator::GenTimeEvent::IDECAYTIME ;
  infoNames["smeareddecaytime"] = TimeDependentGenerator::GenTimeEvent::ISMEAREDDECAYTIME ;
  return infoNames ;
}

TimeDependentGenerator::GenTimeEvent::GenTimeEvent(const IDalitzEvent& evt, const int tag, const double decaytime,
						   const double smeareddecaytime) :
  DalitzEvent(evt)
{
  setTag(tag) ;
  setDecayTime(decaytime) ;
  setSmearedDecayTime(smeareddecaytime) ;
} 

TimeDependentGenerator::GenTimeEvent::GenTimeEvent(const IDalitzEvent& evt) :
  DalitzEvent(evt)
{
  while(getVectorOfValues().size() < 3)
    getVectorOfValues().push_back(-999.) ;
}

int TimeDependentGenerator::GenTimeEvent::getTag() const {
  return getValueFromVector(ITAG) ;
}

void TimeDependentGenerator::GenTimeEvent::setTag(int tag) {
  setValueInVector(ITAG, tag) ;
}

double TimeDependentGenerator::GenTimeEvent::getDecayTime() const {
  return getValueFromVector(IDECAYTIME) ;
}

void TimeDependentGenerator::GenTimeEvent::setDecayTime(double decaytime) {
  setValueInVector(IDECAYTIME, decaytime) ;
}

double TimeDependentGenerator::GenTimeEvent::getSmearedDecayTime() const {
  return getValueFromVector(ISMEAREDDECAYTIME) ;
}

void TimeDependentGenerator::GenTimeEvent::setSmearedDecayTime(double decaytime) {
  setValueInVector(ISMEAREDDECAYTIME, decaytime) ;
}

// Take the CP conjugate of the head of the decay pattern.
DalitzEventPattern TimeDependentGenerator::anti(DalitzEventPattern pat) {
  pat[0].antiThis() ;
  return pat ;
}

TimeDependentGenerator::TimeDependentGenerator(const DalitzEventPattern& pattern, double width, double deltam,
					       double deltagamma, double qoverp, double phi, 
					       TRandom3* rndm, TH1F* h_efficiency, 
                                               float resWidth, bool addExpEffects, Eff3piSymmetric* sEfficiency) :
  m_rndm(rndm),
  m_pattern(pattern),
  m_cppattern(anti(pattern)),
  m_model(new FitAmpSum(pattern)),
  m_cpmodel(new FitAmpSum(m_cppattern)),
  m_bothmodel(*m_model),
  m_generator(m_pattern, &m_bothmodel),
  m_width(width),
  m_deltam(deltam),
  m_deltagamma(deltagamma),
  m_qoverp(polar(qoverp, phi)),
  m_scale(1.),
  m_ngen(0),
  m_naccept(0),
  m_h_efficiency(h_efficiency),
  m_efficiencyFit(),
  m_resWidth(resWidth),
  m_addExpEffects(addExpEffects),
  m_sEfficiency(sEfficiency)
{
  init() ;
}


TimeDependentGenerator::TimeDependentGenerator(MINT::counted_ptr<FitAmpSum> model, 
					       MINT::counted_ptr<FitAmpSum> cpmodel,
					       double width, double deltam,
					       double deltagamma, double qoverp, double phi, 
					       TRandom3* rndm, TH1F* h_efficiency, 
                                               float resWidth, bool addExpEffects, Eff3piSymmetric* sEfficiency) :
  m_rndm(rndm),
  m_pattern(model->getAmpPtr(0)->getTreePattern()),
  m_cppattern(cpmodel->getAmpPtr(0)->getTreePattern()),
  m_model(model),
  m_cpmodel(cpmodel),
  m_bothmodel(*m_model),
  m_generator(m_pattern, &m_bothmodel),
  m_width(width),
  m_deltam(deltam),
  m_deltagamma(deltagamma),
  m_qoverp(polar(qoverp, phi)),
  m_scale(1.),
  m_ngen(0),
  m_naccept(0),
  m_h_efficiency(h_efficiency),
  m_efficiencyFit(),
  m_resWidth(resWidth),
  m_addExpEffects(addExpEffects),
  m_sEfficiency(sEfficiency)
{
  init() ;
}

void TimeDependentGenerator::init() {
  m_bothmodel.addAsList(*m_cpmodel) ;
  if( m_h_efficiency != NULL ){
    m_efficiencyFit = TSpline3(m_h_efficiency) ;
  }
  if( m_sEfficiency != NULL ){
    m_model->setEfficiency(m_sEfficiency);
    m_cpmodel->setEfficiency(m_sEfficiency);
  }
  // Ensure the models are fully initialised.
  DalitzEvent evt(m_pattern, m_rndm) ;
  m_model->getVal(evt) ;
  m_cpmodel->getVal(evt) ;
  m_bothmodel.getVal(evt) ;
}

pair<double, double> TimeDependentGenerator::generate_decay_time() const {
  double decaytime = m_rndm->Exp(1./m_width) ;
  double smeareddecaytime = decaytime ;
  
  if ( (m_h_efficiency != NULL) && (m_addExpEffects) ){
    int i = 0 ;
    float efficiency = 0 ;
    int maxiter = 100000 ;
    
    smeareddecaytime =  decaytime + m_rndm->Gaus(0, m_resWidth);

    while(true){
      if( smeareddecaytime > m_efficiencyFit.GetXmax() ){
        efficiency = m_efficiencyFit.Eval( m_efficiencyFit.GetXmax() ) ;
      } 
      else if( smeareddecaytime > m_efficiencyFit.GetXmin() ){
	efficiency = m_efficiencyFit.Eval( smeareddecaytime ) ;
      }
      else{
        efficiency = m_efficiencyFit.Eval( m_efficiencyFit.GetXmin() ) ;
      }      

      if(m_rndm->Rndm() < efficiency){
        break;
      }
 
      decaytime = m_rndm->Exp(1./m_width) ;
      smeareddecaytime =  decaytime + m_rndm->Gaus(0, m_resWidth);

      i += 1 ;
      if( i > maxiter ){ 
        cout << "WARNING: Decay time generation limit exceeded." << endl ; 
        break ;
      }
    }
  }
  return pair<double, double>(decaytime, smeareddecaytime) ;
}

MINT::counted_ptr<IDalitzEvent> TimeDependentGenerator::generate_event(const double s13, const double s23) {
  int tag(0) ;
  double decaytime(0.) ;
  double smeareddecaytime(0.) ;

  MINT::counted_ptr<IDalitzEvent> evt(0) ;
  
  unsigned i = 0 ;
  unsigned maxattempts = 100000 ;
  // Accept/reject using the incoherent sum in m_bothmodel as an envelope to generate
  // events.
  while(true) {
    tag = m_rndm->Rndm() > 0.5 ? 1 : -1 ;
    tie(decaytime, smeareddecaytime) = generate_decay_time() ;
    evt = MINT::counted_ptr<IDalitzEvent>(new DalitzEvent(tag == 1 ? m_pattern : m_cppattern, s13, s23)) ;
    if(accept_or_reject(tag, decaytime, *evt))
      break ;
    i += 1 ;
    if( i >= maxattempts)
      throw range_error("Failed to generate an event in 100000 attempts!") ;
  }

  return MINT::counted_ptr<IDalitzEvent>(new GenTimeEvent(*evt, tag, decaytime, smeareddecaytime)) ;
}

MINT::counted_ptr<IDalitzEvent> TimeDependentGenerator::generate_event() {
  int tag(0) ;
  double decaytime(0.) ;
  double smeareddecaytime(0.) ;

  MINT::counted_ptr<IDalitzEvent> evt(0) ;
  
  unsigned i = 0 ;
  unsigned maxattempts = 100000 ;
  // Accept/reject using the incoherent sum in m_bothmodel as an envelope to generate
  // events.
  while(true) {
    tag = m_rndm->Rndm() > 0.5 ? 1 : -1 ;
    tie(decaytime, smeareddecaytime) = generate_decay_time() ;
    evt = m_generator.newEvent() ;
    if(accept_or_reject(tag, decaytime, *evt))
      break ;
    i += 1 ;
    if( i >= maxattempts)
      throw range_error("Failed to generate an event in 100000 attempts!") ;
  }

  return MINT::counted_ptr<IDalitzEvent>(new GenTimeEvent(*evt, tag, decaytime, smeareddecaytime)) ;
}

double TimeDependentGenerator::envelope_value(const double decaytime, IDalitzEvent& evt) {
  return m_bothmodel.Prob(evt) * m_scale * exp(-decaytime * m_width) * m_width ;
}

bool TimeDependentGenerator::accept_or_reject(const int tag, const double decaytime, IDalitzEvent& evt) {
  ++m_ngen ;

  double pdfval = pdf_value(tag, decaytime, evt) ;

  double maxval = envelope_value(decaytime, evt) ;
  if(pdfval > maxval){
    cout << "pdfval " << pdfval << " maxval " << maxval << endl ;
    //throw out_of_range("pdfval > maxval. That shouldn't happen!") ;
    m_scale *= pdfval / maxval ;
    cout << "scale " << m_scale << endl ;
    ++m_naccept ;
    return true ;
  }
  if(m_rndm->Rndm() * maxval < pdfval){
    ++m_naccept ;
    return true ;
  }
  return false ;
}

double TimeDependentGenerator::get_scale() const {
  return m_scale ;
}

float TimeDependentGenerator::get_gen_efficiency() const {
  return m_ngen > 0 ? float(m_naccept)/m_ngen : 0. ;
}

double TimeDependentGenerator::pdf_value(int tag, double decaytime, IDalitzEvent& evt) {
  complex<double> Ap = m_model->ComplexVal(evt) ;
  complex<double> Am = m_cpmodel->ComplexVal(evt) ;
  if(tag == -1)
    swap(Ap, Am) ;
  // This is gives the same values as the below using the amplitude coefficients.
  /*double magAp = norm(Ap) ;
  double magAm = norm(Am) ;
  complex<double> crossterm = pow(m_qoverp, tag) * conj(Ap) * Am ;
  double magqoverp = pow(norm(m_qoverp), tag) ;
  double deltamt = m_deltam * decaytime ;
  double halfdgammat = 0.5 * m_deltagamma * decaytime ;
  double pdfval = 0.5 * exp(-m_width * decaytime) 
  * ((magAp + magqoverp * magAm) * cosh(halfdgammat)
  + (magAp - magqoverp * magAm) * cos(deltamt)
  - 2 * crossterm.real() * sinh(halfdgammat)
  + 2 * crossterm.imag() * sin(deltamt)) ;
  return pdfval ;*/
  AmpPair coeffs = amplitude_coefficients(tag, decaytime) ;
  complex<double> Amp = coeffs.first * Ap + coeffs.second * Am ;
  double ampnorm = norm(Amp) ;
  return ampnorm ;
}

double TimeDependentGenerator::pdf_value(IDalitzEvent& evt) {
  return pdf_value(evt.getValueFromVector(GenTimeEvent::ITAG),
		   evt.getValueFromVector(GenTimeEvent::IDECAYTIME),
		   evt) ;
}

TimeDependentGenerator::AmpPair 
TimeDependentGenerator::amplitude_coefficients(const int tag, const double decaytime) {
  double coeff = 0.5 * exp(-decaytime * 0.5 * (m_width + 0.5 * m_deltagamma)) ;
  complex<double> expterm = exp(complex<double>(0.5 * m_deltagamma * decaytime, m_deltam * decaytime)) ;
  complex<double> plusterm = 1. + expterm ;
  complex<double> minusterm = 1. - expterm ;
  complex<double> coeffprod = coeff * plusterm ;
  complex<double> coeffmix = pow(m_qoverp, tag) * coeff * minusterm ;
  return AmpPair(coeffprod, coeffmix) ;
}

TH1F TimeDependentGenerator::draw_pdf_vs_time(IDalitzEvent& evt, unsigned nbins, float tmin, float tmax,
					      const std::string& name) {
  TH1F histo(name.c_str(), "", nbins, tmin, tmax) ;
  int tag = int(evt.getValueFromVector(GenTimeEvent::ITAG)) ;
  for(unsigned i = 1 ; i < nbins + 1 ; ++i){
    double time = histo.GetXaxis()->GetBinCenter(i) ;
    double val = pdf_value(tag, time, evt) ;
    histo.SetBinContent(i, val) ;
    histo.SetBinError(i, 0.) ;
  }
  return histo ;
}

TH1F TimeDependentGenerator::draw_envelope_vs_time(IDalitzEvent& evt, unsigned nbins, float tmin, float tmax,
						   const std::string& name) {
  TH1F histo(name.c_str(), "", nbins, tmin, tmax) ;
  for(unsigned i = 1 ; i < nbins + 1 ; ++i){
    double time = histo.GetXaxis()->GetBinCenter(i) ;
    double val = envelope_value(time, evt) ;
    histo.SetBinContent(i, val) ;
    histo.SetBinError(i, 0.) ;
  }
  return histo ;
}

