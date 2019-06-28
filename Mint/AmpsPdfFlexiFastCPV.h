#ifndef __AMPSFLEXIFASTCPV_H__
#define __AMPSFLEXIFASTCPV_H__

#include <Mint/PdfBase.h>
#include <Mint/IDalitzPdf.h>
#include <Mint/IDalitzEvent.h>
#include <Mint/FitParameter.h>
#include <Mint/AmpsPdfFlexiFast.h>
#include <iostream>
#include <complex>

// Time-dependent PDF
class AmpsPdfFlexiFastCPV : public MINT::PdfBase<IDalitzEvent>
, virtual public IDalitzPdf{
    
protected:
    AmpsPdfFlexiFast* _amps1;
    AmpsPdfFlexiFast* _amps2;
    AmpsPdfFlexiFast* _ampsSum;
    
    MINT::FitParameter& _r;
    MINT::FitParameter& _delta;
    MINT::FitParameter& _gamma;
    
    MINT::FitParameter& _tau;
    MINT::FitParameter& _dGamma;
    MINT::FitParameter& _dm;
    MINT::FitParameter& _eff_tag;
    MINT::FitParameter& _w;
    
    double _intA;
    double _intAbar;
    std::complex<double> _intAAbar;
    
public:
    void parametersChanged(){
        _ampsSum->parametersChanged();
        _intA = _ampsSum->integralForMatchingPatterns(true,1);
        _intAbar = _ampsSum->integralForMatchingPatterns(true,-1);        
        _intAAbar = _ampsSum->ComplexSumForMatchingPatterns(false);
    }
    void beginFit(){
        _ampsSum->beginFit();
        printIntegralVals();
    }
    void endFit(){
        printIntegralVals();
        _ampsSum->endFit();
    }
    
    void printIntegralVals(){
        std::cout << "intSum = " << _ampsSum->getIntegralValue() << std::endl;
        std::cout << "intA = " << _intA << std::endl;
        std::cout << "intAbar = " << _intAbar << std::endl;
        std::cout << "intAAbar = " << _intAAbar.real() << std::endl;
    }
    
    inline double un_normalised_noPs(IDalitzEvent& evt){
        const double t = (double) evt.getValueFromVector(0);
        const double dt = (double) evt.getValueFromVector(1);
        const double q = static_cast<double>((int)evt.getValueFromVector(2)); 
        const double w = (double) evt.getValueFromVector(3);
        const double f = static_cast<double>((int)evt.getValueFromVector(4));
        
        const double e_eff = fabs(q)*_eff_tag + (1.-fabs(q))*(1.-_eff_tag);
        
        double r = (double)_r; // * sqrt(_intA/_intAbar);
        const std::complex<double> phase_diff = std::polar((double)r,((double) _delta -(double)_gamma*f)/360.*2*pi);
        
        const std::complex<double> amp = _amps1->ComplexVal_un_normalised_noPs(evt) ;
        const std::complex<double> amp_bar = _amps2->ComplexVal_un_normalised_noPs(evt) * phase_diff;
        
        const double val =  exp(-fabs(t)/(double)_tau) *
        (
         (2.-fabs(q))*(norm(amp) + norm(amp_bar))*cosh((double)_dGamma/2.*t)
         +f*q*(1.-2.*w)*(norm(amp) - norm(amp_bar)) *cos((double)_dm*t)
         -(2.-fabs(q))*2.0*real(amp_bar*conj(amp))*sinh((double)_dGamma/2.*t)
         -f*2.0*q*(1.-2.*w)*imag(amp_bar*conj(amp))*sin((double)_dm*t)
         )*e_eff;
        
        return val;
    }
    
    virtual double getVal(IDalitzEvent& evt){
        
        const double f = static_cast<double>((int)evt.getValueFromVector(4));
        const double val = un_normalised_noPs(evt);
        
        double r = (double)_r; // * sqrt(_intA/_intAbar);
        const double Gamma = 1./((double) _tau);
        const std::complex<double> phase_diff = std::polar((double)r,((double) _delta -(double)_gamma*f)/360.*2*pi);
        const double int_interference = (phase_diff*_intAAbar).real();
        
        if(_intA == -1 ){
            std::cout << "AmpsPdfFlexiFastCPV:: _norm = -1, should not have happened." << std::endl;
            throw "can't deal with that";
        }
        
        return val/(2.* ((_intA + r* r * _intAbar) * 4.*Gamma - int_interference * 2. * _dGamma )/ (4.*Gamma*Gamma-_dGamma*_dGamma));
    }
    
    virtual double getVal_withPs(IDalitzEvent& evt){return getVal(evt);}
    virtual double getVal_noPs(IDalitzEvent& evt){return getVal(evt);}
    
    virtual double getVal(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal(*evt);
    }
    virtual double getVal_withPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_withPs(*evt);
    }
    virtual double getVal_noPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_noPs(*evt);
    }
    
    virtual DalitzHistoSet histoSet(){return _ampsSum->histoSet();}
    
    void doFinalStatsAndSaveForAmp12(MINT::Minimiser* min=0,const std::string& fname = "FitAmpResults",
				     const std::string& fnameROOT="fitFractions"){
        _amps1->redoIntegrator();
        _amps2->redoIntegrator();
        _amps1->doFinalStatsAndSave(min,((std::string)fname+".txt").c_str(),((std::string)fnameROOT+".root").c_str());
        _amps2->doFinalStatsAndSave(min,((std::string)fname+"_CC.txt").c_str(),((std::string)fnameROOT+"_CC.root").c_str());        
    }
    
    AmpsPdfFlexiFastCPV(AmpsPdfFlexiFast* amps1, AmpsPdfFlexiFast* amps2, AmpsPdfFlexiFast* ampsSum, 
                        MINT::FitParameter& r, MINT::FitParameter& delta,MINT::FitParameter& gamma,
                        MINT::FitParameter& tau, MINT::FitParameter& dGamma, MINT::FitParameter& dm,
			MINT::FitParameter& eff_tag, MINT::FitParameter& w ):
    _amps1(amps1),
      _amps2(amps2),
      _ampsSum(ampsSum),
      _r(r),
      _delta(delta),
      _gamma(gamma),
      _tau(tau),
      _dGamma(dGamma),
      _dm(dm),
      _eff_tag(eff_tag),
      _w(w),
      _intA(-1),
      _intAbar(-1),
      _intAAbar(-1)
    {}
};


#endif
