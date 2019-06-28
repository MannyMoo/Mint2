#ifndef __AMPSPDFFLEXIFAST_H__
#define __AMPSPDFFLEXIFAST_H__

#include <TRandom3.h>
#include <Mint/SignalGenerator.h>
#include <Mint/FromFileGenerator.h>
#include <Mint/IDalitzEvent.h>
#include <Mint/IEventGenerator.h>
#include <Mint/NamedParameter.h>
#include <string>
#include <Mint/DalitzEventPattern.h>
#include <Mint/DalitzPdfBaseFlexiFastInteg.h>

/*
Seems to be a way to top up integrator events for an amplitude ...
*/

class AmpsPdfFlexiFast
: public DalitzPdfBaseFlexiFastInteg
{
protected:
    TRandom* _localRnd;
    SignalGenerator* _sgGen;
    FromFileGenerator* _fileGen;
    MINT::IEventGenerator<IDalitzEvent>* _chosenGen;
    MINT::NamedParameter<std::string> _integratorSource;
    std::string _integratorFileName;
public:
    double un_normalised_noPs(IDalitzEvent& evt){
        double ampSq =  _amps->RealVal(evt);
        return ampSq;// * getEvent()->phaseSpace();
    }
    
    std::complex<double> ComplexVal_un_normalised_noPs(IDalitzEvent& evt){
        return  _amps->ComplexVal(evt);
    }
    
    AmpsPdfFlexiFast(const DalitzEventPattern& pat
                     , IFastAmplitudeIntegrable* amps
                     , MINT::MinuitParameterSet* mps
                     , double precision=1.e-4
                     , std::string method="efficient"
                     , std::string fname =  "SignalIntegrationEvents.root", bool genMoreEvents = false
                     )
    : DalitzPdfBaseFlexiFastInteg(pat, 0, amps, precision, mps)
    , _localRnd(0)
    , _sgGen(0)
    , _fileGen(0)
    , _chosenGen(0)
    , _integratorSource("IntegratorSource", (std::string) "new", (char*) 0)
    , _integratorFileName(fname)
    {
        std::cout << " AmpsPdfFlexiFast with integ method " << method << std::endl;
        bool nonFlat = "efficient" == method;
        bool generateNew = ((std::string)_integratorSource == (std::string)"new");
        if(nonFlat){
            std::cout << "AmpsPdfFlexiFast uses nonFlat integration." << std::endl;
            if(generateNew){
                _sgGen =  new SignalGenerator(pat);
                _sgGen->setWeighted();
                //_sgGen->setSaveEvents(_integratorFileName);
                _chosenGen = _sgGen;
            }else{
                // here, SignalGenerator is used by FromFileGenerator, to fill
                // up missing events in case more are needed than found in the
                // file.  Since we don't know with which random seed the
                // events in the file were generated, we supply a random
                // number generator with randomised seed.
                _localRnd = new TRandom3(time(0));
                _sgGen =  new SignalGenerator(pat, _localRnd);
                _sgGen->setWeighted();
                _sgGen->dontSaveEvents();// saving events is done by FromFileGenerator
                if(genMoreEvents) _fileGen   = new FromFileGenerator(_integratorFileName, _sgGen);
                else{
                    _fileGen = new FromFileGenerator(_integratorFileName, 0, "UPDATE");
                    std::cout << "not going to generate any more events" << std::endl;
                }
                _chosenGen = _fileGen;
            }
            this->setEventGenerator(_chosenGen);
        }else{
            std::cout << "AmpsPdfFlexiFast uses flat integration." << std::endl;
        }
    }
    
    IFastAmplitudeIntegrable* getAmpSum(){ return _amps;}
    
    ~AmpsPdfFlexiFast(){
        if(0 != _fileGen)  delete _fileGen;
        if(0 != _sgGen)    delete _sgGen;
        if(0 != _localRnd) delete _localRnd;
    }
};

#endif
