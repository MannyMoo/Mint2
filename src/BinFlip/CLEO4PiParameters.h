#ifndef __CLEO4PIPARAMETERS_H__
#define __CLEO4PIPARAMETERS_H__

#include <Mint/HadronicParameters.h>
#include <HyperHistogram.h>

class CLEO4PiBinning: public HadronicParameters::PhaseBinningBase {
 public:
  CLEO4PiBinning(const HyperHistogram&);
  virtual ~CLEO4PiBinning() {};
  virtual std::string type() const override;
  //virtual void Print(const std::string&, std::ostream& os = std::cout) const override ;
  //static MINT::counted_ptr<PhaseBinningBase> 
  //  fromConfig(const std::string&, const std::string& fname = "") ;
  
  //virtual bool operator==(const PhaseBinningBase&) const override ;

  virtual HadronicParameters::EventBinInfo binInfo(IDalitzEvent&) const override ;

  HadronicParameters::Bins getBins() const;
  
 protected :
  HyperHistogram m_hist;

};

#endif
