#include "TLorentzVector.h"

#ifndef HHLorentzVector_
#define HHLorentzVector_

namespace HHKinFit2{
  class HHLorentzVector : public TLorentzVector {

  public:
    HHLorentzVector(double x, double y, double z, double t);
    HHLorentzVector();
    void SetEEtaPhiM(double E,double Eta,double Phi,double M);
    void SetEkeepM(double E);
    void SetEkeepBeta(double E);
    void SetMkeepE(double m);
    HHLorentzVector operator+(HHLorentzVector const& rhs) const;
    HHLorentzVector operator-(HHLorentzVector const& rhs) const;
    HHLorentzVector operator-() const;	
    HHLorentzVector &operator=(HHLorentzVector const& other);
  };
}
#endif /* HHLorentzVector_ */
