/*
 * class for fit objects with constant beta (=p/E)
 */

#ifndef HHFitObjectConstBeta_
#define HHFitObjectConstBeta_

#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/HHLorentzVector.h"
#include "HHKinFit2/HHKinFit2/HHFitObjectE.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHLorentzVector.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectE.h"
#endif

#include "TMatrixD.h"

namespace HHKinFit2{
class HHFitObjectEConstBeta : public HHFitObjectE {
 public:
  HHFitObjectEConstBeta(HHLorentzVector const& initial4vector);
  HHLorentzVector constrainEtoMinv(double m, HHLorentzVector const& other4vector) const;
  double calculateEConstrainedToMinv(double m, HHLorentzVector const& other4vector) const;
  HHLorentzVector changeE(double E) const;

  void print() const;

};
}
#endif /* HHFitObjectConstBeta_ */
