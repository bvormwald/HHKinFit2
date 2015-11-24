/*
 * class for fit objects
 */

#ifndef HHFitObjectMET_
#define HHFitObjectMET_

#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/HHFitObject.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitObject.h"
#endif

#include "TVector2.h"
#include "TMatrixD.h"

namespace HHKinFit2{
class HHFitObjectMET : public HHFitObject {
 public:
  HHFitObjectMET(TVector2 const& initialvector,TVector2 const& fitvector=TVector2(0,0));
  
  void setCovMatrix(double xx, double yy, double xy);
  void setCovMatrix(TMatrixD const covmat);

  virtual void print() const;
  virtual void printInitial4Vector() const;
  virtual void printFit4Vector() const;
  virtual void printCovMatrix() const;
};
}
#endif /* HHFitObjectMET_ */
