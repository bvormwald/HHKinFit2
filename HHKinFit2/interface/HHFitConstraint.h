/*
 * class for general fit constraint
 */

#ifndef HHFitConstraint_
#define HHFitConstraint_

#ifdef HHKINFIT2
#include "HHFitObject.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitObject.h"
#endif

namespace HHKinFit2{
class HHFitConstraint {
 public:
  HHFitConstraint(HHFitObject* fitobject);
  virtual ~HHFitConstraint() {};

  virtual void prepare(bool respectLimits=true);
  virtual double getChi2() const = 0;
  virtual double getLikelihood() const = 0;
  virtual void printChi2() const {};

 protected:
  HHFitObject* const m_fitobject;

};
}
#endif /* HHFitConstraint_ */
