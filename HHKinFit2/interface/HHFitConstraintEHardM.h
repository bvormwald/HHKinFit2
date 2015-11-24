/*
 * class for constraint on 4 momentum conservation
 */

#ifndef HHFitConstraintEHardM_
#define HHFitConstraintEHardM_

#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObject.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObject.h"
#endif

#include <vector>

namespace HHKinFit2{
class HHFitConstraintEHardM : public HHFitConstraint {
 public:
  HHFitConstraintEHardM(HHFitObject* fitobject, HHFitObject* constrainedobject, double mass);

  void prepare(bool respectLimits=true);
  double getChi2() const;
  double getLikelihood() const;

 private:
  HHFitObject* const m_constrainedobject;
  double const m_mass;
  
};
}

#endif /* HHFitConstraintEHardM_ */
