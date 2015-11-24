/*
 * class for soft energy limit
 */

#ifndef HHFitConstraintSoftBoundary_
#define HHFitConstraintSoftBoundary_

#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2/HHFitObject.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObject.h"
#endif

#include <vector>

namespace HHKinFit2{
class HHFitConstraintSoftBoundary : public HHFitConstraint {
 public:
  HHFitConstraintSoftBoundary(HHFitObject* object, double percent);
  double getChi2() const;
  double getLikelihood() const;

 private:
  double m_percent;
};
}
#endif /* HHFitConstraintSoftBoundary_ */
