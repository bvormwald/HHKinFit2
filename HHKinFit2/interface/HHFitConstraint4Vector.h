/*
 * class for constraint on 4 momentum conservation
 */

#ifndef HHFitConstraint4Vector_
#define HHFitConstraint4Vector_

#ifdef HHKINFIT2
#include "HHFitConstraint.h"
#include "HHFitObject.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObject.h"
#endif

#include <vector>

namespace HHKinFit2{
class HHFitConstraint4Vector : public HHFitConstraint {
 public:
  HHFitConstraint4Vector(HHFitObject* object, bool px, bool py, bool pz, bool E);
  double getChi2() const;
  double getLikelihood() const;

  void printChi2() const;

 private:
  int m_ncomp;
  TMatrixD m_cov;
  TMatrixD m_invcov;
  std::vector<int> m_indices;
  bool m_components[4];
};
}
#endif /* HHFitConstraint4Vector_ */
