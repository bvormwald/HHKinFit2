#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/HHFitConstraint.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint.h"
#endif

HHKinFit2::HHFitConstraint::HHFitConstraint(HHFitObject* fitobject)
  :m_fitobject(fitobject){

}

void
HHKinFit2::HHFitConstraint::prepare(bool respectLimits){
}
