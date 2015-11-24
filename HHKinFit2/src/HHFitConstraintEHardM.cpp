#ifdef HHKINFIT2
#include "HHFitConstraintEHardM.h"
#include "HHFitObjectE.h"
#include "exceptions/HHEnergyRangeException.h"
#include "exceptions/HHInvMConstraintException.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraintEHardM.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectE.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHEnergyRangeException.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHInvMConstraintException.h"
#endif

HHKinFit2::HHFitConstraintEHardM::HHFitConstraintEHardM(HHFitObject* fitobject, HHFitObject* constrainedobject, double mass)
  : HHFitConstraint(fitobject),
    m_constrainedobject(constrainedobject),
    m_mass(mass){

}

double
HHKinFit2::HHFitConstraintEHardM::getChi2() const{
  return(0);
}

double
HHKinFit2::HHFitConstraintEHardM::getLikelihood() const{
  return(1);
}

void
HHKinFit2::HHFitConstraintEHardM::prepare(bool respectLimits){
  HHFitObjectE* new4momentum2 = static_cast<HHFitObjectE*>(m_constrainedobject);
  try{
    new4momentum2->constrainEtoMinvandSave(m_mass, m_fitobject->getFit4Vector(),respectLimits);
  }
  catch(HHEnergyRangeException const& e){
    throw(HHInvMConstraintException(e.what()));
  }
}
