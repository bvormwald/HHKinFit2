/*
 * class for additional likelihood-based constraints
 */

#ifndef HHFitConstraintLikelihood_
#define HHFitConstraintLikelihood_

#ifdef HHKINFIT2
#include "HHFitConstraint.h"
#include "HHFitObject.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObject.h"
#endif

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSpline.h"


namespace HHKinFit2{
class HHFitConstraintLikelihood : public HHFitConstraint {
 public:
  HHFitConstraintLikelihood(HHFitObject* object1,HHFitObject* object2, TF1* likelihood1, TF1* likelihood2 );
  HHFitConstraintLikelihood(HHFitObject* object1,HHFitObject* object2, TH1D* likelihoodhisto1, TH1D* likelihoodhisto2 );
  HHFitConstraintLikelihood(HHFitObject* object1,HHFitObject* object2, TH2D* likelihoodhisto);
  HHFitConstraintLikelihood(HHFitObject* object1,HHFitObject* object2, TSpline3* likelihoodspline1,  TSpline3* likelihoodspline2);
  Double_t getChi2() const;
  double getLikelihood() const;

 private:
  TF1* const m_likelihood1;
  TF1* const m_likelihood2;
  TH1D* const m_likelihoodhisto1;
  TH1D* const m_likelihoodhisto2;
  TH2D* const m_likelihoodhisto;
  TSpline3* m_likelihoodspline1;
  TSpline3* m_likelihoodspline2;
  HHFitObject* m_object2;
  int mode;

  //for splines
  double m_zmax1;
  double m_zcutleft1;
  double m_zcutright1;
  double m_zmax2;
  double m_zcutleft2;
  double m_zcutright2;
};
}
#endif /* HHFitConstraintLikelihood_ */
