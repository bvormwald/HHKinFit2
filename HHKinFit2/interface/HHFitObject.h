/*
 * class for fit objects
 */

#ifndef HHFitObject_
#define HHFitObject_

#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/interface/HHLorentzVector.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHLorentzVector.h"
#endif

#include "TMatrixD.h"

namespace HHKinFit2{
class HHFitObject {
 public:
  HHFitObject();
  HHFitObject(HHLorentzVector const initial4vector);
  virtual ~HHFitObject(){};

  virtual HHLorentzVector getInitial4Vector() const;
  virtual HHLorentzVector getFit4Vector() const;
  virtual TMatrixD getCovMatrix() const;

  void setFit4Vector(HHLorentzVector const vec);
  virtual void setCovMatrix(TMatrixD const covmatrix);

  virtual void reset();

  virtual void printInitial4Vector() const;
  virtual void printFit4Vector() const;
  virtual void printCovMatrix() const;

  virtual void print() const;

 protected:
  HHLorentzVector m_fit4vector;
  HHLorentzVector const m_initial4vector;
  TMatrixD m_covmatrix;
};
}
#endif /* HHFitObject_ */
