/*
 * class for fit objects
 */

#ifndef HHFitObjectE_
#define HHFitObjectE_

#ifdef HHKINFIT2
#include "HHLorentzVector.h"
#include "HHFitObject.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHLorentzVector.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObject.h"
#endif

#include "TMatrixD.h"

namespace HHKinFit2{
class HHFitObjectE : public HHFitObject {
 public:
  HHFitObjectE(HHLorentzVector const& initial4vector);
  
  double getE() const;
  virtual HHLorentzVector changeE(double E) const = 0;
  void changeEandSave(double E, bool respectLimits=true);
  virtual HHLorentzVector scaleE(double scale) const;
  void scaleEandSave(double scale, bool respectLimits=true);
  virtual HHLorentzVector constrainEtoMinv(double minv, HHLorentzVector const& other4vector) const =0;
  virtual double calculateEConstrainedToMinv(double m, HHLorentzVector const& other4vector) const =0;
  void constrainEtoMinvandSave(double minv, HHLorentzVector const& other4vector, bool respectLimits=true);

  double getUpperFitLimitE() const;
  double getLowerFitLimitE() const;

  void setFitLimitsE(double const lowerlimit, double const upperlimit);
  void setFitLimitsE(HHLorentzVector const& own4vectorMin, double const minv, HHLorentzVector const& other4vectorMin);
  void setUpperFitLimitE(double const upperlimit);
  void setUpperFitLimitE(double const minv, HHLorentzVector const& other4vectorMin);
  void setLowerFitLimitE(double const lowerlimit);
  void setLowerFitLimitE(double const minv, HHLorentzVector const& other4vectorMin);
  void setLowerFitLimitE(HHLorentzVector const& other4vectorMin);

  void setInitStart(double start);
  void setInitPrecision(double prec);
  void setInitDirection(double daN);
  void setInitStepWidth(double h);
  double getInitStart();
  double getInitPrecision();
  double getInitDirection();
  double getInitStepWidth();

  void setCovMatrix(double dE);
  void setCovMatrix(TMatrixD cov);

  void reset();
  void resetLimits();

  virtual void print() const;
  void printLimits() const;

 private:
  double m_upperLimitE;
  double m_lowerLimitE;
  double m_initstart;
  double m_initprec;
  double m_initstep;
  double m_initdirection;
};
}
#endif /* HHFitObjectE_ */
