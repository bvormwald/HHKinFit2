/*
 * class for general fit constraint
 */

#ifndef HHKinFit_
#define HHKinFit_

#include "TGraph.h"
#include <vector>

#ifdef HHKINFIT2
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectE.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitObjectE.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint.h"
#endif

namespace HHKinFit2{
class HHKinFit {
 public:
  HHKinFit();

  void setPrintLevel(int printlevel=0);

  void fit();
  double getChi2(bool respectLimits=true) const;
  double getL(bool respectLimits=true) const;
  std::vector<HHFitObjectE*> getListOfFitObjects() const;
  std::vector<HHFitConstraint*> getListOfConstraints() const;
  
  void addFitObjectE(HHFitObjectE* fitobject);
  void addConstraint(HHFitConstraint* constraint);
  
  TGraph* getChi2Function(int steps);
  TGraph* getLFunction(int steps);
  int getConvergence() const;

  void printChi2() const;


 private:
  std::vector<HHFitObjectE*> m_fitobjects;
  std::vector<HHFitConstraint*> m_constraints;

  double m_chi2;
  int m_convergence;
  int m_printlevel;
  int m_maxloops;
};
}
#endif /* HHKinFit_ */
