#ifndef PSMATH_H_
#define PSMATH_H_

#include <Rtypes.h>

class PSMath {
 public:
  PSMath();
  static int PSfit(int iloop, int &iter, int &method, int &mode,
	      bool &noNewtonShifts, int printlevel,
	      int np, double a[], double astart[], double alimit[][2],
	      double aprec[], double daN[], double h[], double aMemory[][5],
	      double chi2, double chi2iter[], double g[], double H[], double Hinv[] );

  static void PSNewtonLimitShift(int sign, int np, double a[], double alimit[][2], double aprec[],
			  double daN[], double h[], double g[], double H[]);
  
  static double PSNewtonAnalyzer(int np, double a[], double alimit[][2], double aprec[],
			    double daN[], double h[],
			    double g[], double H[], double Hinv[], double chi2, bool noNewtonShifts, int printlevel=1);

  static void PSfitShow(int iloop, int convergence, int iter,
		 int method, int mode,
		 int printlevel, int graphiklevel,
		 int np, double a[], double astart[], double alimit[][2],
		 double aprec [], double daN[], double h[], double chi2,
		 double g[], double H[] );

  static double PSLineSearch(int & mode, double hh, double xlimit[],
			double epsx, double epsf, double x[4], double f[],
			double chi2 , int printlevel);

  static void PSLineLimit(int np, double astart[], double daN[], double alimit[][2], double xlimit[] );

  static double  PSVnorm(double x[], int n);

  static void PSVprint(const char* text, double x[], int n);
  static void PSMprint(const char* text, double A[], int ni, int nj);
  static void PSM2print(const char* text, double A[][2], int ni);

  static double PSMinverse(double H[], double Hinv[], int p);

  static double PSMCholtest();

  static double PSMmultiply(double A[],double B[],double C[],int n1,int n2);

  static double PSMmultiplyMRRT(double A[], int n1, int n2);

  static double PSMmultiplyMT(double A[], double B[], int n1, int n2);

  static double PSMRTrianInvert2(double R[], int n);

  static double PSMRTrianInvert(double R[], double Rinv[],  int n);
  
  static double PSMCholesky(double M[], double R[],  int n);

  static double PSfitCheckLimits(int np, double a[], double h[],
			    double alimit[][2], double aprecision[],
			    double daN[], double g[], double d);
  
  static double PSminIterate(double a[], double daN[], double h[], int p,
			double g[], double H[], double Hinv[], double X0);
  
  static double PSfuncQuadratic(double a[], double amean[], double F0,
			   double g[], double H[], int np);


 
  static int PSderivative(int icall, int np, double a[], double h[],
		     double chi2, double chi2iter[],
		     double g[], double H[]);

  static int PSderivative1(int icall, double a[], double h[],
			       double chi2, double g[], double H[]);

  static double PSfitMinStep(int np, double a[], double h[],
			double chi2iter[],
			double g[], double H[], double Hinv[], double daN[]);

  static int PSfitconstrain0(double F, double g, double H, double Fix, double aix[]);


};

#endif
