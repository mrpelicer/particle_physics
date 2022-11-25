#ifndef MV_Nucleus_H
#define MV_Nucleus_H

//Standard & math libs:
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdio>

//gsl_libraries:
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

//Boost libraries for bessel zeros:
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/bessel.hpp>


class MV_nucleus_class
{
public:
//Constructor and destructor:
  MV_nucleus_class(double fm_, double R_, double mg_, double Qs2_, double A_, double RA_);
  ~MV_nucleus_class(){};
//Main function
  double operator()(double r_, double b_, int iN_) const;
  double calcN0(double r_, double b_) const;
  double calcNe(double r_, double b_) const;
//Public parameters - may be passed to other objects!
  double fm, R, mg, Qs2, A, RA, delta, T_norm;
  double bMin, bMax;

private:

  double xMin, xMax, x_div;
  int ix_Max, ib_Max;
  double d_x, d_b;
  int Maxroots;
  std::vector <double> b, x, n0_i, ne_i, TA, N0A_1, N0A_2, NeA_1;
  int ib_grid;

//Params for integration:
  struct n0_params
  {
    double Rp, mgp, xp;
  };

  struct T_params
  {
    double deltap, Ap, RAp, NAp, bp;
  };

//Calculates N0^A and Ne^A
  void gridNA(std::vector<double> & NA_, int iN_);
  void intNA(std::vector<double> & NA_, int iA_);
  static double gsl_N0A_func1(double x_, void * p);
  static double gsl_N0A_func2(double x_, void * p);
  static double gsl_NeA_func(double x_, void * p);
  double N0A_func1(double x_);
  double N0A_func2(double x_);
  double NeA_func(double x_);

//Calculates thickness and derivatives:
  void gridThickness(std::vector<double> & T_);
  void intThickness(std::vector<double> & T_, double b_, int ib_) const;
  static double TA_function(double z, void * p);
  double itpThickness_NAfuncs(std::vector<double> vec_, double b_) const;
  double derivativeThickness(std::vector<double> TA_, double b_, int order_) const;

//Calculate N0 and Ne for proton:
  void grid_ni(std::vector<double> & ni_, int i_);
  double itp_ni(std::vector<double> ni_, double x_) const;
  void int_n0(std::vector<double> &ni, double x_, int ix_) const;
  void int_ne(std::vector<double> &ni, double x_, int ix_) const;
  static double n0_func(double d_, void * p);
  static double n0_func2(double d_, void * p);
  static double ne_func(double d_, void * p);
  static double ne_func2(double d_, void * p);
};

#endif // MV_Funcs_H.
