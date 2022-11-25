#ifndef MV_Proton_H
#define MV_Proton_H

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
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

//Boost libraries for bessel zeros:
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/bessel.hpp>

class MV_proton_class
{
public:

//Constructor and destructor:
    MV_proton_class(double fm_, double R_, double mg_, double Qs2_);
    ~MV_proton_class(){};

//Main function
    double operator()(double r_, double b_, int iN_) const;

    double fm, R, mg, Qs2;
    double bMin, bMax;

    void freeVecs();

private:

//Main Parameters
    int ib_Max;
    double b_div, d_b;
    int Maxroots;
    std::vector <double> b;
    std::vector <double> n0_i;
    std::vector <double> err_0;
    std::vector <double> ne_i;
    std::vector <double> err_e;

    //Interpolated b_
    std::vector<double> b_itp;
    //Params for integration:
    struct N0_params
    {
      double Rp, mgp, bp;
    };

//Grid integration
    void Grid_Ni(std::vector<double> & ni_, std::vector<double> & err_, int i_);
//Interpolation for the grid above:
    double ni_itp(std::vector<double> ni_, double b_) const;
//Calculate N0 and Ne:
    double calc_N0(double r_, double b_) const;
    double calc_Ne(double r_, double b_) const;
    void intN0(std::vector<double> &ni, std::vector<double> &err_, double b_, int ib_) const;
    void intNe(std::vector<double> &ni, std::vector<double> &err_, double b_, int ib_) const;

    static double N0(double x, void * p);
    static double N0_2(double x, void * p);
    static double Ne(double x, void * p);
    static double Ne_2(double x, void * p);

};

#endif // MV_Funcs_H.
