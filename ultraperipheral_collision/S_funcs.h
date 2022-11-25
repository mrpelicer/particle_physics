#ifndef S_FUNCS_H
#define S_FUNCS_H

#include "MV_nucleus.h"
#include "MV_proton.h"

class Smatrix_class
{

public:

    //Smatrix_class(MV_proton_class & model_, double qMin_, double qMax_, double R_target_);
    Smatrix_class(MV_proton_class & model_, double qMin_, double qMax_, double R_target_);
    ~Smatrix_class(){};

    double operator()(double q_, int iN_);
  //Calculate r_t and b_t integrals for a given delta:
    void calc_integrals(double delta_);
  //Get dipole matrix in momentum space:
    double getSmatrix(double q_, int iN_) const;
    double Delta, qMin, qMax;
    MV_proton_class & MVmodel;
    double fm, R_target;

private:
    int rtbesselroots, btbesselroots;
    double bMin, bMax, rmin, rmaxInt, rmax;
    int ir_max, iq_max;
    double d_r, d_q;
    double epsb, epsr;
    std::vector <double> I0b_rgrid, I0r_qgrid, I1b_rgrid, I1r_qgrid, r, q;
    double r_grid, q_grid;



    //Integration on b:
    void intS_b(std::vector<double> & grid_, int iq_, int iN_);
    //Wrappers:
    static double gsl_S0b0(double x, void * p);
    static double gsl_S0b(double x, void * p);
    static double gsl_S1b0(double x, void * p);
    static double gsl_S1b(double x, void * p);
    //Main functions:
    double S0_b_func0(double x);
    double S0_b_func(double x);
    double S1_b_func0(double x);
    double S1_b_func(double x);
    //Interpolation:
    double Ib_interpolation(std::vector<double> grid_, double r_);
    //Grid:
    void gridIb(std::vector<double> &grid_, int iN_);
    double getIb_r(double r_, int iN_);

    //Integration on r:
    void intS_r(std::vector<double> & grid_, int iq_, int iN_);
    //Wrapper
    static double gsl_S0r0(double x, void * p);
    static double gsl_S0r(double x, void * p);
    static double gsl_S1r0(double x, void * p);
    static double gsl_S1r(double x, void * p);
    //Main functions:
    double S0_r_func0(double x);
    double S0_r_func(double x);
    double S1_r_func0(double x);
    double S1_r_func(double x);
    //Interpolation:
    double Ir_interpolation(std::vector<double> grid_, double q_) const;
    //Grid:
    void gridIr(std::vector<double> &grid_, int iN_);

};


#endif //Smatrix_class_H
