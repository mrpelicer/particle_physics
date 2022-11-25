//----------------------------------------------------------------
// Author: Emmanuel G. de Oliveira
// emmanuel.de.oliveira@ufsc.br
// fragfun.cpp
// Last modified: May 5st, 2017
// fragmentation function
// follows hep-ph/0510032
//----------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <vector>
//----------------------------------------------------------------
class fragfun_class {
public:
  virtual double operator()(double const x) const { return 0e0; }
  virtual double at1() const { return 0e0; }
};

//----------------------------------------------------------------
class fragfunCNO_class : public fragfun_class {
public:
  double norm, c, a, b;
  fragfunCNO_class(double norm, double c, double a, double b)
    : norm(norm), c(c), a(a), b(b) {}

  double operator()(double const x) const override {
    return norm * c/(1e0 + c) * pow(1e0-x,a) * pow(x,b)
        / ( tgamma(1e0 + a) * tgamma(1e0 + b) / tgamma(2e0 + a + b) ); // a, b >-1 !!
  }

  double at1() const override {
    return norm/(1+c);
  }
};

//----------------------------------------------------------------
class Dplusff_class : public fragfun_class {
public:
  fragfunCNO_class Ddirect_ff, Dstar_ff;
  double branch_Dplus_pion, branch_Dplus_gamma;
  double Dstar_mass, D_mass, p_D;
  Dplusff_class() :
//    Ddirect_ff(0.263e0/2+0.270e0/2, 4.6e0, 1.1e0, 7.6e0),
    Ddirect_ff(0.173e0, 4.6e0, 1.1e0, 7.6e0),
    Dstar_ff(0.227e0, 2.46e0, 1.8e0, 11.3e0),
    branch_Dplus_pion(0.307e0),
    branch_Dplus_gamma(0.016e0),
    Dstar_mass(2010e0),
    D_mass(1869.3e0)
  {
    p_D = (Dstar_mass * Dstar_mass - D_mass*D_mass)/(2*Dstar_mass);
  }

  double operator()(double const x) const override {
    double res1 = Ddirect_ff(x);

    double res2 = 0e0;
    double x_resc = x * Dstar_mass/D_mass;
    if (x_resc < 1e0)
      res2 = Dstar_ff(x_resc) * Dstar_mass/D_mass;

    double res3 = 0e0;
    for (int i = 0; i<1000; i++) {
      double y = (i+.5e0)/1000;
      double cos_theta = (x*Dstar_mass/y - D_mass)/p_D;
      //std::cout << y << " " << cos_theta << std::endl;
      if (std::abs(cos_theta)<1e0)
        res3 += (1e0/1000) * Dstar_ff(y) / (2* y * p_D /Dstar_mass);
    }

    double cos_theta_at1 = (x*Dstar_mass - D_mass)/p_D;
    //std::cout << 1 << " " << cos_theta_at1 << std::endl;
    if (std::abs(cos_theta_at1)<1e0)
      res3 += Dstar_ff.at1() / (2 * p_D / Dstar_mass);

    return res1
        + branch_Dplus_pion*res2
        + branch_Dplus_gamma*res3;
  }

  double at1() const override {
    return Ddirect_ff.at1();
  }

  double at1resc() const {
    return branch_Dplus_pion * Dstar_ff.at1();
  }

  double resc() const {
    return D_mass/Dstar_mass;
  }
/*

  double memo(double const x) const {
    static constexpr int N_points = 100001;
    static bool ready = false;
    static std::vector<double> values(N_points, 0e0);
    if (!ready) {
      for (int i = 0; i<N_points; i++) {
        double y = (1e0/N_points) * i;
        values[i] = operator()(y);
      }
      ready = true;
    }

    int j = x*(N_points-1);
    double delta_x = x - (j*1e0)/N_points;
    return values[j] + delta_x * (values[j+1] - values[j]) * (N_points - 1);
  }

  void test() const {
    for (int i = 0; i<100; i++) {
      double x = (i+0.5e0)/100;
      std::cout << x << " " << Dstar_ff(x) << " " << operator()(x) << std::endl;
    }
  }
*/

};


//----------------------------------------------------------------
class Dzeroff_class : public fragfun_class {
public:
  fragfunCNO_class Ddirect_ff, Dstar_ff;
  double branch_Dzero_pion, branch_Dzero_gamma;
  double Dstar_mass, D_mass, p_D;
  Dzeroff_class() :
//    Ddirect_ff(0.609e0/2+0.598e0/2, 4.6e0, 1.1e0, 7.6e0),
    Ddirect_ff(.141e0, 4.6e0, 1.1e0, 7.6e0),
    Dstar_ff(0.253e0, 2.46e0, 1.8e0, 11.3e0),
    branch_Dzero_pion(0.677e0 + 0.619e0),
    branch_Dzero_gamma(0.381e0),
    Dstar_mass(2006.7e0),
    D_mass(1864.5e0)
  {
    p_D = (Dstar_mass * Dstar_mass - D_mass*D_mass)/(2*Dstar_mass);
  }

  double operator()(double const x) const override {
    double res1 = Ddirect_ff(x);

    double res2 = 0e0;
    double x_resc = x * Dstar_mass/D_mass;
    if (x_resc < 1e0)
      res2 = Dstar_ff(x_resc) * Dstar_mass/D_mass;

    double res3 = 0e0;
    for (int i = 0; i<1000; i++) {
      double y = (i+.5e0)/1000;
      double cos_theta = (x*Dstar_mass/y - D_mass)/p_D;
      //std::cout << y << " " << cos_theta << std::endl;
      if (std::abs(cos_theta)<1e0)
        res3 += (1e0/1000) * Dstar_ff(y) / (2* y * p_D /Dstar_mass);
    }

    double cos_theta_at1 = (x*Dstar_mass - D_mass)/p_D;
    //std::cout << 1 << " " << cos_theta_at1 << std::endl;
    if (std::abs(cos_theta_at1)<1e0)
      res3 += Dstar_ff.at1() / (2 * p_D / Dstar_mass);

    return res1
        + branch_Dzero_pion*res2
        + branch_Dzero_gamma*res3;
  }

  double at1() const override {
    return Ddirect_ff.at1();
  }

  double at1resc() const {
    return branch_Dzero_pion * Dstar_ff.at1();
  }

  double resc() const {
    return D_mass/Dstar_mass;
  }
/*
  double memo(double const x) const {
    static constexpr int N_points = 100001;
    static bool ready = false;
    static std::vector<double> values(N_points, 0e0);
    if (!ready) {
      for (int i = 0; i<N_points; i++) {
        double y = (1e0/N_points) * i;
        values[i] = operator()(y);
      }
      ready = true;
    }

    int j = x*(N_points-1);
    double delta_x = x - (j*1e0)/N_points;
    return values[j] + delta_x * (values[j+1] - values[j]) * (N_points - 1);
  }

  void test() const {
    for (int i = 0; i<100; i++) {
      double x = (i+0.5e0)/100;
      std::cout << x << " " << Dstar_ff(x) << " " << operator()(x) << std::endl;
    }
  }
*/

};
