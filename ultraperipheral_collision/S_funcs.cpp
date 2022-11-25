//This object makes a Hankel transform of orders 0 and 2 of the dipole-S matrix in configuration space.
//The McLerran-Venugopalan model is used.
//Integration with 0-th order Bessel function results in the component S_0 of the S-matrix in momentum space
// iN_ = 0 -> S_0
//Integration with 2-th order Bessel function results in the component S_e of the S-matrix in momentum space.
// iN_ = 1 -> S_e

#include "S_funcs.h"

Smatrix_class::Smatrix_class(MV_proton_class & model_, double qMin_, double qMax_, double R_target_):
   qMin(qMin_), qMax(qMax_),
   MVmodel(model_),
   fm(MVmodel.fm), R_target(R_target_),
   rtbesselroots(60), btbesselroots(150), bMin(0.1), bMax(20.*fm), rmin(0.0),
   rmaxInt(7.*fm),
   rmax(std::max(boost::math::cyl_bessel_j_zero(0.0, rtbesselroots), boost::math::cyl_bessel_j_zero(2.0, rtbesselroots) ) ),
   ir_max(1000), iq_max(200),
   d_r((rmax-rmin)/(double)ir_max), d_q((qMax - qMin)/(double)iq_max),
   epsb(1./(R_target*R_target)), epsr(1./(0.5*0.5*fm*fm)),
   I0b_rgrid(ir_max+1), I0r_qgrid(iq_max+1), I1b_rgrid(ir_max+1), I1r_qgrid(iq_max+1),
   r(ir_max+1), q(iq_max+1)
{
  std::cout << "Ranges for S-matrix:" << std::endl
            << "bmax, rmax, qmax = " << bMax << " " << rmax << " " << qMax << std::endl
            << "dr, dq = " << d_r << " " << d_q << std::endl
            << "ir, iq (max)= " << ir_max << " " << iq_max << std::endl
            << "===================================================" << std::endl;

}


double Smatrix_class::operator()(double q_, int iN_)
{
    double res = getSmatrix(q_, iN_);
    return res;
}


void Smatrix_class::calc_integrals(double delta_)
{
  Delta = delta_;

  //Make integral in b_t, griding in r_t:
  gridIb(I0b_rgrid, 0);
  gridIb(I1b_rgrid, 1);

  //Make integral in r_t, griding in q_t:
  gridIr(I0r_qgrid, 0);
  gridIr(I1r_qgrid, 1);

  MVmodel.freeVecs();

  std::ofstream T0matrix, Tematrix;
  std::string nameT0, nameTe;
  nameT0=  "T0_" + std::to_string(Delta) + ".dat";
  nameTe=  "T1_" + std::to_string(Delta) + ".dat";
  T0matrix.open(nameT0);
  Tematrix.open(nameTe);

  for(int iq = 0; iq<= iq_max; iq++)
  {
    T0matrix << qMin+d_q*iq << " " << Delta << " "<< I0r_qgrid[iq] << std::endl;
    Tematrix << qMin+d_q*iq << " " << Delta << " "<< I1r_qgrid[iq] << std::endl;
  }

  T0matrix.close();
  Tematrix.close();

}

void Smatrix_class::gridIb(std::vector<double> &grid_, int iN_)
{
  std::cout << "Integrating S_"<<iN_<<" in b_t -- Griding in r_t." << std::endl;
  for(int ir = 0; ir <= ir_max; ir++)
  {
    if(ir%10 == 0) std::cout << ir/10 << "/" << (int)ir_max/10 << "\t" <<  std::flush;

    r_grid = rmin + d_r*ir;
    r[ir] = r_grid;
    intS_b(grid_, ir, iN_);
  }
  std::cout << std::endl << "===================================================" << std::endl;
}

void Smatrix_class::gridIr(std::vector<double> &grid_, int iN_)
{
  std::cout << "Integrating S_"<<iN_<<" in r_t -- Griding in q_t." << std::endl;
  for(int iq = 0; iq <= iq_max; iq++)
  {
    if(iq%10 == 0) std::cout << iq/10 << "/" << (int)iq_max/10 << "\t" <<  std::flush;

    q_grid = qMin + d_q*iq;
    q[iq] = q_grid;
    intS_r(grid_, iq, iN_);
  }
  std::cout << std::endl << "===================================================" << std::endl;
}

double Smatrix_class::getIb_r(double r_, int iN_)
{
    double res = 0e0;
    if(iN_ == 0) res = Ib_interpolation(I0b_rgrid, r_);
    if(iN_ == 1) res = Ib_interpolation(I1b_rgrid, r_);

    return res;
}

double Smatrix_class::getSmatrix(double q_, int iN_) const
{//Returns value of the integration in r
    double res = 0e0;
    if(iN_ == 0) res = Ir_interpolation(I0r_qgrid, q_);
    if(iN_ == 1) res = Ir_interpolation(I1r_qgrid, q_);
    return res;
}

double Smatrix_class::Ib_interpolation(std::vector<double> grid_, double r_)
{//Makes the interpolation of a vector grid_ in the variable r_.
    double res = 0.0;

    if(r_ > rmax)
    {
        std::cout << "Error!  Smatrix_class::Ib_interpolation calling r > rmax!"  << std::endl;
        exit(EXIT_FAILURE);
    }

    gsl_interp_accel *acc =  gsl_interp_accel_alloc();
    gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_linear, ir_max+1);

    gsl_interp_init(interpolation, r.data(), grid_.data(), ir_max+1); //vector.data() transforms std::vector<double> to array!
    res = gsl_interp_eval(interpolation, r.data(), grid_.data(), r_, acc);

    gsl_interp_free(interpolation);
    gsl_interp_accel_free (acc);
    return res;
}

double Smatrix_class::Ir_interpolation(std::vector<double> grid_, double q_) const
{
    double res = 0.0;

    if(q_ > qMax)
    {
        std::cout << "Error!  Smatrix_class::Ir_interpolation calling q > qMax!"  << std::endl;
        exit(EXIT_FAILURE);
    }

    gsl_interp_accel *acc =  gsl_interp_accel_alloc();
    gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_linear, iq_max+1);

    gsl_interp_init(interpolation, q.data(), grid_.data(), iq_max+1);

    res = gsl_interp_eval(interpolation, q.data(), grid_.data(), q_, acc);

    gsl_interp_free(interpolation);
    gsl_interp_accel_free (acc);

    return res;
}


void Smatrix_class::intS_b(std::vector<double> & grid_, int ir_, int iN_)
{
    //std::cout.precision(std::numeric_limits<double>::digits );

    double result = 0e0;
    double er_res = 0e0;
    double er_tmp = 0e0;
    double lower_limit = bMin;
    double upper_limit;
    double err_abs = 0e0;
    double err_rel = 5e-3;
    double deltaDiv =1.;

    size_t max_steps = 1e8;
    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(max_steps);

    gsl_function G;
    G.params = this;
    //double res_compare = 0.;
    double res_tmp;
    //int itest = 0;

    if(iN_ == 0)
    {
      if(Delta<=deltaDiv)
      {
         G.function = &Smatrix_class::gsl_S0b0;
         gsl_integration_qags(&G,  lower_limit, bMax, err_abs, err_rel, max_steps, w, &result, &er_res);
      }else{
        //itest=1;
         err_rel = 1e-3;
         G.function = &Smatrix_class::gsl_S0b;
         for(int counter = 1; counter <= btbesselroots; counter++)
         {
           upper_limit = boost::math::cyl_bessel_j_zero(0.0, counter);
           //res_compare = result;

            gsl_integration_qags(&G, lower_limit, upper_limit, err_abs, err_rel, max_steps, w, &res_tmp, &er_tmp);

            lower_limit = upper_limit;
            result+= res_tmp;
            er_res+=er_tmp;

        }
        //result = (result+res_compare)/2.;
     }
   }
   if(iN_ == 1)
   {
    if(Delta<=deltaDiv)
      {
         G.function = &Smatrix_class::gsl_S1b0;
         gsl_integration_qags(&G,  lower_limit, bMax, err_abs, err_rel, max_steps, w, &result, &er_res);
      }else{
        //itest=1;
        err_rel = 1e-3;
        G.function = &Smatrix_class::gsl_S1b;
        for(int counter = 1; counter <= btbesselroots; counter++)
        {
          upper_limit = boost::math::cyl_bessel_j_zero(2.0, counter);
          //res_compare = result;
          gsl_integration_qags(&G, lower_limit, upper_limit, err_abs, err_rel, max_steps, w, &res_tmp, &er_tmp);
          lower_limit = upper_limit;
          result+= res_tmp;
          er_res+=er_tmp;
          //std::cout <<"counter: " << counter << "\t Res_tmp & Res_prev & err : " << result << "\t" << res_compare
          //          << "\t" << error << std::endl;
       }
       //result = (result+res_compare)/2.;
  }
 }
    //Result:
/*
    std::cout << "  ---> Ib_ = " << result << "\t"
              //<< " Error = " << er_res << "\t"
              << " Err(%) = "<< 100*fabs(er_res/result) << std::endl;
*/
    //if(itest==1){std::cout << "itest= 1: err= " << 100*fabs(er_res/result) << std::endl;}
    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);
    grid_[ir_] = result;
}

double Smatrix_class::S0_b_func0(double b_)
{
  double n0 = MVmodel(r_grid, b_, 0);
  double ne = MVmodel(r_grid, b_, 1);

//ŧest function
  //double integrand = exp(-b*b)*gsl_sf_bessel_J0(b*Delta)*exp(-r_grid);
//real function
  double integrand = b_*gsl_sf_bessel_J0(b_*Delta)*(1. - exp(-n0)*gsl_sf_bessel_I0(ne));

    return integrand*exp(-epsr*r_grid*r_grid)*exp(-epsb*b_*b_);
}

double Smatrix_class::S0_b_func(double b_)
{
  double bt = b_/Delta;
  double n0 = MVmodel(r_grid, bt, 0);
  double ne = MVmodel(r_grid, bt, 1);

//ŧest function
    //double integrand = exp(-bt*bt)*gsl_sf_bessel_J0(bt*Delta)*exp(-r_grid)/Delta;
//real function
    double integrand = bt*gsl_sf_bessel_J0(b_)*(1. - exp(-n0)*gsl_sf_bessel_I0(ne))/Delta;
    return integrand*exp(-epsr*r_grid*r_grid)*exp(-epsb*bt*bt);
}


double Smatrix_class::S1_b_func0(double b_)
{
  double n0 = MVmodel(r_grid, b_, 0);
  double ne = MVmodel(r_grid, b_, 1);

//ŧest function
    //double integrand = exp(-b*b)*gsl_sf_bessel_J0(b*Delta)*exp(-r_grid)/2.;
//real function
  double integrand = b_*gsl_sf_bessel_Jn(2, b_*Delta)*exp(-n0)*gsl_sf_bessel_I1(ne);
  return integrand*exp(-epsr*r_grid*r_grid)*exp(-epsb*b_*b_);
}

double Smatrix_class::S1_b_func(double b_)
{
  double bt = b_/Delta;
  double n0 = MVmodel(r_grid, bt, 0);
  double ne = MVmodel(r_grid, bt, 1);
;
//ŧest function
  //double integrand = exp(-bt*bt)*gsl_sf_bessel_J0(bt*Delta)*exp(-r_grid)/(2.*Delta);
//real function
  double integrand = bt*gsl_sf_bessel_Jn(2, b_)*exp(-n0)*gsl_sf_bessel_I1(ne)/Delta;
  return integrand*exp(-epsr*r_grid*r_grid)*exp(-epsb*bt*bt);
}

double Smatrix_class::gsl_S0b0(double x, void * p)
{
    return ((Smatrix_class *) p)->S0_b_func0(x);
}


double Smatrix_class::gsl_S1b0(double x, void * p)
{
    return ((Smatrix_class *) p)->S1_b_func0(x);
}


double Smatrix_class::gsl_S0b(double x, void * p)
{
    return ((Smatrix_class *) p)->S0_b_func(x);
}

double Smatrix_class::gsl_S1b(double x, void * p)
{
    return ((Smatrix_class *) p)->S1_b_func(x);
}


void Smatrix_class::intS_r(std::vector<double> & grid, int iq_, int iN_)
{
    //std::cout.precision(std::numeric_limits<double>::digits10);

    double result = 0e0;
    double er_res = 0e0;
    double er_tmp = 0e0;
    double lower_limit = .1;
    double upper_limit;
    double err_abs = 0e0;
    double err_rel = 1e-4;

    size_t max_steps = 2e8;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(max_steps);
    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

    gsl_function B;
    B.params = this;

    if(iN_ == 0)
    {
//      std::cout << "r_t integral in S_0 for q_t = " << q_grid << std::flush;

        if(q_grid<2.)
        {
            B.function = &Smatrix_class::gsl_S0r0;
            gsl_integration_qags(&B,  lower_limit, rmaxInt, err_abs, err_rel, max_steps, w, &result, &er_res);
        }else{

            //double res_compare = 0e0;
            double res_tmp;
            B.function = &Smatrix_class::gsl_S0r;

            for(int counter = 1; counter <= rtbesselroots; counter++)
            {
                upper_limit = boost::math::cyl_bessel_j_zero(0.0, counter);

                gsl_integration_qags(&B, lower_limit, upper_limit, err_abs, err_rel, max_steps, w, &res_tmp, &er_tmp);
                //res_compare = result;
                lower_limit = upper_limit;
                result+= res_tmp;
                er_res+=er_tmp;

            }
        //result = (result+res_compare)/2.;
        }
    }
    if(iN_ == 1)
    {
      //std::cout << "r_t integral in S_e for q_t = " << q_grid << std::flush;

        if(q_grid<2.)
        {

            B.function = &Smatrix_class::gsl_S1r0;

            gsl_integration_qags(&B,  lower_limit, rmaxInt,
                            err_abs, err_rel, max_steps, w, &result, &er_res);

        }else{

            //double res_compare = 0e0;
            double res_tmp;
            B.function = &Smatrix_class::gsl_S1r;

            for(int counter = 1; counter <= rtbesselroots; counter++)
            {
              upper_limit = boost::math::cyl_bessel_j_zero(2.0, counter);

                gsl_integration_qags(&B, lower_limit, upper_limit, err_abs, err_rel, max_steps, w, &res_tmp, &er_tmp);

                lower_limit = upper_limit;
                result+= res_tmp;
                er_res+=er_tmp;


                //std::cout <<"counter: " << counter << "\t Res_tmp & Res_prev & err : " << result << "\t" << res_compare
                //          << "\t" << er_res << std::endl;

            }
      }
    }

    //Result:
/*
    std::cout << "  -----> Ir_ = " << result << "\t"
           //<< " Error = " << er_res << "\t"
           << " Err(%) = "<< 100*fabs(er_res/result) << std::endl;
*/
    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);
    grid[iq_]= result;
}


double Smatrix_class::S0_r_func0(double r_)
{
    double B_int = getIb_r(r_, 0);
//test function
    //double integrand = B_int*exp(-r_)*gsl_sf_bessel_J0(q_grid*r_);
//real:
    double integrand = r_*gsl_sf_bessel_J0(q_grid*r_)*B_int/(4.*M_PI*M_PI);
    return integrand;
}

double Smatrix_class::S0_r_func(double r_)
{
    double rt = r_/q_grid;
    double B_int = getIb_r(rt, 0);
//test
    //double integrand = exp(-r_*r_)*gsl_sf_bessel_J0(r_);
//real:
    double integrand = rt*gsl_sf_bessel_J0(r_)*B_int/(4.*M_PI*M_PI*q_grid);
    return integrand;
}


double Smatrix_class::S1_r_func0(double r_)
{
    double B_int = getIb_r(r_, 1);
//test
    //double integrand = B_int*exp(-r_)*gsl_sf_bessel_J0(q_grid*r_);
//real:
    double integrand = r_*gsl_sf_bessel_Jn(2, q_grid*r_)*B_int/(2.*M_PI*M_PI);
    return integrand;
}

double Smatrix_class::S1_r_func(double r_)
{
    double rt = r_/q_grid;
    double B_int = getIb_r(rt, 1);
//test
    //double integrand = exp(-r_*r_)*gsl_sf_bessel_Jn(2, r_);
//real
    double integrand = rt*gsl_sf_bessel_Jn(2, r_)*B_int/(2.*M_PI*M_PI*q_grid);
    return integrand;
}

double Smatrix_class::gsl_S0r0(double x, void * p)
{
    return ((Smatrix_class *) p)->S0_r_func0(x);
}

double Smatrix_class::gsl_S0r(double x, void * p)
{
    return ((Smatrix_class *) p)->S0_r_func(x);
}


double Smatrix_class::gsl_S1r0(double x, void * p)
{
    return ((Smatrix_class *) p)->S1_r_func0(x);
}

double Smatrix_class::gsl_S1r(double x, void * p)
{
    return ((Smatrix_class *) p)->S1_r_func(x);
}
