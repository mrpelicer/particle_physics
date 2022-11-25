/* TO DO:
-Add interpolation error; with error from N_integrals.
-Check what to do for b_ = 1.
*/

//#include "sigma.h"
#include "MV_proton.h"

//=============================== Constructor ===============================
MV_proton_class::MV_proton_class(double fm_, double R_, double mg_, double Qs2_):
    fm(fm_), R(R_), mg(mg_), Qs2(Qs2_),
    bMin(0.0), bMax(std::max(boost::math::cyl_bessel_j_zero(0.0, 150.), boost::math::cyl_bessel_j_zero(2.0, 150.))),
    ib_Max(5000), b_div(1.5),
    d_b((bMax - bMin)/(ib_Max)), Maxroots(200),
    n0_i(ib_Max+1), err_0(ib_Max+1),
    ne_i(ib_Max+1), err_e(ib_Max+1)
{

    std::cout << "------------------- Constructing McLerran Venugopalan model for the proton! -------------------" << std::endl
              << "bmax, db= " << bMax << " " << d_b << std::endl;

    Grid_Ni(n0_i, err_0, 0);
    Grid_Ni(ne_i, err_e, 1);

    double r_ = 0.5; // divide by /(0.197) to print r in fermi with ctes in GeV
    std::ofstream N0;
    std::ofstream Ne;

    N0.open("N0.dat");
    Ne.open("Ne.dat");

    N0 << "#b \t nx_integral \t error \t"
         <<"R = "<< R << "\t mg= "<< mg
         <<"\t Qs2 =" << Qs2 << "\t r = " << 0.197*r_
         << " (r in fermi)" << std::endl;

    for(int ib = 0; ib < ib_Max; ib++)
    {
        N0 <<  0.197*b[ib] << " " << calc_N0(r_, b[ib]) << std::endl;
        Ne <<  0.197*b[ib] << " " << calc_Ne(r_, b[ib]) << std::endl;
    }
    N0.close();
    Ne.close();
};


//======================== Main function: Returns the value of N0 & Ne given b_ and rt_ ===============================
double MV_proton_class::operator()(double r_, double b_, int iN_) const
{
  //iN_ = 0 - N_0
  //iN_ = 1 - N_e
    double res = 0e0;

    if(iN_ == 0) res= calc_N0(r_, b_);
    if(iN_ == 1) res= calc_Ne(r_, b_);

    return res;

}

void MV_proton_class::freeVecs()
{
  std::vector<double>(n0_i).swap(n0_i);
  std::vector<double>(ne_i).swap(ne_i);


}
//=============================== Calculate N0  ===============================
double MV_proton_class::calc_N0(double r_, double b_) const
{
  double n0;

  if(r_ == 0) n0 = 0;
  else n0 = Qs2*log(exp(1.0) + 1./(r_*r_*mg*mg) )*r_*r_*exp(-b_*b_/(4.*R*R))/4.
              + Qs2*R*R*r_*r_*ni_itp(n0_i, b_);
  return n0;
}


//=============================== Calculate Ne ===============================
double MV_proton_class::calc_Ne(double r_, double b_) const
{
    double ne = Qs2*R*R*r_*r_*ni_itp(ne_i, b_);
    return ne;
}

//=============================== Make interpolation  ===============================
double MV_proton_class::ni_itp(std::vector<double> ni_, double b_) const
{
    double res;

    if(b_ > bMax)
    {
        std::cout << "Error!  MV_proton_class::ni_itp calling b > bMax!"  << std::endl;
        exit(EXIT_FAILURE);
    }

    gsl_interp_accel *acc =  gsl_interp_accel_alloc();
    gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_linear, ib_Max+1);

    gsl_interp_init(interpolation, b.data(), ni_.data(), ib_Max+1); //vector.data() transforms std::vector<double> to array!

    res = gsl_interp_eval(interpolation, b.data(), ni_.data(), b_, acc);

    gsl_interp_free(interpolation);
    gsl_interp_accel_free (acc);

    return res;

}

//=============================== Make grid for N0:  ===============================
void MV_proton_class::Grid_Ni(std::vector<double> & ni_, std::vector<double> & err_, int i_)
{
    double b_tmp;
    std::vector <double> n0tmp(ib_Max+1);
    std::vector <double> errtmp(ib_Max+1);
  std::cout << "Griding N_" << i_<< "in b_t." << std::endl;
    //make grid for ni_[ib]:
    for(int ib = 0; ib<= ib_Max; ib++)
    {
        b_tmp = bMin + ib*d_b;
        b.push_back(b_tmp);

        if(i_ == 0) intN0(n0tmp, errtmp, b_tmp, ib);
        if(i_ == 1) intNe(n0tmp, errtmp, b_tmp, ib);

        ni_[ib] = n0tmp[ib];
        err_[ib] = errtmp[ib];

    }
    std::cout << "Grid done!" << std::endl;
}




//============================= Calculus of the double integral in N_0  ===============================
void MV_proton_class::intN0(std::vector<double> & ni_, std::vector<double>  & err_, double b_, int ib_) const
{
    double result = 0e0;
    double er_tmp = 0e0;
    double er_res = 0e0;
    double lower_limit = 0.0;
    double upper_limit;
    double err_abs = 0e0;
    double err_rel = 1e-5;

    //define parameters in structure:
    N0_params N0p;
    N0p.Rp = R;
    N0p.mgp = mg;
    N0p.bp = b_;

    //Define gsl variables:
    size_t Max_steps = 2e8;
    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(Max_steps);

    gsl_function G;
    G.params = reinterpret_cast<void *>(&N0p);

   std::cout.precision(std::numeric_limits<double>::digits10);
   //std::cout << "Calculating N_0 integral for b = "<< b_ << std::flush;

//Do integration:

   if(b_ <b_div)
   {
      G.function = &N0;
      gsl_integration_qags(&G,  lower_limit, bMax, err_abs, err_rel, Max_steps, w, &result, &er_res);

    } else
    {

       G.function = &N0_2;
        double res_tmp;
        double res_compare = 0.;

       for(int counter = 1; counter <= Maxroots; counter ++)
       {
               upper_limit = boost::math::cyl_bessel_j_zero(0.0, counter);
               res_compare = result;

               gsl_integration_qags(&G, lower_limit, upper_limit, err_abs, err_rel, Max_steps, w, &res_tmp, &er_tmp);

               lower_limit = upper_limit;

               result+= res_tmp;
               er_res+=er_tmp;

           }
       result = (result+res_compare)/2.;

    }

   //Result:
/*
   std::cout << "  -----> N0_int = " << result << "\t"
             //<< " Error = " << error << "\t"
             << " Err(%) = "<< 100*fabs(er_res/result) << std::endl;
*/


    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);

        ni_[ib_] = result;
        err_[ib_] = er_res;
}

//============================= Func. to be integrated in N_0  ===============================
double MV_proton_class::N0(double x, void * p)
{
    N0_params &N0p= *reinterpret_cast<N0_params *>(p);


    double mgluon = N0p.mgp;
    double R_= N0p.Rp;
    double b_= N0p.bp;
    double mgluon2 = mgluon*mgluon;

    double x2 = x*x;
    double srt = sqrt(x2 + 4.*mgluon2);
    double num = x*srt - 2.*(x2 + 2.*mgluon2)*atanh(x/srt);
    double integrand = x*exp(-x*x*R_*R_)*gsl_sf_bessel_J0(x*b_)*num/(2.*x*srt);

    return integrand;
}

double MV_proton_class::N0_2(double x, void * p)
{
    N0_params &N0p= *reinterpret_cast<N0_params *>(p);

    double mgluon = N0p.mgp;
    double R_= N0p.Rp;
    double b_= N0p.bp;
    double b2 = b_*b_;
    double mgluon2 = mgluon*mgluon;

    double x2 = x*x;
    double srt = sqrt(x2 + 4.*mgluon2*b2);
    double num = x*srt - 2.*(x2 + 2.*mgluon2*b2)*atanh(x/srt);
    double integrand = x*exp(-x2*R_*R_/b2)*gsl_sf_bessel_J0(x)*num/(2.*x*srt*b2);

    return integrand;
}


//============================= Calculus of the double integral in N_e  ===============================
void MV_proton_class::intNe(std::vector<double> & ni_, std::vector<double>  & err_, double b_, int ib_) const
{
    double result = 0e0;
    double er_tmp = 0e0;
    double er_res = 0e0;
    double lower_limit = 0.;
    double upper_limit;
    double err_abs = 0e0;
    double err_rel = 1e-5;

    //define parameters in struct:
    N0_params Nep;
    Nep.Rp = R;
    Nep.mgp = mg;
    Nep.bp = b_;


    size_t Max_steps = 1e8;

    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(Max_steps);

    gsl_function G;
    G.params = reinterpret_cast<void *>(&Nep);

   //std::cout.precision(std::numeric_limits<double>::digits10);
   //std::cout << "Calculating N_e integral for b = "<< b_ << std::flush;


    if(b_ < b_div)
    {
      G.function = &Ne;
      gsl_integration_qags(&G,  lower_limit, bMax,  err_abs, err_rel, Max_steps, w, &result, &er_res);

    }else
    {

        G.function = &Ne_2;
        double res_tmp;
        double res_compare = 0.;

        for(int counter = 1; counter <= Maxroots; counter ++)
        {
           upper_limit = boost::math::cyl_bessel_j_zero(2.0, counter);

           //std::cout << "Limits: " << lower_limit << "\t" << upper_limit << std::endl;
           res_compare = result;

           gsl_integration_qags(&G, lower_limit, upper_limit, err_abs, err_rel, Max_steps, w, &res_tmp, &er_tmp);

           lower_limit = upper_limit;

           result+= res_tmp;
           er_res+=er_tmp;
        }
        result = (result+res_compare)/2.;

    }

    //Result:
/*
        std::cout << "  -----> Ne_int = " << result << "\t"
        //<< " Error = " << er_res<< "\t"
        << " Err(%) = "<< 100*fabs(er_res/result) << std::endl;
*/

        gsl_set_error_handler(old_handler);
        gsl_integration_workspace_free(w);


        ni_[ib_] = result;
        err_[ib_] = er_res;
}

//============================= Func. to be integrated in N_0  ===============================
double MV_proton_class::Ne(double x, void * p)
{
    N0_params &Nep= *reinterpret_cast<N0_params *>(p);

    double mgluon = Nep.mgp;
    double R_= Nep.Rp;
    double b_= Nep.bp;
    double mgluon2 = mgluon*mgluon;

    double x2 = x*x;
    double srt = sqrt(x2 + 4.*mgluon2);
    double num = x*srt - 4.*mgluon2*atanh(x/srt);
    double integrand = x*exp(-x*x*R_*R_)*gsl_sf_bessel_Jn(2, x*b_)*num/(2.*x*srt);

    return integrand;
}


double MV_proton_class::Ne_2(double x, void * p)
{
    N0_params &Nep= *reinterpret_cast<N0_params *>(p);

    double mgluon = Nep.mgp;
    double R_= Nep.Rp;
    double b_= Nep.bp;
    double b2 = b_*b_;
    double mgluon2 = mgluon*mgluon;

    double x2 = x*x;
    double srt = sqrt(x2 + 4.*mgluon2*b2);
    double num = x*srt - 4.*mgluon2*b2*atanh(x/srt);
    double integrand = x*exp(-x2*R_*R_/b2)*gsl_sf_bessel_Jn(2, x)*num/(2.*x*b2*srt);

    return integrand;
}
