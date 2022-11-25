#include "MV_nucleus.h"

//=============================== Constructor ===============================
MV_nucleus_class::MV_nucleus_class(double fm_, double R_, double mg_, double Qs2_, double A_, double RA_):
    fm(fm_), R(R_), mg(mg_), Qs2(Qs2_), A(A_), RA(RA_),
    delta(0.54*fm), T_norm(0.0000058627112555), //Pb: A = 208; T_norm= (fm)0.0000058627112555   -- (1/mg)0.0000119817
                                                //-- for other nucleus, use normalization in Eq.(46) of Phys. Rev. D 95, 094003 (2017)
    bMin(0.0), bMax( std::max(boost::math::cyl_bessel_j_zero(0.0, 200.), boost::math::cyl_bessel_j_zero(2.0, 200.)) ), //adjust to needs;
    xMin(0.01), xMax(bMax), x_div(55.),  //x_div changes the integration method to keep errors under control.
    ix_Max(1+50*(int)xMax), ib_Max(10000), //50 000 in test! any differents on sunday?
    d_x((xMax - xMin)/(ix_Max)), d_b((bMax - bMin)/(ib_Max)), Maxroots(1000),
    b(ib_Max+1), x(ix_Max+1),
    n0_i(ix_Max+1), ne_i(ix_Max+1),  TA(ib_Max+1), N0A_1(ib_Max+1), N0A_2(ib_Max+1), NeA_1(ib_Max+1)
{


      std::cout << "------------------- Constructing single digluon amplitude w/ McLerran Venugopalan model for the nucleus! -------------------" << std::endl
                << "Bmax, ibmax, db = " << bMax << " " << ib_Max << " " << d_b << std::endl
                << "===================================================" << std::endl;
    //[x = b , b = B in phys.rev.D.97 - 094003 (2017)]

    //Fill values to be grided:
    for(int ix = 0; ix<= ix_Max; ix++) x[ix] =  xMin + ix*d_x;
    for(int ib = 0; ib<= ib_Max; ib++) b[ib] =  bMin + ib*d_b;

//Grid thickness in b
    gridThickness(TA);
/*
    std::ofstream outTA;
    outTA.open("thickness.dat");

    outTA << "#(b, TA, TA', TA'', TA''-TA'/b )" << std::endl;
    for(int ib= 0; ib< ib_Max; ib++)
    {
      outTA << b[ib] << " " << itpThickness_NAfuncs(TA, b[ib]) << " "
        << derivativeThickness(TA, b[ib], 1) << " " << derivativeThickness(TA, b[ib], 2) << " "
        << derivativeThickness(TA, b[ib], 2) - derivativeThickness(TA, b[ib], 1)/b[ib] << " "
        << derivativeThickness(TA, b[ib], 2) + derivativeThickness(TA, b[ib], 1)/b[ib] << std::endl;
    }
    outTA.close();
*/
/* ONLY IF USING EQ. 53 OF YANCU'S PAPER!
//Grid n0 and ne in x:
    grid_ni(n0_i, 0);
    grid_ni(ne_i, 1);

//Integrals in x on the first and second terms of N0, and on Ne:
    gridNA(N0A_1, 1);
    gridNA(N0A_2, 2);
    gridNA(NeA_1, 3);
*/

/*
std::ofstream outN;
outN.open("ntmp.dat");

outN << "#(b, N0, Ne)" << std::endl;
for(int ib=0; ib< ib_Max; ib++)
{
  outN << b[ib] << " " << calcN0(0.5, b[ib]) << " "
       << calcNe(0.5, b[ib])  << std::endl;
}

outN.close();
*/

};


//======================== Main function: Returns the value of N0 & Ne given b_ and rt_ ===============================
double MV_nucleus_class::operator()(double r_, double b_, int iN_) const
{
  //iN_ = 0 - N_0
  //iN_ = 1 - N_e
    double res = 0e0;

    if(iN_ == 0) res= A*calcN0(r_, b_);
    if(iN_ == 1) res= A*calcNe(r_, b_);

    return res;

}

//=============================== Calculate N0  ===============================
double MV_nucleus_class::calcN0(double r_, double b_) const
{
  double n0 = 0e0;
if(r_ == 0) n0 = 0;
else{
//For Eq.(53) of Iancu's model.

   //n0 = 2*M_PI*Qs2*log(exp(1.0) + 1./(r_*r_*mg*mg) )*r_*r_*itpThickness_NAfuncs(N0A_1, b_)/4.
  //            + 2*M_PI*Qs2*R*R*r_*r_*itpThickness_NAfuncs(N0A_2, b_);

//For Eq.(54) of Iancu's model.
  double plus  = derivativeThickness(TA, b_, 2) + derivativeThickness(TA, b_, 1)/b_;
  n0 = M_PI*R*R*Qs2*r_*r_*log(exp(1.0) + 1./(r_*r_*mg*mg) )*(itpThickness_NAfuncs(TA, b_ ) + R*R*plus )
            + M_PI*R*R*Qs2*r_*r_*plus/(3.*mg*mg);
}
  return n0;
}


//=============================== Calculate Ne ===============================
double MV_nucleus_class::calcNe(double r_, double b_) const
{
//For Eq.(53) of Iancu's model.
  //double ne = M_PI*Qs2*R*R*r_*r_*itpThickness_NAfuncs(NeA_1, b_)/4.;

//For Eq.(54) of Iancu's model.
    double Minus  = derivativeThickness(TA, b_, 2) - derivativeThickness(TA, b_, 1)/b_;
    double ne = M_PI*R*R*Qs2*r_*r_*Minus/(6.*mg*mg);

    return ne;
}


//=============================== Make grid for N0:  ===============================

void MV_nucleus_class::gridNA(std::vector<double> & NA_, int iN_)
{

  if(iN_ == 3) //Elliptic function N_theta
  {
    std::vector <double> Netmp(1);
    ib_grid = 0;
    intNA(Netmp, 3);
    std::cout << "Netmp(0) = " << Netmp[0] << std::endl;

    for(ib_grid = 0; ib_grid<= ib_Max; ib_grid++)
      {NA_[ib_grid] = Netmp[0]*(derivativeThickness(TA, b[ib_grid], 2) - derivativeThickness(TA, b[ib_grid], 1)/b[ib_grid]);}

  }else if(iN_ ==1 || iN_ == 2) //Main contribution (N_0) can be dividided in two integrations:
  {
    for(ib_grid = 0; ib_grid<= ib_Max; ib_grid++) intNA(NA_, iN_);
    std::cout << "Grid in NA" <<iN_<<" done!" << std::endl;
  }else
  {
    std::cout << "Calling indice out of range in gridNA" << std::endl;
    exit(EXIT_FAILURE);
  }

}

void MV_nucleus_class::intNA(std::vector<double> &NA_, int iA_)
{
    double result = 0e0;
    double error;
    //double lower_limit = 0.0;
    double err_abs = 0e0;
    double err_rel = 1e-5;

    //Define gsl variables:
    size_t Max_steps = 2e6;
    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(Max_steps);

    gsl_function G;
    if(iA_ == 1) G.function = &MV_nucleus_class::gsl_N0A_func1;
    if(iA_ == 2) G.function = &MV_nucleus_class::gsl_N0A_func2;
    if(iA_ == 3) G.function = &MV_nucleus_class::gsl_NeA_func;
    G.params = this;

   std::cout.precision(std::numeric_limits<double>::digits10);

   gsl_integration_qags(&G,  xMin, xMax, err_abs, err_rel, Max_steps, w, &result, &error);

   NA_[ib_grid] =result;

  std::cout << "  -----> NA_" <<iA_<<" ( b= "<< b[ib_grid]<< " )= " << result << "\t"
                << " Err(%) = "<< 100*fabs(error/result) << std::endl;

    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);
}

double MV_nucleus_class::gsl_N0A_func1(double x, void * p)
{
    return ((MV_nucleus_class *) p)->N0A_func1(x);
}
double MV_nucleus_class::gsl_N0A_func2(double x, void * p)
{
    return ((MV_nucleus_class *) p)->N0A_func2(x);
}

double MV_nucleus_class::gsl_NeA_func(double x, void * p)
{
    return ((MV_nucleus_class *) p)->NeA_func(x);
}

double MV_nucleus_class::N0A_func1(double x_)
{
  double x2 = x_*x_;
  double thick = itpThickness_NAfuncs(TA, b[ib_grid])
                + x2*(derivativeThickness(TA, b[ib_grid], 2) + derivativeThickness(TA, b[ib_grid], 1)/b[ib_grid])/4.;

  double integrand = x_*exp(-x2/(4.*R*R))*thick;
    return integrand;
}

double MV_nucleus_class::N0A_func2(double x_)
{
  double x2 = x_*x_;
  double thick = itpThickness_NAfuncs(TA, b[ib_grid])
                + x2*(derivativeThickness(TA, b[ib_grid], 2) + derivativeThickness(TA, b[ib_grid], 1)/b[ib_grid])/4.;

  double integrand = x_*itp_ni(n0_i, x_)*thick;
  return integrand;
}

double MV_nucleus_class::NeA_func(double x_)
{
    double integrand = pow(x_, 3)*itp_ni(ne_i, x_);
    return integrand;
}


void MV_nucleus_class::gridThickness(std::vector<double> & T_)
{
  std::vector<double> ta_tmp(ib_Max+1);
  for(int ib = 0; ib<= ib_Max; ib++)
  {
    intThickness(ta_tmp, b[ib], ib);
    T_[ib] = ta_tmp[ib];
  }
}

void MV_nucleus_class::intThickness(std::vector<double> & T_, double b_, int ib_) const
{
    double result = 0e0;
    double error;
    double lower_limit = 0.0;
    double err_abs = 0e0;
    double err_rel = 1e-5;

    T_params Tp;
    Tp.deltap = delta;
    Tp.Ap= A;
    Tp.RAp = RA;
    Tp.NAp = T_norm;
    Tp.bp  = b_;

    //Define gsl variables:
    size_t Max_steps = 2e6;
    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(Max_steps);

    gsl_function G;
    G.function = &TA_function;
    G.params = reinterpret_cast<void *>(&Tp);

   std::cout.precision(std::numeric_limits<double>::digits10);

   gsl_integration_qagiu(&G,  lower_limit, err_abs, err_rel, Max_steps, w, &result, &error);

   T_[ib_] =result;

/*
      std::cout << "  -----> T_A ( b= "<< b_ << " )= " << result << "\t"
                //<< " Error = " << error << "\t"
                << " Err(%) = "<< 100*fabs(error/result) << std::endl;
*/
    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);
}

double MV_nucleus_class::TA_function(double z, void * p)
{
    T_params &Tp= *reinterpret_cast<T_params *>(p);
    double b2 = pow(Tp.bp, 2);
    double z2 = z*z;
    double arg = (sqrt(b2 + z2) - Tp.RAp)/Tp.deltap;
    double integrand = 2.*Tp.NAp/( exp(arg) + 1.);

    return integrand;
}

double MV_nucleus_class::itpThickness_NAfuncs(std::vector<double> vec_, double b_) const
{
    double res;
    if(b_ > bMax || b_ < bMin)
    {
        std::cout << "Error!  MV_nucleus_class::itpThickness_NAfuncs calling b out of range! b= "  << b_ << std::endl;
        exit(EXIT_FAILURE);
    }

    gsl_interp_accel *acc =  gsl_interp_accel_alloc();
    gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_linear, ib_Max+1);

    gsl_interp_init(interpolation, b.data(), vec_.data(), ib_Max+1); //vector.data() transforms std::vector<double> in an array!

    res = gsl_interp_eval(interpolation, b.data(), vec_.data(), b_, acc);

    gsl_interp_free(interpolation);
    gsl_interp_accel_free (acc);

    return res;

}

double MV_nucleus_class::derivativeThickness(std::vector<double> TA_, double b_, int order_) const
{
  double res = 0e0;
  double h = 0e0;

  if(order_ == 1)
  {
    h = 1e-3;
    if((b_+h) > bMax ) res = (itpThickness_NAfuncs(TA_, b_) - itpThickness_NAfuncs(TA_, b_-h))/h;
    else if((b_-h) < bMin ) res = (itpThickness_NAfuncs(TA_, b_+h) - itpThickness_NAfuncs(TA_, b_))/h;
    else res = (itpThickness_NAfuncs(TA_, b_+h) - itpThickness_NAfuncs(TA_, b_-h))/(2.*h);
  } else if(order_ == 2)
  {
    h =5e-1;
    res = 0e0;
    if( (b_+h) > bMax )     res = (itpThickness_NAfuncs(TA_, b_)      - 2.*itpThickness_NAfuncs(TA_, b_-h) + itpThickness_NAfuncs(TA_, b_-2.*h))/(h*h);
    else if((b_-h) <= bMin)  res = (itpThickness_NAfuncs(TA_, b_+2.*h) - 2.*itpThickness_NAfuncs(TA_, b_+h) + itpThickness_NAfuncs(TA_, b_))/(h*h);
    else                    res = (itpThickness_NAfuncs(TA_, b_+h)    - 2.*itpThickness_NAfuncs(TA_, b_) + itpThickness_NAfuncs(TA_, b_-h))/(h*h);
  }else
  {
    std::cout << "Error in order_ of derivativeThickness!" << std::endl;
    exit(EXIT_FAILURE);
  }
  return res;
}

void MV_nucleus_class::grid_ni(std::vector<double> & ni_, int i_)
{
    //double x_tmp;
    std::vector <double> n0tmp(ix_Max+1);

    //make grid for ni_[ib]:
    for(int ix = 0; ix<= ix_Max; ix++)
    {
        if(i_ == 0) int_n0(n0tmp, x[ix], ix);
        if(i_ == 1) int_ne(n0tmp, x[ix], ix);
        ni_[ix] = n0tmp[ix];
    }
    std::cout << "grid_ni done!" << std::endl;
}

void MV_nucleus_class::int_n0(std::vector<double> & ni_, double x_, int ix_) const
{
    double result = 0e0;
    double er_tmp = 0e0;
    double er_res = 0e0;
    double lower_limit = 0.0;
    double upper_limit;
    double err_abs = 0e0;
    double err_rel = 1e-5;

    //define parameters in structure:
    n0_params n0p;
    n0p.Rp = R;
    n0p.mgp = mg;
    n0p.xp = x_;

    //Define gsl variables:
    size_t Max_steps = 2e6;
    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(Max_steps);

    gsl_function G;
    G.params = reinterpret_cast<void *>(&n0p);

   std::cout.precision(std::numeric_limits<double>::digits10);
   std::cout << "Calculating n0 integral for x = "<< x_ << std::flush;

//Do integration:

   if(x_ <x_div)
   {
       G.function = &n0_func;
       gsl_integration_qagiu(&G,  lower_limit, err_abs, err_rel, Max_steps, w, &result, &er_res);
    } else
    {

       G.function = &n0_func2;
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
   std::cout << "  -----> int_n0 = " << result << "\t"
             //<< " Error = " << error << "\t"
             << " Err(%) = "<< 100*fabs(er_res/result) << std::endl;



    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);

        ni_[ix_] = result;
}

//============================= Func. to be integrated in N_0  ===============================
double MV_nucleus_class::n0_func(double d_, void * p)
{
    n0_params &n0p= *reinterpret_cast<n0_params *>(p);


    double mgluon = n0p.mgp;
    double R_= n0p.Rp;
    double x_= n0p.xp;
    double mgluon2 = mgluon*mgluon;

    double d2 = d_*d_;
    double srt = sqrt(d2 + 4.*mgluon2);
    double num = d_*srt - 2.*(d2 + 2.*mgluon2)*atanh(d_/srt);
    double integrand = d_*exp(-d2*R_*R_)*gsl_sf_bessel_J0(d_*x_)*num/(2.*d_*srt);

    return integrand;
}

double MV_nucleus_class::n0_func2(double d_, void * p)
{
    n0_params &n0p= *reinterpret_cast<n0_params *>(p);

    double mgluon = n0p.mgp;
    double R_= n0p.Rp;
    double x_= n0p.xp;
    double x2 = x_*x_;
    double mgluon2 = mgluon*mgluon;

    double d2 = d_*d_;
    double srt = sqrt(d2 + 4.*mgluon2*x2);
    double num = d_*srt - 2.*(d2 + 2.*mgluon2*x2)*atanh(d_/srt);
    double integrand = d_*exp(-d2*R_*R_/x2)*gsl_sf_bessel_J0(d_)*num/(2.*d_*srt*x2);

    return integrand;
}


//============================= Calculus of the double integral in N_e  ===============================
void MV_nucleus_class::int_ne(std::vector<double> & ni_,  double x_, int ix_) const
{
    double result = 0e0;
    double er_tmp = 0e0;
    double er_res = 0e0;
    double lower_limit = 0.;
    double upper_limit;
    double err_abs = 0e0;
    double err_rel = 1e-4;

    //define parameters in struct:
    n0_params nep;
    nep.Rp = R;
    nep.mgp = mg;
    nep.xp = x_;


    size_t Max_steps = 100000;

    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(Max_steps);

    gsl_function G;
    G.params = reinterpret_cast<void *>(&nep);

   //std::cout.precision(std::numeric_limits<double>::digits10);
   std::cout << "Calculating int_ne integral for x = "<< x_ << std::flush;


    if(x_ < x_div)
    {
        G.function = &ne_func;

        gsl_integration_qagiu(&G,  lower_limit, err_abs, err_rel, Max_steps, w, &result, &er_res);

    }else
    {

        G.function = &ne_func2;
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
        std::cout << "  -----> int_ne = " << result << "\t"
        //<< " Error = " << er_res<< "\t"
        << " Err(%) = "<< 100*fabs(er_res/result) << std::endl;


        gsl_set_error_handler(old_handler);
        gsl_integration_workspace_free(w);

        ni_[ix_] = result;
}

//============================= Func. to be integrated in N_0  ===============================
double MV_nucleus_class::ne_func(double d_, void * p)
{
    n0_params &nep= *reinterpret_cast<n0_params *>(p);

    double mgluon = nep.mgp;
    double R_= nep.Rp;
    double x_= nep.xp;
    double mgluon2 = mgluon*mgluon;

    double d2 = d_*d_;
    double srt = sqrt(d2 + 4.*mgluon2);
    double num = d_*srt - 4.*mgluon2*atanh(d_/srt);
    double integrand = d_*exp(-d_*d_*R_*R_)*gsl_sf_bessel_Jn(2, d_*x_)*num/(2.*d_*srt);

    return integrand;
}


double MV_nucleus_class::ne_func2(double d_, void * p)
{
    n0_params &nep= *reinterpret_cast<n0_params *>(p);

    double mgluon = nep.mgp;
    double R_= nep.Rp;
    double x_= nep.xp;
    double x2 = x_*x_;
    double mgluon2 = mgluon*mgluon;

    double d2 = d_*d_;
    double srt = sqrt(d2 + 4.*mgluon2*x2);
    double num = d_*srt - 4.*mgluon2*x2*atanh(d_/srt);
    double integrand = d_*exp(-d2*R_*R_/x2)*gsl_sf_bessel_Jn(2, d_)*num/(2.*d_*x2*srt);

    return integrand;
}


//=============================== Make interpolation  ===============================
double MV_nucleus_class::itp_ni(std::vector<double> ni_, double x_) const
{
    double res;

    if(x_ > xMax)
    {
        std::cout << "Error!  MV_nucleus_class::ni_itp calling x > xMax!"  << std::endl;
        exit(EXIT_FAILURE);
    }

    gsl_interp_accel *acc =  gsl_interp_accel_alloc();
    gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_linear, ix_Max+1);

    gsl_interp_init(interpolation, x.data(), ni_.data(), ix_Max+1); //vector.data() transforms std::vector<double> to array!

    res = gsl_interp_eval(interpolation, x.data(), ni_.data(), x_, acc);

    gsl_interp_free(interpolation);
    gsl_interp_accel_free (acc);

    return res;

}
