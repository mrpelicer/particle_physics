#include "sigma.h"

//Constructor:
sigma_class::sigma_class(Smatrix_class & Smatrix_, std::vector<double> Pt_, int iDtMax_, double R_proj_,
                          double R_target_, double mf_, double sqrt_s_): // in GeV
    S_matrix(Smatrix_), fm(S_matrix.fm),
    qMin(S_matrix.qMin), qMax(S_matrix.qMax),
    iDt(0), iDtMax(iDtMax_),
    PtVec(Pt_), DtVec(iDtMax),
    PtMin(PtVec[0]), PtMax(PtVec[PtVec.size()-1]), iPtMax(PtVec.size()),
    Agrid(iPtMax, std::vector<double>(iDtMax)), Bgrid(iPtMax, std::vector<double>(iDtMax)),
    Cgrid(iPtMax, std::vector<double>(iDtMax) ), Dgrid(iPtMax, std::vector<double>(iDtMax) ),
    Z(82), alphaem(1e0/137e0),
    R_projectile(R_proj_), R_target(R_target_), //in Gev^{-1} // 1 fermi = 1/(200e6 eV) = 1/(0,2 GeV)
    Ecm(sqrt_s_),
    Mp(0.938), gamma(Ecm/(2.*Mp)), Nc(3),
    charge(7),  mf(mf_),
          //g, d,u, s, c, b, t -  PDG
    mf2(mf*mf), iFlavour(2)
{
  charge={0.0, -1./3, 2./3, -1./3, 2./3,-1./3, 2./3};
}

sigma_class::~sigma_class()
{}

double sigma_class::getFuncs(double Pt_, double Dt_, int iFunc)
{
  double res = 0.0;
  if(iFunc == 0) res = itp_funcs(Agrid, Pt_, Dt_);
  if(iFunc == 1) res = itp_funcs(Bgrid, Pt_, Dt_);
  if(iFunc == 2) res = itp_funcs(Cgrid, Pt_, Dt_);
  if(iFunc == 3) res = itp_funcs(Dgrid, Pt_, Dt_);
  return res;
}


double sigma_class::operator()(double y1_, double y2_, double Pt_, double Dt_, double dPhi_, int iSigma_)
{
  return 0;
}


double sigma_class::getSigma(double y1_, double y2_, double Pt_, double Dt_, double dPhi_, int iSigma_)
{
  //iSigma_ == 0 -> full differential (diff) cross section
  //iSigma_ == 1 -> c.s. w/ D=0
  //iSigma_ == 2 -> c.s. w/ C=D=0
  //iSigma_ == 3 -> c.s. w/ B=C=D=0
  //iSigma_ == 4 -> c.s. w/ A=B=0

  Pt2 = Pt_*Pt_;
  Dt2 = Dt_*Dt_;

  double k1t = sqrt(Pt2 + Pt_*Dt_*cos(dPhi_) + Dt2/4.);
  double k2t = sqrt(Pt2 - Pt_*Dt_*cos(dPhi_) + Dt2/4.);

  double m1t = sqrt(k1t*k1t + mf2);
  double m2t = sqrt(k2t*k2t + mf2);
  double z = m1t*exp(y1_)/(m1t*exp(y1_) + m2t*exp(y2_));
  double omega = m1t*exp(y1_) + m2t*exp(y2_)/2.;
  double photon_flux = get_photonflux(omega);
  double zm = 1.-z;
  double term1 = photon_flux*omega*2.*pow((2*M_PI), 2)*Nc*alphaem*charge[iFlavour]*z*zm/Pt2;

  double A = itp_funcs(Agrid, Pt_, Dt_);
  double B = itp_funcs(Bgrid, Pt_, Dt_);
  if(iSigma_ == 3) {B=0;}
  if(iSigma_ == 4) {A=0; B=0;}
  double transverse_massless  = 0;
  transverse_massless = (z*z + zm*zm)*(A*A + 2.*A*B*cos(2.*dPhi_) + B*B*cos(2.*dPhi_)*cos(2.*dPhi_) );

  double transverse_massive = 0.;
  if(mf != 0)
  {
  double C = itp_funcs(Cgrid, Pt_, Dt_);
  double D = itp_funcs(Dgrid, Pt_, Dt_);
  if(iSigma_ == 1) {D=0;}
  if(iSigma_ == 2 || iSigma_ == 3) {C=0; D=0;}
  transverse_massive = (mf2/Pt2)*(C*C + 2.*C*D*cos(2.*dPhi_) +D*D*cos(2.*dPhi_)*cos(2.*dPhi_) );
}


	double res = term1*(transverse_massless + transverse_massive);

	return res;
}

double sigma_class::avgSigma(double y1_, double y2_, double Pt_, double Dt_, int itype, int iRes_)
{

  	int n = 1000; 	//Steps of integration -- LAST CHANGEd!!!!!
  	double min = 0.0, max = 2.*M_PI; //Limits of integration.

  	double h = (max - min)/(double)n;
  	double sum=0.0, result = 0.0;
    double phi = min+h;

    if(itype==0)
    {
      for(int i = 1; i<n; i++)
	    {
  		  sum+=getSigma(y1_, y2_, Pt_, Dt_,phi, iRes_);
        phi+= h;
      }
      result=  2.*M_PI*h*( 2.*sum + getSigma(y1_, y2_, Pt_, Dt_,min, iRes_)
                                  + getSigma(y1_, y2_, Pt_, Dt_,max, iRes_) )/2.0;
  	}else if(itype==1)
    {
      for(int i = 1; i<n; i++)
	    {
        sum+=cos(2.*phi)*getSigma(y1_, y2_, Pt_, Dt_, phi,  iRes_);
        phi+=h;
      }
      result=  2.*M_PI*h*(2.*sum  + cos(2.*min)*getSigma(y1_, y2_, Pt_, Dt_, min, iRes_)
                                  + cos(2.*max)*getSigma(y1_, y2_, Pt_, Dt_, max, iRes_) )/2.0;

    }else{std::cout << " Invalid itype! " << std::endl;}


  	return result;

}


double sigma_class::getSigma_alternative(double y1_, double y2_, double k1_, double k2_, double Dt_, double Pt_, double dPhi_, double iF_, int iSigma_)
{
  //iSigma_ == 0 -> full differential (diff) cross section
  //iSigma_ == 1 -> c.s. w/ D=0
  //iSigma_ == 2 -> c.s. w/ C=D=0
  //iSigma_ == 3 -> c.s. w/ B=C=D=0
  //iSigma_ == 4 -> c.s. w/ A=B=0

  Pt2 = Pt_*Pt_;
  Dt2 = Dt_*Dt_;

  double k1t = k1_;
  double k2t = k2_;
  double m1t = sqrt(k1t*k1t + mf2);
  double m2t = sqrt(k2t*k2t + mf2);
  double z = m1t*exp(y1_)/(m1t*exp(y1_) + m2t*exp(y2_));
  double omega = m1t*exp(y1_) + m2t*exp(y2_)/2.;
  double photon_flux = get_photonflux(omega);
  double zm = 1.-z;
  double term1 = photon_flux*omega*2.*pow((2*M_PI), 2)*Nc*alphaem*pow(charge[iF_], 2)*z*zm/Pt2;

  double A = itp_funcs(Agrid, Pt_, Dt_);
  double B = itp_funcs(Bgrid, Pt_, Dt_);
  if(iSigma_ == 3) {B=0;}
  if(iSigma_ == 4) {A=0; B=0;}
  double transverse_massless  = 0;
  transverse_massless = (z*z + zm*zm)*(A*A + 2.*A*B*cos(2.*dPhi_) + B*B*cos(2.*dPhi_)*cos(2.*dPhi_) );

  double transverse_massive = 0.;
  if(mf != 0)
  {
  double C = itp_funcs(Cgrid, Pt_, Dt_);
  double D = itp_funcs(Dgrid, Pt_, Dt_);
  if(iSigma_ == 1) {D=0;}
  if(iSigma_ == 2 || iSigma_ == 3) {C=0; D=0;}
  transverse_massive = (mf2/Pt2)*(C*C + 2.*C*D*cos(2.*dPhi_) +D*D*cos(2.*dPhi_)*cos(2.*dPhi_) );
}


	double res = term1*(transverse_massless + transverse_massive);

	return res;
}

double sigma_class::get_photonflux(double omega_)
{
  double xi = omega_*(R_projectile + R_target)/gamma;
	double bessel0 = gsl_sf_bessel_K0(xi);
	double bessel1 = gsl_sf_bessel_K1(xi);
	double sum = xi*bessel0*bessel1 - xi*xi*( bessel1*bessel1 - bessel0*bessel0)/2.;

	return 2*Z*Z*alphaem*sum/(M_PI*omega_);
}

void sigma_class::gridFuncs(double delta_)
{
  if(iDt>iDtMax) {std::cout << "In sigma_class::gridFuncs - iDt> iDtMax; iDt= " << iDt << std::endl; }
  DtVec[iDt] = delta_;
  Dt = delta_;
  std::cout << "Griding cross section level functions for Deltat = " << Dt << ", iDt= " << iDt << std::endl;

  for(int iPt = 0; iPt< iPtMax; iPt++)
  {
    Pt = PtVec[iPt];
    Pt2 = Pt*Pt;
    Agrid[iPt][iDt] = intA();
    Bgrid[iPt][iDt] = intB();
    if(mf!=0)
    {
    Cgrid[iPt][iDt] = intC();
    Dgrid[iPt][iDt] = intD();
    }
  }
  iDt++;
  std::cout << "=================================================== " << std::endl;

}

double sigma_class::itp_funcs(std::vector<std::vector<double> > func_, double Pt_, double Dt_ )
{//Makes the interpolation of a vector grid_ in the variable r_.
    double res = 0.0;

    if(Pt_ > PtMax || Dt_ > DtVec[DtVec.size()-1] )
    {
        std::cout << "Error! Argument of sigma_class::itp_funcs out of range! pt, dt: "
                  << Pt_ << " " <<  Dt_ << std::endl;
        exit(EXIT_FAILURE);
    }

    size_t ipt=0, idt=0;
    do{ipt++;} while(Pt_ > PtVec[ipt]);
    do{idt++;} while(Dt_ > DtVec[idt]);

    if(ipt == PtVec.size() && idt == DtVec.size() )
    {
      res = func_[ipt-1][idt-1];
    }else if(ipt == PtVec.size())
    {
      ipt--;
      res = (func_[ipt][idt-1]*(DtVec[idt] - Dt_) + func_[ipt][idt]*(Dt_ - DtVec[idt-1]))/(DtVec[idt] - DtVec[idt-1]);
    }else if(idt == DtVec.size())
    {
      idt--;
      res = (func_[ipt-1][idt]*(PtVec[ipt] - Pt_) + func_[ipt][idt]*(Pt_ - PtVec[ipt-1]))/(PtVec[ipt] - DtVec[ipt-1]);
    }else
    {
  //find square points for the bilinear interpolation;


    if(ipt == PtVec.size()+1)
    {

    }

  	double t = (Pt_ - PtVec[ipt-1])/(PtVec[ipt] - PtVec[ipt-1]);
  	double u = (Dt_ - DtVec[idt-1])/(DtVec[idt] - DtVec[idt-1]);

  	double func11 =  func_[ipt-1][idt-1];
    double func12 =  func_[ipt-1][idt];
  	double func21 =  func_[ipt][idt-1];
  	double func22 =  func_[ipt][idt];

  	res = (1.-t)*(1.-u)*func11 + (1.-t)*u*func12 + t*(1-u)*func21 + t*u*func22;
    }

  	return res;
}

/*-------------------------------------------------------------- */

double sigma_class::intA(void)
{
    double result = 0e0;
    double error = 0e0;
    double lower_limit = qMin;
    double upper_limit = qMax;
    double err_abs = 0e0;
    double err_rel = 1e-4;

    size_t Max_steps = 1e8;
    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(Max_steps);

    gsl_function G;
    G.function = &sigma_class::gsl_A;
    G.params = this;

//     std::cout << "Calculating A for Pt= " << Pt << ", Delta=" << Dt << std::flush;

    gsl_integration_qags(&G, lower_limit, upper_limit, err_abs, err_rel, Max_steps, w, &result, &error);

/*
    //Result:
    std::cout << " ---> A = " << result << "\t"
              //<< " Error = " << error << "\t"
              << " Err(%) = "<< 100*fabs(error/result) << std::endl;
*/
    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);
    return result;
}


double sigma_class::intB(void)
{

    double result = 0e0;
    double error = 0e0;
    double lower_limit = qMin;
    double upper_limit = qMax;
    double err_abs = 0e0;
    double err_rel = 1e-4;

    size_t Max_steps = 3e8;
    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(Max_steps);

//    std::cout << "Calculating B for Pt= " << Pt << ", Delta=" << Dt << std::flush;

    gsl_function G;
    G.function = &sigma_class::gsl_B;
    G.params = this;

/* //Massless case
  double result2 = 0e0;
  double error2 = 0e0;
  double delimiter = Pt;

    gsl_function H;
    H.function = &sigma_class::gsl_B2;
    H.params = this;
    gsl_integration_qags(&G, lower_limit, delimiter, err_abs, err_rel, Max_steps, w, &result, &error);
    gsl_integration_qags(&H, delimiter, upper_limit, err_abs, err_rel, Max_steps, w, &result2, &error2);

    result+=result2;
    error+=error2;
*/

    //Massive case:
    gsl_integration_qags(&G, lower_limit, upper_limit, err_abs, err_rel, Max_steps, w, &result, &error);
/*
    //Result:
    std::cout << "  ---> B = " << result << "\t"
              //<< " Error = " << error << "\t"
              << " Err(%) = "<< 100*fabs(error/result) << std::endl;
*/
    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);
    return result;
}


double sigma_class::intC(void)
{
    double result = 0e0;
    double error = 0e0;
    double lower_limit = qMin;
    double upper_limit = qMax;
    double err_abs = 0e0;
    double err_rel = 1e-4;

    size_t Max_steps = 1e8;
    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(Max_steps);

    gsl_function G;
    G.function = &sigma_class::gsl_C;
    G.params = this;

//     std::cout << "Calculating C for Pt= " << Pt << ", Delta=" << Dt << std::flush;

    gsl_integration_qags(&G, lower_limit, upper_limit, err_abs, err_rel, Max_steps, w, &result, &error);

/*
    //Result:
    std::cout << " ---> C = " << result << "\t"
              //<< " Error = " << error << "\t"
              << " Err(%) = "<< 100*fabs(error/result) << std::endl;
*/
    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);
    return result;
}


double sigma_class::intD(void)
{
    double result = 0e0;
    double error = 0e0;
    double lower_limit = qMin;
    double upper_limit = qMax;
    double err_abs = 0e0;
    double err_rel = 1e-4;

    size_t Max_steps = 3e8;
    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(Max_steps);

    gsl_function G;
    G.function = &sigma_class::gsl_D;
    G.params = this;

//     std::cout << "Calculating D for Pt= " << Pt << ", Delta=" << Dt << std::flush;

    gsl_integration_qags(&G, lower_limit, upper_limit, err_abs, err_rel, Max_steps, w, &result, &error);
/*
    //Result:
    std::cout << " ---> D = " << result << "\t"
              //<< " Error = " << error << "\t"
              << " Err(%) = "<< 100*fabs(error/result) << std::endl;
*/
    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);
    return result;
}


double sigma_class::A_func(double q_)
{
  //Massless:
  //double integrand = q_*S_matrix(q_, 0);

  //Massive quark:
  double q2 = q_*q_;
  double str = sqrt( pow((q2 + Pt2 + mf2), 2) - 4.*Pt2*q2);

  double integrand = q_*Pt2*(1. + (Pt2 + mf2 - q2)/str )*S_matrix(q_, 0)/(q2+Pt2+mf2+str);

  return integrand;
}

double sigma_class::B_func(double q_)
{
    double q2 = q_*q_;

    //Massless case:
    //double integrand = - q_*S_matrix(q_, 1)*q2/Pt2;

    //Massive case:
    double str = sqrt( pow((q2 + Pt2 + mf2), 2) - 4.*Pt2*q2);

    double integrand = q_*S_matrix(q_, 1)*(Pt2 - q2 - mf2)
                      *(mf2*mf2 + Pt2*Pt2 + q2*q2 + 2.*mf2*(q2 + Pt2) - (Pt2 + q2 + mf2)*str)/(2.*Pt2*q2*str);

    return integrand;
}

double sigma_class::B_func2(double q_)
{

  double q2 = q_*q_;

  //Only massless case:
  double integrand = q_*S_matrix(q_, 1)*Pt2/q2;

  return integrand;
}


double sigma_class::C_func(double q_)
{
  //Massive quark:
  double q2 = q_*q_;
  double str = sqrt( pow((q2 + Pt2 + mf2), 2) - 4.*Pt2*q2);

  double integrand = q_*Pt2*S_matrix(q_, 0)/str;

  return integrand;
}


double sigma_class::D_func(double q_)
{
    double q2 = q_*q_;
    double str = sqrt( pow((q2 + Pt2 + mf2), 2) - 4.*Pt2*q2);

    double integrand = q_*S_matrix(q_, 1)*( q2 + Pt2 + mf2 - (pow((q2 + Pt2 + mf2), 2) - 2.*Pt2*q2)/str )/q2;

    return integrand;
}

double sigma_class::gsl_A(double x, void * p)
{
    return ((sigma_class *) p)->A_func(x);
}

double sigma_class::gsl_B(double x, void * p)
{
    return ((sigma_class *) p)->B_func(x);
}

double sigma_class::gsl_B2(double x, void * p)
{
    return ((sigma_class *) p)->B_func2(x);
}


double sigma_class::gsl_C(double x, void * p)
{
    return ((sigma_class *) p)->C_func(x);
}

double sigma_class::gsl_D(double x, void * p)
{
    return ((sigma_class *) p)->D_func(x);
}
