
#ifndef SIGMA_H
#define SIGMA_H

#include "MV_nucleus.h"
#include "MV_proton.h"
#include "S_funcs.h"

class sigma_class
{

public:

	//Constructor:
	sigma_class(Smatrix_class & Smatrix_, std::vector<double> Pt_, int iDt_, double R_proj_,
	                          double R_target_, double mf_, double sqrt_s_);
	~sigma_class();

	//Main function
    double operator()(double y1_, double y2_, double Pt_, double Dt_, double dPhi_, int iSigma_);
		//Function for photon flux:
//Grid A and B for given delta:
		void gridFuncs(double delta_);
//Get A (0) or B (1)
		double getFuncs(double Pt_, double Dt_, int iFunc);
		double getSigma(double y1_, double y2_, double Pt_, double Dt_, double dPhi_, int iSigma_);
		double avgSigma(double y1_, double y2_, double Pt_, double Dt_, int itype, int iRes_);
		double getSigma_alternative(double y1_, double y2_, double k1_, double k2_, double Dt_, double Pt_, double dPhi_, double iF_, int iSigma_);


	//Photon flux for crosss section:
		double get_photonflux(double omega_);

private: //Variables

	Smatrix_class & S_matrix;
	double fm, qMin, qMax;
	int iDt, iDtMax;
	std::vector<double> PtVec, DtVec;
	double PtMin, PtMax;
	int iPtMax;
	std::vector< std::vector <double> > Agrid, Bgrid, Cgrid, Dgrid;

//Constants:
  const double Z; 				//Number of protons in the nucleus
  const double alphaem;		// QED coupling
  const double R_projectile;		//Projectile radius
  const double R_target;		//target radius
  const double Ecm;				//center of mass energy: (sqrt(s_NN))
  const double Mp;				//Proton mass
  const double gamma;			//= Ecm/m_p
	int Nc;
	std::vector<double> charge;
	double mf, mf2; 				//heavy quark charge and mass.
	int iFlavour;
//conjugate transverse Variables and angles:
	double Pt, Dt;
	double Pt2, Dt2;

	//Get A and B values by interpolation in Pt:
	double itp_funcs(std::vector<std::vector<double> > func_, double Pt_, double Dt_);

	//Integrate in q_t:
	double intA(void);
	double intB(void);
	double intC(void);
	double intD(void);

	//Wrapper
	static double gsl_A(double x, void * p);
	static double gsl_B(double x, void * p);
	static double gsl_B2(double x, void * p);
	static double gsl_C(double x, void * p);
	static double gsl_D(double x, void * p);

	//Main functions:
	double A_func(double x);
	double B_func(double x);
	double B_func2(double x);
	double C_func(double x);
	double D_func(double x);
};

#endif // SIGMA_H
