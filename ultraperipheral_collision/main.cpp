/*
Author: Mateus Reinke Pelicer;
Contact: mateusreinke@hotmail.com

-This code calculates the pA cross section  for heavy quark production using the MC-Lerran Venugopalan Model.
-sigma_class calculates the cross section (cs) and integrates all cs level functions in q_t.
-S_funcs calculates the dipole-S matrix nod and elliptical components in momentum space making a Hankel transform of
the dipole-S matrix in config. space, where the MC-Lerran Venugopalan model is used.
-MV-funcs calculates the functions in the  extended Mc-Lerran Venugopalan Model proposed in Phys. Rev. D 95, 094003 (2017)
*/

//#include "MV_funcs.h"


/*To do w/ FF's:
	 Check kinematics;
	 Pion (+, 0, -) FF's for light quarks
	 D-meson FF's for c quarks
	*/
#include "sigma.h"
#include "FF_Charm.cpp"
#include "LHAPDF/LHAPDF.h"

double sigma_hadronic(sigma_class sigma_, const LHAPDF::PDF* ff_pi_, std::vector<double> ptVec_, double y1_, double y2_,
											double x1_, double  x2_, double PtH_, double DtH_, double PhiPD_, double dt_, double P1_, double  P2_);

int main()
{
//Main parameters:
	double mgluon = 0.25, fm = 1/0.197,  //(GeV^-1)
					R = 2., Q0s2 = 1./(R*R),
				 	A = 208, R_projectile = 1.12*fm*pow(A, 1./3.), //R_target = 1.12*fm*pow(A, 1./3.);
				 																								R_target = 0.8*fm;
				 //Pb: A = 208; Z=82
				 //Au: A = 196; Z=79
	double y1=1., y2=1.;
	double sqrt_s= 2.76e3; //or 1e3 GeV =1 TeV  -- 2.76 or 5.02 for  AA and 5.02 and 13 for pA
	double qMin = 0.01, qMax = 30.;
	double dphi = 0, mf;
	double nb = 3.89379e5; //nb=1GeV^-2

	int iF = 0; //0: massless, 1:Charm quark, 2:Bottom quark
	if(iF == 0)			{mf = 0.00;}
	else if(iF == 1){mf = 1.40;}
	else if(iF == 2){mf = 4.70;}
	else {std::cout << "Invalid quark mass!" << std::endl; exit(EXIT_FAILURE);}

	//partonic tranverse variables:
	double PtMin, PtMax, dPt, deltaMin, deltaMax, dDelta;
	int iPtMax, iDeltaMax;

	//double sigma_average2;
	std::vector<double> PtVec, deltaVec;
	std::vector<double> Pt_save, Dt_save;


//Define quark momentum limits:
 	PtMin=1., PtMax=25.2, dPt=0.1;
 	deltaMin=0.1, deltaMax=1., dDelta=0.005;

	iPtMax = (int) std::floor( (PtMax - PtMin)/dPt)+1;
	iDeltaMax = (int) std::floor( (deltaMax - deltaMin)/dDelta)+1;
//Resize all vectors:
	PtVec.resize(iPtMax, 0);
	deltaVec.resize(iDeltaMax, 0);
//Fill vectors:
	for(int i = 0; i<iPtMax; i++) 	PtVec[i] = PtMin+i*dPt;
	for(int i = 0; i<iDeltaMax; i++) deltaVec[i] = deltaMin+i*dDelta;

//Save data for fixed Deltat or fixed Pt? What values:
	Dt_save = {0.2};
	Pt_save = {1., 5., 10.};

	for(size_t i=0; i<Dt_save.size(); i++){if(Dt_save[i] > deltaMax || Dt_save[i] <deltaMin)
																  {std::cout << "Dt_save out of range!" << std::endl; return 1;}}

	//Number of  calculated points:
	std::cout << "PtMin = " << PtVec[0] << " PtMax = " << PtVec[iPtMax-1]
						<< " dPt = " << dPt << " iPtMax= " << iPtMax << std::endl;
	std::cout << "deltaMin = " << deltaVec[0] << " deltaMax = " << deltaVec[iDeltaMax-1]
						<< " dDelta = " << dDelta << " iDeltaMax= " << iDeltaMax  << std::endl;

//Construct model Mc-Lerran Venugopalan model - must alter Smatrix constructor if changing this
	MV_proton_class MV(fm, R, mgluon, Q0s2);
	//MV_nucleus_class MV(fm, R, mgluon, Q0s2, A, R_target);

//Construct S0 and S_elliptic components of S matrix in momentum space:
	Smatrix_class Smatrix(MV, qMin, qMax, R_target);
//Construct A and B functions and get differential cross section:
	sigma_class sigma(Smatrix, PtVec, iDeltaMax, R_projectile, R_target, mf, sqrt_s);


//Loop over dDelta_T to save all needed functions:
for(int iDt = 0; iDt< iDeltaMax; iDt++)
{
	std::cout << "Running code for deltat= " << deltaVec[iDt] << std::endl;
	Smatrix.calc_integrals(deltaVec[iDt]);
//Save all cross-section level functions.
	sigma.gridFuncs(deltaVec[iDt]);
}

//what should be saved?
	bool saveDt=true;
	bool savePt = false;
	bool doFF = true;
	bool doFF2 = false;


//do FF's w/ PtH and DtH:
if(doFF)
{

	std::cout << "Hadronizing dipole: Save as (PtH, Dt, DtH) w/ fixed angles! " << std::endl;
	const LHAPDF::PDF* ff_pi1= LHAPDF::mkPDF("NNFF10_PIp_nlo",  0);
	const LHAPDF::PDF* ff_pi2= LHAPDF::mkPDF("NNFF10_PIm_nlo",  0);

	double ff1, ff2;
	//hadronic transverse variables and angle:
	double PtHMin, PtHMax, dPtH, DtHMin, DtHMax, dDtH, phiDMin, phiDMax, dphiD,  phiPMin, phiPMax, dphiP;
	int iPtHMax, iDtHMax, iphiDMax, iphiPMax;
	//Limits of integration in fragmentation:
	PtHMin = 0., PtHMax=11., dPtH=1., DtHMin=0.1, DtHMax=2., dDtH=0.005;
	phiDMin = 0.; phiDMax = M_PI; dphiD = M_PI/10;
	phiPMin = 0.; phiPMax = M_PI; dphiP = M_PI/10;


	iDtHMax = (int) std::floor( (DtHMax - DtHMin)/dDtH) ;
	std::vector<double> DtHVec(iDtHMax);
	for(int iDtH = 0; iDtH<iDtHMax; iDtH++){DtHVec[iDtH] = DtHMin + iDtH*dDtH;}

	iPtHMax = (int) std::floor( (PtHMax - PtHMin)/dPtH) ;
	iphiDMax = (int) std::floor( (phiDMax - phiDMin)/dphiD);
	iphiPMax = (int) std::floor( (phiPMax - phiPMin)/dphiP);

	std::vector<double> PtHVec(3) , phiDVec(4),  phiPVec(9);
	PtHVec = {2.0, 5.0, 8.0};
	phiDVec= {0.0, M_PI/10., M_PI/4. , M_PI/2.};
	phiPVec= {0.0, M_PI/10., M_PI/5. , M_PI/2., M_PI};
 	for(int i = 0; i<=iphiDMax; i++){phiDVec[i] = phiDMin + i*dphiD;}
 	for(int i = 0; i<=iphiPMax; i++){phiPVec[i] = phiPMin + i*dphiP;}

	double PtH, DtH, DtH2, PtH2, phiDH, phiPH, phiPDH;
	double Y1 = y1;
	double Y2 = y2;
	double P1, P2, phi1, phi2, phi21;
	double k1, k2, x1, x2, pt, Dt, phiPD;
	double res_ff;
	double  m_M1, mt_M1, m_M2, mt_M2;
	m_M1 = 0., m_M2=0;
	double mt1, mt2;


//Loop over angle:
	for(int iphid=0; iphid<= iphiDMax; iphid++)
	{
	for(int iphip=0; iphip<= iphiPMax; iphip++)
	{
			phiDH = phiDVec[iphid];
			phiPH = phiPVec[iphip];
			phiPDH =phiPH - phiDH;


			//Loop over all hadronic momentum
			for(int iPtH = 0; iPtH<iPtHMax; iPtH++)
			{
				PtH = PtHVec[iPtH];
				PtH2= PtH*PtH;
				std::ofstream sigmaPions;
//				sigmaPions.open("cs_pion" + std::to_string(PtHVec[iPtH]) +"_"+std::to_string(ip)+"_"+std::to_string(id) + ".txt");
				sigmaPions.open("cs_pionDt_PtH" + std::to_string(PtHVec[iPtH])
																	+ "_phiDH"+std::to_string(phiDH)+ "_phiPH"+std::to_string(phiPH)+".txt");

				sigmaPions << "#sigma, Dt, DtH --- PtH, phiPH, phiDH, Y1,2= "
									<< PtH << " " << phiPH << " " << phiDH << " " << y1 << " " << y2
									<< std::endl;

				for(int iDtH = 0; iDtH<iDtHMax; iDtH++)
				{
					DtH = DtHVec[iDtH];
					DtH2= DtH*DtH;
					//meson momentum
					P1 = sqrt( PtH2 + DtH2/4. + PtH*DtH*cos(phiPDH));
					P2 = sqrt( PtH2 + DtH2/4. - PtH*DtH*cos(phiPDH));
					phi1 = atan2( -PtH*sin(phiPH)- DtH*sin(phiDH)/2. , -PtH*cos(phiPH)- DtH*cos(phiDH)/2. ) ;
					phi2 = atan2( +PtH*sin(phiPH)- DtH*sin(phiDH)/2. , +PtH*cos(phiPH)- DtH*cos(phiDH)/2. );
					phi21 = phi2 - phi1;
					//transverse meson mass
					mt_M1= sqrt(P1*P1 +m_M1*m_M1);
					mt_M2= sqrt(P2*P2 +m_M2*m_M2);

					for(size_t iDt = 0; iDt< deltaVec.size(); iDt++)
					{
						Dt  = deltaVec[iDt];
						//partonic variables:
						x1 = -P1*sin(phi21)/(Dt*sin(phi2));
						x2 =  P2*sin(phi21)/(Dt*sin(phi1));
						k1 = P1/x1;
						k2 = P2/x2;
						pt = sqrt(k1*k1+k2*k2-2.*k1*k2*cos(phi21))/2.;
						phiPD = acos( (-k2*k2+k1*k1)/(2.*pt*Dt) );
						mt1= sqrt(k1*k1 +mf*mf);
						mt2= sqrt(k2*k2 +mf*mf);
						if(mf!=0){Y1+=log(x1*mt1/mt_M1);}
						if(mf!=0){Y2+=log(x1*mt2/mt_M2);}

						res_ff = 0.;

						if(pt< PtVec[0] || pt> PtVec[iPtMax-1] || std::fabs(cos(phiPD))>1
								|| x1 < 0 || x1 > 1 || x2 < 0 || x2 > 1 || k1<=1. || k2 <=1. || x1 != x1 || x2!=x2)
						{
							//if(Dt)
							//std::cout << sin(phi1) << " " << sin(phi2) << " " << x1 << " " << x2 << " " << pt << " " << Dt << " " << phi21 <<std::endl;
							res_ff=0;
						}else
						{
							for(int iFlavour=1; iFlavour<=3; iFlavour++)//1,2,3=d,u,s
							{
								ff1 =  ff_pi1->xfxQ2( iFlavour, x1, k1*k1);
								ff2 =  ff_pi2->xfxQ2(-iFlavour, x2, k2*k2);
								res_ff+=ff1*ff2*sigma.getSigma_alternative(y1, y2, k1, k2, Dt, pt, phiPD, iFlavour, 0)/(P1*P2*std::fabs(sin(phi21)) );
							}
						}
						sigmaPions << res_ff << " " << Dt << " " << DtH << std::endl;
						//sigma_mesons[iDt][iDtH][iPtH] = res_ff;
					}//end  Dt loop
				}//end DtH loop
				sigmaPions.close();
			}//end PtH loop

		}
	}

}


/*------------------------------------------------------------------------------------------------------------------*/
//do FF's w/ PtH and DtH:
if(doFF2)
{

	std::cout << "Hadronizing dipole! Save as (PtH, PhiDH, DtH) w/ fixed phiPH and Dt" << std::endl;
	const LHAPDF::PDF* ff_pi1= LHAPDF::mkPDF("NNFF10_PIp_nlo",  0);
	const LHAPDF::PDF* ff_pi2= LHAPDF::mkPDF("NNFF10_PIm_nlo",  0);

	double ff1, ff2;
	//hadronic transverse variables and angle:
	double PtHMin, PtHMax, dPtH, DtHMin, DtHMax, dDtH, phiDMin, phiDMax, dphiD,  phiPMin, phiPMax, dphiP;
	int iPtHMax, iDtHMax, iphiDMax, iphiPMax;
	//Limits of integration in fragmentation:
	PtHMin = 0., PtHMax=11., dPtH=1., DtHMin=0.1, DtHMax=2., dDtH=0.05;
	phiDMin = 0.; phiDMax = 2.*M_PI; dphiD = M_PI/100;
	phiPMin = 0.; phiPMax = M_PI; dphiP = M_PI/10;

	iPtHMax = (int) std::floor( (PtHMax - PtHMin)/dPtH) ;
	iDtHMax = (int) std::floor( (DtHMax - DtHMin)/dDtH) ;
	iphiDMax = (int) std::floor( (phiDMax - phiDMin)/dphiD);
	iphiPMax = (int) std::floor( (phiPMax - phiPMin)/dphiP);

 std::vector<double> PtHVec(iPtHMax), DtHVec(iDtHMax), phiDVec(iphiDMax),  phiPVec(iphiPMax);
 for(int iPtH = 0; iPtH<iPtHMax; iPtH++){PtHVec[iPtH] = PtHMin + iPtH*dPtH;}
 for(int iDtH = 0; iDtH<iDtHMax; iDtH++){DtHVec[iDtH] = DtHMin + iDtH*dDtH;}
 for(int i = 0; i<iphiDMax; i++){phiDVec[i] = phiDMin + i*dphiD;}
 for(int i = 0; i<iphiPMax; i++){phiPVec[i] = phiPMin + i*dphiP;}

	double PtH, DtH, DtH2, PtH2, phiDH, phiPH, phiPDH;
	double Y1 = y1;
	double Y2 = y2;
	double P1, P2, phi1, phi2, phi21;
	double k1, k2, x1, x2, pt, Dt, phiPD;
	double res_ff;
	double  m_M1, mt_M1, m_M2, mt_M2;
	m_M1 = 0., m_M2=0;
	double mt1, mt2;


		for(int iphip=0; iphip<iphiPMax; iphip++)
		{
		for(size_t idt=0; idt<Dt_save.size(); idt++)
		{
				phiPH =phiPVec[iphip];
				Dt  = Dt_save[idt];
			//Loop over all hadronic momentum
			for(int iPtH = 0; iPtH<iPtHMax; iPtH++)
			{
				PtH = PtHVec[iPtH];
				PtH2= PtH*PtH;
				std::ofstream sigmaPions;
//				sigmaPions.open("cs_pion" + std::to_string(PtHVec[iPtH]) +"_"+std::to_string(ip)+"_"+std::to_string(id) + ".txt");
				sigmaPions.open("cs_pionPhiDH_PtH" + std::to_string(PtHVec[iPtH])
																	+ "_Dt"+std::to_string(Dt_save[idt])+ "_phiPH"+std::to_string(phiPH)+".txt");

				sigmaPions << "#sigma, phiDH, DtH --- PtH, Dt, phiPH, Y1,2= "
									<< PtH << " " << phiPH<< " " << Dt << " " << y1 << " " << y2
									<< std::endl;

				for(int iDtH = 0; iDtH<iDtHMax; iDtH++)
				{
					DtH = DtHVec[iDtH];
					DtH2= DtH*DtH;

					for(int iphiDH = 0; iphiDH<iphiDMax ; iphiDH++)
					{
						phiDH=phiDVec[iphiDH];
						phiPDH=phiPH - phiDH;

						//meson momenta
						P1 = sqrt( PtH2 + DtH2/4. + PtH*DtH*cos(phiPDH));
						P2 = sqrt( PtH2 + DtH2/4. - PtH*DtH*cos(phiPDH));
						phi1 = atan2( -PtH*sin(phiPH)- DtH*sin(phiDH)/2. , -PtH*cos(phiPH)- DtH*cos(phiDH)/2. ) ;
						phi2 = atan2( +PtH*sin(phiPH)- DtH*sin(phiDH)/2. , +PtH*cos(phiPH)- DtH*cos(phiDH)/2. );
						phi21 = phi2 - phi1;
						//transverse meson mass
						mt_M1= sqrt(P1*P1 +m_M1*m_M1);
						mt_M2= sqrt(P2*P2 +m_M2*m_M2);

						//partonic variables:
						x1 = -P1*sin(phi21)/(Dt*sin(phi2));
						x2 =  P2*sin(phi21)/(Dt*sin(phi1));
						k1 = P1/x1;
						k2 = P2/x2;
						pt = sqrt(k1*k1+k2*k2-2.*k1*k2*cos(phi21))/2.;
						phiPD = acos( (-k2*k2+k1*k1)/(2.*pt*Dt) );
						mt1= sqrt(k1*k1 +mf*mf);
						mt2= sqrt(k2*k2 +mf*mf);
						if(mf!=0){Y1+=log(x1*mt1/mt_M1);}
						if(mf!=0){Y2+=log(x1*mt2/mt_M2);}

						res_ff = 0.;

						if(pt< PtVec[0] || pt> PtVec[iPtMax-1] || std::fabs(cos(phiPD))>1
								|| x1 < 0 || x1 > 1 || x2 < 0 || x2 > 1 || k1<=1. || k2 <=1. || x1 != x1 || x2!=x2)
						{
							//if(Dt)
							//std::cout << sin(phi1) << " " << sin(phi2) << " " << x1 << " " << x2 << " " << pt << " " << Dt << " " << phi21 <<std::endl;
							res_ff=0;
						}else
						{
							for(int iFlavour=1; iFlavour<=3; iFlavour++)//1,2,3=d,u,s
							{
								ff1 =  ff_pi1->xfxQ2( iFlavour, x1, k1*k1);
								ff2 =  ff_pi2->xfxQ2(-iFlavour, x2, k2*k2);
								res_ff+=ff1*ff2*sigma.getSigma_alternative(y1, y2, k1, k2, Dt, pt, phiPD, iFlavour, 0)/(P1*P2*std::fabs(sin(phi21)) );
							}
						}
						sigmaPions << res_ff << " " << phiDH << " " << DtH << std::endl;
						//sigma_mesons[iDt][iDtH][iPtH] = res_ff;
					}//end  Dt loop
				}//end DtH loop
				sigmaPions.close();
			}//end PtH loop

		}
		}

}


/*
//Calculate hadronic cross section w/ P1, P2, y1, y2 dependence (delta_t integrated):
if(doFF2)
{

std::cout << "Calculating hadronic cross section: P1, P2, phi21: " << std::endl;

//define Meson properties:
double m_M1, mt_M1, m_M2, mt_M2;
m_M1 = 0., m_M2=0;
//hadronic transverse variables and angle:
std::vector<double> P1, P2, phi21;

double P1Min, P1Max, P2Min, P2Max, dP1, dP2, phi21Min, phi21Max, dphi21;
int iP1Max, iP2Max, iphi21Max;
//Limits of integration in fragmentation:

P1Min = 1., 	P1Max=10.1,	 		iP1Max = 50;
P2Min = 1., 	P2Max=10.1,			iP2Max = 50;
phi21Min= 0.; phi21Max=M_PI;	iphi21Max= 20;
P1.resize(iP1Max+1);
P2.resize(iP2Max+1);
phi21.resize(iphi21Max+1);

dP1 = (P1Max - P1Min)/iP1Max;
dP2 = (P2Max - P2Min)/iP2Max;
dphi21= (phi21Max-phi21Min)/iphi21Max;

for(int iP1 = 0; iP1<=iP1Max; iP1++){P1[iP1] = P1Min + iP1*dP1;}
for(int iP2 = 0; iP2<=iP2Max; iP2++){P2[iP2] = P2Min + iP2*dP2;}
for(int ip = 0; ip<=iphi21Max; ip++){phi21[ip] = phi21Min + ip*dphi21;}

std::cout << "Boundaries: " << std::endl;
std::cout << "P1: " << P1[0] << " " << P1[iP1Max] << std::endl
					<< "P2: " << P2[0] << " " << P2[iP2Max] << std::endl
					<< "phi21: "<< phi21[0]<< " " << phi21[iphi21Max]<< std::endl;

//momentum fraction definitions:

double x1Min, x1Max, x2Min, x2Max, dx1, dx2;
double Nx1, Nx2, x1, x2;
//lower boundaries defined inside transverse momentum loop;
x1Max = 1., x2Max=1.;
Nx1 = 1000., Nx2 = 1000.;

std::vector< std::vector <std::vector<double> > > sigma_mesons;
sigma_mesons.resize( iP1Max+1, std::vector<std::vector<double>>(iP2Max+1,std::vector<double>(iphi21Max+1,0e0)));

std::cout << "Size of sigma_vector: " << sigma_mesons.size() << " " << sigma_mesons[0].size() << " " <<  sigma_mesons[0][0].size()  << std::endl;

//Define fragmentation function:
const LHAPDF::PDF* ff_pi= LHAPDF::mkPDF("NNFF10_PIp_nlo",  0); //Pion plus!
double ff1, ff2;
double res_ff;


//Start saving cross sections:
for(int ip = 0; ip<=iphi21Max; ip++)
{

std::ofstream cs_pions;
cs_pions.open("cs_pions"+ std::to_string(ip) +".txt");

cs_pions << "P1, P2, sigma. --- phi21 = " << phi21[ip] << " y1, y2= " << y1 << " " << y2 << std::endl;

	for(int i1 = 0; i1<= iP1Max; i1++)
	{
	mt_M1 = sqrt(P1[i1]*P1[i1] + m_M1*m_M1);
	x1Min = (mt_M1/sqrt_s)*exp(y1);
	dx1 = (x1Max - x1Min)/Nx1;

		for(int i2 = 0; i2<=iP2Max; i2++)
		{

		std::cout << "Integrating cross section for: " << P1[i1] << " " << P2[i2] << " " << phi21[ip] << std::endl;
		mt_M2 = sqrt(P2[i2]*P2[i2] + m_M2*m_M2);
		x2Min = (mt_M2/sqrt_s)*exp(y2);
		 //wrong! it should be the photon-proton sqrt_s, since the nucleus effect is fully englobed in the the WW photon dist.
		dx2 = (x2Max - x2Min)/Nx2;

		res_ff=0.;

		//is integration necessary or are the points out of phase space?
		double dtmin = sqrt(P1[i1]*P1[i1] + P2[i2]*P2[i2] + 2.*P1[i1]*P2[i2]*cos(phi21[ip]));
		double dtmax = sqrt(P1[i1]*P1[i1]/(x1Min*x1Min) + P2[i2]*P2[i2]/(x2Min*x2Min) + 2.*P1[i1]*P2[i2]*cos(phi21[ip])/(x1Min*x2Min));
		double ptmin = sqrt(P1[i1]*P1[i1] + P2[i2]*P2[i2] - 2.*P1[i1]*P2[i2]*cos(phi21[ip]))/2.;
		double ptmax = sqrt(P1[i1]*P1[i1]/(x1Min*x1Min) + P2[i2]*P2[i2]/(x2Min*x2Min) - 2.*P1[i1]*P2[i2]*cos(phi21[ip])/(x1Min*x2Min))/2.;

			if(dtmax< deltaVec[0] || dtmin>deltaVec[iDeltaMax-1] || ptmax< PtVec[0] || ptmin> PtVec[iPtMax-1])
			{
				std::cout << "Completely out of phase space! " << dtmin << " " << dtmax << " "
				 																							 << ptmin << " " << ptmax << std::endl;
			}else
			{
				//DO INTEGRATION!
				for(int ix1 = 0; ix1<Nx1; ix1++)
				{
					x1 = x1Min  +(2*ix1+1)*dx1/2.;

					for(int ix2 = 0; ix2<Nx2; ix2++)
					{
						x2 = x2Min + (2*ix2+1)*dx2/2.;

						double k1t  = P1[i1]/x1;
						double k2t  = P2[i2]/x2;
						double pt 	= sqrt(k1t*k1t + k2t*k2t - 2.*k1t*k2t*cos(phi21[ip]) )/2.;
						double dt 	= sqrt(k1t*k1t + k2t*k2t + 2.*k1t*k2t*cos(phi21[ip]) );
						double cosPD=(k1t*k1t - k2t*k2t)/(2.*pt*dt);
						double phiPD= acos(cosPD);

						if(dt< deltaVec[0] || dt>deltaVec[iDeltaMax-1] || pt<= PtVec[0] || pt> PtVec[iPtMax-1] || std::fabs(cosPD)>1
								|| x1Min > x1Max || x2Min > x2Max || pt<=dt || k1t<=1. || k2t <=1.)
						{
						}else
						{
							for(int iFlavour=1; iFlavour<=3; iFlavour++)
							{
								ff1 =  ff_pi->xfxQ2( iFlavour, x1, P1[i1]*P1[i1]);
								ff2 =  ff_pi->xfxQ2(-iFlavour, x2, P2[i2]*P2[i2]);
								res_ff+= sigma.getSigma_alternative(y1, y2, k1t, k2t, dt, pt, phiPD, iFlavour, 0)
												*ff1*ff2/(x1*x1*x2*x2);
							}
						}
					}//end x2 loop
				}//end x1 loop
				res_ff*=nb*dx1*dx2; //pass result to nanobars/GeV^4
			}
				sigma_mesons[i1][i2][ip] = res_ff;

			  cs_pions << P1[i1] << " " << P2[i2] << " " << sigma_mesons[i1][i2][ip] << std::endl;

			}//end P1 loop
	}//end P2 loop

cs_pions.close();
}//end  c12 loop

}//end do FF2

*/

std::string nameA, nameB, nameC, nameD, nameS, nameSavg, nameSavgABC, nameSavgAB, nameSavgA, nameSavgCD;
double sigma_differential, sigma_average, sigmaCos_average;


if(saveDt==true)
{


for(size_t iDt = 0; iDt< Dt_save.size(); iDt++)
{

//Define and open output files:
	std::ofstream outA, outB, outC, outD, outSigma, outSigmaAvg;
	std::ofstream outSigmaAvgABC, outSigmaAvgAB, outSigmaAvgA, outSigmaAvgCD;


	nameA= "A_dlt" + std::to_string(Dt_save[iDt]) + ".dat";
	nameB= "B_dlt" + std::to_string(Dt_save[iDt]) + ".dat";
	nameC= "C_dlt" + std::to_string(Dt_save[iDt]) + ".dat";
	nameD= "D_dlt" + std::to_string(Dt_save[iDt]) + ".dat";
	nameS = "sigma_dlt" + std::to_string(Dt_save[iDt]) + ".dat";
	nameSavg = "sigmaAvg_dlt" + std::to_string(Dt_save[iDt]) + ".dat";
	nameSavgABC  = "sigmaAvgABC_dlt" + std::to_string(Dt_save[iDt]) + ".dat";
	nameSavgAB  = "sigmaAvgAB_dlt" + std::to_string(Dt_save[iDt]) + ".dat";
	nameSavgA  = "sigmaAvgA_dlt" + std::to_string(Dt_save[iDt]) + ".dat";
  nameSavgCD  = "sigmaAvgCD_dlt" + std::to_string(Dt_save[iDt]) + ".dat";

	outA.open(nameA);
	outB.open(nameB);
	outA << "#P_t, Delta_t, A" << std::endl;
	outB << "#P_t, Delta_t, B" << std::endl;
	if(mf!=0)
	{
	outC.open(nameC);
	outD.open(nameD);
	outC << "#P_t, Delta_t, C" << std::endl;
	outD << "#P_t, Delta_t, D" << std::endl;
	}

/*
	outSigma.open(nameS);
	outSigmaAvg.open(nameSavg);
	outSigmaAvgABC.open(nameSavgABC);
	outSigmaAvgAB.open(nameSavgAB);
	outSigmaAvgA.open(nameSavgA);
	outSigmaAvgCD.open(nameSavgCD);

	outSigma << "#P_t, Delta_t, sigma(phi= 0, M_PI/4, M_PI/2, 3*M_PI/4, M_PI) ----  y1, y2 = " << y1 <<" "<< y2 << std::endl;
	outSigmaAvg << "#P_t, Delta_t, <sigma>  <sigma*cos(2phi)>   Ratio(<sgm*cos> /<sgm>)  ----  y1, y2 = " << y1 <<" "<< y2 << std::endl;
	outSigmaAvgABC << "#P_t, Delta_t, <sigma>  <sigma*cos(2phi)>   Ratio(<sgm*cos> /<sgm>)  ----  y1, y2 = " << y1 <<" "<< y2 << std::endl;
	outSigmaAvgAB << "#P_t, Delta_t, <sigma>  <sigma*cos(2phi)>   Ratio(<sgm*cos> /<sgm>)  ----  y1, y2 = " << y1 <<" "<< y2 << std::endl;
	outSigmaAvgA << "#P_t, Delta_t, <sigma>  <sigma*cos(2phi)>   Ratio(<sgm*cos> /<sgm>)  ----  y1, y2 = " << y1 <<" "<< y2 << std::endl;
	outSigmaAvgCD << "#P_t, Delta_t, <sigma>  <sigma*cos(2phi)>   Ratio(<sgm*cos> /<sgm>)  ----  y1, y2 = " << y1 <<" "<< y2 << std::endl;

*/
	std::cout << "Saving partonic grids in delta_t= " << Dt_save[iDt] <<  std::endl;

//Loop over Pt and save all
	for(int iPt = 0; iPt<iPtMax; iPt++)
	{

		outA << PtVec[iPt] << " " << Dt_save[iDt] << " " << sigma.getFuncs(PtVec[iPt], Dt_save[iDt], 0) << std::endl;
		outB << PtVec[iPt] << " " << Dt_save[iDt] << " " << sigma.getFuncs(PtVec[iPt], Dt_save[iDt], 1) << std::endl;
		if(mf!= 0)
		{
		outC << PtVec[iPt] << " " << Dt_save[iDt] << " " << sigma.getFuncs(PtVec[iPt], Dt_save[iDt], 2) << std::endl;
		outD << PtVec[iPt] << " " << Dt_save[iDt] << " " << sigma.getFuncs(PtVec[iPt], Dt_save[iDt], 3) << std::endl;
		}
/*
		sigma_average 	  = nb*sigma.avgSigma(y1, y2, PtVec[iPt], Dt_save[iDt], 0, 0);
		sigmaCos_average  = nb*sigma.avgSigma(y1, y2, PtVec[iPt], Dt_save[iDt], 1, 0);
		outSigmaAvg << PtVec[iPt] << " " << Dt_save[iDt]  << " "
								<< sigma_average  << " " <<  sigmaCos_average << " "
								<< sigmaCos_average/sigma_average
								<< std::endl;

		double sigma_averageABC 	  = nb*sigma.avgSigma(y1, y2, PtVec[iPt], Dt_save[iDt], 0, 1);
		double sigmaCos_averageABC  = nb*sigma.avgSigma(y1, y2, PtVec[iPt], Dt_save[iDt], 1, 1);
		outSigmaAvgABC << PtVec[iPt] << " " << Dt_save[iDt]  << " "
									 << sigma_averageABC  << " " <<  sigmaCos_averageABC << " "
									 << sigmaCos_averageABC/sigma_average
									 << std::endl;

		double sigma_averageAB 	 	 = nb*sigma.avgSigma(y1, y2, PtVec[iPt], Dt_save[iDt], 0, 2);
		double sigmaCos_averageAB  = nb*sigma.avgSigma(y1, y2, PtVec[iPt], Dt_save[iDt], 1, 2);
		outSigmaAvgAB << PtVec[iPt] << " " << Dt_save[iDt]  << " "
									<< sigma_averageAB  << " " <<  sigmaCos_averageAB << " "
									<< sigmaCos_averageAB/sigma_average
									<< std::endl;

		double sigma_averageA 	  = nb*sigma.avgSigma(y1, y2, PtVec[iPt], Dt_save[iDt], 0, 3);
		double sigmaCos_averageA  = nb*sigma.avgSigma(y1, y2, PtVec[iPt], Dt_save[iDt], 1, 3);
		outSigmaAvgA << PtVec[iPt] << " " << Dt_save[iDt]  << " "
								 << sigma_averageA  << " " <<  sigmaCos_averageA << " "
							 	 << sigmaCos_averageA/sigma_average
								 << std::endl;

		double sigma_averageCD 	  = nb*sigma.avgSigma(y1, y2, PtVec[iPt], Dt_save[iDt], 0, 4);
		double sigmaCos_averageCD  = nb*sigma.avgSigma(y1, y2, PtVec[iPt], Dt_save[iDt],1, 4);
		outSigmaAvgCD << PtVec[iPt] << " " << Dt_save[iDt]  << " "
									<< sigma_averageCD  << " " <<  sigmaCos_averageCD << " "
									<< sigmaCos_averageCD/sigma_average
									<< std::endl;


*/
	}
	outA.close();
	outB.close();
//	outSigma.close();
//	outSigmaAvg.close();
	if(mf!=0){outC.close(); outD.close(); }

}
}
//============================== Save in Pt: ====================================

if(savePt==true)
{

	for(size_t iPt = 0; iPt< Pt_save.size(); iPt++)
	{

	//Define and open output files:
		std::ofstream outA, outB, outC, outD, outSigma, outSigmaAvg;
		std::ofstream outSigmaAvgABC, outSigmaAvgAB, outSigmaAvgA, outSigmaAvgCD;

		nameA= "A_pt" + std::to_string(Pt_save[iPt]) + ".dat";
		nameB= "B_pt" + std::to_string(Pt_save[iPt]) + ".dat";
		nameC= "C_pt" + std::to_string(Pt_save[iPt]) + ".dat";
		nameD= "D_pt" + std::to_string(Pt_save[iPt]) + ".dat";
		nameS = "sigma_pt" + std::to_string(Pt_save[iPt]) + ".dat";
		nameSavg = "sigmaAvg_pt" 				+ std::to_string(Pt_save[iPt]) + ".dat";
		nameSavgABC  = "sigmaAvgABC_pt" + std::to_string(Pt_save[iPt]) + ".dat";
		nameSavgAB  = "sigmaAvgAB_pt" 	+ std::to_string(Pt_save[iPt]) + ".dat";
		nameSavgA  = "sigmaAvgA_pt" 		+ std::to_string(Pt_save[iPt]) + ".dat";
	  nameSavgCD  = "sigmaAvgCD_pt" 	+ std::to_string(Pt_save[iPt]) + ".dat";

		outA.open(nameA);
		outB.open(nameB);
		outSigma.open(nameS);
		outSigmaAvg.open(nameSavg);
		outSigmaAvgABC.open(nameSavgABC);
		outSigmaAvgAB.open(nameSavgAB);
		outSigmaAvgA.open(nameSavgA);
		outSigmaAvgCD.open(nameSavgCD);

		outA << "#Delta_t, P_t, A" << std::endl;
		outB << "#Delta_t, P_t, B" << std::endl;
		outSigma << "#Delta_t, P_t, sigma(phi= 0, M_PI/4, M_PI/2, 3*M_PI/4, M_PI) ----  y1, y2 = " << y1 <<" "<< y2 << std::endl;
		outSigmaAvg << "#Delta_t, P_t, <sigma>  <sigma*cos(2phi)>   Ratio(<sgm*cos> /<sgm>)  ----  y1, y2 = " << y1 <<" "<< y2 << std::endl;
		outSigmaAvgABC << "#Delta_t, P_t, <sigma>    <sigma*cos(2phi)>   Ratio(<sgm*cos> /<sgm>)  ----  y1, y2 = " << y1 <<" "<< y2 << std::endl;
		outSigmaAvgAB << "#Delta_t, P_t, <sigma>    <sigma*cos(2phi)>   Ratio(<sgm*cos> /<sgm>)  ----  y1, y2 = " << y1 <<" "<< y2 << std::endl;
		outSigmaAvgA << "#Delta_t, P_t, <sigma>  <sigma*cos(2phi)>   Ratio(<sgm*cos> /<sgm>)  ----  y1, y2 = " << y1 <<" "<< y2 << std::endl;
		outSigmaAvgCD << "#Delta_t, P_t, <sigma>    <sigma*cos(2phi)>   Ratio(<sgm*cos> /<sgm>)  ----  y1, y2 = " << y1 <<" "<< y2 << std::endl;

		if(mf!=0)
		{
		outC.open(nameC);
		outD.open(nameD);
		outC << "#Delta_t, P_t, C" << std::endl;
		outD << "#Delta_t, P_t, D" << std::endl;
		}

		std::cout << "Averaging over angle for P_t= " << Pt_save[iPt] <<  std::endl;

	//Loop over Pt and save all
		for(int iDt = 0; iDt<iDeltaMax; iDt++)
		{

			outA << deltaVec[iDt] << " " << Pt_save[iPt] << " " << sigma.getFuncs(Pt_save[iPt], deltaVec[iDt], 0) << std::endl;
			outB << deltaVec[iDt] << " " << Pt_save[iPt] << " " << sigma.getFuncs(Pt_save[iPt], deltaVec[iDt], 1) << std::endl;
			if(mf!= 0)
			{
			outC << deltaVec[iDt] << " " << Pt_save[iPt] << " " << sigma.getFuncs(Pt_save[iPt], deltaVec[iDt], 2) << std::endl;
			outD << deltaVec[iDt] << " " << Pt_save[iPt] << " " << sigma.getFuncs(Pt_save[iPt], deltaVec[iDt], 3) << std::endl;
			}

	//Loop over angle:
			outSigma << deltaVec[iDt] << " " <<  Pt_save[iPt] << " " << std::flush;
			for(int iPhi = 0; iPhi<=4; iPhi++)
			{
				dphi = iPhi*M_PI/4.;
				sigma_differential = nb*sigma.getSigma(y1, y2, Pt_save[iPt], deltaVec[iDt], dphi, 0);

				outSigma << sigma_differential << " " << std::flush;
			}
			outSigma << std::endl;

			sigma_average	  = nb*sigma.avgSigma(y1, y2, Pt_save[iPt], deltaVec[iDt], 0, 0);
			sigmaCos_average = nb*sigma.avgSigma(y1, y2, Pt_save[iPt], deltaVec[iDt], 1, 0);
			outSigmaAvg << deltaVec[iDt]  << " "<< Pt_save[iPt] << " "
									<< sigma_average  << " "
									<< sigmaCos_average << " " << sigmaCos_average/sigma_average
									<< std::endl;

			double sigma_averageABC	  = nb*sigma.avgSigma(y1, y2, Pt_save[iPt], deltaVec[iDt], 0, 1);
			double sigmaCos_averageABC = nb*sigma.avgSigma(y1, y2, Pt_save[iPt], deltaVec[iDt], 1, 1);
			outSigmaAvgABC << deltaVec[iDt]  << " "<< Pt_save[iPt] << " "
										 << sigma_averageABC  << " " <<  sigmaCos_averageABC << " "
										 << sigmaCos_averageABC/sigma_average << std::endl;

			double sigma_averageAB	  = nb*sigma.avgSigma(y1, y2, Pt_save[iPt], deltaVec[iDt], 0, 2);
			double sigmaCos_averageAB = nb*sigma.avgSigma(y1, y2, Pt_save[iPt], deltaVec[iDt], 1, 2);
			outSigmaAvgAB << deltaVec[iDt]  << " "<< Pt_save[iPt] << " "
									  << sigma_averageAB  << " "<<  sigmaCos_averageAB << " "
										<< sigmaCos_averageAB/sigma_average << std::endl;

			double sigma_averageA	  = nb*sigma.avgSigma(y1, y2, Pt_save[iPt], deltaVec[iDt], 0, 3);
			double sigmaCos_averageA = nb*sigma.avgSigma(y1, y2, Pt_save[iPt], deltaVec[iDt], 1, 3);
			outSigmaAvgA << deltaVec[iDt]  << " "<< Pt_save[iPt] << " "
									 << sigma_averageA  << " "<<  sigmaCos_averageA << " "
									 << sigmaCos_averageA/sigma_average << std::endl;

			double sigma_averageCD	  = nb*sigma.avgSigma(y1, y2, Pt_save[iPt], deltaVec[iDt], 0, 4);
			double sigmaCos_averageCD = nb*sigma.avgSigma(y1, y2, Pt_save[iPt], deltaVec[iDt], 1, 4);
			outSigmaAvgCD << deltaVec[iDt]  << " "<< Pt_save[iPt] << " "
										<< sigma_averageCD  << " "<<  sigmaCos_averageCD << " "
										<< sigmaCos_averageCD/sigma_average << std::endl;

		}
		outA.close();
		outB.close();
		outSigma.close();
		outSigmaAvg.close();
		if(mf!=0){outC.close(); outD.close(); }
	}
}


return 0;
}


double sigma_hadronic(sigma_class sigma_, const LHAPDF::PDF* ff_pi_, std::vector<double> ptVec_, double y1_, double y2_,
											double x1_, double  x2_, double PtH_, double DtH_, double PhiPD_, double dt_, double P1_, double  P2_)
{

	double k1t = P1_/x1_;
	double k2t = P2_/x2_;
	double phi21 = acos((-PtH_*PtH_ + DtH_*DtH_/4.)/(P1_*P2_) );
	double pt = sqrt(k1t*k1t + k2t*k2t -2.*k1t*k2t*cos(phi21))/2.;
	double cosPD=(k1t*k1t-k2t*k2t)/(2.*pt*dt_); //this can be defined purely from hadronic + momentum fraction vars
																							//, so it can be constrained to be IN phase space from beginning!
	double phiPD = acos( cosPD );
	double Dirac_factor, ff1, ff2;
	double res= 0;

	if(pt<ptVec_[0] || pt>ptVec_[ptVec_.size()-1] || cosPD<-1 || cosPD>1){}
	else{

		for(int iFlavour=1; iFlavour<=3; iFlavour++)
		{
			ff1 =  ff_pi_->xfxQ2( iFlavour, x1_, k1t*k1t);
			ff2 =  ff_pi_->xfxQ2(-iFlavour, x2_, k2t*k2t);
			Dirac_factor=std::fabs( (P1_*P1_/x1_ -(PtH_*PtH_ - DtH_*DtH_/4.)/x2_)/
															(x1_*x1_*sqrt( k1t*k1t + k2t*k2t - 2.*(PtH_*PtH_ - DtH_*DtH_/4.)/(x1_*x1_)))  );

	 res+= sigma_.getSigma_alternative(y1_, y2_, k1t, k2t, dt_, pt, phiPD, iFlavour, 0)
									*ff1*ff2/(x1_*x1_*x2_*x2_*Dirac_factor*dt_);
		}
	}
		return res;
}
