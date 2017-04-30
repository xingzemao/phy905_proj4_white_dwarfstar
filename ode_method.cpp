#include<iostream>
#include<math.h>
#include"ode_method.h"
using namespace std;


void RK4_class::operator()(const double &r,const double &h,  const double &rho_old, const double &m_old, double &rho_new, double &m_new)
{
	//first point
	double gamma=gamma_f(rho_old);
	double k1_m = h* m_function( rho_old, r);
	double k1_rho = h*rho_function( -m_old/gamma, r ,rho_old );

	double m1 = m_old+ k1_m/2.0;
	double rho1 = rho_old + k1_rho/2.0;
	//first mid-point
	gamma=gamma_f(rho1);
	double k2_m = h* m_function(rho1,r+0.5*h);
	double k2_rho =  h*rho_function(-m1/gamma,r+0.5*h, rho1);

	double m2 = m_old+ k2_m/2.0;
	double rho2 = rho_old+ k2_rho/2.0;
	//second mid-point
	gamma=gamma_f(rho2);
	double k3_m = h*m_function(rho2, r+0.5*h);
	double k3_rho = h* rho_function(-m2/gamma, r+ 0.5*h, rho2);

	double m3 = m_old+ k3_m;
	double rho3 = rho_old+ k3_rho;
	//last point
	gamma=gamma_f(rho3);
	double k4_m = h*m_function(rho3, r+h);
	double k4_rho = h*rho_function(-m3/gamma, r+h, rho3);
	m_new = m_old + (k1_m + 2.0*k2_m + 2.0*k3_m + k4_m)/6.0;
	rho_new =  rho_old +(k1_rho + 2.0*k2_rho + 2.0*k3_rho + k4_rho)/6.0;
}

double  RK4_class::gamma_f(const double &rho)const
{
	double x= pow(rho,1.0/3.0);
	return x*x/(3*pow(1.0+x*x,0.5));
}


double RK4_class::rho_function(const double &c, const double &r, const double &rho) const
{	const double f = c*rho/(r*r);
	if (r==0){  return 0;}
	else{	return f;}
}


double RK4_class::m_function( const double &c, const double &r) const
{
	const double f = c*r*r;
	return f;
}

double RK4_class::epsilon_f(const double &r)const
{
	return 3.0*(r*(1+2*r*r)*pow(1+r*r,0.5)-log(r+pow(1+r*r,0.5)))/8.0;
}


void RK2_class::operator()( const double&r, const double&h,const double&rho_old,const double&m_old,double&rho_new,double&m_new)
{
	double gamma= gamma_f(rho_old);
	double k1_m= h*m_function(rho_old,r);
	double k1_rho= h*rho_function(-m_old/gamma,rho_old,r);


	gamma = gamma_f(rho_old+k1_rho*0.5);
	double k2_m= h*m_function(rho_old+ k1_rho*0.5, r+0.5*h);
	double k2_rho= h*rho_function(-(m_old+k1_m*0.5)/gamma, rho_old+ 0.5*k1_rho, r+0.5*h);

	rho_new= rho_old+ k2_rho;
	m_new= m_old+k2_m;
}

double RK2_class::rho_function(const double&c, const double &rho, const double&r) const
{
	if(r==0){return 0;}
	else {return c*rho/(r*r);}
}

double RK2_class::m_function(const double &rho, const double &r)const
{
	return rho*r*r;
}

double RK2_class::gamma_f(const double &rho)const
{
	double x=pow(rho,1.0/3.0);
	return x*x/(3*pow(1+x*x,0.5));
}
