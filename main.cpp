#include<iostream>
#include<math.h>
#include<fstream>
#include"ode_method.h"
using namespace std;

int main()
{
	double Ye=0.5;
//	double rho_old=10;
//	double m_new,rho_new;
	class RK4_class dwarf_star;
//	double rest_kinetic_e=0,grav_e=0;
	double rho_array[8]={-1,0,1,2,3,4,5,6};

	ofstream myfile;
	myfile.open("C_r_m");
	for (int i=0;i<8;i++)
	{
	//	h= pow(10,h_array[i]);
		double rest_kinetic_e=0,grav_e=0;
		double rho_old= pow(10,rho_array[i]);
		double r=0,m_old=0,h=1.0e-6,rho_new=0,m_new=0;
	//	cout<<rho_array[i]<<"centra density is: "<<rho_old<<endl;
		while(rho_old>1.0e-6)
		{
			//obtain the new m and rho by Runge-Kutta 4th order method;
			dwarf_star(r,h,rho_old,m_old, rho_new, m_new);
	//		cout<<"r is "<<(r+h)*Ye*7.72/6.95/100<<"sun radius"<<"  rho is:"<<rho_new<<"  m is:"<<m_new*Ye*Ye*5.67/1.99<<"solar mass"<<endl;
			//update parameters	
			rho_old=rho_new;
			m_old=m_new;
			r+=h;

			//calculate kinetic energy
			double x=pow(rho_new,1.0/3.0);
//			cout<<dwarf_star.epsilon_f(x)<<"   "<<x<<endl;
			rest_kinetic_e+=dwarf_star.epsilon_f(x)*r*r*h;

			//calculate the gravitational energy
			grav_e+= rho_new*m_new*r*h;
	
	//		if(i==5){
		//		cout<<"r is "<<(r+h)*Ye*7.72/6.95/100<<"sun radius"<<"  rho is:"<<rho_new<<"  m is:"<<m_new*Ye*Ye*5.67/1.99<<"solar mass"<<"   "<<m_old<<endl;
	
	//		}
		}
		//add unit and magnitude to unitless values
		rest_kinetic_e*=1.7402*pow(Ye,3.0);
		cout<<"The rest kinetic energy is: "<<rest_kinetic_e<<"10^57 MeV"<<endl;
		grav_e*=1.743157*pow(Ye,3.0);
		cout<<"The gravitational energy is: "<<grav_e<<"10^57 MeV"<<endl;
		cout<<"central density"<<rho_array[i]<<" r is "<<"  "<<(r)*Ye*7.72/6.95/100<<"sun radius"<<"  m is:"<<m_old*Ye*Ye*5.67/1.99<<"solar mass"<<endl;
		myfile<< rho_array[i]<<"r:   "<<r*Ye*7.72/6.95/100<<"m:   "<<m_old*Ye*Ye*5.67/1.99<<endl;		
	}
	myfile.close();
	return 0;
}
