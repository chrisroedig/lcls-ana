#include <math.h>
#include "ebeam.hh"

double getFelPh(double ebenergy,double pkcurr){

	
	
	double f=9.496;					//factor for formula
	double kparm   =3.5;			//[unitless] undulator K param
	double pund    =0.03;			//[m] undulator period
	double efel    =ebenergy/1000;		//[GeV] energy of ebeam in undulator (needs ipk correction)
	
	
	double ephoton=f*efel*efel/((1.0+kparm*kparm/2.0)*pund);//jackson p. 692
	return ephoton;
}
