#include <TROOT.h>
#include <TH1F.h>
#include <TProfile.h>
#include <math.h>
#include "etof.hh"

double* transData;
int transDataCount;
double maxTransE;
int initEtofTrans(char* name){
	printf("\nREADING ETOF TRANSMISSION FILE...");
	FILE * scaleFile;
	scaleFile=fopen("indata/etofTrans.dat","r");
	double transRaw[10000];
	float emax=0.0;
	int i=0;int count=0;
	for(i=0;i<10000;i++){
		float t=0.0;float e=0.0;
		if(fscanf(scaleFile,"%f %f",&e,&t)!=EOF){
			count++;emax=e;
		}		
		transRaw[i]=(double)t;
	}
	maxTransE=(double)emax;
	transData=&transRaw[0];
	transDataCount=count;
	printf("DONE found %d values",count);		
	return 0;
}
void etofTrans(double* energies, double* counts,double* countsOut, int bins, double vret ){
	int i=0;
	for(i=0;i<bins;i++){
		int transBin=(int)round(((energies[i]-vret)/vret)*((double)transDataCount/maxTransE));
		if((transBin<0)||(transBin>transDataCount)){
			countsOut[i]=0.0;	
		}else{
			countsOut[i]=(transData[transBin]>0.10)?(counts[i]/transData[transBin]):0.0;
		}
	}	
}

void fillCF(double* t, double* v, unsigned numSamples, float baseline,
           float thresh, TH1* hist) {
  // find the boundaries where the pulse crosses the threshold
  double   peak=0.0;

  double detTime=0.0;
  unsigned start  =0;
  bool     crossed=false;
  bool     rising = thresh > baseline;
  for(unsigned k=0; k<numSamples; k++) {
    double y = v[k];
    bool over = 
      ( rising && y>thresh) ||
      (!rising && y<thresh);
    if (!crossed && over) {
      crossed = true;
      start   = k;
      peak    = y;
    }
    else if (crossed && !over) {
      //  find the edge
      double edge_v = 0.5*(peak+baseline);
      unsigned i=start;
      if (rising) { // leading edge +
    while(v[i] < edge_v)
      i++;
      }
      else {        // leading edge -
    while(v[i] > edge_v)
      i++;
      }
      detTime=(t[i]+t[i-1])/2;
      hist->Fill(detTime,1.0);
      crossed = false;
    }
    else if (( rising && y>peak) ||
      (!rising && y<peak))
      peak = y;
  }
}

double avgCfWindow(
			double* t, double* v, unsigned numSamples,
			float baseline, float thresh,
			double tp,double vret,
			double eMin, double eMax){
  // find the boundaries where the pulse crosses the threshold
  double   peak=0.0;
  double detTime=0.0;
  double detE=0.0;
  unsigned start  =0;
  bool     crossed=false;
  bool     rising = thresh > baseline;
  
  double eAvg=0.0;
  double eCount=0.0;
  
  for(unsigned k=0; k<numSamples; k++) {
    double y = v[k];
    bool over = 
      ( rising && y>thresh) ||
      (!rising && y<thresh);
    if (!crossed && over) {
      crossed = true;
      start   = k;
      peak    = y;
    }
    else if (crossed && !over) {
      //  find the edge
      double edge_v = 0.5*(peak+baseline);
      unsigned i=start;
      if (rising) { // leading edge +
    while(v[i] < edge_v)
      i++;
      }
      else {        // leading edge -
    while(v[i] > edge_v)
      i++;
      }
      detTime=(t[i]+t[i-1])/2;
      detE=t2E( detTime, tp, vret );
      if((detE>eMin)&&(detE<eMax)){
        double weight=detE/pow(detE,1.5);//jacobian correction

		  eAvg+=detE*weight;
		  eCount+=1.0*weight;
		  }
      crossed = false;
    }
    else if (( rising && y>peak) ||
      (!rising && y<peak))
      peak = y;
  }
	if(eCount>0.0){
		return eAvg/eCount;
	}
	return 0.0;
	}
int countCfEnWindow(
			double* t, double* v, unsigned numSamples,
			float baseline, float thresh,
			double tp,double vret,
			double eMin, double eMax){
	 // find the boundaries where the pulse crosses the threshold
  double   peak=0.0;
  double detTime=0.0;
  double detE=0.0;
  unsigned start  =0;
  bool     crossed=false;
  bool     rising = thresh > baseline;
  
  int eCount=0;
  
  for(unsigned k=0; k<numSamples; k++) {
    double y = v[k];
    bool over = 
      ( rising && y>thresh) ||
      (!rising && y<thresh);
    if (!crossed && over) {
      crossed = true;
      start   = k;
      peak    = y;
    }
    else if (crossed && !over) {
      //  find the edge
      double edge_v = 0.5*(peak+baseline);
      unsigned i=start;
      if (rising) { // leading edge +
		while(v[i] < edge_v)
		i++;
      }
      else {        // leading edge -
		while(v[i] > edge_v)
		i++;
	  }
	  
      detTime=(t[i]+t[i-1])/2;
      printf(">>>>>>%e",vret);
      detE=t2E( detTime, tp, vret );
      if((detE>eMin)&&(detE<eMax)){eCount+=1;}
      crossed = false;
    }
    else if (( rising && y>peak) ||
      (!rising && y<peak))
      peak = y;
  }
	
	return eCount;
}
int countCfWindow(
			double* t, double* v, unsigned numSamples,
			float baseline, float thresh,
			double tMin, double tMax){
	 // find the boundaries where the pulse crosses the threshold
  double   peak=0.0;
  double detTime=0.0;
  unsigned start  =0;
  bool     crossed=false;
  bool     rising = thresh > baseline;
  
  int eCount=0;
  
  for(unsigned k=0; k<numSamples; k++) {
    double y = v[k];
    bool over = 
      ( rising && y>thresh) ||
      (!rising && y<thresh);
    if (!crossed && over) {
      crossed = true;
      start   = k;
      peak    = y;
    }
    else if (crossed && !over) {
      //  find the edge
      double edge_v = 0.5*(peak+baseline);
      unsigned i=start;
      if (rising) { // leading edge +
		while(v[i] < edge_v)
		i++;
      }
      else {        // leading edge -
		while(v[i] > edge_v)
		i++;
	  }
	  
      detTime=(t[i]+t[i-1])/2;


      if((detTime>tMin)&&(detTime<tMax)){eCount+=1;}
      crossed = false;
    }
    else if (( rising && y>peak) ||
      (!rising && y<peak))
      peak = y;
  }
	
	return eCount;
}
void fillCal(double* t, double* v, unsigned numSamples, float baseline,
           float thresh,double tp,double vret,double eCorrect, TProfile* profile,double roiStart, double roiEnd){
		   return fillCal(t,v,numSamples,baseline,thresh,tp,vret,eCorrect,1.0,profile,roiStart,roiEnd);
		   }
void fillCal(double* t, double* v, unsigned numSamples, float baseline,
           float thresh,double tp,double vret,double eCorrect,double weight, TProfile* profile,double roiStart, double roiEnd){
	unsigned i =0;			//iterator
	double en=0.0;		//calculated energy
	double amp=0.0;		//corrected amplitude
	bool above=false;	//noise gate state
	//for each sample
	for(i=0;i<numSamples;i++){
		
		//figure out gate state from thresh, base and singal
		if(baseline>thresh){above=(v[i]<thresh);}	//negative
		else{above=(v[i]>=thresh);}					//positive
		
		if(above){						//if above noise
			//get energy for this time then correct by eCorrect
			en=t2E( t[i], tp, vret )-eCorrect;
			amp=fabs(v[i])/pow(en,1.5);	//jacobian correction
			profile->Fill(en,amp*weight);		//fill into profile
		}
	}//end for
}
void fill(double* t, double* v, unsigned numSamples, float baseline,
		float thresh,TProfile* profile){
	unsigned i =0;			//iterator
	bool above=false;	//noise gate state
	//for each sample
	for(i=0;i<numSamples;i++){
		//figure out gate state from thresh, base and singal
		if(baseline>thresh){above=(v[i]<thresh);}	//negative
		else{above=(v[i]>=thresh);}					//positive	
		if(above){						//if above noise	
			profile->Fill(t[i],v[i]);		//fill into profile
		}else{							//below noise
			profile->Fill(t[i],0.0);		//fill zero profile
		}
	}//end for
		   
		   }
void calStats(
	double* t,double* v, unsigned numSamples, 
	double baseline, double thresh, 
	double emin, double emax,
	double tp, double vret,
	double* eMean,double *eSum)
	{
		unsigned i=0;
		bool above=false;
		double en=0.0;
		double amp=0.0;
		double sum=0.0;
		double weights=0.0;
		for(i=0;i<numSamples;i++){//for each sample
			
			//figure out gate state from thresh, base and singal
			if(baseline>thresh){above=(v[i]<thresh);}	//negative
			else{above=(v[i]>=thresh);}					//positive
			
			if(above){
				
				en=t2E( t[i], tp, vret );//get the energy
				//printf("tp:%e----vret:%e\n",tp,vret);//debug
				//printf("t:%e---v:%e--e:%e\n",t[i],v[i],en);//debug
				amp=fabs(v[i])/pow(en,1.5);	//jacobian correction
				sum+=en*amp;//add weighted energy to sum 
				weights+=amp;//add amplitude to sum
				}
			
			}
		*eSum=weights;
		if(weights>0){
			*eMean=sum/weights;
		}else{
			*eMean=0.0;
		}
		return; 
		
	}

int cfCalStats(
			double* t, double* v, unsigned numSamples,
			float baseline, float thresh,
			double tp,double vret,
			double eMin, double eMax,
			double* eMean, double* eStdev,
			int* counts, double* eDetMin, double* eDetMax
			)
	{
	// find the boundaries where the pulse crosses the threshold
  double   peak=0.0;
  double detTime=0.0;
  double detE=0.0;
  unsigned start  =0;
  bool     crossed=false;
  bool     rising = thresh > baseline;

  int eCount=0;
  double hits[100];
  
  
  for(unsigned k=0; k<numSamples; k++) {
	double y = v[k];
    bool over = (rising && y>thresh)||(!rising && y<thresh);
    if (!crossed && over) {
		crossed = true;
		start   = k;
		peak    = y;
    }else if (crossed && !over) {
		//  find the edge
		double edge_v = 0.5*(peak+baseline);
		unsigned i=start;
		if (rising) { // leading edge +
			while(v[i] < edge_v){i++;}
		}
		else {        // leading edge -
			while(v[i] > edge_v){i++;}
		}
		
		detTime=(t[i]+t[i-1])/2;//time of threshold crossing
		//printf("\n>>>>>> t:%e\n",detTime);
		detE=t2E( detTime, tp, vret );//energy at threshold crossing
		//printf(">>>>>>min:%e max:%e det:%e",eMin,eMax,detE);
		
		if((detE>eMin)&&(detE<eMax)){//check if energy in desired range
	
	
	
			if(eCount<100){//chekc if within first 100 detections
				hits[eCount]=detE;
				eCount++;
			}
	
		}else{
		  
		}
		crossed = false;
    }else if (( rising && y>peak) ||
      (!rising && y<peak))
      peak = y;
	}//end detection loop
	*counts=eCount;
	if(eCount<=0){return 0;}
	
	//BEGIN STATS
	int i=0;
	double detEMax=0.0;
	double detEMin=0.0;
	double eAvg=0.0;
	
	for(i=0;i<eCount;i++){// 1st loop: mean, min/max
		if(i==0){//first
			detEMax=hits[i];detEMin=hits[i];
		}else{//others
			detEMax=hits[i]>detEMax?hits[i]:detEMax;
			detEMin=hits[i]<detEMin?hits[i]:detEMin;
		}
		eAvg+=hits[i];
	}
	eAvg=eAvg/eCount;
	*eMean=eAvg;
	*eDetMax=detEMax;
	*eDetMin=detEMin;
	
	double stDevSum=0.0;
	for(i=0;i<eCount;i++){// 2nd loop: stdev 
		stDevSum+=pow(hits[i]-eAvg,2);
	}
	*eStdev=sqrt(stDevSum/eCount);
	
	
	return 0;
	}

double getEtofPrompt(double* t, double* v, unsigned numSamples){
	 // find the boundaries where the pulse crosses the threshold
bool found=false;
  double   peak=0.0;

  double detTime=0.0;
  unsigned start  =0;
  bool     crossed=false;
  bool     rising = false;
  double 	baseline=-0.000;
  double	thresh=-0.005;
  double tmin=6.25e-7;
  double tmax=6.50e-7;
  for(unsigned k=0; k<numSamples; k++) {
    double y = v[k];
    bool over = 
      ( rising && y>thresh) ||
      (!rising && y<thresh);
    if (!crossed && over) {
      crossed = true;
      start   = k;
      peak    = y;
    }
    else if (crossed && !over) {
      //  find the edge
      double edge_v = 0.5*(peak+baseline);
      unsigned i=start;
      if (rising) { // leading edge +
    while(v[i] < edge_v)
      i++;
      }
      else {        // leading edge -
    while(v[i] > edge_v)
      i++;
      }
      detTime=(t[i]+t[i-1])/2;
      if((detTime>tmin)&&(detTime<tmax)){
		if(!found){
			return detTime;
			found=true;
		}
	  }
      crossed = false;
    }
    else if (( rising && y>peak) ||
      (!rising && y<peak))
      peak = y;
	
  }
  return 0.0;

}
void fillCalCF(double* t, double* v, unsigned numSamples, float baseline,
float thresh,double tp,double vret,double eCorrect, TH1* hist){
	return fillCalCF(t,v,numSamples, baseline, thresh,tp,vret,eCorrect,1.0,  hist);
	}
void fillCalCF(double* t, double* v, unsigned numSamples, float baseline,
           float thresh,double tp,double vret,double eCorrect,double weight, TH1* hist) {
  // find the boundaries where the pulse crosses the threshold
  double   peak=0.0;
	
  double detTime=0.0;
  double detE=0.0;
  unsigned start  =0;
  bool     crossed=false;
  bool     rising = thresh > baseline;
              
  for(unsigned k=0; k<numSamples; k++) {
    double y = v[k];
    bool over = 
      ( rising && y>thresh) ||
      (!rising && y<thresh);
    if (!crossed && over) {
      crossed = true;
      start   = k;
      peak    = y;
    }
    else if (crossed && !over) {
      //  find the edge
      double edge_v = 0.5*(peak+baseline);
      unsigned i=start;
      if (rising) { // leading edge +
    while(v[i] < edge_v)
      i++;
      }
      else {        // leading edge -
    while(v[i] > edge_v)
      i++;
      }
      detTime=(t[i]+t[i-1])/2;


      detE=t2E( detTime, tp, vret );//get energy from time
      
	  //printf("---------E:%e\n",detE+eCorrect);
      hist->Fill((detE+eCorrect),weight);
      crossed = false;
    }
    else if (( rising && y>peak) ||
      (!rising && y<peak))
      peak = y;
  }
}

// *********************************
// *     BERTOLD'S CODE            *
// *********************************
//
//
//  T2E - time of flight to kinetic energy conversion for LCLS-AMO ETOFs 
//
//  Version 0.9:  2010-AUG-06
//  Author: Bertold Kraessig
//
//  This code calculates the kinetic energy of an on-axis electron
//  observed at flight time t(sec) relative to the prompt peak for the nominal 
//  dimensions of the AMO ETOFs and a given retardation voltage (in Volts). The 
//  parametrization is based on SIMION simulations of the ETOFs. The case
//  of zero retardation is calculated analytically, with an additional two-
//  parameter correction term to better match the SIMION results. The 
//  difference of the retarded versus the unretarded cases is parametrized
//  with four parameters. The parametrization was based on simultations for
//  1V to 4000V retardations and 1eV to 4000eV excess kinetic energies over
//  the retardation potential. The agreement between the parametrization and
//  the simulations is better than 1% except at kinetic energies of less
//  than 10 eV above the retardation potential. Deviations from the
//  simulations are to be expected when the elecron starting position is
//  displaced from the nominal source point of the ETOF. Also, electron
//  trajectories for a given kinetic energy with non-zero starting angles 
//  within the acceptance angle will have somewhat different flight
//  times.
//  
//
//  
//  
//  ************************************************************************
//  *                                                                      *
//  *   PLEASE NOTE: I claim ownership for the parametrization used in     *
//  *   this code.                                                         *
//  *   You are welcome to use this parametrization in your data analysis. *
//  *   I request to be acknowleged as the source for the tof-to-energy    *
//  *   conversions in any published or otherwise disseminated work        *
//  *   that uses this parametrization and in which I am not an author.    *
//  *   I also request, if you should pass this code or any transcriptions *
//  *   of this code on to others, that you also communicate my request    * 
//  *   and include this note in the code.                                 *
//  *                                                                      *
//  *   For questions send e-mail to: kraessig@anl.gov                     *
//  *                                                                      *
//  *   Thank you.                                                         *
//  *   Bertold Kraessig                                                   *
//  *                                                                      *
//  ************************************************************************
//
//
#include <stdio.h>
#include <math.h>

//
// calculates the first root of a monic cubic polynomial with coefficients a
//

double cuberoot ( double* a ) 
{
//  const double Pi = 3.141592653589793;
  double Q;
  double R;
  double S;
  double T;
  double sqmQ;
  double Rpsqdis;
  double Rmsqdis;
  double dis;
  double th;
  double v1;

  Q = ( 3*a[1] - pow( a[2], 2) )/9.;
  R = ( 9*a[2]*a[1] -27*a[0] - 2*pow( a[2], 3 ) )/54.;
  dis = pow( Q,3 ) + pow( R,2 );

  if ( dis < 0. ) 
  {
    sqmQ = sqrt(-Q);
    th = acos( R / pow( sqmQ,3 ) );
    v1 = 2. * sqmQ * cos( th/3. ) - a[2]/3.;
//    v2 = 2. * sqmQ * cos( (th + 2*Pi)/3. ) - a[2]/3;
//    v3 = 2. * sqmQ * cos( (th + 4*Pi)/3. ) - a[2]/3;
   }
   else 
   {
     Rpsqdis = R + sqrt(dis);
     Rmsqdis = R - sqrt(dis);
	 S = pow( Rpsqdis, 1./3. );
     T = pow( Rmsqdis, 1./3. );
     v1 = -a[2]/3. + (S+T);
//   v2 and v3 are imaginary
//
   }

       return (v1) ;
       
}

// calculates kinetic energy for a given flight time 
// for 0V retardation, t in nanseconds, E in eV, lengths in mm
//

double Eoft( double t ) 
{ 
//
 const double em=0.1758820150;    // e/m
 const double c=299.792458;       // speed of light
 const double mc2=510998.91;      // electron mass in eV
 const double U=2000.;             // acceleration voltage detector
 const double L=598.;              // eTOF length in mm
 const double L2=12.;              // acceleration length
 
 double a[3];
 double Fm;
 double v;
 double gamma;
 double E;

 Fm = em * U/L2;                      // acceleration force divided by mass
 
 a[1] = -(L-L2)*Fm;
 a[2] = Fm/2.*t-L/t;
 a[0] = pow( L-L2 , 2 ) * Fm/2./t;
 v = cuberoot(a);
//
// relativistic correction
//
 gamma=1./sqrt( 1. - pow(v/c,2) );
 E = ( gamma - 1. ) * mc2;
// printf("v=%f  gamma=%f\n",v,gamma);
//printf("Eoft before  E=%f eV\n",E);
// empirical correction to match SIMION: E=E+509.536*t^-2.6185
 E = E + 509.536 * pow(t,-2.6185);
 
//printf("Eoft after E=%f eV\n",E);

 return (E);
}

// The difference between the Energy at time t without retardation
// and the Energy at the same time t with retardation
// there is a single difference curve for the scaled quantities
// T=sqrt(ret)*t and D=d/ret
//

 double Difference( double t, double ret )
 {
  double a[2];
  double b[2];
  double s;
  double u;
  double v;
  double x;
  double D;

  if ( ret != 0. ) 
  {
      a[0] = -10.5738; 
      a[1] = -1.7433;
      b[0] = 2.03875;  
      b[1] = -0.00398517;
      s = 1./sqrt(ret)/t;
      x = log(s);
      u = exp(a[0]+a[1]*x)+exp(b[0]+b[1]*x);
      v = -log(u);
      D = exp(v)*ret;
  }
  else 
  { 
    D=0.;
  }
// printf("Difference D=%f eV\n",D);
  return (D);
 }
//
// conversion routine t2E( tof, tp, ret )
// tof in seconds, prompt time tp in seconds, ret in Volts
//
double t2E ( double tof, double tp, double ret )
{
 const double delta=1.993E-9;  // prompt photon travel time
 double t;
 double Ekin;

 if ( tp == 0. ) tp = delta;
 t = ( tof - tp + delta )*1.e9;                  // in nsec
 Ekin = Eoft( t ) + ret - Difference( t, ret );

 return ( Ekin );
}	
