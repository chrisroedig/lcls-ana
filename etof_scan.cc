/* $Id: myana.cc,v 1.14 2010/07/22 22:26:57 caf Exp $ */
#include <TROOT.h>
#include <TH1F.h>
#include <TProfile.h>
#include <math.h>

#include "myana.hh"
#include "main.hh"
#include "etof.hh"
#include "ebeam.hh"
#include "rootDump.hh"
/*
*
* THIS CODE PRODUCES ETOF SPECTRAL AVERAGES THAT ARE SORTED BY XFEL/IR timing
*
*/

//misc Globals
static int n_event=0;
static bool bad1000=false;

//ETOF DAQ GLOBAL VARS
static int     numChannelsETof = 5; //should be 5 :)
static int     numSamplesETof = 0;
static double  sampleIntervalETof = 0;//total acquisition interval        

//fixed params for each etof
static double tPrompt[5];
static double vRet[5];

//timing gridparameters , units are ps
static double timingStart=-0.50;
static double timingEnd=0.50;
const int timingBins=100;

//static vars to hold timing info

static double motor1=0.0; 		//last recovered valid motor position
static double motorCtr=42.960; //approximate center of most scans

//static double angsft1=0.0;
//static double angsftCtr=1.0653300e+7;

const int energyBins=9000;
const double energyMin=0;
const double energyMax=900;
const int bcWidth=5;

//the is a 2d array of etof histograms (number of histograms)x(number of timing bins)
static TH1F* etofTime[5][timingBins];

//this is a histogram that holds timing statistics
static TH1F* timingHist;


// This function is called once at the beginning of the analysis job,
// You can ask for detector "configuration" information here.

void beginjob() 
{

	int fail = 0;
	int j=0;int i=0;
	char name[32];
	fail = getAcqConfig( AmoETof, numChannelsETof, numSamplesETof, sampleIntervalETof);    
	if ( fail != 0 )
    		printf( "-----ERROR-----> COULD NOT GRAB ETOF CONFIG!!!!!!");
	
	printf("numsamples:%d\n",numSamplesETof);
	printf("sampint:%e\n",sampleIntervalETof);
	
	//correct prompt/t-zero times
	tPrompt[0]=6.42953e-7;  //ETOF#1
	tPrompt[1]=6.41819e-7;  //ETOF#2
	tPrompt[2]=6.40565e-7;  //ETOF#4
	tPrompt[3]=6.40565e-7;  //ETOF#3
	tPrompt[4]=6.42382e-7;  //ETOF#5
	
	//correct vret for each spectrometer
	vRet[0]=700.0;			//ETOF#1
	vRet[1]=700.0;			//ETOF#2
	vRet[2]=0.0;			//ETOF#4
	vRet[3]=700.0;			//ETOF#3
	vRet[4]=700.0;			//ETOF#5
	
	
	for( j=0 ; j<timingBins ; j++ ){//for each timing bin

		for(i=0;i<5;i++){//for each histogram

			//make a name
			sprintf(name,"etof%dtime%d",i,j);
			//set up a ROOT histogram
			etofTime[i][j]=new TH1F(name,name,energyBins,energyMin,energyMax);
		}

	}
	sprintf(name,"timing");
	//set up a histgram for timing stats
	timingHist=new TH1F(name,name,timingBins,timingStart,timingEnd);

}

// This function is called once for each run.  You should check to see
// if detector configuration information has changed.

void beginrun() 
{
  int fail = 0;
  	printf("[ BEGIN RUN ]\n");

  //grab etof config - check for changes sicne begin job  
  int     numChannelsETof2 = 0;
  int     numSamplesETof2 = 0;
  double  sampleIntervalETof2 = 0;    
  fail = getAcqConfig( AmoETof, numChannelsETof2, numSamplesETof2, sampleIntervalETof2);    
  if ( fail != 0 || numChannelsETof2 != numChannelsETof || numSamplesETof2 != numSamplesETof ||
   sampleIntervalETof2 != sampleIntervalETof )
    printf( "-----WARNING---> ETof configuration has been changed between runs!\n" );  


}

void begincalib()
{
}

// This is called once every shot.  You can ask for
// the individual detector shot data here.

void event() 
{
	int i;
	int fail = 0;
	
	
	double laser_delay=0.0;
	double speed_of_light = 299792458.0;	
	
	if(n_event%10000==0){
		printf("\n[ %d events ] \n",n_event);
	}
	if(n_event%1000==0){
		if(bad1000){printf("X");}else{printf("|");}
		bad1000=false;
	}
	else if(n_event%100==0){printf(".");}
	fflush(stdout); 
	
	//get PCAV timing data 
	double pCT1, pCT2, phaseCavityCharge1, phaseCavityCharge2;
	fail = getPhaseCavity( pCT1, pCT2, phaseCavityCharge1, phaseCavityCharge2 );
	if(fail){printf("x");fflush(stdout);return;}//on fail return

	//get Laser delay timing data
	float value;
	fail=0;	
	
	
	// THIS BLCOK USES THE DELAY STAGE	
	//replace this with whatever you currently use to measure FEL/IR timing
	fail = getPvFloat("AMO:LAS:DLS:05:MTR.RBV",value); 
	if(!fail){
		motor1=value;
	}

	laser_delay = ((motor1-motorCtr)*0.001*2/speed_of_light*1e12)+pCT2-.170;	

	// laser_delay should be the FEL / IR timing delay in ps
  	if(n_event%1000==0){
		//every 100 events print the delay
		printf("\n[delay %e]\n",laser_delay);
	}
	//fill the timing value into a histogram
	timingHist->Fill(laser_delay,1.0);

	//turn the measured delay into an array index based on the range and steps declared at the beginning of this file

	int bin = (int) round( (double) timingBins *( laser_delay - timingStart ) / ( timingEnd - timingStart ) ) ;
	
	if(((laser_delay)<timingStart)||((laser_delay)>timingEnd)){
		//delay was outisde the predefined range
		n_event++;
		bad1000=true;
		printf("^");
		return;
		}
	if((bin<0)||(bin>(timingBins-1))){
		//bin index doesn't fit into array
		n_event++;
		printf("*");
		bad1000=true;
		return;
		}
	/*
	*
	* AT THIS POINT:
	* we have measured the timing delay and are ready to grab the etof data
	*
	*/
	//get etof data..
	double* timeETof[5];
	double* voltageETof[5];

	//get the acquird data from each etof
    	for(i=0;i<5;i++){
		fail = getAcqValue( AmoETof, i, timeETof[i], voltageETof[i]);  			
		
		}

	if(fail){printf("E");return;}//on fail return
	
	//fill each etof into the correct timing bin
	for(i=0;i<5;i++){
		//fillCallCF() is defined in etof.cc it does: CFD, time-to-energy conversion and fills into a histogram
		//fillCalCF( time_array, voltage_array, num_of_samples, cfd_baseline, cfd_threshold, t_zero, retarding, energy offset, target histogram )		
		fillCalCF( timeETof[i], voltageETof[i], numSamplesETof, -.001, -.006, tPrompt[i], vRet[i], 0.0, etofTime[i][bin] );
	}
	

	
	n_event++;return;	
}//END EVENT

void endcalib()
{
}

void endrun() 
{
  printf("\n[ END RUN ]\n");
}

void endjob()
{

	printf("\n[ END JOB %d EVENTS ]\n\n",n_event);

	int i=0;int j=0;int k=0;
	char name[32];
	
	
	//set up to use etof transmission file
	sprintf(name,"indata/etofTrans.dat");
	initEtofTrans(name);
	
		
	double eAx[energyBins];//temp array to hold energy axis
	FILE* eaxOut;
	FILE* taxOut;
	FILE* rawOut;
	FILE* corrOut;
	FILE* smoothOut;
	
	
	for(i=0;i<5;i++){
		
		printf("Saving data from etof%d: \n",i);

		
		sprintf(name,"outdata/etof%dTimeScanRaw.dat",i);
		printf("	dumping 2d array to file: ");printf(name);printf("\n");
		rawOut=fopen(name,"w");
		
		sprintf(name,"outdata/etof%dTimeScanCorrected.dat",i);
		printf("	dumping 2d array to file: ");printf(name);printf("\n");
		corrOut=fopen(name,"w");
		
		sprintf(name,"outdata/etof%dTimeScanSmooth%d.dat",i,bcWidth);
		printf("	dumping 2d array to file: ");printf(name);printf("\n");
		smoothOut=fopen(name,"w");
		
		for(j=0;j<timingBins;j++){
			
			double iAx[energyBins];//temp array to hold raw line-out
			double cAx[energyBins];//temp array to hold corrected line out
			double sAx[energyBins];//temp array to hold corrected line out
			
			//fill the raw array from histo		
			arrayTH1F(eAx,iAx,energyBins,etofTime[i][j]);
			
			//correct the line-out for spec transmission
			etofTrans(eAx,iAx,cAx,energyBins,vRet[i]);
			
			//smooth according to desired box car width
			smoothArray(cAx,sAx,energyBins,bcWidth);
			
			//write a line to each file 
			for(k=0;k<energyBins;k++){//write out the values one by one
				fprintf(rawOut," %e ",iAx[k]);	
				fprintf(corrOut," %e ",cAx[k]);	
				fprintf(smoothOut," %e ",sAx[k]);	
			}
			//terminate the lines
			fprintf(rawOut   ,"\n");
			fprintf(corrOut  ,"\n");
			fprintf(smoothOut,"\n");

		}
		fclose(rawOut);
		fclose(corrOut);
		fclose(smoothOut);
					printf("done...\n\n");
	}
	printf("writing energy axis file: ");printf(name);fflush(stdout);
	sprintf(name,"outdata/etofEnergy.dat");
	eaxOut=fopen(name,"w");
	for(i=0;i<energyBins;i++){
		fprintf(eaxOut,"%e ",eAx[i]);
		}
	fclose(eaxOut);
	printf("  done\n\n");fflush(stdout);
	
	
	sprintf(name,"outdata/etofTime.dat");
	printf("writing time axis file: ");printf(name);fflush(stdout);
	taxOut=fopen(name,"w");
	double hbs=(timingEnd-timingStart)/timingBins/2.0;
	for(i=0;i<timingBins;i++){
		fprintf(taxOut,"%e ",(timingStart+hbs+(double)i*2.0*hbs));
		}
	fclose(taxOut);
	printf("  done\n\n");fflush(stdout);
	
	sprintf(name,"outdata/timingHist.dat");
	dumpThFile(name,timingHist,timingBins);
	
	printf("\n\n BYE BYE NOW!\n");
	
}//END:endjob()
