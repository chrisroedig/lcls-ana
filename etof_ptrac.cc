/* $Id: myana.cc,v 1.14 2010/07/22 22:26:57 caf Exp $ */
#include <TROOT.h>
#include <TH1F.h>
#include <TProfile.h>
#include <math.h>

//headers for analysis framework
#include "myana.hh"
#include "main.hh"

//headers for things chris wrote
//etof analysis library
#include "etof.hh"
//ebeam library, calculates expected photon energy
#include "ebeam.hh"
//root dump, writes easy to use ASCII formatted data
#include "rootDump.hh"


//misc Globals
static int n_event=0;
static int good_event=0;

//ETOF DAQ GLOBAL VARS
static int     numChannelsETof = 0; //should be 5 :)
static int     numSamplesETof = 0;
static double  sampleIntervalETof = 0;//total acquisition interval        

//fixed params for each etof
static double tPrompt[5];
static double vRet[5];


//etof energy bins
const int eBinCount	=9000;
const double eBinMin=0.0;
const double eBinMax=900.0;

//streaking bins
const int 	streakBins	=100;
const double 	streakMin	=-10.0;
const double 	streakMax	=10.0;

//timing vars
//this may be irrelevant depending on the method of timing
static double motor1=0.0;      //last recovered valid motor position
static double motorCtr=42.960; //approximate center of most scans

//sorted shit
static TH1F* 		mapHist0[streakBins];
static TH1F* 		mapHist1[streakBins];
static TH1F* 		mapHist2[streakBins];
static TH1F*	 	mapHist3[streakBins];
static TH1F*            mapHist4[streakBins];

//streak stats 
static TH1F*		streakHist;

//full time profile
static TProfile* timeProf;

// This function is called once at the beginning of the analysis job,
// You can ask for detector "configuration" information here.

void beginjob() 
{
	int i=0; //iterator in for loops
	int fail = 0;//error flag
	char name[32]; //placeholder for a name
  
	//grabs the intrument configuration
	fail = getAcqConfig( AmoETof, numChannelsETof, numSamplesETof, sampleIntervalETof);    
	
	if ( fail != 0 )
    		printf( "-----ERROR-----> COULD NOT GRAB ETOF CONFIG!!!!!!");
	
	printf("numsmaples:%d\n",numSamplesETof);
	printf("sampint:%e\n",sampleIntervalETof);
	
	//ENTER THE CORRECT PROMPT TIMES HERE
	//numbering accounts for switch of #3/#4 acq lines for this run 	
	tPrompt[0]=6.42953e-7;  //ETOF#1
	tPrompt[1]=6.41819e-7;  //ETOF#2
	tPrompt[2]=6.40565e-7;  //ETOF#3
	tPrompt[3]=6.40565e-7;  //ETOF#4
	tPrompt[4]=6.42382e-7;  //ETOF#5
	
	//ENTER THE CORRECT RETARDING HERE
	vRet[0]=700.0;			//ETOF#1
	vRet[1]=700.0;			//ETOF#2
	vRet[2]=0.0;			//ETOF#3
	vRet[3]=700.0;			//ETOF#4
	vRet[4]=700.0;			//ETOF#5
	
	
	
	/*
	* this looop iterates over the streaking bins and sets up 
	* histogram instances for each etof
	*/
	for(i=0;i<streakBins;i++){
		/*
		 	sprintf(name,"etofP%dbin%d",i,j);
			etofP[i][j]=new TH1F(name,name,600,700.0,900.0); 
		 */ 
		
		sprintf(name,"histBin0_%d",i); 
		mapHist0[i]= new TH1F(name,name,eBinCount,eBinMin,eBinMax);
               
		sprintf(name,"histBin1_%d",i);
                mapHist1[i]= new TH1F(name,name,eBinCount,eBinMin,eBinMax);
		
		sprintf(name,"histBin2_%d",i); 
		mapHist2[i]= new TH1F(name,name,eBinCount,eBinMin,eBinMax);
		
		sprintf(name,"histBin3_%d",i); 
		mapHist3[i]= new TH1F(name,name,eBinCount,eBinMin,eBinMax);
		
		sprintf(name,"histBin4_%d",i); 
		mapHist4[i]= new TH1F(name,name,eBinCount,eBinMin,eBinMax);
				
	}
	//this is a histogram that tracks occurrance of streaking values
	sprintf(name,"streakHist"); 
	streakHist= new TH1F(name,name,streakBins,streakMin,streakMax);
	
	//half wdith of a bin, needed when setting up profile (trace averaging instances)
	double hbs = sampleIntervalETof/numSamplesETof/2.0;
	//
	sprintf(name,"timeProf"); 
	timeProf= new TProfile(name,name,numSamplesETof,0.0-hbs,sampleIntervalETof-hbs);
	
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
	double laser_delay=0.0;
	double speed_of_light = 299792458.0;	
	
	int i;
	int fail = 0;

	//nice progress indication....	
	if(n_event%10000==0){
		//print this update every 10000 events
		printf("\n[ %d events, %d good ] \n",n_event,good_event);
	}

	n_event++;
	//print a bar every 1000 events
	if(n_event%1000==0){printf("|");}
	//print a dot every 100 events
	else if(n_event%100==0){printf(".");}
	fflush(stdout); 

        /*
        * TIMING/DELAY INFO
        * This area may need to change to be up to date with the current method of tracking timing delay
        */

	//get phase cavity timing data 
	//tracks electron bunch arrival against master synch clock
	double pCT1, pCT2, phaseCavityCharge1, phaseCavityCharge2;
	fail = getPhaseCavity( pCT1, pCT2, phaseCavityCharge1, phaseCavityCharge2 );
	if(fail){printf("x");fflush(stdout);return;}//on fail return

	//get Laser delay timing data
	float value;
	fail=0;	
	
	 // THIS BLOCK USES THE DELAY STAGE
	fail = getPvFloat("AMO:LAS:DLS:05:MTR.RBV",value); 
	if(!fail){
		motor1=value;
	}
	laser_delay = ((motor1-motorCtr)*0.001*2/speed_of_light*1e12)+pCT2-.170;	
	
	if(abs(laser_delay)>0.50){return;}//rejecty if more than 50fs delay
	/*
	 * EBEAM DATA (commonly used for photon energy)
	 * Get data from EBeam (i.e. accelerator)
	 */
	double ebcharge; double ebenergy; double posx; double posy;
	double angx; double angy; double pkcurr;
	fail =  getEBeam(ebcharge, ebenergy, posx, posy, angx, angy, pkcurr);
	//getFelPh() is defined in ebeam.cc, go there to adjust K-param and udulator period (shouldn't have changed)
	double eph=getFelPh(ebenergy,pkcurr);
	
    	//ADJUST THE OPERATION PH ENERGY HERE (it's set to 950)
	// eCorrect is the diff between operating photon energy and actual photon energy for this shot 
	double eCorrect = -1.0*(eph-950.0);//displacement from central value 
	
	//get etof data..
	double* timeETof[5];
	double* voltageETof[5];
	
	//the following loop extracts the raw detector trace from each etof 
    
    	for(i=0;i<5;i++){
		fail = getAcqValue( AmoETof, i, timeETof[i], voltageETof[i]);  			
			if(fail){return;}//on fail return
	}
	
	/*
	* STREAK TRACKING
	* here we pick a spectrometer and attempt to figure out by how much 
	* the primary photoelectron has been streaked
	* this completely depends on where you are measuring the streaking
	* the example here extracts the value in counting mode from one of the etofs
	* odds are you'll be looking at an etof channel but in current mode
	*/

	//these are vars that will hold stats 
	double eMean;double eSum;double eStd;int counts;double eMin; double eMax;
	
	/*
	* this passes the time and amp arrays into a function 
	* that performs CFD, time-to-energy, and averaging
	* it outputs the mean energy: eMean , energy stdev, eStd
	* number of counts: counts, and min/max detected energies
	* this should be replaced 
	* with whatever it takes to properly track the primary streak
	*/

	cfCalStats(
		timeETof[2],voltageETof[2],numSamplesETof,
		0.0,-0.006,tPrompt[2],vRet[2],65.0,100.0,
		&eMean, &eStd,&counts,&eMin,&eMax);
	eSum=(double)counts;
	
	//this line calculates the amount by which the primary electron was streaked
	// this is basically the calculation:
	// avg_energy - photon_energy_setting + primary_ionization_potential + eCorrect 
	//this was for Ne 1s and 950 eV xrays
	double eStreak=eMean-950+868+eCorrect;
	
	//here we reject the shot if it has bad or insufficient data
	//change these in accordance with your requirements 

	if(eSum<2.0){return;}//atleast 2 electrons	
	if(eStd>1.0){return;}//no more than 1 eV apart

	//check that the streak value falls within our histo range
	if((eStreak<streakMin)||(eStreak>streakMax)){return;}//out of range, return
	
	//using root's API we add the streak value to the histogram
	streakHist->Fill(eStreak,eSum);//fill into streak Hist using intergral of trace
	
	//here we turn the streak value in eV into an array index
	int bin = (int) round( ( ( (double) streakBins ) * ( ( eStreak-streakMin ) / ( streakMax-streakMin ) ) ) );

	//this is important to prevent segfault
	//if the bin is outside of our array don't try to access the array
	if( ( bin<0 ) || ( bin>=streakBins ) ){ return; } 
	

	//here we add the detection counts to the histogram that corresponds to the current streaking value
	//fillCalCF performs CFD, time-to-energy and calls root to add the detected counts to a histogram instance
	//fillCalCF( timeArray, voltageArray, number_of_samples, signal_baseline, CFD_threshold , t_zero, final_retarding, weight, histogram_to_fill) 

	fillCalCF( timeETof[0], voltageETof[0], numSamplesETof, 0.0, -0.006, tPrompt[0], vRet[0], 0.0, 1.0, mapHist0[bin] );
	fillCalCF( timeETof[1], voltageETof[1], numSamplesETof, 0.0, -0.006, tPrompt[1], vRet[1], 0.0, 1.0, mapHist1[bin] );
	fillCalCF( timeETof[2], voltageETof[2], numSamplesETof, 0.0, -0.006, tPrompt[2], vRet[2], 0.0, 1.0, mapHist2[bin] );
	fillCalCF( timeETof[3], voltageETof[3], numSamplesETof, 0.0, -0.006, tPrompt[3], vRet[3], 0.0, 1.0, mapHist3[bin] );
	fillCalCF( timeETof[4], voltageETof[4], numSamplesETof, 0.0, -0.006, tPrompt[4], vRet[4], 0.0, 1.0, mapHist4[bin] );
	
	//done
	good_event++;return;	
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
	/*
	* This runs at the very end, we simply just dump out data here
	*/
	printf("\n[ END JOB %d EVENTS ]\n\n",n_event);

	int i=0;char name[32]; 

	/// here we dump the arrays of streak sorted histgrams as 2d ascii files

	sprintf(name,"outdata/mapHist0.dat");
	dumpThMatrixFile(name,mapHist0,eBinCount,streakBins);
	
        sprintf(name,"outdata/mapHist1.dat");   
        dumpThMatrixFile(name,mapHist1,eBinCount,streakBins);

        sprintf(name,"outdata/mapHist2.dat");   
        dumpThMatrixFile(name,mapHist2,eBinCount,streakBins);

	sprintf(name,"outdata/mapHist3.dat");
	dumpThMatrixFile(name,mapHist3,eBinCount,streakBins);
	
	sprintf(name,"outdata/mapHist4.dat");	
	dumpThMatrixFile(name,mapHist4,eBinCount,streakBins);
	
	//here we dump the streaking statistics
	sprintf(name,"outdata/streakHist.dat");
	dumpThFile(name,streakHist,streakBins);

	
	printf("\n\n BYE BYE NOW!\n");
  
}//END:endjob()
