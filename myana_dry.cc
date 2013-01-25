/* $Id: myana.cc,v 1.14 2010/07/22 22:26:57 caf Exp $ */
#include <TROOT.h>
#include <TH1F.h>
#include <TProfile.h>
#include <math.h>

#include "myana.hh"
#include "main.hh"

//misc Globals
static int n_event=0;

//ETOF DAQ GLOBAL VARS
static int     numChannelsETof = 0; //should be 5 :)
static int     numSamplesETof = 0;
static double  sampleIntervalETof = 0;//total acquisition interval        


// This function is called once at the beginning of the analysis job,
// You can ask for detector "configuration" information here.

void beginjob() 
{
	int i=0;
	int fail = 0;
	char name[32]; 
  
	fail = getAcqConfig( AmoETof, numChannelsETof, numSamplesETof, sampleIntervalETof);    
	if ( fail != 0 )
    printf( "-----ERROR-----> COULD NOT GRAB ETOF CONFIG!!!!!!");
	
	
	
	
	
	printf("numsmaples:%d\n",numSamplesETof);
	printf("sampint:%e\n",sampleIntervalETof);

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
	
	
	
	if(n_event%10000==0){
		printf("\n[ %d events, ha ha ha ] \n",n_event);
		
	}
	if(n_event%1000==0){
		
		printf("|");}
	else if(n_event%100==0){printf(".");}
	fflush(stdout); 

	n_event++;

	return;	
	


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
	
	
	
	
	
	



	

	
	
	char outname[20];
	
	
	
	

	
	
	
	
	
	
	
	

	printf("\n\n BYE BYE NOW!\n");
  
}//END:endjob()
