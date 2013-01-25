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


//misc Globals
static int n_event=0;

//ETOF DAQ GLOBAL VARS
static int     numChannelsETof = 0; //should be 5 :)
static int     numSamplesETof = 0;
static double  sampleIntervalETof = 0;//total acquisition interval        

//fixed params for each etof
static double tPrompt[5];
static double vRet[5];


//etof energy bins
const int eBinCount=3000;//DO NOT CHANGE THIS

//ETOF HISTOGRAMS
static TH1F*	etofTotal[5];

//smoothing widths in energy bins
const int bcWidth1=5; //smoothed for side bands
const int bcWidth2=60; //smoothed to obscure sidebands


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
	
	
	//numbering accounts for switch of #3/#4 acq lines for this run 	
	tPrompt[0]=6.42953e-7;  //ETOF#1
	tPrompt[1]=6.41819e-7;  //ETOF#2
	tPrompt[2]=6.40565e-7;  //ETOF#4
	tPrompt[3]=6.40565e-7;  //ETOF#3
	tPrompt[4]=6.42382e-7;  //ETOF#5
	
	vRet[0]=700.0;			//ETOF#1
	vRet[1]=700.0;			//ETOF#2
	vRet[2]=700.0;			//ETOF#4
	vRet[3]=700.0;			//ETOF#3
	vRet[4]=700.0;			//ETOF#5
	
	
	
	printf("\n ETOF PARAMS:\n");
	for(i=0;i<5;i++){
		printf("ETOF %d -- VRET %e TPROMPT %e\n",i,vRet[i],tPrompt[i]);
		fflush(stdout);
		sprintf(name,"etof%dtotal",i);	
		etofTotal[i]=new TH1F(name,name,eBinCount, 700,1000); 	
		
	}


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
		printf("\n[ %d events ] \n",n_event);
	}
	if(n_event%1000==0){printf("|");}
	else if(n_event%100==0){printf(".");}
	fflush(stdout); 

	//get etof data..
	double* timeETof[5];
	double* voltageETof[5];
    
    for(i=0;i<5;i++){
		fail = getAcqValue( AmoETof, i, timeETof[i], voltageETof[i]);  			
		}
  	
	if(fail){printf("E");return;}//on fail return
	
	//count all traces into respective accumulation histograms
	for(i=0;i<5;i++){
		fillCalCF(timeETof[i],voltageETof[i],numSamplesETof,-.00,-.008,tPrompt[i],vRet[i],0.0,etofTotal[i]);
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
	int i;char name[32]; 
	
	//set up to use etof transmission file
	sprintf(name,"indata/etofTrans.dat");
	initEtofTrans(name);

	for(i=0;i<5;i++){
		printf("Post Proc on ETOF%d\n",i);
		double eAx[eBinCount];//energy (x) axis
		double iAx[eBinCount];//intensity (y) axis
		double cAx[eBinCount];//tranmission corrected intesity (y) axis
		double sAx1[eBinCount];//smoothed intesity (y) axis
		double sAx2[eBinCount];//smoothed intesity (y) axis

		//create a double array from histo		
		arrayTH1F(eAx,iAx,eBinCount,etofTotal[i]);

		
		printf("----correcting for transmission...");
		etofTrans(eAx,iAx,cAx,eBinCount,vRet[i]);
	
		printf("done\n");
		

		printf("----smoothing by %d bin running avg...",bcWidth1);
		smoothArray(cAx,sAx1,eBinCount,bcWidth1);
		printf("done\n");

		printf("----smoothing by %d bin running avg...",bcWidth2);
		smoothArray(cAx,sAx2,eBinCount,bcWidth2);
		printf("done\n");

		
		//write to file
		sprintf(name,"outdata/etof%dRaw.dat",i);//name of file
		dumpArrays(name,eAx,iAx,eBinCount);//dump array to file	

		//write to file
		sprintf(name,"outdata/etof%dCorrected.dat",i);//name of file
		dumpArrays(name,eAx,cAx,eBinCount);//dump array to file	

		//write to file
		sprintf(name,"outdata/etof%dSmooth%d.dat",i,bcWidth1);//name of file
		dumpArrays(name,eAx,sAx1,eBinCount);//dump array to file	

		//write to file
		sprintf(name,"outdata/etof%dSmooth%d.dat",i,bcWidth2);//name of file
		dumpArrays(name,eAx,sAx2,eBinCount);//dump array to file	
		printf("done with ETOF%d\n\n",i);
		
	}
	
	printf("\n\n BYE BYE NOW!\n");
  
}//END:endjob()
