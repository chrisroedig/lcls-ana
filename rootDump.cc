#include <TROOT.h>
#include <TH1F.h>
#include <TProfile.h>
#include <math.h>
#include "rootDump.hh"
void dumpArrays( char* name,double* xArr, double* yArr, int length ){
	int i=0;
	FILE* Outfile;
	
	printf("Dumping Array of doubles to File: ");
	printf(name);
	fflush(stdout);

	Outfile=fopen(name,"w");
	for(i=0;i<length;i++){			
		fprintf(Outfile,"%e %e\n",xArr[i],yArr[i]);			
	}
	printf("  done!\n");
	fclose(Outfile);
}
void arrayTH1F(double* xArr,double* yArr,int length,TH1F* hist){
	int i=0;
	double xbin=0.0;
	Float_t* ydat;	
	ydat = hist -> GetArray();
	
	for(i=0;i<length;i++){
		
		xbin = (double)hist -> GetBinCenter(i);
		xArr[i]=xbin;
		yArr[i]=(double)ydat[i];
		
	}	
}

void smoothArray(double* arrayIn,double* arrayOut,int length, int boxWidth){
		/* PERFORMS MOVING AVG

		   			
		*/
		int j=0;
		int halfBox=(int)ceil((double)boxWidth/2.0);
		for(j=0;j<length;j++){
			double avg=0.0;
			if((j>halfBox)&&(j<(length-halfBox))){
				int k=0;
				for(k=0;k<boxWidth;k++){
					avg+=arrayIn[j-halfBox+k];			
				}
			avg=avg/boxWidth;	
			}
			arrayOut[j]=avg;
		}
}
void divArray(double* arrayIn1,double* arrayIn2,double* arrayOut,int length, double limit){
		int j=0;
		
		for(j=0;j<length;j++){
			arrayOut[j]=(arrayIn2[j]>limit)?(arrayIn1[j]/arrayIn2[j]):0.0;
		}
}
void dumpTpFile(char* name,TProfile* data,int numSamples){
	int i=0;
	FILE* Outfile;
	
	printf("Dumping TProfile to File: ");printf(name);fflush(stdout);
	Outfile=fopen(name,"w");
	
	Double_t xbin=0.0;
	Double_t* ydat;	
	ydat = data -> GetArray();
	
	for(i=0;i<numSamples;i++){
		
		xbin = data -> GetBinCenter(i);
		
		fprintf(Outfile,"%e %e\n",xbin,ydat[i]);			
	}
	printf("  done!\n");
	fclose(Outfile);
}

void dumpThFile(char* name,TH1F* data,int numSamples){
	int i=0;
	FILE* Outfile;
	
	printf("Dumping T1H to File: ");printf(name);fflush(stdout);
	Outfile=fopen(name,"w");
	
	Double_t xbin=0.0;
	Float_t* ydat;	
	ydat = data -> GetArray();
	
	for(i=0;i<numSamples;i++){
		
		xbin = data -> GetBinCenter(i);
		
		fprintf(Outfile,"%e %e\n",xbin,ydat[i]);			
	}
	printf("  done!\n");
	fclose(Outfile);
}
void dumpThMatrixFile(char* name,TH1F** data,int numSamples,int steps){
	int i=0;
	int j=0;
	FILE* Out;
	
	printf("Dumping T1H to File: ");printf(name);fflush(stdout);
	Out=fopen(name,"w");
	
	
	Double_t xbin;
	Float_t* ydat;	
	
	for(j=0;j<steps;j++){
		ydat = data[j] -> GetArray();
		for(i=0;i<numSamples;i++){
			xbin = data[j] -> GetBinCenter(i);
			fprintf(Out," %e ",ydat[i]);			
		}
		fprintf(Out,"\n");
		
	}
	printf("  done!\n");
	fclose(Out);
}
void dumpThNormMatrixFile(char* name,TH1F** data,int numSamples,int steps,double* norm){
	int i=0;
	int j=0;
	FILE* Out;
	
	printf("Dumping T1H to File: ");printf(name);fflush(stdout);
	Out=fopen(name,"w");
	
	Double_t xbin;
	Float_t* ydat;	
	
	for(j=0;j<steps;j++){
		ydat = data[j] -> GetArray();
		for(i=0;i<numSamples;i++){
			xbin = data[j] -> GetBinCenter(i);
			double yval = (double)ydat[i]/norm[j];	
			fprintf(Out," %e ",yval);			
		}
		fprintf(Out,"\n");
		
	}
	printf("  done!\n");
	fclose(Out);
}
void dumpTpMatrixFile(char* name,TProfile** data,int numSamples,int steps){
	int i=0;
	int j=0;
	FILE* Out;
	
	printf("Dumping T1H to File: ");printf(name);fflush(stdout);
	Out=fopen(name,"w");
	
	Double_t xbin;
	Double_t* ydat;	
	
	for(j=0;j<steps;j++){
		ydat = data[j] -> GetArray();
		for(i=0;i<numSamples;i++){
			xbin = data[j] -> GetBinCenter(i);
			fprintf(Out," %e ",ydat[i]);			
		}
		fprintf(Out,"\n");
		
	}
	printf("  done!\n");
	fclose(Out);
}

