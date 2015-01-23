#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define TIME 500
#define STEPCOUNT 180000*1.5
#define DAYSECS 86400
#define SWTIME 3*DAYSECS    // Switch after 3 days
#define SIMTIME 17*DAYSECS  // Run for a total of 17 days
#define MAXSTEPS 20000000
#define WRITESTEPMOD 20  
#define _SWITCH_ 1   // Replace _SWITCH_ with switch name configuration

int tri = 0;
int two = 0;
int sw2 = 0;
int sw3 = 0;

/*
	Switching mechanism;
	swState :
		0 ... no switch;
		1 ... MCT(iEC) -> BFP(iPC)
		2 ... BFP(iPC) -> MCT(iEC)
*/

int swState = 1;
int avg_bfp = 0;
int avg_mct = 0;
int sw_bfp = 0;
int sw_mct = 0;


// Random generator Hybrid Taus
unsigned z1, z2, z3, z4;  

unsigned TausStep(unsigned *z, int S1, int S2, int S3, unsigned M)  
{  
  unsigned b=(((*z << S1) ^ *z) >> S2);  
  return *z =  (((*z & M) << S3) ^ b);
}  

unsigned LCGStep(unsigned *z, unsigned A, unsigned C)  
{  
	return (A*(*z)+C);  
}  

float HybridTaus()  
{  
  return 2.3283064365387e-10 * (              // Periods  
    TausStep(&z1, 13, 19, 12, 4294967294UL) ^  // p1=2^31-1  
    TausStep(&z2, 2, 25, 4, 4294967288UL) ^    // p2=2^30-1  
    TausStep(&z3, 3, 11, 17, 4294967280UL) ^   // p3=2^28-1  
    LCGStep(&z4, 1664525, 1013904223UL)        // p4=2^32  
   );  
} 



int main(int argc, char *argv[]){

	// Output file
	FILE * fp;
	fp = fopen(argv[1],"w");
	
	//Switching
	swState = atoi(argv[3]);

    //Parameters - change values for parameter scan
    float vector[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    float na = 1.5;  //particle number
    
    //Parameters
	float ta                = 1.0*vector[0];					//binding probability
	float tda				= 0.0001*vector[1];                 //dissociation probability
	float tb                = 1.0*vector[2];                    //binding probability
	float tdb				= 0.0001*vector[3];                 //dissociation probability
	float tdg				= 0.00001*vector[4];                //degradation probability
	float r_vp16_vp16		= 0.0001*vector[5];                 //fully activated transcription

	float r_vp16_krab		= 0.01*r_vp16_vp16*vector[12];		//competitive binding

	float tk                = 3.0526*vector[14];					//KRAB competition factor; Determined using min RMSE;

	// float rdg				= 0.00001*vector[6];				//reporter degradation
	// float rdg			= 0.000006*vector[6];				//reporter degradation
    float rdg				= 0.000003*vector[6];				//reporter degradation

	float cta				= 0.5*r_vp16_vp16*0.1*vector[10];	//constitutive promoter transcription
    //float cmv				= 0.5*r_vp16_vp16*vector[10];	    //constitutive promoter transcription
    float cmv				= r_vp16_vp16;	                    //constitutive promoter transcription
    float cmv_a             = 10.0;                             //activated cmv
    float cmv_ar            = 2.0;                              //activated/reppressed cmv
    float cmv_leak          = 0.01;
    //float cta				= 1.0;	                            //constitutive promoter transcription

	float p_leak_pir		= 1.0*0.1*vector[7];				//pir promotor leakage
	float p_leak_etr		= 1.0*0.1*vector[8];				//etr promotor leakage
	float tas				= 1.0*vector[11];					//TAL A/B asymmetry
	float pas				= 1.0*vector[13];					//Induction asymmetry
	
	float p_leak_krab 		= 0.01;


        // Separate binding sites
	float r_vp16_krab_eq    = 0.5 * cmv;
	float r_vp16_vp16_eq    = 2.0 * cmv;

	
	//Float
	// float fbfp = 0.001*vector[9];
	// float fmct = 0.001*vector[9];
	float fbfp = 0.0032;  // Determined using min RMSE
	float fmct = 0.0032;  // Determined using min RMSE

	//TAL A/B factors
	float ta_factor = 1.0;
	float tb_factor = 1.0;
		
	//Load vectors
	    #ifdef CMVFEED
		#include "models/model_cmvfeed.h"
	    #endif
	    #ifdef MINFEED
		#include "models/model_minfeed.h"
	    #endif
	    #ifdef TOGGLE
		#include "models/model_toggle.h"
	    #endif
	    #ifdef NOLOOP
		#include "models/model_noloop.h"
	    #endif
	    #ifdef NOCOMP
		#include "models/model_nocomp.h"
	    #endif
	    #ifndef iDC  // Warning: taking care of legacy.
		#define iDC -1
	    #endif

	// Seed random
	srand(time(0));

	//High log
	#define HL_AGM 20
	#define HL_PI 3.14159
	#define HL_LN2 0.69314
	#define HL_MBITS 8
	float hl_f;
	float hl_a;
	float hl_g;
	float hl_at;
	float hl_s;
	float hl_x;
	float hl_y;
	float hl_lnf;

	//Misc
	int i		= 0;
	int j		= 0;
	int mi		= 0;
	int tstep	= 0;
	float tau	= 0;
	float t	= 0;
	
	//Ssa misc
	float r1;
	float r2;
	float wsum = 0;
	float wthresh[REACTIONS];
	float wtcumsum = 0;
	float csum = 0;
	int sw = 0;


	//HybridTaus random generator
	z1 = 131; z2 = 129; z3 = 130; z4 = rand();
	
	for (i = 0; i<REACTIONS; i++){
	    for (j = 0; j<SPECIES; j++){
	        dif[i][j] = -reactants[i][j]+products[i][j];
	    }
	}

    
	clock_t start = clock();
	while(t < SIMTIME && tstep < MAXSTEPS){

        r1 = ((float)rand())/RAND_MAX;
	r2 = ((float)rand())/RAND_MAX;

        #ifdef CMVFEED
            #include "weights/weights_cmvfeed.h"
        #endif
        #ifdef MINFEED
            #include "weights/weights_minfeed.h"
        #endif
        #ifdef TOGGLE
            #include "weights/weights_toggle.h"
        #endif
        #ifdef NOLOOP
            #include "weights/weights_noloop.h"
        #endif
        #ifdef NOCOMP
            #include "weights/weights_nocomp.h"
        #endif

		//Calculate wsum
		//weight thresholds
		wsum=0;
		wtcumsum = 0;
		for (i = 0; i<REACTIONS; i++) 	wsum = wsum + w[i];
		if (wsum <= 0){printf("Break\n"); break;}
		for (i = 0; i<REACTIONS; i++){
			wtcumsum = wtcumsum + w[i]/wsum;
			wthresh[i] = wtcumsum;
		}

		//Determine tau;
		//tau = 1.0/wsum * logf(1.0/(1-r1));

		//High log
		hl_f = 1/(1-r1);
		hl_s = hl_f*pow(2,HL_MBITS);		
		hl_y = 4.0/hl_s;
		hl_x = 1;
		//agm start
		hl_a = 0.5*(hl_x+hl_y);
		hl_g = sqrt(hl_x*hl_y);
		for(j=0; j<HL_AGM; j++){
			hl_at = hl_a;
			hl_a = 0.5*(hl_a+hl_g);
			hl_g = sqrt(hl_g*hl_at);
		}
		//agm end
		hl_lnf = HL_PI / (2*hl_a) - HL_MBITS * HL_LN2;

		//Determine tau.
		tau = 1/wsum * hl_lnf;

		//Determine reaction index.
		mi = -1;
		for (i = 0; i<REACTIONS; i++)	mi = (r2 < wthresh[i] && mi == -1)*i + (r2 >= wthresh[i] || mi != -1)*mi;
		mi = (mi > 0)*mi; // Keep from -1. if nothing meets the threshold, set zero

		//Update species matrix;
		for (i=0; i<SPECIES; i++){
		    species[i] = species[i] + dif[mi][i];
		    if(species[i] < 0) {
		        species[i] = 0;
		        //printf("Minor %d\n", i);
		    }
		}

		t = t+tau;
		tstep++;

	    // Print results to file
	    if (tstep % WRITESTEPMOD == 0 ){
            if(tri){
                fprintf(fp,"%f %f %f %d %f\n", species[0]*fbfp, species[1]*fmct, species[2]*fbfp, mi, t);
            }
            else{
                fprintf(fp,"%f %f %d %f\n", species[0]*fbfp, species[1]*fmct, mi, t);
            }
	    }

		//Switch mechanism;
            if (!tri && t >= SWTIME && !sw){
            sw = 1;

		//Toggle
		species[iPC] = 0;
		species[iEC] = 0;
		
	
		//1 ... MCT(iEC) -> BFP(iPC)
		if(swState == 1){
		//printf("switch 1! %d\n", tstep);
			species[iPC] = INDMAX;
			species[iEC] = 0;
			//species[E_KRAB] += species[iEC_E_KRAB];
			//species[iEC_E_KRAB]  = 0;

		}
		//2 ... BFP(iPC) -> MCT(iEC)
		if(swState == 2){
		//printf("switch 2! %d \n", tstep);
			species[iPC] = 0;
			species[iEC] = INDMAX;
			//species[PIP_KRAB] += species[iPC_PIP_KRAB];
			//species[iPC_PIP_KRAB]  = 0;
		}

	}
	
	avg_bfp += species[BFP];
	avg_mct += species[mCITRINE];

	}
	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;

	avg_bfp = avg_bfp/STEPCOUNT;
	avg_mct = avg_mct/STEPCOUNT;
	
	return 0;
}
