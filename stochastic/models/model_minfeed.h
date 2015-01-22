#define SPECIES 28
#define REACTIONS 46
#define INDMAX 70

int ctmp;
int rindex = 0;
int species[SPECIES]				={0};
int reactants[REACTIONS][SPECIES]		= {0};
int products[REACTIONS][SPECIES]		= {0};
int dif[REACTIONS][SPECIES]			= {0};
float w[REACTIONS]				= {0};
float c[REACTIONS]				= {0};



/*Species*/
#define BFP				0
#define mCITRINE		1
#define TALA_KRAB		2
#define TALA_VP16		3
#define TALB_KRAB		4
#define TALB_VP16		5
#define P0						6
#define P0_TALA_VP16_KRAB		7	
#define P0_TALA_KRAB_KRAB		8
#define P0_TALA_VP16_VP16		9
#define P1						10
#define P1_TALA_VP16_KRAB		11	
#define P1_TALA_KRAB_KRAB		12
#define P1_TALA_VP16_VP16		13
#define P2						14
#define P2_TALB_VP16_KRAB		15
#define P2_TALB_KRAB_KRAB		16
#define P2_TALB_VP16_VP16		17
#define P3						18
#define P3_TALB_VP16_KRAB		19
#define P3_TALB_KRAB_KRAB		20
#define P3_TALB_VP16_VP16		21
#define P4						22
#define P5						23
#define P6						24
#define P7						25
#define iPC						26
#define iEC						27


species[P0] = 15 * na;
species[P1] = 3 * na;
species[P2] = 15 * na;
species[P3] = 3 * na;
species[P4] = 15 * na;
species[P5] = 5 * na;
species[P6] = 15 * na;
species[P7] = 5 * na;
species[iPC] = (atoi(argv[2])==1)*INDMAX;
species[iEC] = (atoi(argv[2])==2)*INDMAX;





//Competitive binding
/*Pmin 0 TAL BINDING 0-5*/
reactants[rindex][P0] = 1; reactants[rindex][TALA_KRAB] = 2; 
products[rindex][P0_TALA_KRAB_KRAB] = 1;
c[rindex] = ta*tk;  rindex++;
reactants[rindex][P0_TALA_KRAB_KRAB] = 1;
products[rindex][P0] = 1; products[rindex][TALA_KRAB] = 2; 
c[rindex] = tda;  rindex++;

reactants[rindex][P0] = 1; reactants[rindex][TALA_VP16] = 2; 
products[rindex][P0_TALA_VP16_VP16] = 1;
c[rindex] = ta;  rindex++;
reactants[rindex][P0_TALA_VP16_VP16] = 1;
products[rindex][P0] = 1; products[rindex][TALA_VP16] = 2; 
c[rindex] = tda;  rindex++;

reactants[rindex][P0] = 1; reactants[rindex][TALA_VP16] = 1; reactants[rindex][TALA_KRAB] = 1; 
products[rindex][P0_TALA_VP16_KRAB] = 1;
c[rindex] = ta*tk;  rindex++;
reactants[rindex][P0_TALA_VP16_KRAB] = 1; 
products[rindex][P0] = 1; products[rindex][TALA_VP16] = 1; products[rindex][TALA_KRAB] = 1; 
c[rindex] = tda;  rindex++;


/*Pmin 1 TAL BINDING 6-11*/
reactants[rindex][P1] = 1; reactants[rindex][TALA_KRAB] = 2; 
products[rindex][P1_TALA_KRAB_KRAB] = 1;
c[rindex] = ta*tk;  rindex++;
reactants[rindex][P1_TALA_KRAB_KRAB] = 1;
products[rindex][P1] = 1; products[rindex][TALA_KRAB] = 2; 
c[rindex] = tda;  rindex++;

reactants[rindex][P1] = 1; reactants[rindex][TALA_VP16] = 2; 
products[rindex][P1_TALA_VP16_VP16] = 1;
c[rindex] = ta;  rindex++;
reactants[rindex][P1_TALA_VP16_VP16] = 1;
products[rindex][P1] = 1; products[rindex][TALA_VP16] = 2; 
c[rindex] = tda;  rindex++;

reactants[rindex][P1] = 1; reactants[rindex][TALA_VP16] = 1; reactants[rindex][TALA_KRAB] = 1; 
products[rindex][P1_TALA_VP16_KRAB] = 1;
c[rindex] = ta*tk;  rindex++;
reactants[rindex][P1_TALA_VP16_KRAB] = 1;
products[rindex][P1] = 1; products[rindex][TALA_VP16] = 1; products[rindex][TALA_KRAB] = 1; 
c[rindex] = tda;  rindex++;

/*Pmin 2 TAL BINDING 12-17*/
reactants[rindex][P2] = 1; reactants[rindex][TALB_KRAB] = 2; 
products[rindex][P2_TALB_KRAB_KRAB] = 1;
c[rindex] = tb*tk;  rindex++;
reactants[rindex][P2_TALB_KRAB_KRAB] = 1;
products[rindex][P2] = 1; products[rindex][TALB_KRAB] = 2;
c[rindex] = tdb;  rindex++;

reactants[rindex][P2] = 1; reactants[rindex][TALB_VP16] = 2; 
products[rindex][P2_TALB_VP16_VP16] = 1;
c[rindex] = tb;  rindex++;
reactants[rindex][P2_TALB_VP16_VP16] = 1;
products[rindex][P2] = 1; products[rindex][TALB_VP16] = 2; 
c[rindex] = tdb;  rindex++;

reactants[rindex][P2] = 1; reactants[rindex][TALB_VP16] = 1; reactants[rindex][TALB_KRAB] = 1; 
products[rindex][P2_TALB_VP16_KRAB] = 1;
c[rindex] = tb*tk;  rindex++;
reactants[rindex][P2_TALB_VP16_KRAB] = 1;
products[rindex][P2] = 1; products[rindex][TALB_VP16] = 1; products[rindex][TALB_KRAB] = 1; 
c[rindex] = tdb;  rindex++;

/*Pmin 3 TAL BINDING 18-23*/
reactants[rindex][P3] = 1; reactants[rindex][TALB_KRAB] = 2; 
products[rindex][P3_TALB_KRAB_KRAB] = 1;
c[rindex] = tb*tk;  rindex++;
reactants[rindex][P3_TALB_KRAB_KRAB] = 1;
products[rindex][P3] = 1; products[rindex][TALB_KRAB] = 2;
c[rindex] = tdb;  rindex++;

reactants[rindex][P3] = 1; reactants[rindex][TALB_VP16] = 2; 
products[rindex][P3_TALB_VP16_VP16] = 1;
c[rindex] = tb;  rindex++;
reactants[rindex][P3_TALB_VP16_VP16] = 1;
products[rindex][P3] = 1; products[rindex][TALB_VP16] = 2; 
c[rindex] = tdb;  rindex++;

reactants[rindex][P3] = 1; reactants[rindex][TALB_VP16] = 1; reactants[rindex][TALB_KRAB] = 1; 
products[rindex][P3_TALB_VP16_KRAB] = 1;
c[rindex] = tb*tk;  rindex++;
reactants[rindex][P3_TALB_VP16_KRAB] = 1;
products[rindex][P3] = 1; products[rindex][TALB_VP16] = 1; products[rindex][TALB_KRAB] = 1;
c[rindex] = tdb;  rindex++;



/*P0 Transcription 24-25*/
reactants[rindex][P0_TALA_VP16_KRAB] = 1; products[rindex][P0_TALA_VP16_KRAB] = 1;
products[rindex][TALB_KRAB] = 1; products[rindex][BFP] = 1;
c[rindex] = r_vp16_krab*1.0/ta_factor;  rindex++;

reactants[rindex][P0_TALA_VP16_VP16] = 1; products[rindex][P0_TALA_VP16_VP16] = 1;
products[rindex][TALB_KRAB] = 1; products[rindex][BFP] = 1;
c[rindex] = r_vp16_vp16*ta_factor;  rindex++;


/*P1 Transcription 26-27*/
reactants[rindex][P1_TALA_VP16_KRAB] = 1; products[rindex][P1_TALA_VP16_KRAB] = 1;
products[rindex][TALA_VP16] = 1; 
c[rindex] = r_vp16_krab*1.0/ta_factor;  rindex++;

reactants[rindex][P1_TALA_VP16_VP16] = 1; products[rindex][P1_TALA_VP16_VP16] = 1;
products[rindex][TALA_VP16] = 1; 
c[rindex] = r_vp16_vp16*ta_factor;  rindex++;

/*P2 Transcription 28-29*/
reactants[rindex][P2_TALB_VP16_KRAB] = 1; products[rindex][P2_TALB_VP16_KRAB] = 1;
products[rindex][TALA_KRAB] = 1; products[rindex][mCITRINE] = 1;
c[rindex] = r_vp16_krab*1.0/tb_factor;  rindex++;

reactants[rindex][P2_TALB_VP16_VP16] = 1; products[rindex][P2_TALB_VP16_VP16] = 1;
products[rindex][TALA_KRAB] = 1; products[rindex][mCITRINE] = 1;
c[rindex] = r_vp16_vp16*tb_factor;  rindex++;

/*P3 Transcription 30-31*/
reactants[rindex][P3_TALB_VP16_KRAB] = 1; products[rindex][P3_TALB_VP16_KRAB] = 1;
products[rindex][TALB_VP16] = 1; 
c[rindex] = r_vp16_krab*1.0/tb_factor;  rindex++;

reactants[rindex][P3_TALB_VP16_VP16] = 1; products[rindex][P3_TALB_VP16_VP16] = 1;
products[rindex][TALB_VP16] = 1; 
c[rindex] = r_vp16_vp16*tb_factor;  rindex++;


/*Induction 32-35*/
reactants[rindex][P4] = 1; reactants[rindex][iPC] = 1; 
products[rindex][P4] = 1; products[rindex][iPC] = 1; products[rindex][TALB_KRAB] = 1;
c[rindex] = cta;  rindex++;

reactants[rindex][P5] = 1; reactants[rindex][iPC] = 1; 
products[rindex][P5] = 1; products[rindex][iPC] = 1; products[rindex][TALA_VP16] = 1;
c[rindex] = cta;  rindex++;

reactants[rindex][P6] = 1; reactants[rindex][iEC] = 1; 
products[rindex][P6] = 1; products[rindex][iEC] = 1; products[rindex][TALA_KRAB] = 1;
c[rindex] = cta;  rindex++;

reactants[rindex][P7] = 1; reactants[rindex][iEC] = 1; 
products[rindex][P7] = 1; products[rindex][iEC] = 1; products[rindex][TALB_VP16] = 1;
c[rindex] = cta;  rindex++;

/*Leakage 36-39*/
reactants[rindex][P4] = 1;
products[rindex][P4] = 1;  products[rindex][TALB_KRAB] = 1;
c[rindex] = cta*p_leak_pir;  rindex++;

reactants[rindex][P5] = 1; 
products[rindex][P5] = 1; products[rindex][TALA_VP16] = 1;
c[rindex] = cta*p_leak_pir;  rindex++;

reactants[rindex][P6] = 1;  
products[rindex][P6] = 1; products[rindex][TALA_KRAB] = 1;
c[rindex] = cta*p_leak_etr;  rindex++;

reactants[rindex][P7] = 1; 
products[rindex][P7] = 1; products[rindex][TALB_VP16] = 1;
c[rindex] = cta*p_leak_etr;  rindex++;

/*Degradation 39-44*/
reactants[rindex][BFP] = 1; 
c[rindex] = rdg;  rindex++;
reactants[rindex][mCITRINE] = 1; 
c[rindex] = rdg;  rindex++;
reactants[rindex][TALA_KRAB] = 1; 
c[rindex] = tdg;  rindex++;
reactants[rindex][TALA_VP16] = 1; 
c[rindex] = tdg;  rindex++;
reactants[rindex][TALB_KRAB] = 1; 
c[rindex] = tdg;  rindex++;
reactants[rindex][TALB_VP16] = 1; 
c[rindex] = tdg;  rindex++;


//Weights;
if(argc == 1){
    int ii;
    int jj;
    int rr;
    FILE * fw;
    fw = fopen("weights_minfeed.h", "w");
    for(ii=0; ii<REACTIONS; ii++){
        fprintf(fw, "w[%d] = c[%d] ",ii,ii,ii);
        for(jj=0; jj<SPECIES; jj++){
            rr = reactants[ii][jj];
            while(rr > 0) {
                //printf("* (species[%d]-%d)", jj, rr-1);
                fprintf(fw, "* (species[%d]-%d)", jj, rr-1);
                rr--;
            }
        }
        fprintf(fw,";\n");
    }
    fclose(fw);
    exit(0);
}
