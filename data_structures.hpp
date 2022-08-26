#ifndef _GKaiC_H_
#define _GKaiC_H_

/*---------------------------INCLUDE'S---------------------------------------*/
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>

/* SIMULATION WIDE FUNCTIONS */
#define F_COMP(a,b) (fabs ( a-b ) < 1e-10)
#define MIN(a,b)    ( a < b ? a : b )
#define MAX(a,b)    ( a > b ? a : b )

/* SIMULATON WIDE CONSTANTS */
#define EPSILON     1e-9
#define SCALE       1 - EPSILON
#define PI          3.141592653589793
#define INF_POS     1.0/0.0

#define CANCELLED   -1
#define FALSE 0
#define TRUE  1

#define DATA_NBR 3

#define HIST_RES 1

/* Structure definitions */
class Monomer;
class Hexamer;
class PropensityContainer;

/* Array of integers or doubles */
typedef struct int_array 
{
  int len,  ///< The length of the array
      *val; ///< A pointer to the first element
} IntArray;

typedef struct double_array 
{
  int len;  ///< The length of the array
  double *val; ///< A pointer to the first element
} DoubleArray;


/* Structure for all reaction constants. */
typedef struct reaction_constants
{ 
	//CI-nucleotide exchange
	double kCIhyd;
	double kCIADPoff0;
	double kCIADPoffA;

	double kAkIDoff;
	double kconf0;
	double kBswitch_f;
	double kBswitch_r;
	double ddGTDconf;
	double ddGTDAbind;  
	int    nIADPAref;
	int    nIADPconfref;  
	double dGIbindB;


	//CI-KaiA/KaiB sequestration.
	double kACIAon;
	double kACIAoff;
	double kICIAon;
	double kICIAoff;
	double kICIBGon;
	double kICIBFSon;
	double kICIBoff;  
	double kACIBGon;
	double kACIBFSon;
	double kACIBoff;  
	double nAseq;
	double nBseq;  

	//Energy differences in KaiA affinity for the CII domain due to its phosphorylation state.
	double dgACIIU;
	double dgACIIT;
	double dgACIID;
	double dgACIIS;      

	double dGICII;

	//Activation energy of ADP off rate in CI domain, Active state
	double dgACIActU;
	double dgACIActT;
	double dgACIActD;
	double dgACIActS;

	//Activation energy of ADP off rate in CI domain, Inactive state
	double dgICIActU;
	double dgICIActT;
	double dgICIActD;
	double dgICIActS;

	//CII-nucleotide exchange
	double kCIIAon;
	double kCIIAoff;
	double kCIIhyd0;
	double kCIIhydA;  
	double kCIInucloff0;
	double kCIInucloffA;
	double KATPoKADP;

	//CII Phosphotransfer rates.
	double kUT;
	double kTU;
	double kTD;
	double kDT;
	double kSD;
	double kDS;
	double kUS;
	double kSU;

	double kKidAon;
	double kKidAoff;	

} ReactionConstants;

/* Data container for output data */
typedef struct hexamer_avr
{
  double t;
  double p;
  double Afree;
  double ACI;
  double BCI;  
  double ACII;
  double Atot;
  int Btot;
  double CIATP;
  double CIIATP;
  double Stot;
  double Ttot;  
  double pU;
  double pT;
  double pS;
  double pD;
  double CIATPcons;
  double CIIATPcons;  
  double dCIATPcons;
  double dCIIATPcons;  
  double active;

  double dGACIIbind;
  double dGICIbind;
  double dGconfADP;
  double kADPoff;
  double CIKaiA_bound;
  double CIKaiB_bound;
  double CIKidA_bound;
  double BCI_any;
  int B_active;
  int B_inactive;
  int KaiBKidA;
  int n_B_rebind;
  long int n_CIKaiB_on;
  double prop_CIBon;
  double prop_CIAon;
  int n_max_CIKaiB_bound;
  
} HexamerAvr;

/* Global state variables */
typedef struct system_variables
{
  //Static vars.
  double ATPfrac;
  double KaiA0; //Initial conditions in uM.
  double KaiC0;
  double KaiB0;
  double fraction_fs = 0;
  double KidA0;
  double volume; //Volume is 1 cubic micron -> 600.
  
  double tend;
  double tsim;
  double tequ;
  double tincu;
  double t_hex_equil;
  
  bool start_phosphorylated;
  
  //Number of hexamers in system.
  int N_hexamers;
  
  //File names
  char output_filename[128];
  char param_filename[128];
  
  //Dynamic vars
  int Afree; //Free KaiA dimers;
  int B_active; //fsKaiB
  int B_inactive; //gsKaiB
  int B1_active; //fsKaiB tagged to track rebinding
  int KidA_free;
  int KaiBKidA; //KaiB-KidA complexes
  int KaiB1KidA; //complexes tagged to track rebinding
  double cAfree; //Free KaiA concentration 
  int CIATPcons; //Number of consumed ATP molecules in CI domain.
  int CIIATPcons; //Number of consumed ATP molecules in CII domain.
      
  //Counters
  int step_cntr;
  double t_sample_incr;
  int sample_cnt; //Max sample number.
  int sample_cntr; //Current sample number.
  
  int rnd_seed;
  long int n_CIKaiB_on;
  int n_B_rebind;
    
  //State flips
  uint FlipCntrFw;
  uint FlipCntrBw;
  
  /* Pointer to data struct for output data, 
     split to active and inactive hexamers*/
  HexamerAvr *Aoutput_data;
  HexamerAvr *Ioutput_data;
} SystemVariables;

struct simulation_parameters {

	double KaiA0 = -1;
	double KaiC0 = -1;
	double ATPfrac = -1;
	double ATPfrac_lo = -1;
	double amp_guess = -1;
	double weight_corr = 1;
	double weight_baseline = 1;
	double weight_amp = 1;
	double weight_ATP = 1;
	double t_equ = 0;
};

#endif
