/*
  Program to simulate the in-vitro Kai system using the Kinetic Montecarlo Algorithm.
  The model describes a set of individually simulated hexamers, which, in turn,
  consist of six, individuall described, monomers each.
  
  USAGE:
  
  > ./KMCKaiC <parameter_filename> (Assumes default.par if nothing is given)
  
  OUTPUT: Timetraces of variables relating to the hexamer states, averaged over all hexamers.
  Program creates three files; spliot in active, inactive and total hexamers.
*/

#include <sstream>
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <random>

#include "data_structures.hpp"
#include "Hexamer.hpp"
#include "PropensityContainer.hpp"
#include "propagate.hpp"

using namespace std;

/* Local function definitions */
ReactionConstants initialize_reaction_consts(SystemVariables *sys);
HexamerAvr calc_hex_averages( Hexamer *hexamers, SystemVariables *sys, ReactionConstants *react_const );
void initialize_system_vars(SystemVariables *sys);
void write_outputfile(SystemVariables *sys);
void allocate_data_array(SystemVariables *sys);
string delUnnecessary(string &str);

/* Let's roll */
int main(int argc, char *argv[])
{   
// First set name of parameter file (default.par if no argument is given).
	SystemVariables sys;
	random_device rd;
	mt19937_64 engine(rd());
	uniform_real_distribution<double> u01(0, 1);

	if(argc > 1)
	{
		strncpy(sys.param_filename, argv[1], 128);
	}
	else
	{
		strcpy(sys.param_filename,"default.par");
	}
	cerr << "Parameter file set to: " << sys.param_filename << endl; 

	//Read system variables from parameter file and initialize.
	initialize_system_vars(&sys);

	//Read reaction rates from parameter file.
	ReactionConstants reaction_consts( initialize_reaction_consts(&sys) );

	//Initialize random number generator with (time-dependent, 0) seed.

	//Create proensity container, Hexamers and Monomer objects.
	PropensityContainer prop_cont(sys.N_hexamers); 
	Hexamer *hexamers = new Hexamer[sys.N_hexamers];
	//Monomer *monomers = new Monomer[6*sys.N_hexamers];

	allocate_data_array(&sys);

//Initialize hexamer and monomer states and initialize their propensities.
	for(int i(0); i<sys.N_hexamers; i++)
	{
		hexamers[i].set_prop_container(&prop_cont);
		hexamers[i].set_reaction_consts(&reaction_consts); 
		hexamers[i].set_sysvars(&sys);
		hexamers[i].initialize_state(i,1,0,0,0,6,6,0,0,0); 
		//hexamers[i].set_sextet(monomers, u01, engine); 
		hexamers[i].set_propensities();
	}

	double t_sample(sys.tequ);
	cerr << "System initialized with " << sys.Afree << " KaiA dimers and " << sys.N_hexamers 
	<< " KaiC hexamers. " << endl;
	std::cerr << "fsKaiB: " << sys.B_active << " gsKaiB:" << sys.B_inactive << endl;
	std::cerr << "KidA: " << sys.KidA_free << endl;
	cerr << "Equilibration time: " << sys.tequ << ", Max simulation time: " << sys.tend << endl;   

	/*int j = 0;
	std::ifstream infile;
	infile.open(argv[2]);
	std::vector<double> pulse_value;
	double val;

	while(infile >> val) {

		pulse_value.push_back(val);
	}

	infile.close();
	infile.open(argv[3]);
	std::vector<double> pulse_time;

	while(infile >> val) {

		pulse_time.push_back(val);
	}

	infile.close();*/

	/*engine.seed(0);
	int KidA_high = sys.KidA_free;
	sys.KidA_free = 0;
	bool KidA_changed = false;*/

	//Propagate simulation until tend is reached.
	int step_count = 0;
	while( propagate(&sys, hexamers, &prop_cont, &reaction_consts, u01, engine) ) 
	{ 
		/*if(!KidA_changed && sys.tsim > 44) {

			KidA_changed = true;
			sys.KidA_free = KidA_high;
		}*/
		step_count += 1;
		// Record samples for output data 
		if(sys.tsim > t_sample)
		{   
			//Calculate quantities of interest.
			HexamerAvr hex_avr( calc_hex_averages( hexamers, &sys, &reaction_consts ) );
  
			//Print to screen.
			//if(sys.sample_cntr%100 == 0)
			//{     
			//cerr << endl << "### SIMULATOR STATE ###" << endl;
			//cerr << "t: " << sys.tsim << ", p(t): " << hex_avr.p << ", Act(t): " << hex_avr.active << endl;
			//cerr << "Afree: " << hex_avr.Afree << ", Atot: " << hex_avr.Atot 
			//     << ", ACII: " << hex_avr.ACII << ", ACI: " << hex_avr.ACI << endl;
			//cerr << "CIATP: " << hex_avr.CIATP << ", CIIATP: " << hex_avr.CIIATP << endl;
			//cerr << "pU: " << hex_avr.pU << ", pT: " << hex_avr.pT 
			//     << ", pD: " << hex_avr.pD << ", pS: " << hex_avr.pS << endl;
			//}

			t_sample += sys.t_sample_incr;
			sys.sample_cntr++;
		}

		/*if(j < pulse_value.size() && sys.tsim > pulse_time[j]) {

			sys.ATPfrac = pulse_value[j];

			for(int i = 0; i < sys.N_hexamers; i += 1) {

				hexamers[i].set_propensities();
			}
			j += 1;
		}*/
	}

	// ### For normal output ###
	write_outputfile(&sys);  
	
	std::cout << step_count << "\n";
	cerr << "System Ready" << endl;
}
  
  
//Loops over all hexamers and calculates avergaes which are store in hex_avr_data;
HexamerAvr calc_hex_averages( Hexamer *hexamers, SystemVariables *sys, ReactionConstants *reaction_consts )
{
  HexamerAvr Ahex_avr_data, Ihex_avr_data;
  
  //State dependent variables.
  int ICIBA(0), ACIBA(0);
  int ACIB(0), ICIB(0); 
  int IACII(0), AACII(0);
  int ICIATP(0), ACIATP(0), ICIATPstd(0), ACIATPstd(0);
  int ICIIATP(0), ACIIATP(0);
  int IStot(0), AStot(0);
  int ITtot(0), ATtot(0);
  int IpU(0), ApU(0);
  int IpT(0), ApT(0);
  int IpS(0), ApS(0);
  int IpD(0), ApD(0);
  int Ip(0), Ap(0);
  int ACIB_any(0);
  int ICIB_any(0);
  int ACIKidA_bound = 0;
  int ICIKidA_bound = 0;

  //State independent variables;
  int prev_CIATPcons(0), prev_CIIATPcons(0);
  int active(0), inactive(0);
  int Atot = 0;
  
  //Energies
  double AdGconfADP(0.), IdGconfADP(0.), 
         IdGACIIbind(0.), AdGACIIbind(0.), 
         AADPoff(0.), IADPoff(0.), 
	 ACIKaiA_bound(0), ICIKaiA_bound(0),
	 ACIKaiB_bound(0), ICIKaiB_bound(0);
  double A_prop_CIBon = 0;
  double I_prop_CIBon = 0;
  double A_prop_CIAon = 0;
  double I_prop_CIAon = 0;
  int A_n_max_CIKaiB_bound = 0;
  int I_n_max_CIKaiB_bound = 0;
  
  double N_hexamers(sys->N_hexamers);
  //double N_hexamers(1);  
  
  if(sys->sample_cntr > 0)
  {
    prev_CIATPcons = sys->Aoutput_data[sys->sample_cntr-1].CIATPcons * (6*N_hexamers);
    prev_CIIATPcons = sys->Aoutput_data[sys->sample_cntr-1].CIIATPcons * (6*N_hexamers);    
  }
      
  for(int i(0); i<sys->N_hexamers; i++)
  {  
    if( hexamers[i].get_active() )
    {
      Ap += hexamers[i].get_hex_pT() + hexamers[i].get_hex_pS() + hexamers[i].get_hex_pD();
      if( hexamers[i].get_CIKaiB_bound() == reaction_consts->nBseq && 
          hexamers[i].get_CIKaiA_bound() == 0 ) ACIB++;
      if( hexamers[i].get_CIKaiA_bound() > 0 ) ACIBA += 1;
      if( hexamers[i].get_CIKaiB_bound() > 0) ACIB_any += 1;
      AACII += hexamers[i].get_CIIKaiA_bound();
      ACIATP += hexamers[i].get_hex_CIATP();
      ACIATPstd += hexamers[i].get_hex_CIATP()*hexamers[i].get_hex_CIATP();     
      ACIIATP += hexamers[i].get_hex_CIIATP();
      ATtot += hexamers[i].get_hex_T();    
      AStot += hexamers[i].get_hex_S();
      ApU += hexamers[i].get_hex_pU();
      ApT += hexamers[i].get_hex_pT();
      ApD += hexamers[i].get_hex_pD();
      ApS += hexamers[i].get_hex_pS();
      AdGconfADP += hexamers[i].kconf_f();
      AdGACIIbind += hexamers[i].kCIIAon();
      AADPoff += hexamers[i].kCIADPoff();      
      ACIKaiA_bound += hexamers[i].get_CIKaiA_bound();
      ACIKaiB_bound += hexamers[i].get_CIKaiB_bound();
      ACIKidA_bound += hexamers[i].get_CIKidA_bound();
      A_prop_CIBon += hexamers[i].get_prop(2);
      A_prop_CIBon += hexamers[i].get_prop(3);
      A_prop_CIAon += hexamers[i].get_prop(10);
      if( hexamers[i].get_CIKaiB_bound() == reaction_consts->nBseq) A_n_max_CIKaiB_bound += 1;
      active++;
    }
    else
    {  
      Ip += hexamers[i].get_hex_pT() + hexamers[i].get_hex_pS() + hexamers[i].get_hex_pD();     
      if( hexamers[i].get_CIKaiB_bound() == reaction_consts->nBseq && 
          hexamers[i].get_CIKaiA_bound() == 0 ) ICIB++;
      if( hexamers[i].get_CIKaiA_bound() > 0 ) ICIBA += 1;
      if( hexamers[i].get_CIKaiB_bound() > 0) ICIB_any += 1; 
      IACII += hexamers[i].get_CIIKaiA_bound();
      ICIATP += hexamers[i].get_hex_CIATP();
      ICIATPstd += hexamers[i].get_hex_CIATP()*hexamers[i].get_hex_CIATP();
      ICIIATP += hexamers[i].get_hex_CIIATP();
      ITtot += hexamers[i].get_hex_T();    
      IStot += hexamers[i].get_hex_S();
      IpU += hexamers[i].get_hex_pU();
      IpT += hexamers[i].get_hex_pT();
      IpD += hexamers[i].get_hex_pD();
      IpS += hexamers[i].get_hex_pS();
      IdGconfADP += hexamers[i].kconf_b();            
      IADPoff += hexamers[i].kCIADPoff();
      IdGACIIbind += hexamers[i].kCIIAon();
      ICIKaiA_bound += hexamers[i].get_CIKaiA_bound();
      ICIKaiB_bound += hexamers[i].get_CIKaiB_bound();
      ICIKidA_bound += hexamers[i].get_CIKidA_bound();
      I_prop_CIBon += hexamers[i].get_prop(2);
      I_prop_CIBon += hexamers[i].get_prop(3);
      I_prop_CIAon += hexamers[i].get_prop(10);
      if( hexamers[i].get_CIKaiB_bound() == reaction_consts->nBseq) I_n_max_CIKaiB_bound += 1;
    }
    Atot += hexamers[i].get_CIKaiA_bound() + hexamers[i].get_CIIKaiA_bound();      
  }
  Atot += sys->Afree;
  
  inactive = sys->N_hexamers - active;  
 
  double KaiA0(sys->KaiA0);
  double KaiC0(sys->KaiC0);
  if(KaiA0 == 0)
    KaiA0 = 1;
  if(KaiC0 == 0)
    KaiC0 = 1;

  Ahex_avr_data.CIATPcons  = (double) sys->CIATPcons/(6*N_hexamers);
  Ahex_avr_data.dCIATPcons = 
    (double) 24. * (sys->CIATPcons - prev_CIATPcons)/(6*N_hexamers*sys->t_sample_incr);
  Ahex_avr_data.CIIATPcons  = (double) sys->CIIATPcons/(6*N_hexamers);
  Ahex_avr_data.dCIIATPcons = 
    (double) 24. * (sys->CIIATPcons - prev_CIIATPcons)/(6*N_hexamers*sys->t_sample_incr);

  Ahex_avr_data.active   = (double) active;  
  Ihex_avr_data.active   = (double) inactive;

  Ahex_avr_data.Atot     = (double) Atot / (KaiA0 * sys->volume);    
  Ahex_avr_data.t        = (double) sys->tsim;
  Ahex_avr_data.Afree    = (double) sys->Afree / (KaiA0 * sys->volume);  
  Ihex_avr_data.Atot     = (double) Atot / (KaiA0 * sys->volume);    
  Ihex_avr_data.t        = (double) sys->tsim;
  Ihex_avr_data.Afree    = (double) sys->Afree / (KaiA0 * sys->volume);
    
  Ahex_avr_data.p        = (double) Ap/(6*N_hexamers);
  Ahex_avr_data.ACI      = (double) ACIBA/(KaiA0 * sys->volume);
  Ahex_avr_data.BCI      = (double) ACIB/(KaiC0 * sys->volume);    
  Ahex_avr_data.ACII     = (double) AACII/(KaiA0 * sys->volume);
  Ahex_avr_data.CIATP    = (double) ACIATP/N_hexamers;
  Ahex_avr_data.CIIATP   = (double) ACIIATP/N_hexamers;
  Ahex_avr_data.Ttot     = (double) ATtot/N_hexamers;
  Ahex_avr_data.Stot     = (double) AStot/N_hexamers;  
  Ahex_avr_data.pU       = (double) ApU/N_hexamers;
  Ahex_avr_data.pT       = (double) ApT/N_hexamers;
  Ahex_avr_data.pD       = (double) ApD/N_hexamers;
  Ahex_avr_data.pS       = (double) ApS/N_hexamers;

  Ahex_avr_data.dGconfADP = (double) AdGconfADP/N_hexamers;
  Ahex_avr_data.kADPoff = (double) AADPoff/N_hexamers;  
  Ahex_avr_data.dGACIIbind = (double) AdGACIIbind/N_hexamers;
  Ahex_avr_data.CIKaiA_bound = (double) ACIKaiA_bound / N_hexamers;
  Ahex_avr_data.CIKaiB_bound = (double) ACIKaiB_bound / N_hexamers;
  Ahex_avr_data.CIKidA_bound = (double) ACIKidA_bound / N_hexamers;
  Ahex_avr_data.BCI_any = (double) ACIB_any / N_hexamers;
  Ahex_avr_data.B_active = sys->B_active + sys->B1_active;
  Ahex_avr_data.KaiBKidA = sys->KaiBKidA + sys->KaiB1KidA;
  Ahex_avr_data.n_CIKaiB_on = sys->n_CIKaiB_on;
  Ahex_avr_data.prop_CIBon = A_prop_CIBon;
  Ahex_avr_data.prop_CIAon = A_prop_CIAon;
  Ahex_avr_data.n_max_CIKaiB_bound = A_n_max_CIKaiB_bound;
  Ahex_avr_data.n_B_rebind = sys->n_B_rebind;

  Ihex_avr_data.p        = (double) Ip/(6*N_hexamers);
  Ihex_avr_data.ACI      = (double) ICIBA/(KaiA0 * sys->volume);
  Ihex_avr_data.BCI      = (double) ICIB/(KaiC0 * sys->volume);
  Ihex_avr_data.ACII     = (double) IACII/(KaiA0 * sys->volume);
  Ihex_avr_data.CIATP    = (double) ICIATP/N_hexamers;
  Ihex_avr_data.CIIATP   = (double) ICIIATP/N_hexamers;
  Ihex_avr_data.Ttot     = (double) ITtot/N_hexamers;
  Ihex_avr_data.Stot     = (double) IStot/N_hexamers;  
  Ihex_avr_data.pU       = (double) IpU/N_hexamers;
  Ihex_avr_data.pT       = (double) IpT/N_hexamers;
  Ihex_avr_data.pD       = (double) IpD/N_hexamers;
  Ihex_avr_data.pS       = (double) IpS/N_hexamers;
    
  Ihex_avr_data.dGconfADP = (double) IdGconfADP/N_hexamers;   
  Ihex_avr_data.kADPoff = (double) IADPoff/N_hexamers;
  Ihex_avr_data.dGACIIbind = (double) IdGACIIbind/N_hexamers;
  Ihex_avr_data.CIKaiA_bound = (double) ICIKaiA_bound / N_hexamers;
  Ihex_avr_data.CIKaiB_bound = (double) ICIKaiB_bound / N_hexamers;
  Ihex_avr_data.CIKidA_bound = (double) ICIKidA_bound / N_hexamers;
  Ihex_avr_data.BCI_any = (double) ICIB_any / N_hexamers;
  Ihex_avr_data.B_active = sys->B_active + sys->B1_active;
  Ihex_avr_data.KaiBKidA = sys->KaiBKidA + sys->KaiB1KidA;
  Ihex_avr_data.n_CIKaiB_on = sys->n_CIKaiB_on;
  Ihex_avr_data.prop_CIBon = I_prop_CIBon;
  Ihex_avr_data.prop_CIAon = I_prop_CIAon;
  Ihex_avr_data.n_max_CIKaiB_bound = I_n_max_CIKaiB_bound;
  Ihex_avr_data.n_B_rebind = sys->n_B_rebind;


   

  //Put data in array for writing to file.
  sys->Aoutput_data[sys->sample_cntr] = Ahex_avr_data;
  sys->Ioutput_data[sys->sample_cntr] = Ihex_avr_data;  

  //Calculate some cummulatives for printing to screen.
  HexamerAvr hex_avr;

  hex_avr.Afree = (double) sys->Afree / (KaiA0 * sys->volume);
  hex_avr.Atot = Atot / (KaiA0 * sys->volume);

  hex_avr.active = (double)Ahex_avr_data.active / N_hexamers;
  hex_avr.p = (Ahex_avr_data.p + Ihex_avr_data.p);
  hex_avr.ACI = Ahex_avr_data.ACI + Ihex_avr_data.ACI;
  hex_avr.ACII = Ahex_avr_data.ACII + Ihex_avr_data.ACII;
  hex_avr.CIATP = Ahex_avr_data.CIATP + Ihex_avr_data.CIATP;
  hex_avr.CIIATP = Ahex_avr_data.CIIATP + Ihex_avr_data.CIIATP;
  hex_avr.pU = Ahex_avr_data.pU + Ihex_avr_data.pU;
  hex_avr.pT = Ahex_avr_data.pT + Ihex_avr_data.pT;
  hex_avr.pD = Ahex_avr_data.pD + Ihex_avr_data.pD;
  hex_avr.pS = Ahex_avr_data.pS + Ihex_avr_data.pS;
  hex_avr.CIKaiA_bound = Ahex_avr_data.CIKaiA_bound + Ihex_avr_data.CIKaiA_bound;
  hex_avr.CIKaiB_bound = Ahex_avr_data.CIKaiB_bound + Ihex_avr_data.CIKaiB_bound;
  
  /* Show #Conformational changes */
  //cout << sys->tsim << " " << sys->FlipCntrFw << " " << sys->FlipCntrBw << endl;
  /* Calculate standard deviations in CI-ATP
  cout << sys->tsim << " " 
       << sqrt((double)(ACIATPstd - ACIATP*ACIATP/active)/active) << " "
       << sqrt((double)(ICIATPstd - ICIATP*ICIATP/inactive)/inactive) << "  " 
       << sqrt(((ACIATPstd + ICIATPstd) - (ACIATP+ICIATP)*(ACIATP+ICIATP)/N_hexamers)/N_hexamers) 
       << endl; 
   */
  
    
  /* Count CI-ATP states in ensamble, destinguishing Active and Inactive states 
  int ACIATPlist[7], ICIATPlist[7];
  for(int j(0); j<7; j++)
  {
    ACIATPlist[j] = 0;
    ICIATPlist[j] = 0;
  }
  
  for(int i(0); i<sys->N_hexamers; i++)
  {
    if( hexamers[i].get_active() )
    {
      ACIATPlist[hexamers[i].get_hex_CIATP()] += 1;
    }
    else
    {  
      ICIATPlist[hexamers[i].get_hex_CIATP()] += 1;
    }
  }

  cout << sys->tsim << " "; 
  for(int j(0); j<7; j++)
  {
    cout << (double) (ACIATPlist[j]+ICIATPlist[j])/sys->N_hexamers << " ";
  }  
  for(int j(0); j<7; j++)
  {
    cout << (double) ACIATPlist[j]/sys->N_hexamers << " ";
  }
  for(int j(0); j<7; j++)
  {
    cout << (double) ICIATPlist[j]/sys->N_hexamers << " ";
  }
  
  cout << endl;   
  */

  /*
  if(sys->tsim > 47.3 && sys->tsim < 49.6)
  {
    cerr << "acces main power grid" << endl;
    for(int i(0); i<sys->N_hexamers; i++)
    {
      cout << hexamers[i].get_hex_pS() << endl;
    }
  }
  */
  
  if(inactive == 0) inactive = 1;
  if(active == 0) active = 1;
  
  //FILE *Hfp = fopen( "Hexamers.dat", "a" );
  //fprintf( Hfp, "\n");
  /*Print each hexamer state as integer. 
  0 Active - No CII*A
  1 Active - CII*A
  2 Inactive - CII*A - CI*B < nBseq
  3 Inactive - No CII*A - CI*B < nBseq
  4 Inactive - CI*B = nBseq
  5 Inactive - CI*B*A 
  */

  
  //Print hexamer sequestration state
  /*
  for(int i(0); i<sys->N_hexamers; i++)
  {
    fprintf( Hfp, "%e\t", sys->tsim );
    if( hexamers[i].get_active() && !hexamers[i].get_CIIKaiA_bound() )
      fprintf( Hfp, "%d\t%d\t", hexamers[i].get_index(), 0);
    if( hexamers[i].get_active() && hexamers[i].get_CIIKaiA_bound() )
      fprintf( Hfp, "%d\t%d\t", hexamers[i].get_index(), 1);
    if( !hexamers[i].get_active() && hexamers[i].get_CIKaiB_bound() >= 0 && hexamers[i].get_CIKaiB_bound() < reaction_consts->nBseq && hexamers[i].get_CIIKaiA_bound() )
      fprintf( Hfp, "%d\t%d\t", hexamers[i].get_index(), 2);    
    if( !hexamers[i].get_active() && hexamers[i].get_CIKaiB_bound() >= 0 && hexamers[i].get_CIKaiB_bound() < reaction_consts->nBseq && !hexamers[i].get_CIIKaiA_bound() )
      fprintf( Hfp, "%d\t%d\t", hexamers[i].get_index(), 3);
    if( !hexamers[i].get_active() && hexamers[i].get_CIKaiB_bound() == reaction_consts->nBseq && hexamers[i].get_CIKaiA_bound() == 0)
      fprintf( Hfp, "%d\t%d\t", hexamers[i].get_index(), 4);
    if( !hexamers[i].get_active() && hexamers[i].get_CIKaiB_bound() == reaction_consts->nBseq && hexamers[i].get_CIKaiA_bound() > 0)
      fprintf( Hfp, "%d\t%d\t", hexamers[i].get_index(), 5);
    fprintf( Hfp, "\n");
  }
  */
  
  
  //Print hexamer CI-ATP number.
  /*
  for(int i(0); i<sys->N_hexamers; i++)
  {
    fprintf(Hfp, "%e\t%d\t%d\n", sys->tsim, hexamers[i].get_index(), hexamers[i].get_hex_CIATP());
  }
  */
  
  
  // Get Hexamer S-T phosphorylation state.
  /*
  for(int i(0); i<sys->N_hexamers; i++)
  {
    fprintf(Hfp, "%e\t%d\t%d\n", sys->tsim, hexamers[i].get_index(), hexamers[i].get_hex_S()-hexamers[i].get_hex_T());
  }
  */

  // Get Hexamer B sequestration.
  /*
  for(int i(0); i<sys->N_hexamers; i++)
  {
    if( hexamers[i].get_CIKaiB_bound() == 6 && hexamers[i].get_CIKaiA_bound() > 0 )
      fprintf(Hfp, "%e\t%d\t%d\n", sys->tsim, hexamers[i].get_index(), 10);
    else
      fprintf(Hfp, "%e\t%d\t%d\n", sys->tsim, hexamers[i].get_index(), hexamers[i].get_CIKaiB_bound());
  }  
  */
  
  
  //Print out all hexamer state variables.
  /*
  for(int i(0); i<sys->N_hexamers; i++)
  {
    fprintf( Hfp, "%e\t", sys->tsim ); //1
    fprintf( Hfp, "%d\t", hexamers[i].get_active()); //2
    fprintf( Hfp, "%d\t", hexamers[i].get_CIIKaiA_bound() ); //3
    fprintf( Hfp, "%d\t", hexamers[i].get_CIKaiA_bound() ); //4
    fprintf( Hfp, "%d\t", hexamers[i].get_CIKaiB_bound() ); //5
    fprintf( Hfp, "%d\t", hexamers[i].get_hex_CIATP() ); //6
    fprintf( Hfp, "%d\t", hexamers[i].get_hex_CIIATP() ); //7
    fprintf( Hfp, "%d\t", hexamers[i].get_hex_pU() ); //8
    fprintf( Hfp, "%d\t", hexamers[i].get_hex_pT() ); //9
    fprintf( Hfp, "%d\t", hexamers[i].get_hex_pD() ); //10
    fprintf( Hfp, "%d\t", hexamers[i].get_hex_pS() ); //11   
    fprintf( Hfp, "\n");
  } */
  //fclose(Hfp);
  return hex_avr;  
}


//Allocate data for writing datparam_name.
void allocate_data_array(SystemVariables *sys)
{
  /* set-up data arrays for timetraces */
  sys->sample_cnt = (int) ((sys->tend - sys->tequ) / sys->t_sample_incr + 0.5);
  //Number of timetraces;
  sys->Aoutput_data = (HexamerAvr*) calloc( sys->sample_cnt ,sizeof(HexamerAvr) );
  sys->Ioutput_data = (HexamerAvr*) calloc( sys->sample_cnt ,sizeof(HexamerAvr) );    
}


/* Write timetraces to file
   3 files are produced: Averages for Active state hexamers
                         Averages for Inactive state hexamers
                         Averages for All hexamers 
*/
void write_outputfile(SystemVariables *sys)
{

  
  char Tfilename[128]; memset(&Tfilename[0], 0, sizeof(Tfilename));
  char Afilename[128]; memset(&Afilename[0], 0, sizeof(Afilename));
  char Ifilename[128]; memset(&Ifilename[0], 0, sizeof(Ifilename));
 
  strcat(Tfilename, sys->output_filename);  
  strcat(Tfilename, ".dat");  

  strcpy(Afilename,"A");  
  strcat(Afilename, Tfilename);
  strcpy(Ifilename,"I");
  strcat(Ifilename, Tfilename);
  
  FILE *Afp = fopen( Afilename, "w" );
  FILE *Ifp = fopen( Ifilename, "w" );
  FILE *Tfp = fopen( Tfilename, "w" );      
  
  cerr << "Writing output of Active hexamers to file: " << Afilename << endl;

  // Write file for Active hexamers.
  fprintf(Afp,"#time\tphos_frac\tAfree\tACI\tACII\tAtot\tCIATP\tCIIATP\tTtot\tStot\tpU\tpT\tpD\tpS\tCIATPcons\tdCIATPcons\tCIIATPcons\tdCIIATPcons\tdGconfADP\tdGACIIbind\tkADPoff\tactive\tCIBB\n");  

  // iterate through all samples
  for(int j(0); j < sys->sample_cnt; j++ ) 
  {
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].t );       //1
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].p );       //2
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].Afree );   //3
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].ACI );     //4
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].ACII );    //5
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].Atot );    //6
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].CIATP );   //7                      
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].CIIATP );  //8
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].Ttot );    //9
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].Stot );    //10
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].pU );      //11
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].pT );      //12
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].pD );      //13
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].pS );      //14
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].CIATPcons ); //15
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].dCIATPcons );//16                      
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].CIIATPcons ); //17
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].dCIIATPcons );//18    
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].dGconfADP ); //19     
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].dGACIIbind ); //20
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].kADPoff ); //21   
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].active/sys->N_hexamers );  //22
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].BCI );  //23
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].CIKaiA_bound); //24
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].CIKaiB_bound); //25
    fprintf( Afp, "%e\n", sys->Aoutput_data[j].BCI_any);
  }
    
  // close file
  fclose(Afp);

  cerr << "Writing output of Inactive hexamers to file: " << Ifilename << endl;
 
  // Write file for Inactive hexamers.
  fprintf(Ifp,"#time\tphos_frac\tAfree\tACI\tACII\tAtot\tCIATP\tCIIATP\tTtot\tStot\tpU\tpT\tpD\tpS\tCIATPcons\tdCIATPcons\tCIIATPcons\tdCIIATPcons\tdGconfADP\tdGACIIbind\tkADPoff\tinactive\tCIBB\n");

  // iterate through all samples
  for(int j(0); j < sys->sample_cnt; j++ ) 
  {
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].t );       //1
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].p );       //2
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].Afree );   //3
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].ACI );     //4
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].ACII );    //5
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].Atot );    //6
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].CIATP );   //7                      
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].CIIATP );  //8
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].Ttot );    //9
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].Stot );    //10
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].pU );      //11
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].pT );      //12
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].pD );      //13
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].pS );      //14
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].CIATPcons ); //15
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].dCIATPcons );//16                      
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].CIIATPcons ); //17
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].dCIIATPcons );//18                     
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].dGconfADP ); //19     
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].dGACIIbind ); //20
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].kADPoff ); //21
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].active/sys->N_hexamers ); //22            
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].BCI );  //23    
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].CIKaiA_bound); //24
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].CIKaiB_bound); //25
    fprintf( Ifp, "%e\n", sys->Ioutput_data[j].BCI_any);
  }
    
  // close file
  fclose(Ifp);
  
  cerr << "Writing output of all hexamers to file: " << Tfilename << endl;
  
  // Write file for All hexamers.
    fprintf(Tfp,"#time\tphos_frac\tAfree\tACI\tACII\tAtot\tCIATP\tCIIATP\tTtot\tStot\tpU\tpT\tpD\tpS\tCIATPcons\tdCIATPcons\tCIIATPcons\tdCIIATPcons\tdGconfADP\tdGACIIbind\tkADPoff\tactive\tCIBB\n");  

  // iterate through all samples
  for(int j(0); j < sys->sample_cnt; j++ ) 
  { 
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].t );       //1
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].p + sys->Ioutput_data[j].p );       //2
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].Afree );   //3
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].ACI + sys->Ioutput_data[j].ACI );     //4
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].ACII + sys->Ioutput_data[j].ACII );    //5
    fprintf( Tfp, "%e\t", sys->KaiA0 );    //6
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].CIATP + sys->Ioutput_data[j].CIATP ) ; //7                      
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].CIIATP + sys->Ioutput_data[j].CIIATP ); //8
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].Ttot + sys->Ioutput_data[j].Ttot );    //9
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].Stot + sys->Ioutput_data[j].Stot );    //10
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].pU + sys->Ioutput_data[j].pU );      //11
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].pT + sys->Ioutput_data[j].pT );      //12
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].pD + sys->Ioutput_data[j].pD );      //13
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].pS + sys->Ioutput_data[j].pS );      //14
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].CIATPcons ); //15
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].dCIATPcons );//16                      
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].CIIATPcons ); //17
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].dCIIATPcons );//18                     
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].dGconfADP + sys->Ioutput_data[j].dGconfADP ); //19     
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].dGACIIbind + sys->Ioutput_data[j].dGACIIbind ); //20
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].kADPoff + sys->Ioutput_data[j].kADPoff ); //21
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].active/sys->N_hexamers );  //22   
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].BCI + sys->Ioutput_data[j].BCI ); //23
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].CIKaiA_bound + sys->Ioutput_data[j].CIKaiA_bound); //24
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].CIKaiB_bound + sys->Ioutput_data[j].CIKaiB_bound); //25
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].BCI_any + sys->Ioutput_data[j].BCI_any);
    fprintf( Tfp, "%d\t", sys->Aoutput_data[j].B_active);
    fprintf( Tfp, "%d\t", sys->Aoutput_data[j].KaiBKidA);
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].CIKidA_bound + sys->Ioutput_data[j].CIKidA_bound);
    fprintf( Tfp, "%d\t", sys->Aoutput_data[j].n_CIKaiB_on);
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].prop_CIBon + sys->Ioutput_data[j].prop_CIBon);
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].prop_CIAon + sys->Ioutput_data[j].prop_CIAon);
    fprintf( Tfp, "%d\t", sys->Ioutput_data[j].n_B_rebind);
    fprintf( Tfp, "%d\n", sys->Aoutput_data[j].n_max_CIKaiB_bound + sys->Ioutput_data[j].n_max_CIKaiB_bound);
  }  
  
  // close file
  fclose(Tfp);  
}


/* Load system variables from file defined in command line argument. */
void initialize_system_vars(SystemVariables *sys)
{ 
  ifstream infile(sys->param_filename);
  string line, param_name, param_string, temp;
  double param_value;
  bool paramset(0);
  
  //Set parameter default values (can be overwritten by .par file)
  sys->volume                    = 1200;
  sys->KaiA0                     = 0.6;
  sys->KaiC0                     = 0.6;
  sys->tend                      = 120;
  sys->tequ                      = 0;
  sys->ATPfrac                   = 1.0;
  sys->t_sample_incr             = 0.1;
  sys->tincu = 0;
  sys->t_hex_equil = 0;
  sys->start_phosphorylated      = 0;
  sys->rnd_seed                  = 42;
  
  //cerr << endl << "SYSTEM VARIABLES USED IN THIS SIMULATION: " << endl;
  while (getline(infile, line))
  {
    //remove unnecessary spaces in line.
    line = delUnnecessary( line ); 
    istringstream iss(line);
    paramset = 0;
    
    // Check if line has formatting: '<string> <number>'
    if (iss >> param_name >> param_value) 
    {        
      /* Loading system variables */   
      if( param_name.compare("volume") == 0 )
      {
        //Convert cubic micron to particles per micromolar.
        //Roughly 600 particles go in one cubcic micron of 1 micromolar concentration.
        sys->volume = param_value * 600;
        paramset = 1;
      }  

      if( param_name.compare("KaiA0") == 0 )
      {
        sys->KaiA0 = param_value;
        paramset = 1;
      }  
     
      if( param_name.compare("KaiC0") == 0 )
      {
        sys->KaiC0 = param_value;
        paramset = 1;
      } 

      if( param_name.compare("KaiB0") == 0 )
      {
 	sys->KaiB0 = param_value;
	paramset = 1;
      }

      if( param_name.compare("fraction_fs") == 0)
      {
	sys->fraction_fs = param_value;
	paramset = 1;
      }
      
      if( param_name.compare("KidA0") == 0 )
      {
	sys->KidA0 = param_value;
	paramset = 1;
      } 
      if( param_name.compare("tend") == 0 )
      {
        sys->tend = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("tequ") == 0 )
      {
        sys->tequ = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("tincu") == 0 )
      {
        sys->tincu = param_value;
        paramset = 1;
      }  
     
      if( param_name.compare("t_hex_equil") ==0 )
      {
        sys->t_hex_equil = param_value;
	paramset = 1;
      }

      if( param_name.compare("ATPfrac") == 0 )
      {
        sys->ATPfrac = param_value;
        paramset = 1;
      }        
      
      if( param_name.compare("t_sample_incr") == 0 )
      {
        sys->t_sample_incr = param_value;
        paramset = 1;
      }       
      
      if( param_name.compare("start_phosphorylated") == 0 )
      {
        sys->start_phosphorylated = param_value;
        paramset = 1;
      }                  
      
      if( param_name.compare("rnd_seed") == 0 )
      {
        sys->rnd_seed = param_value;
        paramset = 1;
      }               
             
      if(!paramset)
      {
        //cerr << param_name << " - is not a valid parameter name." << endl;
      }
      else
      {
        //Uncomment to see system values at program start.
        cerr << param_name << " = " << param_value << endl;
        //fprintf( Pfp, "%s\t%e\n", param_name.c_str(), param_value );     
      }
          
    }
    else
    {
      // Check if line has formatting: '<string> <string>'
      iss.clear();
      iss.str( line );
      if (iss >> param_name >> param_string)
      {          
         if( param_name.compare("output_filename") == 0 )
         {
           strcpy(sys->output_filename, param_string.c_str());
           paramset = 1;
         }
                     
         if(!paramset)
         {
           //cerr << param_name << " - is not a valid parameter name." << endl;
         }
         else
         {
           //uncomment to see system values at program start.
           cerr << param_name << " = " << param_string << endl;
           //fprintf( Pfp, "%s\t%s\n", param_name.c_str(), param_string.c_str() );     
         }  
      }
      else
      {       
        //Check for comment and empty lines
        if(!(param_name.compare(0,1,"#") == 0) && !line.empty())
        {
           cerr << "Error - Line with incorrect formatting: -" << param_name << 
           "-" << param_string << "." << endl; 
           exit(1); 
        }   
        continue;
      }
    } //if not variable_name variable_value
  }//while 

  //Set system variables
  int B_total = sys->KaiB0 * sys->volume;
  sys->tsim = 0.0;
  sys->N_hexamers = sys->KaiC0 * sys->volume;
  sys->Afree = sys->KaiA0 * sys->volume;
  sys->B_active = B_total * sys->fraction_fs;
  sys->B_inactive = B_total - sys->B_active;
  sys->B1_active = 0;
  sys->KaiBKidA = 0;
  sys->KaiB1KidA = 0;
  sys->KidA_free = sys->KidA0 * sys->volume;
  sys->cAfree = sys->Afree / sys->volume;
  sys->CIATPcons = 0; sys->CIIATPcons = 0;
  
  sys->step_cntr = 0;
  sys->sample_cntr = 0;
}


/* Function reads in reaction constants from file. */
ReactionConstants initialize_reaction_consts(SystemVariables *sys)
{
  ReactionConstants reaction_consts;

  ifstream infile(sys->param_filename);
  string line, param_name;
  double param_value;
  string param_string;
  bool paramset(0);

  //cerr << endl << "PARAMETER VALUES USED IN THIS SIMULATION: " << endl;
  while (getline(infile, line))
  {
    //remove unnecessary spaces in line.
    line = delUnnecessary( line );
  
    istringstream iss(line);
    paramset = 0;
    
    // Check if line has formatting: '<string> <number>'
    if (iss >> param_name >> param_value) 
    {
      /* Loading model parameters */              
      if( param_name.compare("kCIhyd") == 0 )
      {
        reaction_consts.kCIhyd = param_value; 
        paramset = 1;
      }            

      if( param_name.compare("kCIADPoff0") == 0 )
      {
        reaction_consts.kCIADPoff0 = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kCIADPoffA") == 0 )
      {
        reaction_consts.kCIADPoffA = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kAkIDoff") == 0 )
      {
        reaction_consts.kAkIDoff = param_value; 
        paramset = 1;
      }

      if( param_name.compare("ddGTDconf") == 0 )
      {
        reaction_consts.ddGTDconf = param_value; 
        paramset = 1;
      }
            
      if( param_name.compare("ddGTDAbind") == 0 )
      {
        reaction_consts.ddGTDAbind = param_value; 
        paramset = 1;
      }
            
            
      if( param_name.compare("kconf0") == 0 )
      {
        reaction_consts.kconf0 = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kBswitch_f") == 0 )
      {
        reaction_consts.kBswitch_f = param_value; 
        paramset = 1;
      }

      if( param_name.compare("kBswitch_r") == 0 )
      {
        reaction_consts.kBswitch_r = param_value; 
        paramset = 1;
      }

      if( param_name.compare("nIADPconfref") == 0 )
      {
        reaction_consts.nIADPconfref = param_value; 
        paramset = 1;
      }

      if( param_name.compare("nIADPAref") == 0 )
      {
        reaction_consts.nIADPAref = param_value; 
        paramset = 1;
      }

      if( param_name.compare("dGIbindB") == 0 )
      {
        reaction_consts.dGIbindB = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("dgACIIU") == 0 )
      {
        reaction_consts.dgACIIU = param_value; 
        paramset = 1;
      }

      if( param_name.compare("dgACIIT") == 0 )
      {
        reaction_consts.dgACIIT = param_value; 
        paramset = 1;
      }

      if( param_name.compare("dgACIID") == 0 )
      {
        reaction_consts.dgACIID = param_value; 
        paramset = 1;
      }

      if( param_name.compare("dgACIIS") == 0 )
      {
        reaction_consts.dgACIIS = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("dGICII") == 0 )
      {
        reaction_consts.dGICII = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("dgACIActU") == 0 )
      {
        reaction_consts.dgACIActU = param_value; 
        paramset = 1;
      }   
      
      if( param_name.compare("dgACIActT") == 0 )
      {
        reaction_consts.dgACIActT = param_value; 
        paramset = 1;
      }  
      
      if( param_name.compare("dgACIActD") == 0 )
      {
        reaction_consts.dgACIActD = param_value; 
        paramset = 1;
      }  
      
      if( param_name.compare("dgACIActS") == 0 )
      {
        reaction_consts.dgACIActS = param_value; 
        paramset = 1;
      }  
      
      if( param_name.compare("dgICIActU") == 0 )
      {
        reaction_consts.dgICIActU = param_value; 
        paramset = 1;
      }   
      
      if( param_name.compare("dgICIActT") == 0 )
      {
        reaction_consts.dgICIActT = param_value; 
        paramset = 1;
      }  
      
      if( param_name.compare("dgICIActD") == 0 )
      {
        reaction_consts.dgICIActD = param_value; 
        paramset = 1;
      }  
      
      if( param_name.compare("dgICIActS") == 0 )
      {
        reaction_consts.dgICIActS = param_value; 
        paramset = 1;
      }                                    

      if( param_name.compare("kICIBGon") == 0 )
      {
        reaction_consts.kICIBGon = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kICIBFSon") == 0 )
      {
        reaction_consts.kICIBFSon = param_value; 
        paramset = 1;
      }

      if( param_name.compare("kICIBoff") == 0 )
      {
        reaction_consts.kICIBoff = param_value; 
        paramset = 1;
      }         

      if( param_name.compare("kACIBGon") == 0 )
      {
        reaction_consts.kACIBGon = param_value; 
        paramset = 1;
      }

      if( param_name.compare("kACIBFSon") == 0 )
      {
        reaction_consts.kACIBFSon = param_value; 
        paramset = 1;
      }

      if( param_name.compare("kACIBoff") == 0 )
      {
        reaction_consts.kACIBoff = param_value; 
        paramset = 1;
      }
                
      if( param_name.compare("kACIAon") == 0 )
      {
        reaction_consts.kACIAon = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kACIAoff") == 0 )
      {
        reaction_consts.kACIAoff = param_value; 
        paramset = 1;
      }

      if( param_name.compare("kICIAon") == 0 )
      {
        reaction_consts.kICIAon = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kICIAoff") == 0 )
      {
        reaction_consts.kICIAoff = param_value; 
        paramset = 1;
      }

      if( param_name.compare("nAseq") == 0 )
      {
        reaction_consts.nAseq = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("nBseq") == 0 )
      {
        reaction_consts.nBseq = param_value; 
        paramset = 1;
      }      
      
      if( param_name.compare("kCIIAon") == 0 )
      {
        reaction_consts.kCIIAon = param_value; 
        paramset = 1;
      }                                                      

      if( param_name.compare("kCIIAoff") == 0 )
      {
        reaction_consts.kCIIAoff = param_value; 
        paramset = 1;
      }                                                      

      if( param_name.compare("kCIIhyd0") == 0 )
      {
        reaction_consts.kCIIhyd0 = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kCIIhydA") == 0 )
      {
        reaction_consts.kCIIhydA = param_value; 
        paramset = 1;
      }                                                             

      if( param_name.compare("kCIInucloff0") == 0 )
      {
        reaction_consts.kCIInucloff0 = param_value; 
        paramset = 1;
      }                                                      

      if( param_name.compare("kCIInucloffA") == 0 )
      {
        reaction_consts.kCIInucloffA = param_value; 
        paramset = 1;
      }  
      
      if( param_name.compare("KATPoKADP") == 0 )
      {
        reaction_consts.KATPoKADP = param_value; 
        paramset = 1;
      }      

      if( param_name.compare("kUT") == 0 )
      {
        reaction_consts.kUT = param_value; 
        paramset = 1;
      }  

      if( param_name.compare("kTU") == 0 )
      {
        reaction_consts.kTU = param_value;
        paramset = 1;
      }  

      if( param_name.compare("kTD") == 0 )
      {
        reaction_consts.kTD = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("kDT") == 0 )
      {
        reaction_consts.kDT = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("kSD") == 0 )
      {
        reaction_consts.kSD = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("kDS") == 0 )
      {
        reaction_consts.kDS = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("kUS") == 0 )
      {
        reaction_consts.kUS = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("kSU") == 0 )
      {
        reaction_consts.kSU = param_value;
        paramset = 1;
      }
      if( param_name.compare("kKidAon") == 0 )
      {
        reaction_consts.kKidAon = param_value;
        paramset = 1;
      }
      if( param_name.compare("kKidAoff") == 0 )
      {
        reaction_consts.kKidAoff = param_value;
        paramset = 1;
      }
      
      if(!paramset)
      {
        //cerr << param_name << " - is not a valid parameter name." << endl;
      }
      else
      {
        //Uncomment to see parameter values at program start.
        cerr << param_name << " = " << param_value << endl;
        //fprintf( Pfp, "%s\t%e\n", param_name.c_str(), param_value );             
      }
          
    }
    else
    {
      // Check if line has formatting: '<string> <string>'
      iss.clear();
      iss.str( line );
      if (iss >> param_name >> param_string)
      {                
         if(!paramset)
         {
           //cerr << param_name << " - is not a valid parameter name." << endl;
         }
         else
         {
           //Uncomment to see parameter values at program start.
           cerr << param_name << " = " << param_string << endl;
           //fprintf( Pfp, "%s\t%s\n", param_name.c_str(), param_string.c_str() );                
         }              
      }
      else
      { 
        //Check for comment and empty lines
        if(!(param_name.compare(0,1,"#") == 0) && !line.empty())
        {
           cerr << "Error - Line with incorrect formatting: " << param_name << endl; 
           exit(1); 
        }   
        continue;
      }
      

    }
  }//while 

  //fclose(Pfp);    
  return reaction_consts;  
}

string delUnnecessary(string &str)
{
    int size = str.length();
    for(int j = 0; j<=size; j++)
    {
        for(int i = 0; i <=j; i++)
        {
            if(str[i] == ' ' && str[i+1] == ' ')
            {
                str.erase(str.begin() + i);
            }
            else if(str[0]== ' ')
            {
                str.erase(str.begin());
            }
            else if(str[i] == '\0' && str[i-1]== ' ')
            {
                str.erase(str.end() - 1);
            }
        }
    }
    return str;
}
