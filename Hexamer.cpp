#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>
#include <iomanip>
#include <random>

#include "Hexamer.hpp"
#include "PropensityContainer.hpp"

using namespace std;

/* Set state variables of hexamer */
void Hexamer::initialize_state( int index, bool active, bool CIIKaiA_bound, int CIKaiA_bound, int CIKaiB_bound, int CIATP_bound, int CIIATP_bound, int pT, int pD, int pS)
{
	this->index = index;
	this->active = active;
	this->CIIKaiA_bound = CIIKaiA_bound;
	this->CIKaiA_bound = CIKaiA_bound;
	this->CIKaiB_bound = CIKaiB_bound;
	this->CIATP_bound = CIATP_bound;
	this->CIIATP_bound = CIIATP_bound;
	this->pT = pT;
	this->pD = pD;
	this->pS = pS;
}

/* Set which monomers are in the hexamer, and vice versa.
   and initialize the monomer states. */
/*void Hexamer::set_sextet(Monomer *monomers, uniform_real_distribution<double> u01, mt19937_64 engine)
{
  bool bernoullip = false;
  if(sys->ATPfrac > u01(engine)) {
	bernoullip = true;
  } 

  for(int i(6*index); i < (index + 1) * 6; i++)
  {
    sextet[i-6*index] = &monomers[i];
    monomers[i].set_hexamer( this );
    monomers[i].set_prop_container( &(this->mon_prop_cont) );
    monomers[i].set_reaction_consts(this->reaction_consts);
    monomers[i].set_sysvars(this->sys);
    
    if(this->sys->start_phosphorylated)
    {   
      double rnd( u01(engine) );
      bool S(0),T(0);      
      
      //different starting fractions of monomer phosphorylation states.
      //double Tfrac(0.12), Sfrac(0.13), Dfrac(0.08);
      //double Tfrac(0.06), Sfrac(0.25), Dfrac(0.10);      
      //double Tfrac(0.07), Sfrac(0.33), Dfrac(0.18);
      //double Tfrac(0.00), Sfrac(0.0), Dfrac(1.0);
      double Tfrac(0.10), Sfrac(0.15), Dfrac(0.40);     
      
      if(rnd < Dfrac) 
      {
        T=1; S=1; //D-state
      }    
      if(rnd > Dfrac && rnd < Dfrac + Tfrac) 
      {
        T=1; S=0; //T-state
      }
      if(rnd > Dfrac + Tfrac && rnd < Dfrac + Tfrac + Sfrac) 
      {
        T=0; S=1; //S-state
      }
           
      monomers[i].initialize_state(i - 6*index,T,S,1, bernoullip);
    }
    else
    {
      monomers[i].initialize_state(i - 6*index,0,0,1, bernoullip);
    }
  }
  
  // Monomer propensities depend on hexamer state,
  //   so they can only be calculated after all monomer states are set.
  
  for(int i(0); i < 6; i++)
  {
    sextet[i]->set_propensities();
  }
}*/

/* Propagate hexamer, returns TRUE if reaction is fired
   that induces a change in external variables */
int Hexamer::propagate(uniform_real_distribution<double> &u01, mt19937_64 &engine)
{
  double rnd( u01(engine) );

  int reaction_idx( this->choose_reaction(rnd) );   
  int ext_change( this->fire_reaction( reaction_idx, u01, engine ) );
   
  /* Update propensities due to state change */
  set_propensities();
   
  //Reactions 0,1,3,4,5 change Afree.   
  return ext_change;
}

/* Give random number rnd [0,1), function returns index of fired reaction.
   Index to reaction mapping is defined in calc_propensities.
*/
int Hexamer::choose_reaction(double rnd)
{
  double qacc(this->prop_list[0]);
  int i(0);
  
  rnd*=SCALE;
  rnd*=hex_prop_cont->get_qtot_el(this->index);
  while( rnd > qacc && i < HEXAMER_N_REACTS)
  {
    qacc += this->prop_list[++i];
  }
    
  return i;
}

/* Given fired reaction channel index, execute state change in hexamer or monomer.
   Returns TRUE if a change in a monomers alters the hexamer propensities. */
int Hexamer::fire_reaction(int reaction_channel, uniform_real_distribution<double> &u01, mt19937_64 &engine)
{
	int ext_change(0);
	/*if(reaction_channel < 8)
	{
	printf("%d\n", reaction_channel);
	}*/
	switch (reaction_channel)
	{ 
		case 0: 
			CIIKaiA_bound = 1; 
			sys->Afree--; 
			sys->cAfree = (double) sys->Afree / sys->volume;
			ext_change = 1;
			break;
		case 1: 
			CIIKaiA_bound = 0; 
			sys->Afree++;
			sys->cAfree = (double) sys->Afree / sys->volume;   
			ext_change = 1;     
			break;
		case 2: 
			CIKaiB_bound++;
			sys->B_active -= 1;
			ext_change = 2;
			break;
		case 3:
			CIKaiB_bound++;
			CIKidA_bound += 1;
			sys->KaiBKidA -= 1;
			ext_change = 2;
			break;
		case 4:
			CIKaiB_bound++;
			sys->B_inactive -= 1;
			ext_change = 2;
			break;
		case 5: 
			CIKaiB_bound--;
			sys->B_active += 1;
			ext_change = 2;
			break;
		case 6:
			CIKaiB_bound -= 1;
			sys->B_inactive += 1;
			ext_change = 2;
			break;
		case 7:
			CIKaiB_bound -= 1;
			CIKidA_bound -= 1;
			sys->KaiBKidA += 1;
			ext_change = 2;
			break;
		case 8:
			CIKidA_bound += 1;
			sys->KidA_free -= 1;
			ext_change = 3;
			break;
		case 9:
			CIKidA_bound -= 1;
			sys->KidA_free += 1;
			ext_change = 3;
			break;

		case 10:
			CIKaiA_bound++; 
			sys->Afree--;
			sys->cAfree = (double) sys->Afree / sys->volume;
			ext_change = 1;    
			break;
		case 11: 
			CIKaiA_bound--; 
			sys->Afree++;
			sys->cAfree = (double) sys->Afree / sys->volume; 
			ext_change = 1;       
			break;
		case 12: 
			active = 0;
			break;
		case 13: 
			active = 1;
			break;
		case 14:
			sys->CIATPcons += 1;
			CIATP_bound -= 1;
			break;
		case 15:
			CIATP_bound += 1;
			break;
		case 16:
			sys->CIIATPcons += 1;
			CIIATP_bound -= 1;
			break;
		case 17:
			CIIATP_bound += 1;
			break;
		case 18:
			pT += 1;
			CIIATP_bound -= 1;
			break;
		case 19:
			pT -= 1;
			CIIATP_bound += 1;
			break;
		case 20:
			pT -=1;
			pD += 1;
			CIIATP_bound -= 1;
			break;
		case 21:
			pT += 1;
			pD -= 1;
			CIIATP_bound += 1;
			break;
		case 22:
			pD -= 1;
			pS += 1;
			CIIATP_bound += 1;
			break;
		case 23:
			pD += 1;
			pS -= 1;
			CIIATP_bound -= 1;
			break;
		case 24:
			pS -= 1;
			CIIATP_bound += 1;
			break;
		case 25:
			pS += 1;
			CIIATP_bound -= 1;
			break;
	}

	return ext_change;
}

/* Add propensities of all reaction channels, and return total. */
void Hexamer::update_prop_cont()
{
  double hex_qint( prop_list[1] + prop_list[5] + prop_list[6] + prop_list[7] + prop_list[9]);
  double hex_qA( prop_list[0] + prop_list[10]);
  double hex_qB(prop_list[2] + prop_list[3] + prop_list[4]);
  double hex_qKidA(prop_list[8]);
  
  for(int i(11); i<HEXAMER_N_REACTS; i++)
  {
    hex_qint += prop_list[i];
  }

  hex_prop_cont->update_qAll_el( hex_qint, hex_qA, hex_qB, hex_qKidA, this->index);
}

/* Calculate all reaction propensities in this hexamer
   that depend on the external state variable Afree. */
double Hexamer::update_KaiA_prop()
{
  // CII: CII + KaiA -> CII.KaiA
  prop_list[0] = prefactor_r0 * sys->cAfree;
   
  // CI: CI.KaiB + KaiA -> CI.KaiB.KaiA
  prop_list[10] = prefactor_r4 * sys->cAfree; 

  return prop_list[0] + prop_list[10];
}

double Hexamer::update_KaiB_prop()
{
  // CI: CI + KaiB_active -> CI.KaiB
	prop_list[2] = ( CIKaiB_bound < reaction_consts->nBseq && sys->tsim > sys->tincu + sys->t_hex_equil) ? kCIBFSon() * sys->B_active / sys->volume: 0.0;
	prop_list[3] = ( CIKaiB_bound < reaction_consts->nBseq && sys->tsim > sys->tincu + sys->t_hex_equil) ? kCIBFSon() * sys->KaiBKidA / sys->volume: 0.0;
	prop_list[4] = ( CIKaiB_bound < reaction_consts->nBseq && sys->tsim > sys->tincu + sys->t_hex_equil) ? kCIBGon() * sys->B_inactive / sys->volume: 0.0;

	return prop_list[2] + prop_list[3] + prop_list[4];
}

double Hexamer::update_KidA_prop()
{
	prop_list[8] = max((double) CIKaiB_bound - CIKidA_bound, 0.) * sys->KidA_free * reaction_consts->kKidAon;
	return prop_list[8];
}

/* Calculate all reaction propensities possible in this Hexamer.
   Hexamer has HEXAMER_N_REACTS + 6 reaction channels, 
   which propensities are stored in prop_list. */
void Hexamer::set_propensities()
{  
	double kCIIoffnucl = kCIInucloff();
	int pU = 6 - pT - pD - pS;
	double kCIIhyd;

	if(CIIKaiA_bound) {

		kCIIhyd = reaction_consts->kCIIhydA;
	}
	else {

		kCIIhyd = reaction_consts->kCIIhyd0;
	}

	/* Calculate reactionchannel propensities */   
	// CII: CII + A <-> ACII
	prefactor_r0 = (1-CIIKaiA_bound) * kCIIAon();
	prop_list[0] = prefactor_r0 * sys->cAfree;
	prop_list[1] = CIIKaiA_bound * kCIIAoff();

	// CI: CI + B <-> CI*B ([KaiB] dependent)
	prop_list[2] = ( CIKaiB_bound < reaction_consts->nBseq && sys->tsim > sys->tincu + sys->t_hex_equil) ? kCIBFSon() * sys->B_active / sys->volume: 0.0;
	prop_list[3] = ( CIKaiB_bound < reaction_consts->nBseq && sys->tsim > sys->tincu + sys->t_hex_equil) ? kCIBFSon() * sys->KaiBKidA / sys->volume: 0.0;
	prop_list[4] = ( CIKaiB_bound < reaction_consts->nBseq && sys->tsim > sys->tincu + sys->t_hex_equil) ? kCIBGon() * sys->B_inactive / sys->volume: 0.0;
	prop_list[5] = (CIKaiB_bound - CIKidA_bound) * kCIBFSoff();
	prop_list[6] = (CIKaiB_bound - CIKidA_bound) * kCIBGoff();
	prop_list[7] = CIKidA_bound * kCIBFSoff();

	//CI: CI*B + KidA <-> CI*B*KidA
	prop_list[8] = max((double) CIKaiB_bound - CIKidA_bound, 0.) * sys->KidA_free * reaction_consts->kKidAon;
	prop_list[9] = CIKidA_bound * reaction_consts->kKidAoff;

	// CI: CI + A <-> CI*A
	prefactor_r4 = ( CIKaiB_bound == reaction_consts->nBseq ) ? max(reaction_consts->nAseq - CIKaiA_bound, 0.) * kCIAon() : 0.0;
	prop_list[10] = prefactor_r4 * sys->cAfree;
	prop_list[11] = CIKaiA_bound * kCIAoff();

	// CI: active <-> inactive
	prop_list[12] = (double) active * kconf_f();
	prop_list[13] = (double) (1-active) * kconf_b();

	// CI: CI-ATP <-> CI-ADP
	prop_list[14] = CIATP_bound * reaction_consts->kCIhyd;
	prop_list[15] = (6 - CIATP_bound) * kCIADPoff();

	// CII: ATP <-> ADP
	prop_list[16] = CIIATP_bound * (kCIIhyd + (1 - sys->ATPfrac) * kCIIoffnucl);
	prop_list[17] = (6 - CIIATP_bound) * sys->ATPfrac * kCIIoffnucl / reaction_consts->KATPoKADP;
	
	double dgUTbind( reaction_consts->dgACIIT - reaction_consts->dgACIIU );
	double dgTDbind( reaction_consts->dgACIID - reaction_consts->dgACIIT );
	double dgUSbind( reaction_consts->dgACIIS - reaction_consts->dgACIIU );
	double dgSDbind( reaction_consts->dgACIID - reaction_consts->dgACIIS );

	// CII: U <-> T
	prop_list[18] = pU * CIIATP_bound * reaction_consts->kUT * exp(-CIIKaiA_bound * dgUTbind / 2) / 6;
	prop_list[19] = pT * (6 - CIIATP_bound) * reaction_consts->kTU * exp(CIIKaiA_bound * dgUTbind / 2) / 6;

	// CII: T <-> D
	prop_list[20] = pT * CIIATP_bound * reaction_consts->kTD * exp(-CIIKaiA_bound * dgTDbind / 2) / 6;
	prop_list[21] = pD * (6 - CIIATP_bound) * reaction_consts->kDT * exp(CIIKaiA_bound * dgTDbind / 2) / 6;

	// CII: D <-> S
	prop_list[22] = pD * (6 - CIIATP_bound) * reaction_consts->kDS * exp(CIIKaiA_bound * dgSDbind / 2) / 6;
	prop_list[23] = pS * CIIATP_bound * reaction_consts->kSD * exp(-CIIKaiA_bound * dgSDbind / 2) / 6;

	// CII: S <-> U
	prop_list[24] = pS * (6 - CIIATP_bound) * reaction_consts->kSU * exp(CIIKaiA_bound * dgUSbind / 2) / 6;
	prop_list[25] = pU * CIIATP_bound * reaction_consts->kUS * exp(-CIIKaiA_bound * dgUSbind / 2) / 6;

	update_prop_cont();
}

/*CII nucleotide exchange rate, depends on KaiA bound state.*/
double Hexamer::kCIInucloff()
{
  double kCIInucloff0( reaction_consts->kCIInucloff0 );
  double kCIInucloffA( reaction_consts->kCIInucloffA );
  double CIIAbound( get_CIIKaiA_bound() );
  
  return (double) kCIInucloff0 * (1. - CIIAbound) + kCIInucloffA * CIIAbound;
}


/* Conformation switch rate from Active to Inactive state. */
double Hexamer::kconf_f()
{
  double beta_energy( dGconfADP() + dGconfKaiA() + dGconfKaiB() ); 

  return reaction_consts->kconf0 * exp(-0.5 * beta_energy);
}

double Hexamer::kconf_b()
{
  double beta_energy( dGconfADP() + dGconfKaiA() + dGconfKaiB() ); 

  return reaction_consts->kconf0 * exp(0.5 * beta_energy);
}


/* Free energy difference between Active and Inactive state,
   dependent on ADP in CI domain. */
double Hexamer::dGconfADP()
{
  int    nIADPtot( 6 - this->get_hex_CIATP() );
   
  return reaction_consts->ddGTDconf * (reaction_consts->nIADPconfref - nIADPtot);
}


/* Free energy difference between Active and Inactive state,
   dependent on KaiB bound to the CI domain. */
double Hexamer::dGconfKaiB()
{
  const double KdAct( reaction_consts->kACIBoff/reaction_consts->kACIBFSon );
  const double KdIna( reaction_consts->kICIBoff/reaction_consts->kICIBFSon );  
    
  double beta_energy( (double)CIKaiB_bound * log( KdIna/KdAct ) );
  
  return beta_energy;
  //return 0;
}


/* Free energy difference between Active and Inactive state,
   dependent on KaiA bound to the CI and/or CI domain. */
double Hexamer::dGconfKaiA()
{
  const double KdAct( reaction_consts->kACIAoff/reaction_consts->kACIAon );
  const double KdIna( reaction_consts->kICIAoff/reaction_consts->kICIAon );  
    
  double beta_energy( (double)CIKaiA_bound * log( KdIna/KdAct )  );

  beta_energy += CIIKaiA_bound * reaction_consts->dGICII;
  
  return beta_energy;
}


/* Binding free energy difference between CII and CII.A, in Active and Inactive state */
double Hexamer::dGACIIbind()
{
  double dgnU((double) get_hex_pU() * reaction_consts->dgACIIU);
  double dgnT((double) get_hex_pT() * reaction_consts->dgACIIT);
  double dgnD((double) get_hex_pD() * reaction_consts->dgACIID);
  double dgnS((double) get_hex_pS() * reaction_consts->dgACIIS);
  
  double beta_energy( dgnU + dgnT + dgnD + dgnS + (1-active) * reaction_consts->dGICII );

  return beta_energy;
}


/* Activation energy between ADP bound state ADP unbound state,
   in Active state dependent, on CII phos state of whole hexamer. */
double Hexamer::dGACIActADPoff()
{
  double dgnU((double) get_hex_pU() * reaction_consts->dgACIActU);
  double dgnT((double) get_hex_pT() * reaction_consts->dgACIActT);
  double dgnD((double) get_hex_pD() * reaction_consts->dgACIActD);
  double dgnS((double) get_hex_pS() * reaction_consts->dgACIActS);
  
  return dgnU + dgnT + dgnD + dgnS;
}


/* Activation energy between ADP bound state ADP unbound state,
   in Inactive state, dependent on CII phos state of whole hexamer. */
double Hexamer::dGICIActADPoff()
{
  double dgnU((double) get_hex_pU() * reaction_consts->dgICIActU);
  double dgnT((double) get_hex_pT() * reaction_consts->dgICIActT);
  double dgnD((double) get_hex_pD() * reaction_consts->dgICIActD);
  double dgnS((double) get_hex_pS() * reaction_consts->dgICIActS);
  
  return dgnU + dgnT + dgnD + dgnS;
}


double Hexamer::kCIBFSon()
{  
  return (double) active * reaction_consts->kACIBFSon + (1-active) * reaction_consts->kICIBFSon;
}

double Hexamer::kCIBGon()
{  
  return (double) active * reaction_consts->kACIBGon + (1-active) * reaction_consts->kICIBGon;
}

double Hexamer::kCIBFSoff()
{
  return (double) active * reaction_consts->kACIBoff + (1-active) * reaction_consts->kICIBoff;
}

double Hexamer::kCIBGoff()
{
	return kCIBGon() * kCIBFSoff() * reaction_consts->kBswitch_r / (reaction_consts->kBswitch_f * kCIBFSon());
}

/* KaiA on-rate to CI domain, Active/Inactive -conformation dependent. */
double Hexamer::kCIAon()
{
  return (double) active * reaction_consts->kACIAon + (1-active) * reaction_consts->kICIAon;
}


/* KaiA off-rate to CI domain, Active/Inactive -conformation dependent. */
double Hexamer::kCIAoff()
{
  return (get_CIKaiB_bound() == reaction_consts->nBseq) ? reaction_consts->kICIAoff : reaction_consts->kACIAoff;
}


/*  KaiA on-rate to CII domain */
double Hexamer::kCIIAon()
{
  return reaction_consts->kCIIAon * exp( -0.5*dGACIIbind() );
}


/*  KaiA off-rate to CII domain */
double Hexamer::kCIIAoff()
{
  return reaction_consts->kCIIAoff * exp( 0.5*dGACIIbind() );
}


/* ADP off rate in CI domain */
double Hexamer::kCIADPoff()
{
  double kCIADPoff(reaction_consts->kCIADPoff0), 
         beta_actenergy(0);

  //Energy of transition state.
  if(get_CIKaiB_bound() == reaction_consts->nBseq && reaction_consts->nBseq > 0)
    beta_actenergy += this->dGICIActADPoff();
  else
    beta_actenergy += this->dGACIActADPoff();    

  // KaiB lowers ADP off rate by stabilization of CI-ADP state
  beta_actenergy += get_CIKaiB_bound() * reaction_consts->kAkIDoff;
  
  return kCIADPoff * exp( -beta_actenergy );
}


//Print reaction propensities of this Hexamer.
void Hexamer::print_hex_props()
{
  cout << endl;
  cout << "### Propensities of Hexamer - idx " << index << endl;
  cout << "Hex state - active: " << active << ", CIIA: " << CIIKaiA_bound 
       << ", CIB " << CIKaiB_bound << ", CIA " << CIKaiA_bound << endl;
  cout << "Mon state - S/T tot: " << get_hex_S() << "/" << get_hex_T() 
       << ", CI/CII-ADP tot: " << 6-get_hex_CIATP() << "/" << 6-get_hex_CIIATP() << endl;
  cout << "CII + A -> CIIA    (0) " << prop_list[0] << endl;
  cout << "CIIA -> A + CII    (1) " << prop_list[1] << endl;
  cout << "CIB(m) + B -> CIB(m+1) (2) " << prop_list[2] << endl;
  cout << "CIB(m) -> B + CIB(m-1) (3) " << prop_list[3] << endl;  
  cout << "CI + A -> CIA      (4) " << prop_list[4] << endl; 
  cout << "CIA -> A + CI      (5) " << prop_list[5] << endl;  
  cout << "Active -> Inactive (6) " << prop_list[6] << endl;
  cout << "Inactive -> Active (7) " << prop_list[7] << endl;  
  //cout << "Monomer qtot 1     (8) " << prop_list[8] << endl;
  //cout << "Monomer qtot 2     (9) " << prop_list[9] << endl;
  //cout << "Monomer qtot 3     (10) " << prop_list[10] << endl;
  //cout << "Monomer qtot 4     (11) " << prop_list[11] << endl;
  //cout << "Monomer qtot 5     (12) " << prop_list[12] << endl;      
  //cout << "Monomer qtot 6     (13) " << prop_list[13] << endl;
  cout << "Props - tot: " << hex_prop_cont->get_qtot_el(index) 
       << ", int: " << hex_prop_cont->get_qint_el(index) 
       << ", A: " << hex_prop_cont->get_qA_el(index) << endl;
}

//Total #nucleotide binding pockets with ATP in CI domain.
int Hexamer::get_hex_CIATP()
{
  return CIATP_bound;
}

int Hexamer::get_hex_CIIATP()
{
  return CIIATP_bound;
}

// Return #monomers in U-state.
int Hexamer::get_hex_pU()
{
  return 6 - pT - pD - pS;
}

// Return #monomers in T-state.
int Hexamer::get_hex_pT()
{
  return pT;
}

// Return #monomers in S-state.
int Hexamer::get_hex_pS()
{
  return pS;
}

// Return #monomers in D-state.
int Hexamer::get_hex_pD()
{
  return pD;
}

// Return #monomers with T site phosphorylated.
int Hexamer::get_hex_T()
{
  return pT + pD;
}

// Return #monomers with S site phosphorylated.
int Hexamer::get_hex_S()
{
  return pD + pS;
}


