/* propagates the simulator one event,
   Returns TRUE if simulation has not reached required end time yet.
*/

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>
#include <random>

#include "main.hpp"
#include "Hexamer.hpp"
#include "PropensityContainer.hpp"
#include "propagate.hpp"

using namespace std;

bool propagate(SystemVariables *sys, Hexamer *hexamers, PropensityContainer *prop_cont, ReactionConstants *reaction_consts, uniform_real_distribution<double> &u01, mt19937_64 &engine)
{
	//Draw two random numbers.
	double rnd1( u01(engine) ), rnd2( u01(engine)), rnd3( u01(engine));
	double prop_B_f = reaction_consts->kBswitch_f * sys->B_inactive;
	double prop_B_r = reaction_consts->kBswitch_r * sys->B_active;
	double prop_KidA_bind = reaction_consts->kKidAon * sys->B_active * sys->KidA_free / sys->volume;
	double prop_KidA_unbind = reaction_consts->kKidAoff * sys->KaiBKidA;
	double prop_total = prop_cont->get_qtot() + prop_B_f + prop_B_r + prop_KidA_bind + prop_KidA_unbind;
	//Update tsim. 
	sys->tsim -= log(rnd1) / prop_total;
	double choice = rnd2 * prop_total;
	//printf("%f ", (*sys).tsim);
	if(choice < prop_cont->get_qtot()) {
		//Choose hexamer to propagate with rnd2.
		int fired_hex_idx( prop_cont->choose_index(rnd3) ); 

		/* When KaiA associates with or disociates from the CI or CII domains of a hexamer,
		the KaiA concentration in solution changes and propensities of all hexamers change. */
		int ext_changed = hexamers[fired_hex_idx].propagate(u01, engine);
		if(ext_changed == 1)
		{   
		//Change Afree dependent reactions in all hexamers.
			prop_cont->set_qA_all(hexamers);
		}
		else if(ext_changed == 2) 
		{
			prop_cont->set_qB_all(hexamers);
		}
		else if(ext_changed == 3) {

			prop_cont->set_qKidA_all(hexamers);
		}
	}
	else {
		if(choice < prop_cont->get_qtot() + prop_B_f) {
			
			sys->B_inactive -= 1;
			sys->B_active += 1;
		}
		else if(choice < prop_cont->get_qtot() + prop_B_f + prop_B_r) {
			
			sys->B_inactive += 1;
			sys->B_active -= 1;
		}
		else if(choice < prop_cont->get_qtot() + prop_B_f + prop_B_r + prop_KidA_bind) {

			sys->B_active -= 1;
			sys->KidA_free -= 1;
			sys->KaiBKidA += 1;
		}
		else {

			sys->B_active += 1;
			sys->KaiBKidA -= 1;
			sys->KidA_free += 1;
		}
		prop_cont->set_qB_all(hexamers);
	}
	sys->step_cntr++;
	return sys->tsim < sys->tend;
}
