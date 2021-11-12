/* 
  Propagate function propagates the system one reaction.
*/

#ifndef _PROPAGATE_H_
#define _PROPAGATE_H_

#include <random>
#include "data_structures.hpp"
#include "PropensityContainer.hpp"

bool propagate(SystemVariables *sys, Hexamer *hexamers, PropensityContainer *prop_cont, ReactionConstants *reaction_consts, std::uniform_real_distribution<double> &u01, std::mt19937_64 &engine);

#endif
