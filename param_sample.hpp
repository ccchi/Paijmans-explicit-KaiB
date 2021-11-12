#ifndef PARAM_SAMPLE
#define PARAM_SAMPLE

#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <algorithm>
#include <vector>
#include <tuple>
#include <math.h>
#include <string>
#include <cstring>
#include <cassert>
#include "io.hpp"
#include "sampler.hpp"
#include "moves.hpp"
#include "data_structures.hpp"
#include "likelihood.hpp"

const int TAG_LIKELIHOODS = 0;
const int TAG_STRETCH_MOVE_SIZE = 1;
const int TAG_STRETCH_MOVE_PARTNER = 2;
const int TAG_STRETCH_MOVE_ACCEPTANCE = 3;

std::vector<ReactionConstants> parse_reaction_constants(const char* filename);
std::string delUnnecessary(std::string &str);
ReactionConstants initialize_reaction_consts(SystemVariables *sys);
void initialize_system_vars(SystemVariables *sys);

#endif
