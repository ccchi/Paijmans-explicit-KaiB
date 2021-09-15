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
#include <mpi.h>
#include <cassert>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

const int ROOT_PROCESS = 0;
const int TAG_LIKELIHOODS = 0;
const int TAG_STRETCH_MOVE_SIZE = 1;
const int TAG_STRETCH_MOVE_PARTNER = 2;
const int TAG_STRETCH_MOVE_ACCEPTANCE = 3;

std::vector<ReactionConstants> parse_reaction_constants(const char* filename);

ReactionConstants operator+(ReactionConstants x, ReactionConstants const &y);
ReactionConstants& operator+=(ReactionConstants &x, ReactionConstants const &y);
ReactionConstants operator*(double a, ReactionConstants x);
ReactionConstants operator*(ReactionConstants x, double a);
ReactionConstants& operator*=(ReactionConstants &x, double a);
ReactionConstants operator-(ReactionConstants x, ReactionConstants const &y);
ReactionConstants& operator-=(ReactionConstants &x, ReactionConstants const &y);
ReactionConstants operator/(ReactionConstants const &x, double a);
ReactionConstants& operator/=(ReactionConstants &x, double a);
ReactionConstants& set_zero(ReactionConstants &x);
std::string delUnnecessary(std::string &str);
ReactionConstants initialize_reaction_consts(SystemVariables *sys);
void initialize_system_vars(SystemVariables *sys);
std::vector<std::vector<double>> parse_sampler_coordinates(std::string &filename);
Eigen::VectorXd cross_correlation(const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y);
double check_oscillation(const Eigen::VectorXd autocorrelation, int min_crossings, double corr_threshold);
