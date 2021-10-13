#include "Hexamer.hpp"
#include "PropensityContainer.hpp"
#include "propagate.hpp"
#include "sampler.hpp"
#include "main.hpp"

const int NUM_REACTION_CONSTANTS = 10;

using namespace std;

int main(int argc, char *argv[]){
	//input format: <simulation parameters> <reaction constants> <data file> <bootstrap standard devation> <bootstrap ensemble size> <goodman-weare ensemble size> <goodman-weare alpha value> <maximum iterations> 
	
	/* MPI tags used
	 * 0: likelihoods
	 * 1: random number z generated for stretch move
	 * 2: index j of partner walker in stretch move
	 * 3: flag for whether a proposed change was accepted
	 */

	MPI_Init(NULL, NULL);
	int world_size;
	int world_id;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_id);

	string reaction_constants_filename;
	string system_variables_filename;
	string output_filename;
	int maxiter = 0;
	double alpha = 0;
	int max_sim_iter = 0;
	
	if(argc == 1) {

		cout << "Usage:\n";
		cout << "-sys\t System variables filename.\n";
		cout << "-rc\t Reaction constants filename. This specifies the starting state of sampler. Each line should specify the reaction coordinates for one walker in the ensemble sampler.\n";
		cout << "-a\t Alpha value for ensemble sampler. Must be greater than 1.\n";
		cout << "-max\t Maximum number of iterations to run the sampler for.\n";
		cout << "-out\t Output filename prefix. The average sampler state over all walkers and all iterations will be written to <prefix>_average.txt and the state of the sampler at the end of the last iteration will be written to <prefix>_sampler_state.txt \n";

		return 1;
	}

	for(int i = 1; i < argc - 1; i += 1) {
		
		if(strcmp(argv[i], "-sys") == 0) {
		
			system_variables_filename = argv[i + 1];
			i += 1;
		} 
		else if(strcmp(argv[i], "-rc") == 0 && reaction_constants_filename.empty()) {
			
			reaction_constants_filename = argv[i + 1];
			i += 1;
		}
		else if(strcmp(argv[i], "-a") == 0 && alpha < 1) {
			
			alpha = stod(argv[i + 1]);
			i += 1;
		}
		else if(strcmp(argv[i], "-max") == 0 && maxiter == 0) {
			
			maxiter = stoi(argv[i + 1]);
			i += 1;
		}
		else if(strcmp(argv[i], "-out") == 0 && output_filename.empty()) {
			
			output_filename = argv[i + 1];
			i += 1;
		}
		else if(strcmp(argv[i], "-simiter") == 0 && max_sim_iter == 0) {

			max_sim_iter = stoi(argv[i + 1]);
			i += 1;
		}
	}
	

	bool bad_input = false;

	if(output_filename.empty()) {

		cerr << "Missing output filename\n";
		bad_input = true;
	}
	if(system_variables_filename.size() == 0) {

		cerr << "Missing system variables input file\n";
		bad_input = true;
	}
	if(reaction_constants_filename.empty()) {

		cerr << "Missing reaction constants input file\n";
		bad_input = true;
	}
	if(alpha <= 1) {

		cerr << "Sampler alpha value must be a real number greater than 1\n";
		bad_input = true;
	}
	if(maxiter < 1) {
		
		cerr << "Sampler must be run for at least one step\n";
		bad_input = true;
	}
	if(max_sim_iter < 1) {

		cerr << "Simulations must be run for at least one step\n";
		bad_input = true;
	}
	if(bad_input) {

		return 1;
	}

	std::vector<std::vector<double>> starting_position = parse_sampler_coordinates(reaction_constants_filename);
	assert(starting_position.size() > 0);
	int NUM_REACTION_CONSTANTS = starting_position[0].size();
	std::vector<SystemVariables> sys(starting_position.size());
	SystemVariables sys_base;
	strcpy(sys_base.param_filename, system_variables_filename.c_str());
	initialize_system_vars(&sys_base);
	std::vector<ReactionConstants> reaction_constants(starting_position.size());

	for(int i = 0; i < starting_position.size(); i += 1) {

		reaction_constants[i] = initialize_reaction_consts(&sys_base);
		reaction_constants[i].kACIBFSon = starting_position[i][0];
		reaction_constants[i].kACIBGon = starting_position[i][1];
		reaction_constants[i].kACIBoff = starting_position[i][2];
		reaction_constants[i].kICIBFSon = starting_position[i][3];
		reaction_constants[i].kICIBGon = starting_position[i][4];
		reaction_constants[i].kICIBoff = starting_position[i][5];
		reaction_constants[i].kBswitch_f = starting_position[i][6];
		reaction_constants[i].kBswitch_r = starting_position[i][7];
		reaction_constants[i].kKidAon = starting_position[i][8];
		reaction_constants[i].kKidAoff = starting_position[i][9];
	}

	vector<ReactionConstants> reaction_constants_propose(reaction_constants.size());
	vector<double> z(reaction_constants.size());
	vector<int> partner(reaction_constants.size());
	vector<int> accept(reaction_constants.size());
	vector<double> log_likelihoods(reaction_constants.size(), -INFINITY);
	std::vector<double> log_likelihoods_propose(reaction_constants.size());

	random_device rd;
	uniform_real_distribution<double> u01(0, 1);
	mt19937_64 engine(rd());

	if(world_id == ROOT_PROCESS) {

		for(int i = 0; i < reaction_constants.size(); i += 1) {

			SystemVariables sys = sys_base;
			PropensityContainer prop_cont(sys.N_hexamers);
			std::vector<Hexamer> hexamers(sys.N_hexamers);
			double t_sample = sys.tequ;
			std::vector<double> data;

			for(int j = 0; j < sys.N_hexamers; j += 1) {

				hexamers[j].set_prop_container(&prop_cont);
				hexamers[j].set_reaction_consts(&reaction_constants[i]);
				hexamers[j].set_sysvars(&sys);
				hexamers[j].initialize_state(j,1,0,0,0,6,6,0,0,6);
				hexamers[j].set_propensities();
			}

			int step_count = 0;

			while(propagate(&sys, &hexamers[0], &prop_cont, &reaction_constants[i], u01, engine) && step_count < max_sim_iter) {

				step_count += 1;

				if(sys.tsim > t_sample) {

					double CIKaiB_bound;

					for(auto& hexamer: hexamers) {

						CIKaiB_bound += hexamer.get_CIKaiB_bound();
					}

					CIKaiB_bound /= sys.N_hexamers;
					data.push_back(CIKaiB_bound);
					t_sample += sys.t_sample_incr;
				}
			}

			if(step_count < max_sim_iter) {

				Eigen::Map<Eigen::VectorXd> x(&data[0], data.size());
				double mean = x.sum() / x.size();
				x.array() -= mean;
				Eigen::VectorXd autocorrelation = cross_correlation(x, x);
				log_likelihoods[i] = check_oscillation(autocorrelation, 6, sys.t_sample_incr, 23);
			}
		}
	}

	for(int step = 0; step < maxiter; step += 1) { 
		
		//Process 0 proposes z and j values for stretch move and sends to all other processes
		if(world_id == ROOT_PROCESS) {

			for(int i = 0; i < reaction_constants.size(); i += 1) {
				
				z[i] = pow((alpha - 1) * u01(engine) + 1, 2) / alpha;
				partner[i] = (int) ((reaction_constants.size() - 1) * u01(engine));
				partner[i] += (partner[i] >= i);
			}

			printf("Stretch move sizes:\n");

			for(int i = 0; i < z.size(); i += 1) {

				printf("%f\t", z[i]);
			}
			printf("\n");
			printf("Stretch move partners:\n");
			for(int i = 0; i < partner.size(); i += 1) {

				printf("%d\t", partner[i]);
			}
			printf("\n");
			printf("Sending stretch move parameters to processes:\t");
			for(int i = 1; i < world_size; i += 1) {

				MPI_Send(&z[0], z.size(), MPI_DOUBLE, i, TAG_STRETCH_MOVE_SIZE, MPI_COMM_WORLD);
				MPI_Send(&partner[0], partner.size(), MPI_INT, i, TAG_STRETCH_MOVE_PARTNER, MPI_COMM_WORLD);
				printf("%d\t", i);
			}
			printf("\n");
		}
		//All other processes receieve z and j values from process 0
		else {
			MPI_Recv(&z[0], z.size(), MPI_DOUBLE, ROOT_PROCESS, TAG_STRETCH_MOVE_SIZE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&partner[0], partner.size(), MPI_INT, ROOT_PROCESS, TAG_STRETCH_MOVE_PARTNER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		
		for(int i = 0; i < reaction_constants_propose.size(); i += 1) {
		
			reaction_constants_propose[i] = reaction_constants[partner[i]] + z[i] * (reaction_constants[i] - reaction_constants[partner[i]]);
		}
		
		for(auto& logp : log_likelihoods_propose) {

			logp = 0;
		}

		if(world_id == ROOT_PROCESS) {

			for(const auto& elem : reaction_constants_propose) {

				std::cout << elem.kACIBFSon << "\t";
				std::cout << elem.kACIBGon << "\t";
				std::cout << elem.kACIBoff << "\t";
				std::cout << elem.kICIBFSon << "\t";
				std::cout << elem.kICIBGon << "\t";
				std::cout << elem.kICIBoff << "\t";
				std::cout << elem.kBswitch_f << "\t";
				std::cout << elem.kBswitch_r << "\t";
				std::cout << elem.kKidAon << "\t";
				std::cout << elem.kKidAoff << "\n";
			}
		}

		//Each process computes partial likelihoods for each walker in ensemble
		for(int i = world_id; i < reaction_constants.size(); i += world_size) {

			SystemVariables sys = sys_base;
			PropensityContainer prop_cont(sys.N_hexamers);
			std::vector<Hexamer> hexamers(sys.N_hexamers);
			double t_sample = sys.tequ;
			std::vector<double> data;

			for(int j = 0; j < sys.N_hexamers; j += 1) {

				hexamers[j].set_prop_container(&prop_cont);
				hexamers[j].set_reaction_consts(&reaction_constants_propose[i]);
				hexamers[j].set_sysvars(&sys);
				hexamers[j].initialize_state(j,1,0,0,0,6,6,0,0,6);
				hexamers[j].set_propensities();
			}

			int step_count = 0;

			while(propagate(&sys, &hexamers[0], &prop_cont, &reaction_constants_propose[i], u01, engine) && step_count < max_sim_iter) {

				int step_count= 0;

				if(sys.tsim > t_sample) {

					double CIKaiB_bound;

					for(auto& hexamer: hexamers) {

						CIKaiB_bound += hexamer.get_CIKaiB_bound();
					}

					CIKaiB_bound /= sys.N_hexamers;
					data.push_back(CIKaiB_bound);
					t_sample += sys.t_sample_incr;
				}
			}

			if(step_count < max_sim_iter) {

				Eigen::Map<Eigen::VectorXd> x(&data[0], data.size());
				double mean = x.sum() / x.size();
				x.array() -= mean;
				Eigen::VectorXd autocorrelation = cross_correlation(x, x);
				log_likelihoods_propose[i] = check_oscillation(autocorrelation, 6, sys.t_sample_incr, 23);
			}
		}	

		if(world_id == ROOT_PROCESS) {
			
			//Process 0 receives and combines partial likelihoods
			/*printf("log likelihoods:\t");
			for(int i = 0; i < log_likelihoods.size(); i += 1) {


				printf("%f\t", log_likelihoods[i]);
			}
			printf("\n");*/
			vector<double> log_likelihoods_propose_recv(log_likelihoods_propose.size());
			
			printf("Receiving proposal log likelihoods from process:\n");
			for(int i = 1; i < world_size; i += 1) {		

				MPI_Recv(&log_likelihoods_propose_recv[0], log_likelihoods_propose.size(), MPI_DOUBLE, i, TAG_LIKELIHOODS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				printf("%d:\t", i);
				for(int j = 0; j < log_likelihoods_propose_recv.size(); j += 1) {

					printf("%f\t", log_likelihoods_propose_recv[j]);
				}
				printf("\n");
				
				for(int j = 0; j < log_likelihoods_propose.size(); j += 1) {
					
					log_likelihoods_propose[j] += log_likelihoods_propose_recv[j];
				}
			}
			printf("Proposal log likelihoods:\n");
			for(int i = 0; i < log_likelihoods_propose.size(); i += 1) {

				printf("%f\t", log_likelihoods_propose[i]);
			}
			printf("\n");
			//Process 0 computes whether proposed move for each walker is accepted or not
			for(int i = 0; i < reaction_constants.size(); i += 1) {
				
				double log_likelihood_ratio = log_likelihoods_propose[i] - log_likelihoods[i];

				if(((NUM_REACTION_CONSTANTS - 1) * log(z[i]) + log_likelihood_ratio) > log(u01(engine))) {

					accept[i] = 1;
					log_likelihoods[i] = log_likelihoods_propose[i];
				}
				else {

					accept[i] = 0;
				}	
			}
			//Process 0 sends accept/reject information to other processes
			printf("Acceptance:\n");
			for(int i = 0; i < accept.size(); i += 1) {

				printf("%d\t", accept[i]);
			}
			printf("\n");
			printf("Sending acceptance information to processes:\n");
			for(int i = 1; i < world_size; i += 1) {
				
				printf("%d\t", i);
				MPI_Send(&accept[0], accept.size(), MPI_INT, i, TAG_STRETCH_MOVE_ACCEPTANCE, MPI_COMM_WORLD);
				MPI_Send(&log_likelihoods[0], log_likelihoods.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
			printf("\n");
		}
		else {
			//All other processes send process 0 their computed partial likelihoods
			MPI_Send(&log_likelihoods_propose[0], log_likelihoods_propose.size(), MPI_DOUBLE, ROOT_PROCESS, TAG_LIKELIHOODS, MPI_COMM_WORLD);
			//Receive update information from process 0
			MPI_Recv(&accept[0], accept.size(), MPI_INT, ROOT_PROCESS, TAG_STRETCH_MOVE_ACCEPTANCE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&log_likelihoods[0], log_likelihoods.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		for(int i = 0; i < reaction_constants.size(); i += 1) {

			if(accept[i]) {				

				reaction_constants[i] = reaction_constants_propose[i];
			}
		}
	}
	if(world_id == ROOT_PROCESS) {

		ofstream outfile;
		outfile.open(output_filename + "_sampler_state.txt");

		if(outfile.is_open()) {

			for(int i = 0; i < reaction_constants.size(); i += 1) {

				outfile << reaction_constants[i].kACIBFSon << "\t";
				outfile << reaction_constants[i].kACIBGon << "\t";
				outfile << reaction_constants[i].kACIBoff << "\t";
				outfile << reaction_constants[i].kICIBFSon << "\t";
				outfile << reaction_constants[i].kICIBGon << "\t";
				outfile << reaction_constants[i].kICIBoff << "\t";
				outfile << reaction_constants[i].kBswitch_f << "\t";
				outfile << reaction_constants[i].kBswitch_r << "\t";
				outfile << reaction_constants[i].kKidAon << "\t";
				outfile << reaction_constants[i].kKidAoff;
				outfile << endl;
			}

			outfile.close();
		}
	}
	MPI_Finalize();
}

std::vector<std::vector<double>> parse_sampler_coordinates(std::string &filename) {

	std::vector<std::vector<double>> coordinates;
	std::ifstream infile(filename);
	std::string line;

	while(getline(infile, line)) {

		std::vector<double> walker;
		std::string num = "";

		for(const auto& ch : line) {

			if(!isspace(ch)) {

				assert(isdigit(ch) || ch == '.' || ch == 'e' || ch == '-' || ch == '+');
				num += ch;
			}
			else {

				walker.push_back(stod(num));
				num = "";
			}
		}

		if(num.size() > 0) {

			walker.push_back(stod(num));
		}

		coordinates.push_back(walker);
	}

	for(const auto& elem : coordinates) {

		assert(elem.size() == coordinates[0].size());
	}

	return coordinates;
}

ReactionConstants operator+(ReactionConstants x, ReactionConstants const &y) {

	return x += y;
}

ReactionConstants& operator+=(ReactionConstants &x, ReactionConstants const &y) {

	x.kCIhyd *= y.kCIhyd;
	x.kCIADPoff0 *= y.kCIADPoff0;
	x.kCIADPoffA *= y.kCIADPoffA;
	x.kAkIDoff += y.kAkIDoff;
	x.kconf0 *= y.kconf0;
	x.ddGTDconf += y.ddGTDconf;
	x.ddGTDAbind += y.ddGTDAbind;
	x.nIADPAref += y.nIADPAref;
	x.nIADPconfref += y.nIADPconfref;
	x.dGIbindB += y.dGIbindB;
	x.kACIAon *= y.kACIAon;
	x.kACIAoff *= y.kACIAoff;
	x.kICIAon *= y.kICIAon;
	x.kICIAoff *= y.kICIAoff;
	x.nAseq += y.nAseq;
	x.nBseq += y.nBseq;
	x.dgACIIU += y.dgACIIU;
	x.dgACIIT += y.dgACIIT;
	x.dgACIID += y.dgACIID;
	x.dgACIIS += y.dgACIIS;
	x.dGICII += y.dGICII;
	x.dgACIActU += y.dgACIActU;
	x.dgACIActT += y.dgACIActT;
	x.dgACIActD += y.dgACIActD;
	x.dgACIActS += y.dgACIActS;
	x.dgICIActU += y.dgICIActU;
	x.dgICIActT += y.dgICIActT;
	x.dgICIActD += y.dgICIActD;
	x.dgICIActS += y.dgICIActS;
	x.kCIIAon *= y.kCIIAon;
	x.kCIIAoff *= y.kCIIAoff;
	x.kCIIhyd0 *= y.kCIIhyd0;
	x.kCIIhydA *= y.kCIIhydA;
	x.kCIInucloff0 *= y.kCIInucloff0;
	x.kCIInucloffA *= y.kCIInucloffA;
	x.KATPoKADP *= y.KATPoKADP;
	x.kUT *= y.kUT;
	x.kTU *= y.kTU;
	x.kTD *= y.kTD;
	x.kDT *= y.kDT;
	x.kSD *= y.kSD;
	x.kDS *= y.kDS;
	x.kUS *= y.kUS;
	x.kSU *= y.kSU;
	x.kACIBFSon *= y.kACIBFSon;
	x.kACIBGon *= y.kACIBGon;
	x.kACIBoff *= y.kACIBoff;
	x.kICIBFSon *= y.kICIBFSon;
	x.kICIBGon *= y.kICIBGon;
	x.kICIBoff *= y.kICIBoff;
	x.kBswitch_f *= y.kBswitch_f;
	x.kBswitch_r *= y.kBswitch_r;
	x.kKidAon *= y.kKidAon;
	x.kKidAoff *= y.kKidAoff;

	return x;

}

ReactionConstants operator*(double a, ReactionConstants x) {
	
	return x *= a;
}

ReactionConstants operator*(ReactionConstants x, double a) {
	
	return x *= a;
}

ReactionConstants& operator*=(ReactionConstants &x, double a) {
	
	x.kCIhyd = pow(x.kCIhyd, a);
	x.kCIADPoff0 = pow(x.kCIADPoff0, a);
	x.kCIADPoffA = pow(x.kCIADPoffA, a);
	x.kAkIDoff *= a;
	x.kconf0 = pow(x.kconf0, a);
	x.ddGTDconf *= a;
	x.ddGTDAbind *= a;
	x.nIADPAref = (int) rint(a * x.nIADPAref);
	x.nIADPconfref = (int) rint(a * x.nIADPconfref);
	x.dGIbindB *= a;
	x.kACIAon = pow(x.kACIAon, a);
	x.kACIAoff = pow(x.kACIAoff, a);
	x.kICIAon = pow(x.kICIAon, a);
	x.kICIAoff = pow(x.kICIAoff, a);
	x.nAseq = (int) rint(a * x.nAseq);
	x.nBseq = (int) rint(a * x.nBseq);
	x.dgACIIU *= a;
	x.dgACIIT *= a;
	x.dgACIID *= a;
	x.dgACIIS *= a;
	x.dGICII *= a;
	x.dgACIActU *= a;
	x.dgACIActT *= a;
	x.dgACIActD *= a;
	x.dgACIActS *= a;
	x.dgICIActU *= a;
	x.dgICIActT *= a;
	x.dgICIActD *= a;
	x.dgICIActS *= a;
	x.kCIIAon = pow(x.kCIIAon, a);
	x.kCIIAoff = pow(x.kCIIAoff, a);
	x.kCIIhyd0 = pow(x.kCIIhyd0, a);
	x.kCIIhydA = pow(x.kCIIhydA, a);
	x.kCIInucloff0 = pow(x.kCIInucloff0, a);
	x.kCIInucloffA = pow(x.kCIInucloffA, a);
	x.KATPoKADP = pow(x.KATPoKADP, a);
	x.kUT = pow(x.kUT, a);
	x.kTU = pow(x.kTU, a);
	x.kTD = pow(x.kTD, a);
	x.kDT = pow(x.kDT, a);
	x.kSD = pow(x.kSD, a);
	x.kDS = pow(x.kDS, a);
	x.kUS = pow(x.kUS, a);
	x.kSU = pow(x.kSU, a);
	x.kACIBFSon = pow(x.kACIBFSon, a);
	x.kACIBGon = pow(x.kACIBGon, a);
	x.kACIBoff = pow(x.kACIBoff, a);
	x.kICIBFSon = pow(x.kICIBFSon, a);
	x.kICIBGon = pow(x.kICIBGon, a);
	x.kICIBoff = pow(x.kICIBoff, a);
	x.kBswitch_f = pow(x.kBswitch_f, a);
	x.kBswitch_r = pow(x.kBswitch_r, a);
	x.kKidAon = pow(x.kKidAon, a);
	x.kKidAoff = pow(x.kKidAoff, a);

	return x;
}

ReactionConstants operator-(ReactionConstants x, ReactionConstants const &y) {
	
	return x -= y;
}

ReactionConstants& operator-=(ReactionConstants &x, ReactionConstants const &y) {
	
	x.kCIhyd /= y.kCIhyd;
	x.kCIADPoff0 /= y.kCIADPoff0;
	x.kCIADPoffA /= y.kCIADPoffA;
	x.kAkIDoff -= y.kAkIDoff;
	x.kconf0 /= y.kconf0;
	x.ddGTDconf -= y.ddGTDconf;
	x.ddGTDAbind -= y.ddGTDAbind;
	x.nIADPAref -= y.nIADPAref;
	x.nIADPconfref -= y.nIADPconfref;
	x.dGIbindB -= y.dGIbindB;
	x.kACIAon /= y.kACIAon;
	x.kACIAoff /= y.kACIAoff;
	x.kICIAon /= y.kICIAon;
	x.kICIAoff /= y.kICIAoff;
	x.nAseq -= y.nAseq;
	x.nBseq -= y.nBseq;
	x.dgACIIU -= y.dgACIIU;
	x.dgACIIT -= y.dgACIIT;
	x.dgACIID -= y.dgACIID;
	x.dgACIIS -= y.dgACIIS;
	x.dGICII -= y.dGICII;
	x.dgACIActU -= y.dgACIActU;
	x.dgACIActT -= y.dgACIActT;
	x.dgACIActD -= y.dgACIActD;
	x.dgACIActS -= y.dgACIActS;
	x.dgICIActU -= y.dgICIActU;
	x.dgICIActT -= y.dgICIActT;
	x.dgICIActD -= y.dgICIActD;
	x.dgICIActS -= y.dgICIActS;
	x.kCIIAon /= y.kCIIAon;
	x.kCIIAoff /= y.kCIIAoff;
	x.kCIIhyd0 /= y.kCIIhyd0;
	x.kCIIhydA /= y.kCIIhydA;
	x.kCIInucloff0 /= y.kCIInucloff0;
	x.kCIInucloffA /= y.kCIInucloffA;
	x.KATPoKADP /= y.KATPoKADP;
	x.kUT /= y.kUT;
	x.kTU /= y.kTU;
	x.kTD /= y.kTD;
	x.kDT /= y.kDT;
	x.kSD /= y.kSD;
	x.kDS /= y.kDS;
	x.kUS /= y.kUS;
	x.kSU /= y.kSU;
	x.kACIBFSon /= y.kACIBFSon;
	x.kACIBGon /= y.kACIBGon;
	x.kACIBoff /= y.kACIBoff;
	x.kICIBFSon /= y.kICIBFSon;
	x.kICIBGon /= y.kICIBGon;
	x.kICIBoff /= y.kICIBoff;
	x.kBswitch_f /= y.kBswitch_f;
	x.kBswitch_r /= y.kBswitch_r;
	x.kKidAon /= y.kKidAon;
	x.kKidAoff /= y.kKidAoff;

	return x;

}

ReactionConstants operator/(ReactionConstants const &x, double a) {
	
	return x * (1 / a);
}

ReactionConstants& operator/=(ReactionConstants &x, double a) {
	
	return x *= 1 / a;
}

ReactionConstants& set_zero(ReactionConstants &x) {
	
	x.kCIhyd -= 0;
	x.kCIADPoff0 = 0;
	x.kCIADPoffA = 0;
	x.kAkIDoff = 0;
	x.kconf0 = 0;
	x.ddGTDconf = 0;
	x.ddGTDAbind = 0;
	x.nIADPAref = 0;
	x.nIADPconfref = 0;
	x.dGIbindB = 0;
	x.kACIAon = 0;
	x.kACIAoff = 0;
	x.kICIAon = 0;
	x.kICIAoff = 0;
	x.nAseq = 0;
	x.nBseq = 0;
	x.dgACIIU = 0;
	x.dgACIIT = 0;
	x.dgACIID = 0;
	x.dgACIIS = 0;
	x.dGICII = 0;
	x.dgACIActU = 0;
	x.dgACIActT = 0;
	x.dgACIActD = 0;
	x.dgACIActS = 0;
	x.dgICIActU = 0;
	x.dgICIActT = 0;
	x.dgICIActD = 0;
	x.dgICIActS = 0;
	x.kCIIAon = 0;
	x.kCIIAoff = 0;
	x.kCIIhyd0 = 0;
	x.kCIIhydA = 0;
	x.kCIInucloff0 = 0;
	x.kCIInucloffA = 0;
	x.KATPoKADP = 0;
	x.kUT = 0;
	x.kTU = 0;
	x.kTD = 0;
	x.kDT = 0;
	x.kSD = 0;
	x.kDS = 0;
	x.kUS = 0;
	x.kSU = 0;
	x.kACIBFSon = 0;
	x.kACIBGon = 0;
	x.kACIBoff = 0;
	x.kICIBFSon = 0;
	x.kICIBGon = 0;
	x.kICIBoff = 0;
	x.kBswitch_f = 0;
	x.kBswitch_r = 0;
	x.kKidAon = 0;
	x.kKidAoff = 0;
	
	return x;
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
        //cerr << param_name << " = " << param_value << endl;
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
           //cerr << param_name << " = " << param_string << endl;
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
  sys->tsim = 0.0;
  sys->N_hexamers = sys->KaiC0 * sys->volume;
  sys->Afree = sys->KaiA0 * sys->volume;
  sys->B_inactive = sys->KaiB0 * sys->volume;
  sys->B_active = 0;
  sys->KaiBKidA = 0;
  sys->KidA_free = sys->KidA0 * sys->volume;
  sys->cAfree = sys->Afree / sys->volume;
  sys->CIATPcons = 0; sys->CIIATPcons = 0;
  
  sys->step_cntr = 0;
  sys->sample_cntr = 0;
}

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
        //cerr << param_name << " = " << param_value << endl;
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
           //cerr << param_name << " = " << param_string << endl;
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

double check_oscillation(const Eigen::VectorXd autocorrelation, int min_crossings, double t_sample_incr, double expected_period) {

	Eigen::VectorXd deriv = autocorrelation.tail(autocorrelation.size() - 1) - autocorrelation.head(autocorrelation.size() - 1);
	Eigen::VectorXi rect = (deriv.array() <= 0).cast<int>();
	Eigen::VectorXi peaks = ((rect.tail(rect.size() - 1) - rect.head(rect.size() - 1)).array() > 0).cast<int>();
	std::vector<int> peak_indices;

	if(peaks.count() < min_crossings) {

		return -INFINITY;
	}

	for(int i = 0; i < peaks.size(); i += 1) {

		if(peaks[i] > 0) {

			peak_indices.push_back(i);
		}
	}

	Eigen::Map<Eigen::VectorXi> indices_map(peak_indices.data(), peak_indices.size());
	double period = t_sample_incr * (indices_map.tail(indices_map.size() - 1) - indices_map.head(indices_map.size() - 1)).sum() / (indices_map.size() - 1);

	return -pow(period - expected_period, 2);
}

Eigen::VectorXd cross_correlation(const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y) {


	Eigen::VectorXd x_pad;
	Eigen::VectorXd y_pad;

	if(x.size() > y.size()) {

		x_pad = Eigen::VectorXd::Zero(x.size() * 2);
		x_pad.head(x.size()) = x;
		y_pad = Eigen::VectorXd::Zero(x.size() * 2);
		y_pad.head(y.size()) = y;
	}
	else {

		x_pad = Eigen::VectorXd::Zero(y.size() * 2);
		x_pad.head(x.size()) = x;
		y_pad = Eigen::VectorXd::Zero(y.size() * 2);
		y_pad.head(y.size()) = y;
	}

	Eigen::FFT<double> fft;
	Eigen::VectorXcd F_x(x_pad.size());
	Eigen::VectorXcd F_y(y_pad.size());
	fft.fwd(F_x, x_pad);
	fft.fwd(F_y, y_pad);
	Eigen::VectorXcd F_corr = F_x.conjugate().array() * F_y.array();
	Eigen::VectorXcd corr(F_corr.size());
	fft.inv(corr, F_corr);
	return corr.real().head(corr.size() / 2) / corr.real()[0];
}
