#include "param_sample.hpp"
#include "data_structures.hpp"

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

	Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> walkers = loadtxt(reaction_constants_filename); 
	assert(walkers.size() > 0);
	int NUM_REACTION_CONSTANTS = walkers.cols();
	SystemVariables sys_base;
	strcpy(sys_base.param_filename, system_variables_filename.c_str());
	ReactionConstants reaction_consts_base = initialize_reaction_consts(&sys_base);
	initialize_system_vars(&sys_base);

	std::random_device rd;
	unsigned int seed;

	if(world_id == ROOT_PROCESS) {
	
		seed = rd();
	}

	MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
	uniform_real_distribution<double> u01_sampler(0, 1);
	uniform_real_distribution<double> u01_ssa(0, 1);
	mt19937_64 engine_sampler(seed);

	if(world_id == ROOT_PROCESS) {

		seed = rd();
	}

	MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
	mt19937_64 engine_ssa(seed);
	Eigen::Array<double, Eigen::Dynamic, 2, Eigen::RowMajor> log_posteriors = compute_initial_log_posteriors(walkers, log_posterior_function, sys_base, reaction_consts_base, max_sim_iter, engine_ssa, u01_ssa);

	for(int step = 0; step < maxiter; step += 1) { 
		
		int accepted_step;
		int failed_step;
		std::tie(walkers, log_posteriors, accepted_step, failed_step) = affine_invariant_ensemble_mpi(walkers, log_posteriors, alpha, engine_sampler, u01_sampler, walk_move_DE<double>, log_posterior_function, sys_base, reaction_consts_base, max_sim_iter, engine_ssa, u01_ssa);

		if(world_id == ROOT_PROCESS) {

			std::cout << "log posterior:\n";
			std::cout << log_posteriors.rowwise().sum() << "\n";
		}
	}
	if(world_id == ROOT_PROCESS) {

		savetxt(output_filename, walkers);
	}
	MPI_Finalize();
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
