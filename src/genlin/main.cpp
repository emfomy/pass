////////////////////////////////////////////////////////////////////////////////
// Particle Swarm Stepwise (PaSS) Algorithm                                   //
//                                                                            //
// main.cpp                                                                   //
// The main funcitons                                                         //
//                                                                            //
// Author: emfo<emfomy@gmail.com>                                             //
////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <numeric>
#include <getopt.h>
#include <mkl.h>
#include <mpi.h>
#include <omp.h>
#include "pass.hpp"

using namespace pass;

// Default arguments
const unsigned int kTest = 100;
const char *kDataRoot    = "genlin.dat";

// Global variables
bool *J0;        // vector, 1 by p, the chosen indices (solution)
char *dataname;  // string, the name of data
int mpi_size;    // the size of MPI communicator
int mpi_rank;    // the rank of MPI processes

// Functions
void PassHelp( const char *cmd );
void PassLoad( const char *fileroot );

////////////////////////////////////////////////////////////////////////////////
// Display help messages                                                      //
////////////////////////////////////////////////////////////////////////////////
void PassHelp( const char *cmd ) {
  Parameter temp;
  printf("Usage: %s [options] ...\n", cmd);
  printf("\n%-48s%-48s%s\n", "Option", "Detail", "Defalut Value");
  printf("\nMain Options:\n");
  printf("  %-46s%-48s%s\n",
         "-f <file>, --file <file>",
         "load data from <file>",
         kDataRoot);
  printf("  %-46s%-48s%d\n",
         "-i ###, --iteration ###",
         "the number of iterations",
         temp.num_iteration);
  printf("  %-46s%-48s%d\n",
         "-p ###, --particle ###",
         "the number of particles per thread",
         temp.num_particle_thread);
  printf("  %-46s%-48s%d\n",
         "-t ###, --test ###",
         "the number of tests",
         kTest);
  printf("  %-46s%-40s\n",
         "-h, --help", "display help messages");

  printf("\nProbability Options:\n");
  printf("  --prob <pfg> <pfl> <pfr> <pbl> <pbr>\n");
  printf("    %-44s%-48s%.2f\n",
         "<pfg>",
         "the probabilities of forward step: global",
         temp.prob_forward_global);
  printf("    %-44s%-48s%.2f\n",
         "<pfl>",
         "the probabilities of forward step: local",
         temp.prob_forward_local);
  printf("    %-44s%-48s%.2f\n",
         "<pfr>",
         "the probabilities of forward step: random",
         temp.prob_forward_random);
  printf("    %-44s%-48s%.2f\n",
         "<pbl>",
         "the probabilities of backward step: local",
         temp.prob_backward_local);
  printf("    %-44s%-48s%.2f\n",
         "<pbr>",
         "the probabilities of backward step: random",
         temp.prob_backward_random);

  printf("\nCriterion Options:\n");
  printf("  %-46s%-48s\n",
         "--AIC",
         "Akaike information criterion");
  printf("  %-46s%-48s\n",
         "--BIC",
         "Bayesian information criterion");
  printf("  %-46s%-48s\n",
         "--EBIC=<gamma>",
         "Extended Bayesian information criterion");
  printf("    %-44s%-48s%.2f\n",
         "<gamma> (optional)",
         "the parameter of EBIC",
         temp.ebic_gamma);
  printf("  %-46s%-48s\n",
         "--HDBIC (default)",
         "High-dimensional Bayesian information criterion");
  printf("  %-46s%-48s\n",
         "--HQC",
         "Hannan-Quinn information criterion");
  printf("  %-46s%-48s\n",
         "--HDHQC",
         "High-dimensional Hannan-Quinn information criterion");
}

////////////////////////////////////////////////////////////////////////////////
// Main function                                                              //
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char **argv ) {
  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  // Initialize random seed
  srand(time(NULL) ^ mpi_rank);
  srand(rand());

  // Initialize arguments
  auto num_test = kTest;
  auto dataroot = kDataRoot;

  // Load arguments
  int optidx = 0;
  char opts[] = "f:i:p:t:\1\2\3\4::\5\6\7h", c;
  opterr = (mpi_rank == 0);
  option long_opts[] = {
    {"file",      required_argument, nullptr, 'f'},
    {"iteration", required_argument, nullptr, 'm'},
    {"particle",  required_argument, nullptr, 't'},
    {"test",      required_argument, nullptr, 'o'},
    {"prob",      no_argument,       nullptr, '\1'},
    {"AIC",       no_argument,       nullptr, '\2'},
    {"BIC",       no_argument,       nullptr, '\3'},
    {"EBIC",      optional_argument, nullptr, '\4'},
    {"HDBIC",     no_argument,       nullptr, '\5'},
    {"HQC",       no_argument,       nullptr, '\6'},
    {"HDHQC",     no_argument,       nullptr, '\7'},
    {"help",      no_argument,       nullptr, 'h'}
  };
  while ( (c = getopt_long(argc, argv, opts, long_opts, &optidx)) != -1 ) {
    switch ( c ) {
      case 'f': {
        dataroot = optarg;
        break;
      }
      case 'i': {
        parameter.num_iteration = atoi(optarg);
        break;
      }
      case 'p': {
        parameter.num_particle_thread = atoi(optarg);
        break;
      }
      case 't': {
        num_test = atoi(optarg);
        break;
      }
      case '\1': {
        if ( argv[optind] == nullptr ) {
          if ( mpi_rank == 0 ) {
            fprintf(stderr, "%s: option '--%s' requires requires 5 arguments\n",
                   argv[0], long_opts[optidx].name);
            PassHelp(argv[0]);
          }
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        parameter.prob_forward_global = atof(argv[optind]);
        optind++;
        if ( argv[optind] == nullptr ) {
          if ( mpi_rank == 0 ) {
            fprintf(stderr, "%s: option '--%s' requires requires 5 arguments\n",
                   argv[0], long_opts[optidx].name);
            PassHelp(argv[0]);
          }
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        parameter.prob_forward_local = atof(argv[optind]);
        optind++;
        if ( argv[optind] == nullptr ) {
          if ( mpi_rank == 0 ) {
            fprintf(stderr, "%s: option '--%s' requires requires 5 arguments\n",
                   argv[0], long_opts[optidx].name);
            PassHelp(argv[0]);
          }
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        parameter.prob_forward_random = atof(argv[optind]);
        optind++;
        if ( argv[optind] == nullptr ) {
          if ( mpi_rank == 0 ) {
            fprintf(stderr, "%s: option '--%s' requires requires 5 arguments\n",
                   argv[0], long_opts[optidx].name);
            PassHelp(argv[0]);
          }
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        parameter.prob_backward_local = atof(argv[optind]);
        optind++;
        if ( argv[optind] == nullptr ) {
          if ( mpi_rank == 0 ) {
            fprintf(stderr, "%s: option '--%s' requires requires 5 arguments\n",
                   argv[0], long_opts[optidx].name);
            PassHelp(argv[0]);
          }
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        parameter.prob_backward_random = atof(argv[optind]);
        break;
      }
      case '\2': {
        parameter.criterion = AIC;
        break;
      }
      case '\3': {
        parameter.criterion = BIC;
        break;
      }
      case '\4': {
        parameter.criterion = EBIC;
        if ( optarg ) {
          parameter.ebic_gamma = atof(optarg);
        }
        break;
      }
      case '\5': {
        parameter.criterion = HDBIC;
        break;
      }
      case '\6': {
        parameter.criterion = HQC;
        break;
      }
      case '\7': {
        parameter.criterion = HDHQC;
        break;
      }
      case 'h': {
        if ( mpi_rank == 0 ) {
          PassHelp(argv[0]);
        }
        MPI_Abort(MPI_COMM_WORLD, 0);
      }
      default: {
        if ( mpi_rank == 0 ) {
          PassHelp(argv[0]);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
  }

  // Check parameters
  auto num_thread = omp_get_max_threads();
  auto num_proc = omp_get_num_procs();
  if ( num_thread > num_proc ) {
    omp_set_num_threads(num_proc);
    num_thread = num_proc;
  }
  auto num_particle = mpi_size * num_thread * parameter.num_particle_thread;

  ////////////////////////////////////////////////////////////////////////////
  // Load data                                                              //
  ////////////////////////////////////////////////////////////////////////////

  if ( mpi_rank == 0 ) {
    printf("================================================================"
        "================================================================\n");
  }

  // Load data and allocate memory
  PassLoad(dataroot);
  I0 = new bool[p];

  // Display parameters
  if ( mpi_rank == 0 ) {
    if ( parameter.criterion == EBIC ) {
      printf("%s: n=%d, p=%d, #Node=%d, #Thread=%d, "
             "#Particle=%d, #Iteration=%d, #Test=%d, Criterion=%s%.1f\n",
             dataname, n, p, mpi_size, num_thread,
             num_particle, parameter.num_iteration, num_test,
             Criterion2String(parameter.criterion), parameter.ebic_gamma);
    } else {
      printf("%s: n=%d, p=%d, #Node=%d, #Thread=%d, "
             "#Particle=%d, #Iteration=%d, #Test=%d, Criterion=%s\n",
             dataname, n, p, mpi_size, num_thread,
             num_particle, parameter.num_iteration, num_test,
             Criterion2String(parameter.criterion));
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // Centralize and normalize the original data                             //
  ////////////////////////////////////////////////////////////////////////////

  if ( mpi_rank == 0 ) {
    printf("Normalizing data... ");
  }

  // Centralize and normalize X0
  for ( auto j = 0; j < p; ++j ) {
    float stemp = 0.0f;
    for ( auto i = 0; i < n; ++i ) {
      stemp += X0[i+j*n];
    }
    stemp /= n;
    vsLinearFrac(n, X0+j*n, X0+j*n, 1.0f, -stemp, 0.0f, 1.0f, X0+j*n);
    cblas_sscal(n, (1.0f/cblas_snrm2(n, X0+j*n, 1)), X0+j*n, 1);
  }

  // Centralize and normalize Y0
  {
    float stemp = 0.0f;
    for ( auto i = 0; i < n; ++i ) {
      stemp += Y0[i];
    }
    stemp /= n;
    vsLinearFrac(n, Y0, Y0, 1.0f, -stemp, 0.0f, 1.0f, Y0);
    cblas_sscal(n, (1.0f/cblas_snrm2(n, Y0, 1)), Y0, 1);
  }

  parameter.is_normalized = true;

  if ( mpi_rank == 0 ) {
    printf("Done.\n");
  }

  ////////////////////////////////////////////////////////////////////////////
  // Run PaSS                                                               //
  ////////////////////////////////////////////////////////////////////////////

  // Declare variables
  int num_correct, num_incorrect, num_true_selection, num_test_selection;
  double start_time = 0.0, total_time = 0.0;
  float *rate_positive_selection = nullptr, *rate_false_discovery = nullptr;

  if ( mpi_rank == 0 ) {
    printf("================================================================"
        "================================================================\n");

    // Initialize statistic data
    rate_positive_selection = new float[num_test];
    rate_false_discovery    = new float[num_test];

    // Create solution model
    Particle solution;
    bool btemp = true;
    for ( auto i = 0; i < p; i++ ) {
      if ( J0[i] ) {
        if ( btemp ) {
          solution.InitializeModel(i);
          btemp = false;
        } else {
          solution.UpdateModel(i);
        }
      }
    }
    solution.ComputeCriterion();

    // Display solution model
    num_true_selection = 0;
    auto isize = static_cast<int>(log10(p))+1;
    printf("True(%02d):\t%12.6f; ", mpi_rank, solution.phi);
    for ( auto i = 0; i < p; i++ ) {
      if ( J0[i] ) {
        printf("%-*d ", isize, i);
        num_true_selection++;
      }
    }
    printf("\n\n");
  } else {
    num_true_selection = 0;
    for ( auto i = 0; i < p; i++ ) {
      if ( J0[i] ) {
        num_true_selection++;
      }
    }
  }

  for ( auto t = 0u; t < num_test; ++t ) {
    // Record beginning time
    MPI_Barrier(MPI_COMM_WORLD);
    if ( mpi_rank == 0 ) {
      start_time = omp_get_wtime();
    }

    // Run PaSS
    GenLin();

    // Find best model
    struct { float value; int rank; } send, recv;
    send.value = phi0;
    send.rank = mpi_rank;
    MPI_Allreduce(&send, &recv, 1, MPI_FLOAT_INT, MPI_MINLOC, MPI_COMM_WORLD);

    // Record ending time
    MPI_Barrier(MPI_COMM_WORLD);
    if ( mpi_rank == 0 ) {
      total_time += omp_get_wtime() - start_time;
    }

    if ( mpi_rank == recv.rank ) {
      // Display model
      num_correct = 0;
      num_incorrect = 0;
      num_test_selection = 0;
      auto isize = static_cast<int>(log10(p))+1;
      printf("%4d(%02d):\t%12.6f; ", t, mpi_rank, phi0);
      for ( auto i = 0; i < p; i++ ) {
        if ( I0[i] ) {
          printf("%-*d ", isize, i);
          num_test_selection++;
          if ( J0[i] ) {
            num_correct++;
          } else {
            num_incorrect++;
          }
        }
      }
      printf("\n");

      // Compute accuracy rate
      if ( mpi_rank == 0 ) {
        rate_positive_selection[t] =
            static_cast<float>(num_correct) / num_true_selection;
        rate_false_discovery[t] =
            static_cast<float>(num_incorrect) / num_test_selection;
      } else {
        float rate_temp[2];
        rate_temp[0] = static_cast<float>(num_correct) / num_true_selection;
        rate_temp[1] = static_cast<float>(num_incorrect) / num_test_selection;
        MPI_Send(rate_temp, 2, MPI_FLOAT, 0, t, MPI_COMM_WORLD);
      }
    } else if ( mpi_rank == 0 ) {
      float rate_temp[2];
      MPI_Recv(rate_temp, 2, MPI_FLOAT, recv.rank, t,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      rate_positive_selection[t] = rate_temp[0];
      rate_false_discovery[t] = rate_temp[1];
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // Display statistic report                                               //
  ////////////////////////////////////////////////////////////////////////////
  if ( mpi_rank == 0 ) {
    printf("================================================================"
        "================================================================\n");

    // Compute means and standard deviations
    auto rate_positive_selection_mean =
        std::accumulate(rate_positive_selection,
                        rate_positive_selection+num_test,
                        0.0f) / num_test;
    auto rate_positive_selection_meansq =
        rate_positive_selection_mean * rate_positive_selection_mean;
    auto rate_positive_selection_sqmean =
        cblas_sdot(num_test, rate_positive_selection, 1,
                             rate_positive_selection, 1) / num_test;
    auto rate_positive_selection_sd =
        std::sqrt(rate_positive_selection_sqmean -
                  rate_positive_selection_meansq);
    auto rate_false_discovery_mean =
        std::accumulate(rate_false_discovery,
                        rate_false_discovery+num_test,
                        0.0f) / num_test;
    auto rate_false_discovery_meansq =
        rate_false_discovery_mean * rate_false_discovery_mean;
    auto rate_false_discovery_sqmean =
        cblas_sdot(num_test, rate_false_discovery, 1,
                             rate_false_discovery, 1) / num_test;
    auto rate_false_discovery_sd =
        std::sqrt(rate_false_discovery_sqmean -
                  rate_false_discovery_meansq);

    // Display statistic report
    printf("%s\n", dataname);
    printf("#Node      = %d\n", mpi_size);
    printf("#Thread    = %d\n", num_thread);
    printf("#Particle  = %d\n", num_particle);
    printf("#Iteration = %d\n", parameter.num_iteration);
    printf("pfg        = %.2f\n", parameter.prob_forward_global);
    printf("pfl        = %.2f\n", parameter.prob_forward_local);
    printf("pfr        = %.2f\n", parameter.prob_forward_random);
    printf("pbl        = %.2f\n", parameter.prob_backward_local);
    printf("pbr        = %.2f\n", parameter.prob_backward_random);
    if ( parameter.criterion == EBIC ) {
      printf("Criterion  = %s%.1f\n",
             Criterion2String(parameter.criterion),
             parameter.ebic_gamma);
    } else {
      printf("Criterion  = %s\n",
             Criterion2String(parameter.criterion));
    }
    printf("#Test      = %d\n", num_test);
    printf("PSR(mean)  = %.6f\n", rate_positive_selection_mean);
    printf("FDR(mean)  = %.6f\n", rate_false_discovery_mean);
    printf("PSR(sd)    = %.6f\n", rate_positive_selection_sd);
    printf("FDR(sd)    = %.6f\n", rate_false_discovery_sd);
    printf("Time       = %.6lf sec\n", total_time / num_test);
    printf("================================================================"
        "================================================================\n");
  }

  ////////////////////////////////////////////////////////////////////////////

  // Delete memory
  delete[] X0;
  delete[] Y0;
  delete[] J0;

  // Finalize MPI
  MPI_Finalize();

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Load data from file                                                        //
//                                                                            //
// Parameters:                                                                //
// fileroot: the root of data file                                            //
////////////////////////////////////////////////////////////////////////////////
void PassLoad( const char *fileroot ) {
  if ( mpi_rank == 0 ) {
    printf("Loading model from '%s'... ", fileroot);
  }

  // Open file
  auto file = fopen(fileroot, "rb");
  if ( !file ) {
    if ( mpi_rank == 0 ) {
      printf("Failed!\n");
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Read data
  int size;
  fread(&size, sizeof(int), 1, file);
  dataname = new char[size];
  fread(dataname, sizeof(char), size, file);
  fread(&n, sizeof(int), 1, file);
  fread(&p, sizeof(int), 1, file);
  X0 = new float[n*p];
  Y0 = new float[n];
  J0 = new bool[p];
  fread(X0, sizeof(float), n*p, file);
  fread(Y0, sizeof(float), n, file);
  fread(J0, sizeof(bool), p, file);

  // Close file
  fclose(file);

  if ( mpi_rank == 0 ) {
    printf("Done.\n");
  }
}
