////////////////////////////////////////////////////////////////////////////////
// Particle Swarm Stepwise (PaSS) Algorithm                                   //
//                                                                            //
// main.cpp                                                                   //
// The main functions of PaSS for general linear regression                   //
//                                                                            //
// Author: Mu Yang <emfomy@gmail.com>                                         //
////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cstring>
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
  printf("  %-46s%-48s\n",
         "--brief",
         "switch to brief mode (default)");
  printf("  %-46s%-48s\n",
         "--verbose",
         "switch to verbose mode");
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

  ////////////////////////////////////////////////////////////////////////////
  //  Load arguments                                                        //
  ////////////////////////////////////////////////////////////////////////////

  // Initialize arguments
  auto num_test = kTest;
  auto dataroot = kDataRoot;
  int verbose   = 0;

  // Set arguments
  int optidx = 0;
  char opts[] = "f:i:p:t:\1\2\3\4::\5\6\7h", c;
  opterr = (mpi_rank == 0);
  option long_opts[] = {
    {"file",      required_argument, nullptr, 'f'},
    {"iteration", required_argument, nullptr, 'm'},
    {"particle",  required_argument, nullptr, 't'},
    {"test",      required_argument, nullptr, 'o'},
    {"brief",     no_argument,       &verbose, 0},
    {"verbose",   no_argument,       &verbose, 1},
    {"help",      no_argument,       nullptr, 'h'},
    {"prob",      no_argument,       nullptr, '\1'},
    {"AIC",       no_argument,       nullptr, '\2'},
    {"BIC",       no_argument,       nullptr, '\3'},
    {"EBIC",      optional_argument, nullptr, '\4'},
    {"HDBIC",     no_argument,       nullptr, '\5'},
    {"HQC",       no_argument,       nullptr, '\6'},
    {"HDHQC",     no_argument,       nullptr, '\7'}
  };

  // Load arguments
  while ( (c = getopt_long(argc, argv, opts, long_opts, &optidx)) != -1 ) {
    switch ( c ) {
      case 0: {
        break;
      }
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

  // Check arguments
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
    fflush(stdout);
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
  int num_true_selection = 0;
  double start_time = 0.0, total_time = 0.0;
  float *rate_positive_selection = nullptr, *rate_false_discovery = nullptr;
  Particle particle;

  if ( mpi_rank == 0 ) {
    printf("================================================================"
        "================================================================\n");

    // Initialize statistic data
    rate_positive_selection = new float[num_test];
    rate_false_discovery    = new float[num_test];

    // Build solution model
    bool btemp = true;
    for ( auto i = 0; i < p; i++ ) {
      if ( J0[i] ) {
        if ( btemp ) {
          particle.InitializeModel(i);
          btemp = false;
        } else {
          particle.UpdateModel(i);
        }
      }
    }
    if ( btemp ) {
      particle.k = 0;
      cblas_scopy(n, Y0, 1, particle.R, 1);
    }
    particle.ComputeCriterion();

    // Display solution model
    auto isize = static_cast<int>(log10(p))+1;
    printf("True(**):\t%12.6f; ", particle.phi);
    for ( auto i = 0; i < p; i++ ) {
      if ( J0[i] ) {
        printf("%-*d ", isize, i);
        num_true_selection++;
      } else if (verbose) {
        printf("%*s ", isize, "");
      }
    }
    printf("\n\n");
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

    // Transfer best model
    if ( recv.rank != 0 ) {
      if ( mpi_rank == recv.rank ) {
        MPI_Send(&phi0, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(I0, p, MPI_BYTE, 0, 1, MPI_COMM_WORLD);
      } else if ( mpi_rank == 0 ) {
        MPI_Recv(&phi0, 1, MPI_FLOAT, recv.rank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(I0, p, MPI_BYTE, recv.rank, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }

    if ( mpi_rank == 0 ) {
      // Build best model
      bool btemp = true;
      for ( auto i = 0; i < p; i++ ) {
        if ( I0[i] ) {
          if ( btemp ) {
            particle.InitializeModel(i);
            btemp = false;
          } else {
            particle.UpdateModel(i);
          }
        }
      }
      if ( btemp ) {
        particle.k = 0;
        cblas_scopy(n, Y0, 1, particle.R, 1);
      }
      particle.ComputeCriterion();

      // Display best model
      auto num_correct = 0;
      auto num_incorrect = 0;
      auto num_test_selection = 0;
      auto isize = static_cast<int>(log10(p))+1;
      printf("%4d(%02d):\t%12.6f; ", t, recv.rank, particle.phi);
      for ( auto i = 0; i < p; i++ ) {
        if ( I0[i] ) {
          printf("%-*d ", isize, i);
          num_test_selection++;
          if ( J0[i] ) {
            num_correct++;
          } else {
            num_incorrect++;
          }
        } else if (verbose) {
          printf("%*s ", isize, "");
        }
      }
      printf("\n");

      // Compute accuracy rate
      rate_positive_selection[t] =
          static_cast<float>(num_correct) / num_true_selection;
      rate_false_discovery[t] =
          static_cast<float>(num_incorrect) / num_test_selection;
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
    vsLinearFrac(num_test, rate_positive_selection, rate_positive_selection,
                 1.0f, -rate_positive_selection_mean, 0.0f, 1.0f,
                 rate_positive_selection);
    auto rate_positive_selection_sd =
        cblas_snrm2(num_test, rate_positive_selection, 1);
    auto rate_false_discovery_mean =
        std::accumulate(rate_false_discovery,
                        rate_false_discovery+num_test,
                        0.0f) / num_test;
    vsLinearFrac(num_test, rate_false_discovery, rate_false_discovery,
                 1.0f, -rate_false_discovery_mean, 0.0f, 1.0f,
                 rate_false_discovery);
    auto rate_false_discovery_sd =
        cblas_snrm2(num_test, rate_false_discovery, 1);

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
    fflush(stdout);
  }

  // Open file
  auto file = fopen(fileroot, "r");
  if ( !file ) {
    if ( mpi_rank == 0 ) {
      printf("Failed!\n");
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Skip comments
  const int kLineSize = 4096;
  char line[kLineSize];
  do {
    fgets(line, kLineSize, file);
  } while (line[0] == '#');

  // Read data name
  auto size = strlen(line);
  dataname = new char[size];
  memcpy(dataname, line, size-1);
  dataname[size-1] = '\0';

  // Read data size
  fscanf(file, "%d%d\n", &n, &p);

  // Alloocate memory
  X0 = new float[n*p];
  Y0 = new float[n];
  J0 = new bool[p];

  // Read J0
  fscanf(file, "%*s\n");
  for ( auto j = 0; j < p; ++j ) {
    int itemp;
    fscanf(file, "%d", &itemp);
    J0[j] = itemp;
  }

  for ( auto i = 0; i < n; ++i) {
    fscanf(file, "%f", Y0+i);
    for ( auto j = 0; j < p; ++j ) {
      fscanf(file, "%f", X0+i+j*n);
    }
  }

  // Close file
  fclose(file);

  if ( mpi_rank == 0 ) {
    printf("Done.\n");
  }
}
