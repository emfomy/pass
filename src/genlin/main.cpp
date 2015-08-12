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
#include <cstring>
#include <ctime>
#include <cmath>
#include <essl.h>
#include <mpi.h>
#include <omp.h>
#include "pass.hpp"

using namespace pass;

// Global variables
bool *J0;        // vector, 1 by p, the chosen indices (solution)
int num_test;    // scalar, the number of tests
char *dataname;  // string, the name of data
int world_size;  // the size of MPI communicator
int world_rank;  // the rank of MPI process

// Functions
void PassConfig( const char* fileroot );
void PassLoad( const char* fileroot );

////////////////////////////////////////////////////////////////////////////////
// Main function                                                              //
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char **argv ) {
  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Initialize random seed
  srand(time(NULL) ^ world_rank);
  srand(rand());

  // Load arguments
  auto cfgroot  = (argc > 1) ? argv[1] : "genlin.cfg";
  auto dataroot = (argc > 2) ? argv[2] : "genlin.dat";

  ////////////////////////////////////////////////////////////////////////////
  // Load parameters and data                                               //
  ////////////////////////////////////////////////////////////////////////////

  if ( world_rank == 0 ) {
    printf("================================================================\n");
  }

  // Load parameters
  num_test = 10;
  PassConfig(cfgroot);

  // Check parameters
  auto num_thread = omp_get_max_threads();
  auto num_proc = omp_get_num_procs();
  if ( num_thread > num_proc ) {
    omp_set_num_threads(num_proc);
    num_thread = num_proc;
  }

  // Load data and allocate memory
  PassLoad(dataroot);
  I0 = new bool[p];

  // Display parameters
  if ( world_rank == 0 ) {
    if ( parameter.criterion == EBIC ) {
      printf("%s: n=%d, p=%d, numNode=%d, numThread=%d, "
             "numParticle=%d, numIteration=%d, numTest=%d, criterion=%s%.1f\n",
             dataname, n, p, world_size, num_thread,
             world_size*num_thread, parameter.num_iteration, num_test,
             Criterion2String(parameter.criterion), parameter.ebic_gamma);
    } else {
      printf("%s: n=%d, p=%d, numNode=%d, numThread=%d, "
             "numParticle=%d, numIteration=%d, numTest=%d, criterion=%s\n",
             dataname, n, p, world_size, num_thread,
             world_size*num_thread, parameter.num_iteration, num_test,
             Criterion2String(parameter.criterion));
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // Centralize and normalize the original data                             //
  ////////////////////////////////////////////////////////////////////////////

  if ( world_rank == 0 ) {
    printf("Normalizing data... ");
  }

  // Centralize and normalize X0
  for ( auto i = 0; i < p; ++i ) {
    float stemp = 1.0f;
    stemp = sdot(n, X0+i*n, 1, &stemp, 0) / n;
    sves(n, X0+i*n, 1, &stemp, 0, X0+i*n, 1);
    sscal(n, (1.0f/snorm2(n, X0+i*n, 1)), X0+i*n, 1);
  }

  // Centralize and normalize Y0
  {
    float stemp = 1.0f;
    stemp = sdot(n, Y0, 1, &stemp, 0) / n;
    sves(n, Y0, 1, &stemp, 0, Y0, 1);
    sscal(n, (1.0f/snorm2(n, Y0, 1)), Y0, 1);
  }

  parameter.is_normalized = true;

  if ( world_rank == 0 ) {
    printf("Done.\n");
  }

  ////////////////////////////////////////////////////////////////////////////
  // Run PaSS                                                               //
  ////////////////////////////////////////////////////////////////////////////

  // Declare variables
  int num_correct, num_incorrect, num_true_selection, num_test_selection;
  double start_time, total_time = 0.0;
  float *rate_positive_selection, *rate_false_discovery;

  if ( world_rank == 0 ) {
    printf("================================================================\n");

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
    printf("True(%02d):\t%12.6f; ", world_rank, solution.phi);
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

  for ( auto t = 0; t < num_test; ++t ) {
    // Record beginning time
    MPI_Barrier(MPI_COMM_WORLD);
    if ( world_rank == 0 ) {
      start_time = omp_get_wtime();
    }

    // Run PaSS
    GenLin();

    // Find best model
    struct { float value; int rank; } send, recv;
    send.value = phi0;
    send.rank = world_rank;
    MPI_Allreduce(&send, &recv, 1, MPI_FLOAT_INT, MPI_MINLOC, MPI_COMM_WORLD);

    // Record ending time
    MPI_Barrier(MPI_COMM_WORLD);
    if ( world_rank == 0 ) {
      total_time += omp_get_wtime() - start_time;
    }

    if ( world_rank == recv.rank ) {
      // Display model
      num_correct = 0;
      num_incorrect = 0;
      num_test_selection = 0;
      auto isize = static_cast<int>(log10(p))+1;
      printf("%4d(%02d):\t%12.6f; ", t, world_rank, phi0);
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
      if ( world_rank == 0 ) {
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
    } else if ( world_rank == 0 ) {
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
  if ( world_rank == 0 ) {
    printf("================================================================\n");

    // Compute means and standard deviations
    float stemp = 1.0f;
    auto rate_positive_selection_mean =
        sdot(num_test, rate_positive_selection, 1, &stemp, 0)
        / num_test;
    auto rate_positive_selection_meansq =
        rate_positive_selection_mean * rate_positive_selection_mean;
    auto rate_positive_selection_sqmean =
        sdot(num_test, rate_positive_selection, 1, rate_positive_selection, 1)
        / num_test;
    auto rate_positive_selection_sd =
        sqrt(rate_positive_selection_sqmean - rate_positive_selection_meansq);
    auto rate_false_discovery_mean =
        sdot(num_test, rate_false_discovery, 1, &stemp, 0)
        / num_test;
    auto rate_false_discovery_meansq =
        rate_false_discovery_mean * rate_false_discovery_mean;
    auto rate_false_discovery_sqmean =
        sdot(num_test, rate_false_discovery, 1, rate_false_discovery, 1)
        / num_test;
    auto rate_false_discovery_sd =
        sqrt(rate_false_discovery_sqmean - rate_false_discovery_meansq);

    // Display statistic report
    printf("%s\n", dataname);
    printf("nNode     = %d\n", world_size);
    printf("nThread   = %d\n", num_thread);
    printf("nP        = %d\n", num_thread * world_size);
    printf("nI        = %d\n", parameter.num_iteration);
    printf("pfg       = %.2f\n", parameter.prob_forward_global);
    printf("pfl       = %.2f\n", parameter.prob_forward_local);
    printf("pfr       = %.2f\n", parameter.prob_forward_random);
    printf("pbl       = %.2f\n", parameter.prob_backward_local);
    printf("pbr       = %.2f\n", parameter.prob_backward_random);
    if ( parameter.criterion == EBIC ) {
      printf("cri       = %s%.1f\n",
             Criterion2String(parameter.criterion),
             parameter.ebic_gamma);
    } else {
      printf("cri       = %s\n",
             Criterion2String(parameter.criterion));
    }
    printf("nTest     = %d\n", num_test);
    printf("PSR(mean) = %.6f\n", rate_positive_selection_mean);
    printf("FDR(mean) = %.6f\n", rate_false_discovery_mean);
    printf("PSR(sd)   = %.6f\n", rate_positive_selection_sd);
    printf("FDR(sd)   = %.6f\n", rate_false_discovery_sd);
    printf("Time      = %.6lf sec\n", total_time / num_test);
    printf("================================================================\n");
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
// Load parameters from config file                                           //
//                                                                            //
// Input Parameters:                                                          //
// fileroot:  the root of config file                                         //
// parameter: the default parameters of the PaSS algorithm                    //
// num_test:  the default number of tests                                     //
//                                                                            //
// Output Parameters:                                                         //
// parameter: the parameters of the PaSS algorithm                            //
// num_test:  the number of tests                                             //
////////////////////////////////////////////////////////////////////////////////
void PassConfig( const char* fileroot ) {
  const int kBufferSize = 1024;

  if ( world_rank == 0 ) {
    printf("Loading config from '%s'... ", fileroot);
  }

  // Open file
  auto file = fopen(fileroot, "r");

  // Check if file exists
  if ( file ) {
    // Read data
    char line[kBufferSize];
    fgets(line, kBufferSize, file);
    sscanf(line, "%*s %d",
           &parameter.num_iteration);
    fgets(line, kBufferSize, file);
    sscanf(line, "%*s %f %f %f %f %f",
           &parameter.prob_forward_global,
           &parameter.prob_forward_local,
           &parameter.prob_forward_random,
           &parameter.prob_backward_local,
           &parameter.prob_backward_random);

    char cristr[kBufferSize];
    fgets(line, kBufferSize, file);
    sscanf(line, "%*s %s", cristr);
    if ( !strcmp(cristr, "AIC") ) {
      parameter.criterion = AIC;
    }
    else if ( !strcmp(cristr, "BIC") ) {
      parameter.criterion = BIC;
    }
    else if ( !strcmp(cristr, "EBIC") ) {
      parameter.criterion = EBIC;
      parameter.ebic_gamma = 1.0f;
    }
    else if ( cristr[0] == 'E' && cristr[1] == 'B' &&
             cristr[2] == 'I' && cristr[3] == 'C' ) {
      parameter.criterion = EBIC;
      parameter.ebic_gamma = atof(cristr+4);
    }
    else if ( !strcmp(cristr, "HDBIC") ) {
      parameter.criterion = HDBIC;
    }
    else if ( !strcmp(cristr, "HQC") ) {
      parameter.criterion = HQC;
    }
    else if ( !strcmp(cristr, "HDHQC") ) {
      parameter.criterion = HDHQC;
    }
    else {
      if ( world_rank == 0 ) {
        printf("Failed!\nThere is no criterion named '%s'!\n", cristr);
      }
      exit(1);
    }

    fgets(line, kBufferSize, file);
    sscanf(line, "%*s %d", &num_test);

    // Close file
    fclose(file);

    if ( world_rank == 0 ) {
      printf("Done.\n");
    }
  } else if ( world_rank == 0 ) {
    printf("Failed!\nCreating config file '%s'... ", fileroot);

    // Open file
    file = fopen( fileroot, "w" );
    if ( !file ) {
        printf("Failed!\n");
      exit(1);
    }

    // Write data
    fprintf(file, "nI    %d\n",
            parameter.num_iteration);
    fprintf(file, "prob  %.1f %.1f %.1f %.1f %.1f\n",
            parameter.prob_forward_global,
            parameter.prob_forward_local,
            parameter.prob_forward_random,
            parameter.prob_backward_local,
            parameter.prob_backward_random);
    if ( parameter.criterion == EBIC && parameter.ebic_gamma != 1.0f ) {
      fprintf(file, "cri   %s%.1f\n",
              Criterion2String(parameter.criterion),
              parameter.ebic_gamma);
    } else {
      fprintf(file, "cri   %s\n",
              Criterion2String(parameter.criterion));
    }
    fprintf(file, "nTest %d\n", num_test);

    fprintf(file, "\n\nNote:\n");
    fprintf(file, "<nI>:   the number of iterations.\n");
    fprintf(file, "<prob>: <pfg> <pfl> <pfr> <pbl> <pbr>\n");
    fprintf(file, "<pfg>:  the probability of forward step: global\n");
    fprintf(file, "<pfl>:  the probability of forward step: local\n");
    fprintf(file, "<pfr>:  the probability of forward step: random\n");
    fprintf(file, "<pbl>:  the probability of backward step: local\n");
    fprintf(file, "<pbr>:  the probability of backward step: random\n");
    fprintf(file, "<cri>:  the criterion.\n");
    fprintf(file, "        AIC:         "
                  "Akaike information criterion.\n");
    fprintf(file, "        BIC:         "
                  "Bayesian information criterion.\n");
    fprintf(file, "        EBIC:        "
                  "Extended Bayesian information criterion.\n");
    fprintf(file, "        EBIC<gamma>: "
                  "EBIC with parameter gamma.\n");
    fprintf(file, "        HDBIC:       "
                  "High-dimensional Bayesian information criterion.\n");
    fprintf(file, "        HQC:         "
                  "Hannan-Quinn information criterion.\n");
    fprintf(file, "        HDHQC:       "
                  "High-dimensional Hannan-Quinn information criterion.\n");
    fprintf(file, "<nTest> the number of tests.\n");

    // Close file
    fclose(file);

    printf("Done.\nUses default config.\n");
  }
}

////////////////////////////////////////////////////////////////////////////////
// Load data from file                                                        //
//                                                                            //
// Parameters:                                                                //
// fileroot: the root of data file                                            //
////////////////////////////////////////////////////////////////////////////////
void PassLoad( const char* fileroot ) {
  if ( world_rank == 0 ) {
    printf("Loading model from '%s'... ", fileroot);
  }

  // Open file
  auto file = fopen(fileroot, "rb");
  if ( !file ) {
    if ( world_rank == 0 ) {
      printf("Failed!\n");
    }
    exit(1);
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

  if ( world_rank == 0 ) {
    printf("Done.\n");
  }
}
