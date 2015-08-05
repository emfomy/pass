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
#include <omp.h>
#include "pass.hpp"

using namespace pass;

// Global variables
bool *J0;        // vector, 1 by p, the chosen indices (solution)
int num_test;    // scalar, the number of tests
char *dataname;  // string, the name of data

// The statistic structure
struct Statistic {
  int num_correct;                // the total number of correct selection
  int num_incorrect;              // the total number of incorrect selection
  int num_true_selection;         // the total selection of solution model
  int num_test_selection;         // the total selection of test model
  float rate_positive_selection;  // the positive selection rate
  float rate_false_discovery;     // the false discovery rate
};

// Functions
void PassConfig( const char* fileroot );
void PassLoad( const char* fileroot );

////////////////////////////////////////////////////////////////////////////////
// Main function                                                              //
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char **argv ) {
  auto cfgroot  = (argc > 1) ? argv[1] : "pass_genlin.cfg";
  auto dataroot = (argc > 2) ? argv[2] : "pass_genlin.dat";

  ////////////////////////////////////////////////////////////////////////////
  // Load parameters and data                                               //
  ////////////////////////////////////////////////////////////////////////////
  printf("================================================================\n");

  // Load parameters
  num_test = 10;
  PassConfig(cfgroot);

  // Check parameters
  auto num_thread = omp_get_max_threads();
  auto num_proc = omp_get_num_procs();
  if ( num_thread > num_proc ) {
    printf("No enough threads!\n");
    omp_set_num_threads(num_proc);
    num_thread = num_proc;
  }
  if ( parameter.num_particle <= 0 ) {
    parameter.num_particle = num_thread;
  }
  if ( num_test < 0 ) {
    printf("nT must be positive or zero!\n");
    exit(1);
  }

  // Load data and allocate memory
  PassLoad(dataroot);
  I0   = new bool[p];

  // Display parameters
  if ( parameter.criterion == EBIC ) {
    printf("%s: n=%d, p=%d, nP=%d, nI=%d, nthr=%d, nT=%d, cri=%s%.1f\n",
           dataname, n, p, parameter.num_particle, parameter.num_iteration,
           num_thread, num_test, Criterion2String(parameter.criterion),
           parameter.ebic_gamma);
  } else {
    printf("%s: n=%d, p=%d, nP=%d, nI=%d, nthr=%d, nT=%d, cri=%s\n",
           dataname, n, p, parameter.num_particle, parameter.num_iteration,
           num_thread, num_test, Criterion2String(parameter.criterion));
  }

  ////////////////////////////////////////////////////////////////////////////
  // Centralize and normalize the original data                             //
  ////////////////////////////////////////////////////////////////////////////

  printf("Normalizing data... ");

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

  printf("Done.\n");

  ////////////////////////////////////////////////////////////////////////////
  // Run PaSS                                                               //
  ////////////////////////////////////////////////////////////////////////////

  printf("================================================================\n");

  // Initialize statistic data
  int isize = static_cast<int>(log10(p))+1;
  double start_time, total_time = 0.0;
  Statistic statistic;
  statistic.rate_positive_selection = 0.0f;
  statistic.rate_false_discovery    = 0.0f;

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
  statistic.num_true_selection = 0;
  printf("True:\t%12.6f; ", solution.phi);
  for ( auto i = 0; i < p; i++ ) {
    if ( J0[i] ) {
      printf("%-*d ", isize, i);
      statistic.num_true_selection++;
    }
  }
  printf("\n\n");

  for ( auto t = 0; t < num_test; ++t ) {
    // Run PaSS
    start_time = omp_get_wtime();
    GenLin();
    total_time += omp_get_wtime() - start_time;

    // Display model
    statistic.num_correct = 0;
    statistic.num_incorrect = 0;
    statistic.num_test_selection = 0;
    printf("%4d:\t%12.6f; ", t, phi0);
    for ( auto i = 0; i < p; i++ ) {
      if ( I0[i] ) {
        printf("%-*d ", isize, i);
        statistic.num_test_selection++;
        if ( J0[i] ) {
          statistic.num_correct++;
        } else {
          statistic.num_incorrect++;
        }
      }
    }
    printf("\n");

    // Compute accuracy rate
    statistic.rate_positive_selection +=
        static_cast<float>(statistic.num_correct) /
        statistic.num_true_selection;
    statistic.rate_false_discovery +=
        static_cast<float>(statistic.num_incorrect) /
        statistic.num_test_selection;
  }

  printf("================================================================\n");

  ////////////////////////////////////////////////////////////////////////////
  // Display statistic report                                               //
  ////////////////////////////////////////////////////////////////////////////

  // Display statistic report
  printf("%s\n", dataname);
  printf("nthr  = %d\n", num_thread);
  printf("nP    = %d\n", parameter.num_particle);
  printf("nI    = %d\n", parameter.num_iteration);
  printf("pfg   = %.2f\n", parameter.prob_forward_global);
  printf("pfl   = %.2f\n", parameter.prob_forward_local);
  printf("pfr   = %.2f\n", parameter.prob_forward_random);
  printf("pbl   = %.2f\n", parameter.prob_backward_local);
  printf("pbr   = %.2f\n", parameter.prob_backward_random);
  if ( parameter.criterion == EBIC ) {
    printf("cri   = %s%.1f\n",
           Criterion2String(parameter.criterion),
           parameter.ebic_gamma);
  } else {
    printf("cri   = %s\n",
           Criterion2String(parameter.criterion));
  }
  printf("nT    = %d\n", num_test);
  printf("PSR   = %.6f\n", statistic.rate_positive_selection / num_test);
  printf("FDR   = %.6f\n", statistic.rate_false_discovery / num_test);
  printf("Time  = %.6lf sec\n", total_time / num_test);
  printf("================================================================\n");

  ////////////////////////////////////////////////////////////////////////////

  delete[] X0;
  delete[] Y0;
  delete[] J0;

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

  printf("Loading config from '%s'... ", fileroot);

  // Open file
  auto file = fopen(fileroot, "r");

  // Check if file exists
  if ( file ) {
    // Read data
    char line[kBufferSize];
    fgets(line, kBufferSize, file);
    sscanf(line, "%*s %d %d",
           &parameter.num_particle,
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
      printf("Failed!\nThere is no criterion named '%s'!\n", cristr);
      exit(1);
    }

    fgets(line, kBufferSize, file);
    sscanf(line, "%*s %d", &num_test);

    // Close file
    fclose(file);

    printf("Done.\n");
  }
  else {
    printf("Failed!\n");
    printf("Creating config file '%s'... ", fileroot);

    // Open file
    file = fopen( fileroot, "w" );
    if ( !file ) {
      printf("Failed!\n");
      exit(1);
    }

    // Write data
    fprintf(file, "nP/nI %d %d\n",
            parameter.num_particle,
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
    fprintf(file, "nT    %d\n", num_test);

    fprintf(file, "\n\nNote:\n");
    fprintf(file, "<nP>:   the number of particles.\n");
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
    fprintf(file, "        HQC:       "
                  "Hannan-Quinn information criterion.\n");
    fprintf(file, "        HDHQC:       "
                  "High-dimensional Hannan-Quinn information criterion.\n");
    fprintf(file, "<nT>    the number of tests.\n");

    // Close file
    fclose(file);

    printf("Done.\n");
    printf("Uses default config.\n");
  }
}

////////////////////////////////////////////////////////////////////////////////
// Load data from file                                                        //
//                                                                            //
// Parameters:                                                                //
// fileroot: the root of data file                                            //
////////////////////////////////////////////////////////////////////////////////
void PassLoad( const char* fileroot ) {
  FILE *file;
  int size;

  printf("Loading model from '%s'... ", fileroot);

  // Open file
  file = fopen(fileroot, "rb");
  if ( !file ) {
    printf("Failed!\n");
    exit(1);
  }

  // Write data
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

  printf("Done.\n");
}
