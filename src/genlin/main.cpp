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
#include "pass.hpp"

using namespace pass;

// Global variables
bool *J0;        // vector, 1 by p, the chosen indices (solution)
int num_test;    // scalar, the number of tests
char *dataname;  // string, the name of data

using namespace pass;

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
  if ( num_test < 0 ) {
    printf("nT must be positive or zero!\n");
    exit(1);
  }

  // Load data and allocate memory
  PassLoad(dataroot);
  I0   = new bool[p];

  // Display parameters
  if ( parameter.criterion == EBIC ) {
    printf("%s: n=%d, p=%d, nP=%d, nI=%d, nT=%d, cri=%s, gamma=%.1lf\n",
           dataname, n, p, parameter.num_particle, parameter.num_iteration,
           num_test, Criterion2String(parameter.criterion),
           parameter.ebic_gamma);
  } else {
    printf("%s: n=%d, p=%d, nP=%d, nI=%d, nT=%d, cri=%s\n",
           dataname, n, p, parameter.num_particle, parameter.num_iteration,
           num_test, Criterion2String(parameter.criterion));
  }

  ////////////////////////////////////////////////////////////////////////////
  // Run PaSS                                                               //
  ////////////////////////////////////////////////////////////////////////////

  printf("================================================================\n");

  for ( auto i = 0; i < num_test; ++i ) {
    // Run PaSS
    GenLin();

    // Display model
    auto isize = static_cast<int>(log10(p))+1;
    for( i = 0; i < p; i++ ) {
      if( I0[i] ) {
        printf( "%-*d ", isize, i );
      }
    }
    printf( "\n" );
  }

  ////////////////////////////////////////////////////////////////////////////

  delete[] X0;
  delete[] Y0;
  delete[] J0;

  printf("================================================================\n");

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

  printf("Loading config from '%s'...\n", fileroot);

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
    else if ( !strcmp(cristr, "HDHQ") ) {
      parameter.criterion = HDHQ;
    }
    else {
      printf("There is no criterion named '%s'!\n", cristr);
      exit(1);
    }

    fgets(line, kBufferSize, file);
    sscanf(line, "%*s %d", &num_test);

    // Close file
    fclose(file);

    printf("Loaded config from '%s'.\n", fileroot);
  }
  else {
    printf("Unable to open file '%s'!\n", fileroot);

    // Open file
    file = fopen( fileroot, "w" );
    if ( !file ) {
      printf("Unable to create file '%s'!\n", fileroot);
      exit(1);
    }
    printf("Creating config file '%s'...\n", fileroot);

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
    fprintf(file, "        HDHQ:        "
                  "High-dimensional Hannan-Quinn information criterion.\n");
    fprintf(file, "<nT>    the number of tests.\n");

    // Close file
    fclose(file);

    printf("Created config file '%s'.\n", fileroot);
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

  printf("Loading model from '%s'...\n", fileroot);

  // Open file
  file = fopen(fileroot, "rb");
  if ( !file ) {
    printf("Unable to open file '%s'!\n", fileroot);
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

  printf("Loaded model into '%s'.\n", fileroot);
}
