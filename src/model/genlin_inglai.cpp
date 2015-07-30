////////////////////////////////////////////////////////////////////////////////
// Particle Swarm Stepwise (PaSS) Algorithm                                   //
//                                                                            //
// genlin_inglai.cpp                                                          //
// Create a linear data using Ing and Lai's method                            //
//                                                                            //
// Author: emfo<emfomy@gmail.com>                                             //
////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstring>
#include <ctime>
#include <random>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace boost;

// Global variables
int n;           // scalar, the number of statistical units
int p;           // scalar, the number of total effects
int r;           // scalar, the number of given effects
float *X;        // matrix, n by p, the regressors
float *Y;        // vector, n by 1, the regressand
float *Beta;     // vector, r by 1, the effects
bool *J;         // vector, 1 by p, the chosen indices
char *dataname;  // string, the name of data

// Functions
void IngLaiConfig( const char* fileroot );
void IngLaiSave( const char* fileroot );

////////////////////////////////////////////////////////////////////////////////
// Main function                                                              //
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char **argv ) {

  ////////////////////////////////////////////////////////////////////////////
  // Initialize variables                                                   //
  ////////////////////////////////////////////////////////////////////////////

  // Initialize random generator
  mt19937 *rng = new mt19937();
  rng->seed(time(NULL));

  normal_distribution<> distribution(1.0f, 1.0f);
  variate_generator<mt19937, normal_distribution<>> dist(*rng, distribution);

  // Initialize variables
  n = 400;
  p = 4000;
  r = 10;

  auto cfgroot  = (argc > 1) ? argv[1] : "pass_genlin_inglai.cfg";
  auto dataroot = (argc > 2) ? argv[2] : "pass_genlin.dat";

  auto strtemp = "GenLin_IngLai";
  dataname = new char[strlen(strtemp)+1];
  memcpy(dataname, strtemp, sizeof(char) * (strlen(strtemp)+1));

  if ( argc > 1 ) cfgroot  = argv[1];
  if ( argc > 2 ) dataroot = argv[2];

  ////////////////////////////////////////////////////////////////////////////
  // Load parameters                                                        //
  ////////////////////////////////////////////////////////////////////////////

  printf("================================================================\n");

  // Load parameters
  IngLaiConfig(cfgroot);

  // Display parameters
  printf( "\nn = %d, p = %d, r = %d\n", n, p, r );
  printf( "Beta: " );
  for ( auto i = 0; i < r; i++ ) {
    printf( "%8.3f", Beta[i] );
  }
  printf( "\n\n" );

  ////////////////////////////////////////////////////////////////////////////
  // Create data                                                            //
  ////////////////////////////////////////////////////////////////////////////

  printf("Creating a linear data using Ing and Lai's method...\n");

  // Allocate memory
  X = new float[n*p];
  Y = new float[n];
  J = new bool[p];
  auto S = new float[n];

  // Generate X & Y using normal random
  for ( auto i = 0; i < n*p; ++i ) {
    X[i] = dist();
  }
  for ( auto i = 0; i < n; ++i ) {
    Y[i] = dist();
  }

  // S[i] := sum( X[i, 0~r] )
  for ( auto i = 0; i < n; ++i ) {
    for ( auto j = 0; j < r; ++j ) {
      S[i] += X[i+j*n];
    }
  }

  // X[i col] := sqrt(.75/r) * S + .5 * X[i col], i >= r
  // Y += X[0~r cols] * Beta
  float dtemp = sqrt(0.75f/r);
  for ( auto i = r; i < p; ++i ) {
    for ( auto j = 0; j < n; ++j ) {
      Y[j] += Beta[i] * (dtemp * S[j] * 0.5f * X[i*n+j]);
    }
  }

  // Generate J
  memset(J, false, sizeof(bool) * p);
  for ( auto i = 0; i < r; ++i ) {
    J[i] = true;
  }

  ////////////////////////////////////////////////////////////////////////////

  // Save data
  IngLaiSave(dataroot);

  // Free memory
  delete[] X;
  delete[] Y;
  delete[] Beta;
  delete[] J;

  printf("================================================================\n");

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Load parameters from config file                                           //
//                                                                            //
// Parameters:                                                                //
// fileroot: the root of config file                                          //
////////////////////////////////////////////////////////////////////////////////
void IngLaiConfig( const char* fileroot ) {
  const int kBufferSize = 1024;

  printf("Loading config from '%s'...\n", fileroot);

  // Open file
  auto file = fopen(fileroot, "r");

  // Check if file exists
  if ( file ) {
    int offset;
    char line[kBufferSize], *linetemp;

    // Read data
    fgets(line, kBufferSize, file);
    sscanf(line, "%*s %d %d %d", &n, &p, &r);
    Beta = new float[r];
    fgets(line, kBufferSize, file);
    sscanf(line, "%*s %n", &offset);
    linetemp = line;
    for ( auto i = 0; i < r; i++ ) {
      linetemp += offset;
      sscanf(linetemp, "%f %n", &Beta[i], &offset);
    }

    // Close file
    fclose(file);

    printf("Loaded config from '%s'.\n", fileroot);
  } else {
    printf("Unable to open file '%s'!\n", fileroot);

    // Open file
    file = fopen(fileroot, "w");
    if ( !file ) {
      printf("Unable to create file '%s'!\n", fileroot);
      exit(1);
    }
    printf("Creating config file '%s'...\n", fileroot);
  
    // Generate Beta
    Beta = new float[r];
    for ( auto i = 0; i < r; i++ ) {
      Beta[i] = 3.0f + .75f*i;
    }

    // Write data
    fprintf(file, "n/p/r %d %d %d\n", n, p, r);
    fprintf(file, "Beta ");
    for ( auto i = 0; i < r; i++ ) {
      fprintf(file, " %.3f", Beta[i]);
    }
    fprintf(file, "\n");

    // Close file
    fclose(file);

    printf("Created config file '%s'.\n", fileroot);
    printf("Uses default config.\n");
  }
}

////////////////////////////////////////////////////////////////////////////////
// Save data into file                                                        //
//                                                                            //
// Parameters:                                                                //
// fileroot: the root of data file                                            //
////////////////////////////////////////////////////////////////////////////////
void IngLaiSave( const char* fileroot ) {
  FILE *file;
  int size = strlen(dataname)+1;

  printf("Saving model into '%s'...\n", fileroot);

  // Open file
  file = fopen(fileroot, "wb");
  if ( !file ) {
    printf("Unable to open file '%s'!\n", fileroot);
    exit(1);
  }

  // Write data
  fwrite(&size, sizeof(int), 1, file);
  fwrite(dataname, sizeof(char), size, file);
  fwrite(&n, sizeof(int), 1, file);
  fwrite(&p, sizeof(int), 1, file);
  fwrite(X, sizeof(float), n * p, file);
  fwrite(Y, sizeof(float), n, file);
  fwrite(J, sizeof(int), p, file);

  // Close file
  fclose(file);

  printf("Saved model into '%s'.\n", fileroot);
}
