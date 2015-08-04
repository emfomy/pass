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
#include <cblas.h>
#include <lapacke.h>

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
  srand(time(NULL));
  int iseed[4] = {rand()%4096, rand()%4096, rand()%4096, (rand()%2048)*2+1};

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

  printf("Creating a linear data using Ing and Lai's method... ");

  // Allocate memory
  X = new float[n*p];
  Y = new float[n];
  J = new bool[p]();
  auto S = new float[n];
  float stemp;

  // Generate X & Y using normal random
  LAPACKE_slarnv(3, iseed, n*p, X);
  LAPACKE_slarnv(3, iseed, n, Y);

  // S[i] := sum( X[i, 0~r] )
  stemp = 1.0f;
  for ( auto i = 0; i < n; ++i ) {
    S[i] = cblas_sdot(r, X+i*n, n, &stemp, 0);
  }

  // X[i col] := sqrt(.75/r) * S + .5 * X[i col], i >= r
  stemp = sqrt(0.75f/r);
  for ( auto i = 0; i < r; ++i ) {
    cblas_sscal(n, 0.5f, X+i*n, 1);
    cblas_saxpy(n, stemp, S, 1, X+i*n, 1);
  }

  // Y += X[0~r cols] * Beta
  cblas_sgemv(CblasColMajor, CblasNoTrans,
              n, r, 1.0f, X, n, Beta, 1, 1.0f, Y, 1);

  // Generate J
  memset(J, true, sizeof(bool) * r);

  printf("Done.\n");

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

  printf("Loading config from '%s'... ", fileroot);

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

    printf("Done.\n");
  } else {
    printf("Failed!\n");
    printf("Creating config file '%s'... ", fileroot);

    // Open file
    file = fopen(fileroot, "w");
    if ( !file ) {
      printf("Failed!\n");
      exit(1);
    }
  
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

    printf("Done.\n");
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

  printf("Saving model into '%s'... ", fileroot);

  // Open file
  file = fopen(fileroot, "wb");
  if ( !file ) {
    printf("Failed!\n");
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

  printf("Done.\n");
}
