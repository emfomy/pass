////////////////////////////////////////////////////////////////////////////////
// Particle Swarm Stepwise (PaSS) Algorithm                                   //
//                                                                            //
// genlin_inglai.cpp                                                          //
// Create a linear data using Ing and Lai's method                            //
//                                                                            //
// Author: emfo<emfomy@gmail.com>                                             //
//                                                                            //
// Reference:                                                                 //
// Ing, C.-K., & Lai, T. L. (2011, October). A stepwise regression method and //
//   consistent model selection for high-dimensional sparse linear models.    //
//   http://doi.org/10.5705/ss.2010.081                                       //
////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <unistd.h>
#include <cblas.h>
#include <lapacke.h>

// Global variables
int n;                 // scalar, the number of statistical units
int p;                 // scalar, the number of total effects
int r;                 // scalar, the number of given effects
float *X;              // matrix, n by p, the regressors
float *Y;              // vector, n by 1, the regressand
float *Beta;           // vector, r by 1, the effects
bool *J;               // vector, 1 by p, the chosen indices
const char *dataname;  // string, the name of data

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

  // Initialize arguments
  auto cfgroot  = "genlin_inglai.cfg";
  auto dataroot = "genlin.dat";
  dataname      = "GenLin_IngLai";

  // Load arguments
  char c;
  bool input_error = false;
  opterr = false;
  while ( (c = getopt(argc, argv, "c:d:h")) != static_cast<char>(EOF) ) {
    switch ( c ) {
      case 'c': {
        cfgroot = optarg;
        break;
      }
      case 'd': {
        dataroot = optarg;
        break;
      }
      case 'h': {
        input_error = true;
        break;
      }
      default: {
        printf("invalid option -- '%c'\n", optopt);
        input_error = true;
        break;
      }
    }
  }
  if ( input_error ) {
    printf("Usage: %s [options] ...\n", argv[0]);
    printf("-c <file>                       Read config from <file>.\n");
    printf("-d <file>                       Save data into <file>.\n");
    return 0;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Load parameters                                                        //
  ////////////////////////////////////////////////////////////////////////////

  printf("================================================================\n");

  // Load parameters
  IngLaiConfig(cfgroot);

  // Display parameters
  printf("\nn = %d, p = %d, r = %d\n", n, p, r);
  printf("Beta: ");
  for ( auto i = 0; i < r; i++ ) {
    printf("%8.3f", Beta[i]);
  }
  printf("\n\n");

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
  for ( auto i = 0; i < r; ++i ) {
    S[i] = cblas_sdot(r, X+i*n, n, &stemp, 0);
  }

  // X[i col] := sqrt(.75/r) * S + .5 * X[i col], i >= r
  stemp = sqrt(0.75f/r);
  for ( auto i = r; i < p; ++i ) {
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
  delete[] S;

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
    printf("Failed!\nCreating config file '%s'... ", fileroot);

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

    printf("Done.\nUses default config.\n");
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

  printf("Saving data into '%s'... ", fileroot);

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
  fwrite(J, sizeof(bool), p, file);

  // Close file
  fclose(file);

  printf("Done.\n");
}
