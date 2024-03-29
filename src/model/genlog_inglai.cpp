////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @file    model/genlog_inglai.cpp
/// @brief   Create a general logistic regression data using Ing and Lai's method
///
/// @author  Mu Yang <emfomy@gmail.com>

// ===========================================================================================================================
// Reference:
// Ing, C.-K., & Lai, T. L. (2011).
//   A stepwise regression method and consistent model selection for high-dimensional sparse linear models.
//   http://doi.org/10.5705/ss.2010.081
//

#include <iostream>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <getopt.h>
#include <mkl.h>

// Default arguments
const int kN          = 400;                        ///< the default value of n
const int kP          = 4000;                       ///< the default value of p
const int kR          = 10;                         ///< the default value of r
const char *kDataRoot = "genlog.dat";               ///< the default data file root
const char *kDataName = "General_Logistic_IngLai";  ///< the default data name

// Global variables
int n;                                              ///< scalar, the number of statistical units
int p;                                              ///< scalar, the number of total effects
int r;                                              ///< scalar, the number of given effects
float *X;                                           ///< matrix, n by p, the regressors
float *Y;                                           ///< vector, n by 1, the regressand
float *Beta;                                        ///< vector, r by 1, the effects
bool *J;                                            ///< vector, 1 by p, the chosen indices
const char *dataroot;                               ///< string, the root of the data file
const char *dataname;                               ///< string, the name of the data

// Functions
void GenLogIngLaiHelp( const char *cmd );
void GenLogIngLaiSave( const char *fileroot );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Main function
///
int main( int argc, char **argv ) {

  // ======== Initialize variables ======================================================================================== //

  // Initialize random generator
  srand(time(NULL));
  int iseed[4] = {rand()%4096, rand()%4096, rand()%4096, (rand()%2048)*2+1};

  // Initialize arguments
  n        = kN;
  p        = kP;
  r        = kR;
  dataroot = kDataRoot;
  dataname = kDataName;

  // Load arguments
  int optidx = 0;
  bool bflag = false;
  char opts[] = "f:m:n:p:r:bh", c;
  option long_opts[] = {
    {"file", required_argument, nullptr, 'f'},
    {"name", required_argument, nullptr, 'm'},
    {"beta", no_argument,       nullptr, 'b'},
    {"help", no_argument,       nullptr, 'h'}
  };
  while ( (c = getopt_long(argc, argv, opts, long_opts, &optidx)) != -1 ) {
    switch ( c ) {
      case 0: {
        break;
      }
      case 'f': {
        dataroot = optarg;
        break;
      }
      case 'm': {
        dataname = optarg;
        break;
      }
      case 'n': {
        n = atoi(optarg);
        if ( n <= 0 ) {
          fprintf(stderr, "%s: invalid option -- <n> must be a positive integer!\n", argv[0]);
          GenLogIngLaiHelp(argv[0]);
          abort();
        }
        break;
      }
      case 'p': {
        p = atoi(optarg);
        if ( p <= 0 ) {
          fprintf(stderr, "%s: invalid option -- <p> must be a positive integer!\n", argv[0]);
          GenLogIngLaiHelp(argv[0]);
          abort();
        }
        break;
      }
      case 'r': {
        r = atoi(optarg);
        if ( r < 0 ) {
          fprintf(stderr, "%s: invalid option -- <r> must be a non-negative integer!\n", argv[0]);
          GenLogIngLaiHelp(argv[0]);
          abort();
        }
        break;
      }
      case 'b': {
        bflag = true;
        break;
      }
      case 'h': {
        GenLogIngLaiHelp(argv[0]);
        exit(0);
      }
      default: {
        GenLogIngLaiHelp(argv[0]);
        abort();
      }
    }
  }

  // Create Beta
  if ( bflag ) {
    r = argc - optind;
    Beta = static_cast<float*>(mkl_malloc(r * sizeof(float), 64));
    for ( auto i = 0; i < r; ++i ) {
      Beta[i] = atof(argv[i+optind]);
    }
  } else {
    Beta = static_cast<float*>(mkl_malloc(r * sizeof(float), 64));
    for ( auto i = 0; i < r; ++i ) {
      Beta[i] = 0.3f + 0.075f*i;
    }
  }

  // Display arguments
  printf("================================================================"
         "================================================================\n");
  printf("n = %d, p = %d, r = %d\n", n, p, r);
  printf("Beta: ");
  for ( auto i = 0; i < r; ++i ) {
    printf("%8.3f", Beta[i]);
  }
  printf("\n\n");

  // ======== Create data ================================================================================================= //

  printf("Creating a linear data using Ing and Lai's method... ");
  fflush(stdout);

  // Allocate memory
  float *S;
  X = static_cast<float*>(mkl_malloc(n * p * sizeof(float), 64));
  Y = static_cast<float*>(mkl_malloc(n     * sizeof(float), 64));
  S = static_cast<float*>(mkl_calloc(n     , sizeof(float), 64));
  J = static_cast<bool* >(mkl_malloc(p     * sizeof(bool),  64));

  // Generate X & Y using normal random
  LAPACKE_slarnv(3, iseed, n*p, X);
  LAPACKE_slarnv(3, iseed, n, Y);

  // S := sum( X[0~r cols] )
  for ( auto i = 0; i < r; ++i ) {
    cblas_saxpy(n, 1.0f, X+i*n, 1, S, 1);
  }

  // X[i col] := sqrt(0.75/r) * S + 0.5 * X[i col], i >= r
  cblas_sscal(n, sqrt(0.75f/r), S, 1);
  cblas_sscal(n*(p-r), 0.5f, X+r*n, 1);
  for ( auto i = r; i < p; ++i ) {
    cblas_saxpy(n, 1.0f, S, 1, X+i*n, 1);
  }

  // S = X[0~r cols] * Beta
  cblas_sgemv(CblasColMajor, CblasNoTrans, n, r, 1.0f, X, n, Beta, 1, 0.0f, S, 1);

  // S := exp(S)
  vsExp(n, S, S);

  // S := S ./ (1+S)
  vsLinearFrac(n, S, S, 1.0, 0.0, 1.0, 1.0, S);

  // Compute Y
  for( auto i = 0; i < n; i++ ) {
    Y[i] = ( Y[i] < S[i] );
  }

  // Generate J
  memset(J, false, sizeof(bool) * p);
  memset(J, true, sizeof(bool) * r);

  printf("Done.\n");

  // ====================================================================================================================== //

  // Save data
  GenLogIngLaiSave(dataroot);

  // Free memory
  mkl_free(X);
  mkl_free(Y);
  mkl_free(Beta);
  mkl_free(S);
  mkl_free(J);

  printf("================================================================"
         "================================================================\n");

  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Display help messages
///
/// @param  cmd  the command name
///
void GenLogIngLaiHelp( const char *cmd ) {
  printf("Usage: %s [options] ...\n", cmd);
  printf("\n%-32s%-40s%s\n\n", "Option",                   "Detail",                          "Defalut Value");
  printf("%-32s%-40s%s\n",     "-f <file>, --file <file>", "save data into <file>",           kDataRoot);
  printf("%-32s%-40s%s\n",     "-m <name>, --name <name>", "set the data name as <name>",     kDataName);
  printf("%-32s%-40s\n",       "-b <beta>, --beta <beta>", "set the effects as <beta>s");
  printf("%-32s%-40s%d\n",     "-n ###",                   "the number of statistical units", kN);
  printf("%-32s%-40s%d\n",     "-p ###",                   "the number of total effects",     kP);
  printf("%-32s%-40s%d\n",     "-r ###",                   "the number of given effects",     kR);
  printf("%-32s%-40s\n",       "",                         "ignored if '-b' is set");
  printf("%-32s%-40s\n",       "-h, --help",               "display help messages");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Save data info file
///
/// @param  fileroot  the root of data file
///
void GenLogIngLaiSave( const char *fileroot ) {
  FILE *file;

  printf("Saving data into '%s'... ", fileroot);
  fflush(stdout);

  // Open file
  file = fopen(fileroot, "w");
  if ( !file ) {
    printf("Failed!\n");
    abort();
  }

  // Write data
  fprintf(file, "# 1st  line:  data name\n");
  fprintf(file, "# 2st  line:  n p\n");
  fprintf(file, "# 3rd  line:  * J\n");
  fprintf(file, "# rest lines: Y X\n");
  fprintf(file, "# \n");
  fprintf(file, "# X: matrix, n by p, the regressors\n");
  fprintf(file, "# Y: vector, n by 1, the regressand\n");
  fprintf(file, "# J: vector, 1 by p, the chosen indices\n");
  fprintf(file, "# \n");

  fprintf(file, "%s\n", dataname);
  fprintf(file, "%d %d\n", n, p);

  fprintf(file, "%-16c\t", '*');
  for ( auto j = 0; j < p; ++j ) {
    fprintf(file, "%-16d", J[j]);
  }
  fprintf(file, "\n");

  for ( auto i = 0; i < n; ++i) {
    fprintf(file, "%-+16.6e", Y[i]);
    for ( auto j = 0; j < p; ++j ) {
      fprintf(file, "%-+16.6e", X[i+j*n]);
    }
    fprintf(file, "\n");
  }

  // Close file
  fclose(file);

  printf("Done.\n");
}
