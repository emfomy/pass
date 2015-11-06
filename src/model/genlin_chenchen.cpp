////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @file    model/genlin_chenchen.cpp
/// @brief   Create a general linear regression data using Chen and Chen's method
///
/// @author  Mu Yang <emfomy@gmail.com>

// ===========================================================================================================================
// Reference:
// Chen, J., & Chen, Z. (2008). Extended Bayesian information criteria for model selection with large model spaces.
//   Biometrika, 95(3), 759–771. http://www.stat.ubc.ca/~jhchen/paper/Bio08.pdf
//

#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <getopt.h>
#include <mkl.h>

// Default arguments
const int kN          = 200;                        ///< the default value of n
const int kP          = 50;                         ///< the default value of p
const int kR          = 8;                          ///< the default value of r
const int kType       = 3;                          ///< the default type
const float kRho      = 0.2f;                       ///< the default value of rho
const char *kDataRoot = "genlin.dat";               ///< the default data file root
const char *kDataName = "General_Linear_ChenChen";  ///< the default data name

// Global variables
int n;                                              ///< scalar, the number of statistical units
int p;                                              ///< scalar, the number of total effects
int r;                                              ///< scalar, the number of given effects
int type;                                           ///< scalar, the type of covariance structure, 1~3
float *X;                                           ///< matrix, n by p, the regressors
float *Y;                                           ///< vector, n by 1, the regressand
float *Beta;                                        ///< vector, r by 1, the effects
float rho;                                          ///< scaler, the covariance parameter
bool *J;                                            ///< vector, 1 by p, the chosen indices
const char *dataroot;                               ///< string, the root of the data file
const char *dataname;                               ///< string, the name of data
const char *suffix;                                 ///< string, the suffix of the name of data

// Functions
void GenLinChenChenHelp( const char *cmd );
void GenLinChenChenSave( const char *fileroot );

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
  type     = kType;
  rho      = kRho;
  dataroot = kDataRoot;
  dataname = kDataName;

  // Load arguments
  int optidx = 0;
  bool bflag = false;
  char opts[] = "f:m:n:p:r:t:o:bh", c;
  option long_opts[] = {
    {"file", required_argument, nullptr, 'f'},
    {"name", required_argument, nullptr, 'm'},
    {"type", required_argument, nullptr, 't'},
    {"cov",  required_argument, nullptr, 'c'},
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
          GenLinChenChenHelp(argv[0]);
          exit(1);
        }
        break;
      }
      case 'p': {
        p = atoi(optarg);
        if ( p <= 0 ) {
          fprintf(stderr, "%s: invalid option -- <p> must be a positive integer!\n", argv[0]);
          GenLinChenChenHelp(argv[0]);
          exit(1);
        }
        break;
      }
      case 'r': {
        r = atoi(optarg);
        if ( r < 0 ) {
          fprintf(stderr, "%s: invalid option -- <r> must be a non-negative integer!\n", argv[0]);
          GenLinChenChenHelp(argv[0]);
          exit(1);
        }
        break;
      }
      case 't': {
        type = atoi(optarg);
        if ( type < 1 || type > 3 ) {
          fprintf(stderr, "%s: invalid option -- <type> must be 1, 2, or 3!\n", argv[0]);
          GenLinChenChenHelp(argv[0]);
          exit(1);
        }
        break;
      }
      case 'c': {
        rho = atof(optarg);
        if ( rho < 0 || rho > 1 ) {
          fprintf(stderr, "%s: invalid option -- <rho> must be in range [0, 1]!\n", argv[0]);
          GenLinChenChenHelp(argv[0]);
          exit(1);
        }
        break;
      }
      case 'b': {
        bflag = true;
        break;
      }
      case 'h': {
        GenLinChenChenHelp(argv[0]);
        exit(0);
      }
      default: {
        GenLinChenChenHelp(argv[0]);
        exit(1);
      }
    }
  }

  // Create Beta
  if ( bflag ) {
    r = argc - optind;
    Beta = new float[r];
    for ( auto i = 0; i < r; ++i ) {
      Beta[i] = atof(argv[i+optind]);
    }
  } else {
    Beta = new float[r];
    float Beta_temp[8] = {0.7f, 0.9f, 0.4f, 0.3f, 1.0f, 0.2f, 0.2f, 0.1f};
    for ( auto i = 0; i < r; i++ ) {
      Beta[i] = Beta_temp[i%8];
    }
  }

  // Display arguments
  printf("================================================================"
         "================================================================\n");
  printf("n = %d, p = %d, r = %d, type = %d, rho = %.3f\n", n, p, r, type, rho);
  printf("Beta: ");
  for ( auto i = 0; i < r; i++ ) {
    printf("%8.3f", Beta[i]);
  }
  printf("\n\n");

  // ======== Create data ================================================================================================= //

  printf("Creating a linear data using Chen and Chen's method... ");
  fflush(stdout);

  // Allocate memory
  X = new float[n*p];
  Y = new float[n];
  J = new bool[p];
  auto L = new float[p*p];

  // Generate X & Y using normal random
  LAPACKE_slarnv(3, iseed, n*p, X);
  LAPACKE_slarnv(3, iseed, n, Y);

  // Create covariance matrix
  for ( auto i = 0; i < p; ++i ) {
    L[i+i*p] = 1.0f;
  }
  switch ( type ) {
    case 1: {
      for ( auto j = 0; j < p; ++j ) {
        for ( auto i = j+1; i < p; ++i ) {
          L[i+j*p] = rho;
        }
      }
      suffix = "_1";
      break;
    }
    case 2: {
      for ( auto i = 1; i < p; ++i ) {
        L[i+(i-1)*p] = rho;
      }
      suffix = "_2";
      break;
    }
    case 3: {
      for ( auto j = 0; j < p; ++j ) {
        for ( auto i = j+1; i < p; ++i ) {
          L[i+j*p] = pow(rho, i-j);
        }
      }
      suffix = "_3";
      break;
    }
  }

  // Compute the Cholesky factorization of L
  auto info = LAPACKE_spotrf(LAPACK_COL_MAJOR, 'L', p, L, p);
  if ( info ) {
    printf("Failed.\nThe covariance matrix is illegal.\n");
    abort();
  }

  // X := X * L'
  cblas_strmm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, n, p, 1.0f, L, p, X, n);

  // Permutate X
  if ( type != 1 ) {
    for ( auto i = 1; i < p-1; ++i ) {
      auto j = rand() % i;
      cblas_sswap(n, X+i*n, 1, X+j*n, 1);
    }
  }

  // Y += X[0~r cols] * Beta
  cblas_sgemv(CblasColMajor, CblasNoTrans, n, r, 1.0f, X, n, Beta, 1, 1.0f, Y, 1);

  // Generate J
  memset(J, false, sizeof(bool) * p);
  memset(J, true, sizeof(bool) * r);

  printf("Done.\n");

  // ====================================================================================================================== //

  // Save data
  GenLinChenChenSave(dataroot);

  // Free memory
  delete[] X;
  delete[] Y;
  delete[] Beta;
  delete[] J;
  delete[] L;

  printf("================================================================"
         "================================================================\n");

  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Display help messages
///
/// @param  cmd  the command name
///
void GenLinChenChenHelp( const char *cmd ) {
  printf("Usage: %s [options] ...\n", cmd);
  printf("\n%-32s%-40s%s\n\n", "Option",                   "Detail",                                 "Defalut Value");
  printf("%-32s%-40s%s\n",     "-f <file>, --file <file>", "save data into <file>",                  kDataRoot);
  printf("%-32s%-40s%s\n",     "-m <name>, --name <name>", "set the data name as <name>",            kDataName);
  printf("%-32s%-40s\n",       "-b [beta], --beta [beta]", "set the effects as [beta]s");
  printf("%-32s%-40s%d\n",     "-n ###",                   "the number of statistical units",        kN);
  printf("%-32s%-40s%d\n",     "-p ###",                   "the number of total effects" ,           kP);
  printf("%-32s%-40s%d\n",     "-r ###",                   "the number of given effects",            kR);
  printf("%-32s%-40s\n",       "",                         "ignored if '-b' is set");
  printf("%-32s%-40s%d\n",     "-t ###, --type ###",       "the type of covariance structure (1~3)", kType);
  printf("%-32s%-40s%.2f\n",   "-c ###, --cov ###",        "the covariance parameter",               kRho);
  printf("%-32s%-40s\n",       "-h, --help",               "display help messages");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Save data info file
///
/// @param  fileroot  the root of data file
///
void GenLinChenChenSave( const char *fileroot ) {
  FILE *file;

  printf("Saving data into '%s'... ", fileroot);
  fflush(stdout);

  // Open file
  file = fopen(fileroot, "wb");
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

  fprintf(file, "%s%s\n", dataname, suffix);
  fprintf(file, "%d %d\n", n, p);

  fprintf(file, "%-16c", '*');
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
