////////////////////////////////////////////////////////////////////////////////
// Particle Swarm Stepwise (PaSS) Algorithm                                   //
//                                                                            //
// genlin_chenchen.cpp                                                        //
// Create a general linear regression data using Chen and Chen's method       //
//                                                                            //
// Author: emfo<emfomy@gmail.com>                                             //
//                                                                            //
// Reference:                                                                 //
// Chen, J., & Chen, Z. (2008). Extended Bayesian information criteria for    //
//   model selection with large model spaces. Biometrika, 95(3), 759–771.     //
//   http://www.stat.ubc.ca/~jhchen/paper/Bio08.pdf                           //
////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <getopt.h>
#include <mkl.h>

// Default arguments
const int kN          = 200;
const int kP          = 50;
const int kR          = 8;
const int kType       = 3;
const float kRho      = 0.2f;
const char *kDataRoot = "genlin.dat";
const char *kDataName = "GenLin_ChenChen";

// Global variables
int n;                 // scalar, the number of statistical units
int p;                 // scalar, the number of total effects
int r;                 // scalar, the number of given effects
int type;              // scalar, the type of covariance structure, 1~3
float *X;              // matrix, n by p, the regressors
float *Y;              // vector, n by 1, the regressand
float *Beta;           // vector, r by 1, the effects
float rho;             // scaler, the covariance parameter
bool *J;               // vector, 1 by p, the chosen indices
const char *dataroot;  // string, the root of the data file
const char *dataname;  // string, the name of data
const char *suffix;    // string, the suffix of the name of data

// Functions
void ChenChenHelp( const char *cmd );
void ChenChenSave( const char *fileroot );

////////////////////////////////////////////////////////////////////////////////
// Display help messages                                                      //
////////////////////////////////////////////////////////////////////////////////
void IngLaiHelp( const char *cmd ) {
  printf("Usage: %s [options] ...\n", cmd);
  printf("\n%-32s%-40s%s\n\n", "Option", "Detail", "Defalut Value");
  printf("%-32s%-40s%s\n",
         "-f <file>, --file <file>", "save data into <file>", kDataRoot);
  printf("%-32s%-40s%s\n",
         "-m <name>, --name <name>", "set the data name as <name>", kDataName);
  printf("%-32s%-40s\n",
         "-b [beta], --beta [beta]", "set the effects as [beta]s");
  printf("%-32s%-40s%d\n",
         "-n ###", "the number of statistical units" , kN);
  printf("%-32s%-40s%d\n",
         "-p ###", "the number of total effects" , kP);
  printf("%-32s%-40s%d\n",
         "-r ###", "the number of given effects", kR);
  printf("%-32s%-40s\n",
         "", "ignored if '-b' is set");
  printf("%-32s%-40s%d\n",
         "-t ###, --type ###",
         "the type of covariance structure (1~3)", kType);
  printf("%-32s%-40s%.2f\n",
         "-c ###, --cov ###", "the covariance parameter", kRho);
  printf("%-32s%-40s\n",
         "-h, --help", "display help messages");
}

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
          fprintf(stderr, "%s: invalid option -- "
                 "<n> must be a positive integer!\n", argv[0]);
          IngLaiHelp(argv[0]);
          exit(1);
        }
        break;
      }
      case 'p': {
        p = atoi(optarg);
        if ( p <= 0 ) {
          fprintf(stderr, "%s: invalid option -- "
                 "<p> must be a positive integer!\n", argv[0]);
          IngLaiHelp(argv[0]);
          exit(1);
        }
        break;
      }
      case 'r': {
        r = atoi(optarg);
        if ( r < 0 ) {
          fprintf(stderr, "%s: invalid option -- "
                 "<r> must be a non-negative integer!\n", argv[0]);
          IngLaiHelp(argv[0]);
          exit(1);
        }
        break;
      }
      case 't': {
        type = atoi(optarg);
        if ( type < 1 || type > 3 ) {
          fprintf(stderr, "%s: invalid option -- "
                 "<type> must be 1, 2, or 3!\n", argv[0]);
          IngLaiHelp(argv[0]);
          exit(1);
        }
        break;
      }
      case 'c': {
        rho = atof(optarg);
        if ( rho < 0 || rho > 1 ) {
          fprintf(stderr, "%s: invalid option -- "
                 "<rho> must be in range [0, 1]!\n", argv[0]);
          IngLaiHelp(argv[0]);
          exit(1);
        }
        break;
      }
      case 'b': {
        bflag = true;
        break;
      }
      case 'h': {
        IngLaiHelp(argv[0]);
        exit(0);
      }
      default: {
        IngLaiHelp(argv[0]);
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

  ////////////////////////////////////////////////////////////////////////////
  // Create data                                                            //
  ////////////////////////////////////////////////////////////////////////////

  printf("Creating a linear data using Chen and Chen's method... ");

  // Allocate memory
  X = new float[n*p];
  Y = new float[n];
  J = new bool[p]();
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
  cblas_strmm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit,
              n, p, 1.0f, L, p, X, n);

  // Permutate X
  if ( type != 1 ) {
    for ( auto i = 1; i < p-1; ++i ) {
      auto j = rand() % i;
      cblas_sswap(n, X+i*n, 1, X+j*n, 1);
    }
  }

  // Y += X[0~r cols] * Beta
  cblas_sgemv(CblasColMajor, CblasNoTrans,
              n, r, 1.0f, X, n, Beta, 1, 1.0f, Y, 1);

  // Generate J
  memset(J, true, sizeof(bool) * r);

  printf("Done.\n");

  ////////////////////////////////////////////////////////////////////////////

  // Save data
  ChenChenSave(dataroot);

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

////////////////////////////////////////////////////////////////////////////////
// Load parameters from config file                                           //
//                                                                            //
// Parameters:                                                                //
// fileroot: the root of config file                                          //
////////////////////////////////////////////////////////////////////////////////
void ChenChenConfig( const char *fileroot ) {
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
    fgets(line, kBufferSize, file);
    sscanf(line, "%*s %d", &type);
    fgets(line, kBufferSize, file);
    sscanf(line, "%*s %f", &rho);

    // Close file
    fclose(file);

    printf("Done.\n");
  } else {
    printf("Failed!\nCreating config file '%s'... ", fileroot);

    // Open file
    file = fopen(fileroot, "w");
    if ( !file ) {
      printf("Failed!\n");
      abort();
    }
  
    // Generate Beta
    Beta = new float[r];
    float Beta_temp[8] = {0.7f, 0.9f, 0.4f, 0.3f, 1.0f, 0.2f, 0.2f, 0.1f};
    for ( auto i = 0; i < r; i++ ) {
      Beta[i] = Beta_temp[i%8];
    }

    // Write data
    fprintf(file, "n/p/r %d %d %d\n", n, p, r);
    fprintf(file, "Beta ");
    for ( auto i = 0; i < r; i++ ) {
      fprintf(file, " %.3f", Beta[i]);
    }
    fprintf(file, "\n");
    fprintf(file, "type  %d\n", type);
    fprintf(file, "rho   %.3f\n", rho);

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
void ChenChenSave( const char *fileroot ) {
  FILE *file;
  int size0 = strlen(dataname);
  int size1 = strlen(suffix)+1;
  int size = size0+size1;

  printf("Saving data into '%s'... ", fileroot);

  // Open file
  file = fopen(fileroot, "wb");
  if ( !file ) {
    printf("Failed!\n");
    abort();
  }

  // Write data
  fwrite(&size, sizeof(int), 1, file);
  fwrite(dataname, sizeof(char), size0, file);
  fwrite(suffix, sizeof(char), size1, file);
  fwrite(&n, sizeof(int), 1, file);
  fwrite(&p, sizeof(int), 1, file);
  fwrite(X, sizeof(float), n * p, file);
  fwrite(Y, sizeof(float), n, file);
  fwrite(J, sizeof(bool), p, file);

  // Close file
  fclose(file);

  printf("Done.\n");
}
