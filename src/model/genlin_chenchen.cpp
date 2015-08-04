////////////////////////////////////////////////////////////////////////////////
// Particle Swarm Stepwise (PaSS) Algorithm                                   //
//                                                                            //
// genlin_chenchen.cpp                                                        //
// Create a linear data using Chen and Chen's method                          //
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
#include <cblas.h>
#include <lapacke.h>

// Global variables
int n;           // scalar, the number of statistical units
int p;           // scalar, the number of total effects
int r;           // scalar, the number of given effects
int type;        // scalar, the type of covariance structures, 1~3
float *X;        // matrix, n by p, the regressors
float *Y;        // vector, n by 1, the regressand
float *Beta;     // vector, r by 1, the effects
float rho;       // scaler, the covariance parameter
bool *J;         // vector, 1 by p, the chosen indices
char *dataname;  // string, the name of data

// Functions
void ChenChenConfig( const char* fileroot );
void ChenChenSave( const char* fileroot );

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
  n = 200;
  p = 50;
  r = 8;
  type = 3;
  rho = 0.2;

  auto cfgroot  = (argc > 1) ? argv[1] : "pass_genlin_chenchen.cfg";
  auto dataroot = (argc > 2) ? argv[2] : "pass_genlin.dat";

  auto strtemp = "GenLin_ChenChen";
  dataname = new char[strlen(strtemp)+3];
  memcpy(dataname, strtemp, sizeof(char) * (strlen(strtemp)+1));

  if ( argc > 1 ) cfgroot  = argv[1];
  if ( argc > 2 ) dataroot = argv[2];

  ////////////////////////////////////////////////////////////////////////////
  // Load parameters                                                        //
  ////////////////////////////////////////////////////////////////////////////

  printf("================================================================\n");

  // Load parameters
  ChenChenConfig(cfgroot);

  // Check parameters
  if ( type > 3 || type < 1 ) {
    printf("There is no type %d!\n", type);
    exit(1);
  }

  // Display parameters
  printf("\nn = %d, p = %d, r = %d, type = %d, rho = %.3f\n",
         n, p, r, type, rho);
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
  auto R = new float[p*p];

  // Generate X & Y using normal random
  LAPACKE_slarnv(3, iseed, n*p, X);
  LAPACKE_slarnv(3, iseed, n, Y);

  // Create covariance matrix
  for ( auto i = 0; i < p; ++i ) {
    R[i+i*p] = 1.0f;
  }
  switch ( type ) {
    case 1: {
      for ( auto j = 0; j < p; ++j ) {
        for ( auto i = j+1; i < p; ++i ) {
          R[i+j*p] = rho;
        }
      }
      strcat(dataname, "_1");
      break;
    }
    case 2: {
      for ( auto i = 1; i < p; ++i ) {
        R[i+(i-1)*p] = rho;
      }
      strcat(dataname, "_2");
      break;
    }
    case 3: {
      for ( auto j = 0; j < p; ++j ) {
        for ( auto i = j+1; i < p; ++i ) {
          R[i+j*p] = pow(rho, i-j);
        }
      }
      strcat(dataname, "_3");
      break;
    }
  }

  // Compute the Cholesky factorization of R
  auto info = LAPACKE_spotrf(LAPACK_COL_MAJOR, 'L', p, R, p);
  if ( info ) {
    printf("Failed.\nThe covariance matrix is illegal.\n");
    exit(1);
  }

  // X := X * R
  cblas_strmm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit,
              n, p, 1.0f, R, p, X, n);

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
  delete[] R;

  printf("================================================================\n");

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Load parameters from config file                                           //
//                                                                            //
// Parameters:                                                                //
// fileroot: the root of config file                                          //
////////////////////////////////////////////////////////////////////////////////
void ChenChenConfig( const char* fileroot ) {
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
void ChenChenSave( const char* fileroot ) {
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
