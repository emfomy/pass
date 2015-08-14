////////////////////////////////////////////////////////////////////////////////
// Particle Swarm Stepwise (PaSS) Algorithm                                   //
//                                                                            //
// pass.cpp                                                                   //
// The PaSS algorithm for Linear Regression                                   //
//                                                                            //
// Author: emfo<emfomy@gmail.com>                                             //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Notation:                                                                  //
// X:    the regressors                                                       //
// Y:    the regressand                                                       //
// Beta: the effects                                                          //
// R:    the residual                                                         //
//                                                                            //
// Linear model:                                                              //
// Y = X * Beta + error                                                       //
// R = Y - X * Beta                                                           //
//                                                                            //
// Select index in forward step:                                              //
// idx = argmax_{i not in I} abs( X[i col]' * R  )                            //
//                                                                            //
// Select index in backward step:                                             //
// idx = argmin_{i in I} norm( R_{I exclude i} )                              //
//     = argmin_{i in I} norm( R + Beta[i] * X[i col] )                       //
//                                                                            //
// References:                                                                //
// Chen, R.-B., Huang, C.-C., & Wang, W. (2013).                              //
//   Particle Swarm Stepwise (PaSS) Algorithm for Variable Selection.         //
////////////////////////////////////////////////////////////////////////////////

#include "pass.hpp"
#include <cstdlib>
#include <essl.h>
#include <omp.h>

// The log-binomial function
static inline lbinom( const int n, const int k ) {
  int i;
  return (lgammaf_r(n+1, &i) - lgammaf_r(n-k+1, &i) - lgammaf_r(k+1, &i));
}

////////////////////////////////////////////////////////////////////////////////
// The namespace pass                                                         //
////////////////////////////////////////////////////////////////////////////////
namespace pass {
int n;                // scalar, the number of statistical units
int p;                // scalar, the number of total effects
float* X0;            // matrix, n by p, the regressors
float* Y0;            // vector, n by 1, the regressand
bool* I0;             // vector, 1 by p, the chosen indices
float phi0;           // scalar, the value given by criterion
Parameter parameter;  // the parameters
int world_size;       // the size of MPI communicator
int world_rank;       // the rank of MPI process

////////////////////////////////////////////////////////////////////////////////
// The PaSS algorithm for Linear Regression                                   //
//                                                                            //
// Input Parameters:                                                          //
// n:         scalar, the number of statistical units                         //
// p:         scalar, the number of total effects                             //
// X0:        matrix, n by p, the regressors                                  //
// Y0:        vector, n by 1, the regressand                                  //
// parameter: the parameters                                                  //
//                                                                            //
// Output Global Variables:                                                   //
// I0:        vector, 1 by p, the chosen indices                              //
// phi0:      scalar, the value given by criterion                            //
//                                                                            //
// Note:                                                                      //
// Please call srand before using this routine.                               //
////////////////////////////////////////////////////////////////////////////////
void GenLin() {
  // Check parameters
  auto num_thread = omp_get_max_threads();

  ////////////////////////////////////////////////////////////////////////////
  // Centralize and normalize the original data                             //
  ////////////////////////////////////////////////////////////////////////////

  if ( !parameter.is_normalized ) {
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
  }

  ////////////////////////////////////////////////////////////////////////////
  // Run PaSS                                                               //
  ////////////////////////////////////////////////////////////////////////////

  // Allocate particles
  auto particle = new Particle[num_thread];
  phi0 = INFINITY;

  // Use openMP parallel
  #pragma omp parallel
  {
    auto tid = omp_get_thread_num();

    // Initialize particles
    particle[tid].InitializeModel();
    particle[tid].ComputeCriterion();
    particle[tid].phi_old = particle[tid].phi;

    // Copy best model
    if ( phi0 > particle[tid].phi ) {
      phi0 = particle[tid].phi;
      memcpy(I0, particle[tid].I, sizeof(bool) * p);
    }

    #pragma omp barrier

    // Find best model
    for ( auto i = 1; i < parameter.num_iteration; ++i ) {
      // Update model
      int idx;
      particle[tid].SelectIndex(idx);
      particle[tid].UpdateModel(idx);
      particle[tid].ComputeCriterion();

      // Check singularity
      if ( isnan(particle[tid].phi) ) {
        particle[tid].InitializeModel();
        particle[tid].ComputeCriterion();
      }

      // Change status
      if ( particle[tid].phi > particle[tid].phi_old ) {
        particle[tid].status = !particle[tid].status;
      }
      if ( particle[tid].k <= 1 ) {
        particle[tid].status = true;
      }
      if ( particle[tid].k >= n || particle[tid].k >= p-4 ) {
        particle[tid].status = false;
      }

      particle[tid].phi_old = particle[tid].phi;

      // Copy best model
      if ( phi0 > particle[tid].phi ) {
        phi0 = particle[tid].phi;
        memcpy(I0, particle[tid].I, sizeof(bool) * p);
      }
    }
  }

  // Copy best model
  memcpy(I0, I0, sizeof(bool) * p);

  ////////////////////////////////////////////////////////////////////////////

  delete[] particle;
}

////////////////////////////////////////////////////////////////////////////////
// The constructor of Particle                                                //
////////////////////////////////////////////////////////////////////////////////
Particle::Particle() {
  X        = new float[n*n];
  Y        = new float[n];
  Beta     = new float[n];
  Theta    = new float[n];
  M        = new float[n*n];
  R        = new float[n];
  B        = new float[n];
  D        = new float[n];

  Idx_lf   = new int[n];
  Idx_fl   = new int[p];
  Idx_temp = new int[p];
  I        = new bool[p];

  iseed    = rand();
}

////////////////////////////////////////////////////////////////////////////////
// The destructor of Particle                                                 //
////////////////////////////////////////////////////////////////////////////////
Particle::~Particle() {
  delete[] X;
  delete[] Y;
  delete[] Beta;
  delete[] Theta;
  delete[] M;
  delete[] B;
  delete[] D;
  delete[] R;

  delete[] Idx_lf;
  delete[] Idx_fl;
  delete[] Idx_temp;
  delete[] I;
}

////////////////////////////////////////////////////////////////////////////////
// Initialize the model of Particle randomly                                  //
////////////////////////////////////////////////////////////////////////////////
void Particle::InitializeModel() {
  InitializeModel(rand_r(&iseed) % p);
}

////////////////////////////////////////////////////////////////////////////////
// Initialize the model of Particle                                           //
//                                                                            //
// Parameters:                                                                //
// idx: the index of the effect                                               //
////////////////////////////////////////////////////////////////////////////////
void Particle::InitializeModel( const int idx ) {
  // Initialize size
  k = 1;

  // Initialize index
  memset(I, false, sizeof(bool) * p);

  // Update index
  I[idx] = true;
  Idx_lf[0] = idx;
  Idx_fl[idx] = 0;

  // X := X0[idx col]
  scopy(n, X0+idx*n, 1, X, 1);

  // Y := Y0
  scopy(n, Y0, 1, Y, 1);

  // M := 1 / (X' * X)
  M[0] = 1.0f;

  // Theta := X' * Y
  Theta[0] = sdot(n, X, 1, Y, 1);

  // Beta := M * Theta
  Beta[0] = M[0] * Theta[0];

  // R = Y - X * Beta
  szaxpy(n, -Beta[0], X, 1, Y, 1, R, 1);

  // Set status
  status = true;
}

////////////////////////////////////////////////////////////////////////////////
// Update the model of Particle                                               //
//                                                                            //
// Parameters:                                                                //
// idx: the index of the effect                                               //
////////////////////////////////////////////////////////////////////////////////
void Particle::UpdateModel( const int idx ) {
  if ( status ) {  // forward step
    // Update size
    auto km = k;
    k++;

    // Update index
    I[idx] = true;
    Idx_lf[km] = idx;
    Idx_fl[idx] = km;

    // Set Xnew
    auto Xnew = X+km*n;

    // Insert new row of X
    scopy(n, X0+idx*n, 1, Xnew, 1);

    ////////////////////////////////////////////////////////////////////////
    // Solve Y = X * Beta for Beta                                        //
    ////////////////////////////////////////////////////////////////////////

    // B := X' * Xnew
    sgemv("T", n, km, 1, X, n, Xnew, 1, 0.0f, B, 1);

    // D := M * B
    ssymv("U", km, 1.0f, M, n, B, 1, 0.0f, D, 1);

    // a := 1 / (Xnew' * Xnew - B' * D)
    auto a = 1.0f / (1.0f - sdot(km, B, 1, D, 1));

    // insert D by -1.0
    D[km] = -1.0f;

    // insert M by zeros
    for ( auto i = 0; i < k; i++ ) {
      M[km*n+i] = 0.0f;
    }

    // M += a * D * D'
    ssyr("U", k, a, D, 1, M, n);

    // insert Theta by Xnew' * Y
    Theta[km] = sdot(n, Xnew, 1, Y, 1);

    // insert Beta by zero
    Beta[km] = 0.0f;

    // Beta += a * (D' * Theta) * D
    saxpy(k, a*sdot(k, D, 1, Theta, 1), D, 1, Beta, 1);

    ////////////////////////////////////////////////////////////////////////

  } else {  // backward step
    // Update size
    k--;

    // Update index
    I[idx] = false;

    // Find index
    auto j = Idx_fl[idx];

    // Copy index end to index j
    // a := M[end, end], b := Beta[end], D := M[end col]
    float a, b;
    if ( j != k ) {
      a = M[j*n+j];
      b = Beta[j];

      scopy(n, X+k*n, 1, X+j*n, 1);
      Beta[j] = Beta[k];
      Theta[j] = Theta[k];
      Idx_lf[j] = Idx_lf[k];
      Idx_fl[Idx_lf[j]] = j;

      scopy(j, M+j*n, 1, D, 1);
      scopy(k-j-1, M+j*n+n+j, n, D+j+1, 1);
      D[j] = M[k*n+j];

      scopy(j, M+k*n, 1, M+j*n, 1);
      scopy(k-j-1, M+k*n+j+1, 1, M+j*n+n+j, n);
      M[j*n+j] = M[k*n+k];
    } else {
      a = M[k*n+k];
      b = Beta[k];
      scopy(k, M+k*n, 1, D, 1);
    }

    ////////////////////////////////////////////////////////////////////////
    // Solve Y = X * Beta for Beta                                        //
    ////////////////////////////////////////////////////////////////////////

    // M -= 1/a * D * D'
    ssyr("U", k, -1.0f/a, D, 1, M, n);

    // Beta -= b/a * D
    saxpy(k, -b/a, D, 1, Beta, 1);

    ////////////////////////////////////////////////////////////////////////
  }
  
  // R = Y - X * Beta
  scopy(n, Y, 1, R, 1);
  sgemv("N", n, k, -1.0f, X, n, Beta, 1, 1.0f, R, 1);
}

////////////////////////////////////////////////////////////////////////////////
// Select index of the effect to update                                       //
//                                                                            //
// Output Parameters:                                                         //
// idx: the index of the effect                                               //
////////////////////////////////////////////////////////////////////////////////
void Particle::SelectIndex( int& idx ) {
  auto frand = static_cast<float>(rand_r(&iseed)) / RAND_MAX;

  if ( status ) {  // Forward step
    // Idx_temp[0~itemp] := I0 exclude I
    // Idx_temp[0~(p-k)] := complement of I
    int itemp = 0;
    for ( auto i = 0, j = p-k; i < p; ++i ) {
      if ( !I[i] ) {
        if ( I0[i] ) {
          Idx_temp[itemp] = i;
          itemp++;
        } else {
          j--;
          Idx_temp[j] = i;
        }
      }
    }

    int choose;
    if ( itemp ) {
      choose = (frand < parameter.prob_forward_global)
             + (frand < parameter.prob_forward_global+
                        parameter.prob_forward_local);
    } else {
      choose = frand < parameter.prob_forward_local
                     / (parameter.prob_forward_local+
                        parameter.prob_forward_random);
    }

    switch( choose ) {
      case 2: {  // Global best
        idx = Idx_temp[rand_r(&iseed) % itemp];
        break;
      }
      case 1: {  // Local best
        auto phi_temp = -INFINITY;
        for ( auto i = 0; i < p; ++i ) {
          if ( !I[i] ) {
            // ftemp := abs( X0[i col]' * R )
            auto ftemp = fabs(sdot(n, X0+i*n, 1, R, 1));

            // Check if this value is maximum
            if ( phi_temp < ftemp ) {
              phi_temp = ftemp;
              idx = i;
            }
          }
        }
        break;
      }
      case 0: {  // Random
        idx = Idx_temp[rand_r(&iseed) % (p-k)];
        break;
      }
    }
  } else {  // backward step
    if ( frand < parameter.prob_backward_local ) {  // Local best
      auto phi_temp = INFINITY;
      for ( auto i = 0; i < k; ++i ) {
        // B := R + Beta[i] * X[i col]
        szaxpy(n, Beta[i], X+i*n, 1, R, 1, B, 1);

        // ftemp = norm(B)
        auto ftemp = snorm2(n, B, 1);

        // Check if this value is minimal
        if ( phi_temp > ftemp ) {
          phi_temp = ftemp;
          idx = Idx_lf[i];
        }
      }
    } else {  // Random
      idx = Idx_lf[rand_r(&iseed) % k];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Compute the value given by criterion                                       //
////////////////////////////////////////////////////////////////////////////////
void Particle::ComputeCriterion() {
  // e := norm(R)
  e = snorm2(n, R, 1);

  // Compute criterion
  switch(parameter.criterion) {
    case AIC: {    // phi := n*log(e^2/n) + 2k
      phi = n*logf(e*e/n) + 2.0f*k;
      break;
    }
    case BIC: {    // phi := n*log(e^2/n) + k*log(n)
      phi = n*logf(e*e/n) + k*logf(n);
      break;
    }
    case EBIC: {   // phi := n*log(e^2/n) + k*log(n) + 2gamma*log(p choose k)
      phi = n*logf(e*e/n) + k*logf(n) + 2.0f*parameter.ebic_gamma*lbinom(p, k);
      break;
    }
    case HDBIC: {  // phi := n*log(e^2/n) + k*log(n)*log(p)
      phi = n*logf(e*e/n) + k*logf(n)*logf(p);
      break;
    }
    case HQC: {    // phi := n*log(e^2/n) + 2k*log(log(n))
      phi = n*logf(e*e/n) + 2.0f*k*logf(logf(n));
      break;
    }
    case HDHQC: {  // phi := n*log(e^2/n) + 2k*log(log(n))*log(p)
      phi = n*logf(e*e/n) + 2.0f*k*logf(logf(n))*logf(p);
      break;
    }
  }
}

}
